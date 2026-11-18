#include "bcg_quantiles.hpp"
#include <queue>

namespace bcg {
namespace {
struct AbsRecord { double value = 0; std::int32_t region = 0, padding = 0; };
static_assert(sizeof(AbsRecord) == 16, "Magnitude sort records must use sixteen bytes");

// ==========================================
// Class: QuantileHeads
// Method: Derive every exact omitted-region quantile from one sorted magnitude scan.
// ==========================================
class QuantileHeads {
    struct Head { std::int64_t key = 0, dn = 0; int slot = 0, boundary = 0; };
    std::vector<Head> heap_;
    std::vector<int> position_;
    InitialQuantiles &result_;
    int boundaries_;
    std::int64_t ordinal_ = 0;
    // ==========================================
    // Function: exchange
    // Method: Keep replica-to-heap-position indices synchronized.
    // ==========================================
    void exchange(std::size_t a, std::size_t b) {
        std::swap(heap_[a], heap_[b]);
        position_[heap_[a].slot] = static_cast<int>(a);
        position_[heap_[b].slot] = static_cast<int>(b);
    }
    // ==========================================
    // Function: down
    // Method: Restore the indexed minimum heap after a target index increases.
    // ==========================================
    void down(std::size_t i) {
        for (;;) {
            const auto left = 2 * i + 1;
            if (left >= heap_.size()) break;
            auto best = left;
            if (left + 1 < heap_.size() && heap_[left + 1].key < heap_[left].key) best = left + 1;
            if (heap_[i].key <= heap_[best].key) break;
            exchange(i, best); i = best;
        }
    }
public:
    // ==========================================
    // Function: QuantileHeads constructor
    // Method: Initialize floor(n/(B+1)) order statistics for all eligible replicas.
    // ==========================================
    QuantileHeads(const PairStore &store, int boundaries, int regions, InitialQuantiles &result)
        : position_(static_cast<std::size_t>(regions) + 1, -1), result_(result), boundaries_(boundaries) {
        result.values.resize(static_cast<std::size_t>(regions) + 1);
        for (int slot = 0; slot <= regions; ++slot) {
            const auto n = store.count - (slot ? store.region_count[slot - 1] : 0);
            if (n < static_cast<std::int64_t>(boundaries) + 1) continue;
            result.values[slot].assign(boundaries, 0);
            const auto dn = n / (static_cast<std::int64_t>(boundaries) + 1);
            position_[slot] = static_cast<int>(heap_.size());
            heap_.push_back({dn - 1, dn, slot, 0});
        }
        for (std::size_t i = heap_.size() / 2; i > 0; --i) down(i - 1);
    }
    // ==========================================
    // Function: consume
    // Method: Shift only the replica excluding this record, then resolve due targets.
    // ==========================================
    void consume(const AbsRecord &record) {
        const auto slot = static_cast<std::size_t>(record.region) + 1;
        const int index = position_.at(slot);
        if (index >= 0 && heap_[index].boundary < boundaries_) {
            ++heap_[index].key; down(static_cast<std::size_t>(index));
        }
        while (!heap_.empty() && heap_.front().key == ordinal_) {
            auto &h = heap_.front();
            result_.values[h.slot][h.boundary++] = record.value;
            if (h.boundary == boundaries_) h.key = std::numeric_limits<std::int64_t>::max();
            else h.key += h.dn;
            down(0);
        }
        ++ordinal_;
    }
    // ==========================================
    // Function: finish
    // Method: Reject missing quantile targets or an incomplete merged stream.
    // ==========================================
    void finish(std::int64_t expected) const {
        if (ordinal_ != expected) throw std::runtime_error("Magnitude-sort record count changed");
        for (const auto &h : heap_) if (h.boundary != boundaries_)
            throw std::runtime_error("Incomplete exact PDF quantile extraction");
    }
};

// ==========================================
// Function: merge_runs
// Method: Merge at most 32 sorted files with bounded per-file read buffers.
// ==========================================
void merge_runs(const std::vector<std::string> &files,
                const std::function<void(const AbsRecord &)> &visit) {
    struct Item { AbsRecord value; std::size_t file = 0; };
    struct Greater {
        bool operator()(const Item &a, const Item &b) const { return a.value.value > b.value.value; }
    };
    std::vector<std::unique_ptr<BinaryReader<AbsRecord>>> readers;
    std::priority_queue<Item, std::vector<Item>, Greater> heap;
    for (const auto &file : files) {
        readers.emplace_back(new BinaryReader<AbsRecord>(file, 64 * 1024));
        AbsRecord value;
        if (readers.back()->next(value)) heap.push({value, readers.size() - 1});
    }
    while (!heap.empty()) {
        const auto item = heap.top(); heap.pop(); visit(item.value);
        AbsRecord next;
        if (readers[item.file]->next(next)) heap.push({next, item.file});
    }
}
} // namespace

// ==========================================
// Function: prepare_quantiles
// Method: Bounded exact external sort; reuse one ordering for all jackknife replicas.
// ==========================================
InitialQuantiles prepare_quantiles(const PairStore &store, const std::vector<int> &labels,
                                  int component, int boundaries, int regions,
                                  const ScratchSpace &scratch, std::size_t sort_bytes) {
    InitialQuantiles result;
    QuantileHeads heads(store, boundaries, regions, result);
    if (!store.count) return result;
    const std::size_t capacity = std::max<std::size_t>(sort_bytes / sizeof(AbsRecord), 1);
    std::vector<AbsRecord> buffer; buffer.reserve(capacity);
    std::vector<std::vector<std::string>> levels;
    std::uint64_t serial = 0;
    const auto name = [&] { return scratch.path("magnitude_" + std::to_string(serial++) + ".bin"); };
    // ==========================================
    // Function: add_run
    // Method: Collapse full merge levels so path metadata grows only logarithmically.
    // ==========================================
    const auto add_run = [&](std::string file) {
        std::size_t level = 0;
        for (;;) {
            if (level == levels.size()) levels.emplace_back();
            levels[level].push_back(file);
            if (levels[level].size() < 32) break;
            file = name(); BinaryWriter<AbsRecord> out(file, 128 * 1024);
            merge_runs(levels[level], [&](const AbsRecord &v) { out.append(v); }); out.close();
            for (const auto &old : levels[level]) fs::remove(old);
            levels[level].clear(); ++level;
        }
    };
    // ==========================================
    // Function: flush_run
    // Method: Sort one bounded magnitude block and spill it before accepting more rows.
    // ==========================================
    const auto flush_run = [&] {
        std::sort(buffer.begin(), buffer.end(), [](const AbsRecord &a, const AbsRecord &b) {
            return a.value < b.value;
        });
        const auto file = name(); BinaryWriter<AbsRecord> out(file, 128 * 1024);
        for (const auto &value : buffer) out.append(value);
        out.close(); buffer.clear(); add_run(file);
    };
    store.scan([&](const PairRecord &pair) {
        const double value = std::abs(component ? pair.gx : pair.gt);
        if (!std::isfinite(value)) throw std::runtime_error("Nonfinite rotated shear in PDF input");
        buffer.push_back({value, labels.at(static_cast<std::size_t>(pair.lens)), 0});
        if (buffer.size() == capacity) flush_run();
    });
    if (levels.empty()) {
        std::sort(buffer.begin(), buffer.end(), [](const AbsRecord &a, const AbsRecord &b) {
            return a.value < b.value;
        });
        for (const auto &record : buffer) heads.consume(record);
    } else {
        if (!buffer.empty()) flush_run();
        std::vector<AbsRecord>().swap(buffer);
        std::vector<std::string> remaining;
        for (const auto &level : levels) remaining.insert(remaining.end(), level.begin(), level.end());
        while (remaining.size() > 32) {
            std::vector<std::string> next;
            for (std::size_t first = 0; first < remaining.size(); first += 32) {
                const auto last = std::min(first + 32, remaining.size());
                std::vector<std::string> group(remaining.begin() + first, remaining.begin() + last);
                const auto file = name(); BinaryWriter<AbsRecord> out(file, 128 * 1024);
                merge_runs(group, [&](const AbsRecord &v) { out.append(v); }); out.close();
                for (const auto &old : group) fs::remove(old);
                next.push_back(file);
            }
            remaining = std::move(next);
        }
        merge_runs(remaining, [&](const AbsRecord &value) { heads.consume(value); });
        for (const auto &file : remaining) fs::remove(file);
    }
    heads.finish(store.count); return result;
}
} // namespace bcg
