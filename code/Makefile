CXX = mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
PYTHON_CONFIG ?= python3-config
CXX_RUNTIME ?= auto
PATCHELF ?= auto

TARGET = esd_bcg_mpi
BUILD_DIR = build
FIT_OBJECT = $(BUILD_DIR)/esd_fit_backend.o
FIT_LIB = $(BUILD_DIR)/libesd_fit.a
FIT_MODULES = esd_nfw esd_two_halo esd_archive esd_profiles esd_plot esd_fit
FIT_SOURCES = $(addprefix $(BUILD_DIR)/,$(addsuffix _source.inc,$(FIT_MODULES)))
SOURCES = $(filter-out src/esd_fit_backend.cpp,$(wildcard src/*.cpp))
OBJECTS = $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
PYTHON_BUILD_SETTINGS = $(BUILD_DIR)/python-build.conf
CXX_RUNTIME_SCRIPT = script/select_cxx_runtime.sh

PYTHON_PREFIX = $(shell $(PYTHON_CONFIG) --exec-prefix)
PYTHON_INCLUDES = $(shell $(PYTHON_CONFIG) --includes)
PYTHON_LDFLAGS = $(shell $(PYTHON_CONFIG) --ldflags --embed) -Wl,-rpath,$(PYTHON_PREFIX)/lib

all: $(TARGET)

$(TARGET): $(OBJECTS) $(FIT_LIB) $(PYTHON_BUILD_SETTINGS) $(CXX_RUNTIME_SCRIPT)
	@set -eu; \
	  linked=$$(mktemp "$@.tmp.XXXXXX"); \
	  trap 'rm -f -- "$$linked"' EXIT; \
	  trap 'exit 1' HUP INT TERM; \
	  chmod +rw "$$linked"; \
	  $(CXX) $(CXXFLAGS) $(LDFLAGS) -o "$$linked" $(OBJECTS) $(FIT_LIB) $(PYTHON_LDFLAGS) $(LDLIBS); \
	  sh $(CXX_RUNTIME_SCRIPT) "$(CXX_RUNTIME)" "$(PYTHON_PREFIX)" "$(PATCHELF)" "$$linked"; \
	  mv -f -- "$$linked" "$@"

# Track Make overrides because changing the Python or runtime flags can require a rebuild.
$(PYTHON_BUILD_SETTINGS): FORCE $(CXX_RUNTIME_SCRIPT) | $(BUILD_DIR)
	@set -eu; \
	  settings=$$(mktemp "$@.tmp.XXXXXX"); \
	  trap 'rm -f -- "$$settings"' EXIT; \
	  trap 'exit 1' HUP INT TERM; \
	  { printf '%s\n' "$(PYTHON_PREFIX)" "$(PYTHON_INCLUDES)" "$(PYTHON_LDFLAGS)" "$(CXX_RUNTIME)" "$(PATCHELF)"; \
	    sh $(CXX_RUNTIME_SCRIPT) "$(CXX_RUNTIME)" "$(PYTHON_PREFIX)" "$(PATCHELF)"; \
	  } > "$$settings"; \
	  if ! cmp -s "$$settings" "$@"; then mv -f -- "$$settings" "$@"; fi

$(FIT_LIB): $(FIT_OBJECT)
	$(AR) rcs $@ $<

$(FIT_OBJECT): src/esd_fit_backend.cpp $(FIT_SOURCES) Makefile $(PYTHON_BUILD_SETTINGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude -I$(BUILD_DIR) $(PYTHON_INCLUDES) -DFQ_ESD_PYTHON='"$(PYTHON_PREFIX)/bin/python3"' -MMD -MP -c -o $@ $<

$(BUILD_DIR)/%_source.inc: script/%.py Makefile | $(BUILD_DIR)
	{ printf '%s\n' 'static constexpr char $*_source[] = R"FQPY('; cat $<; printf '%s\n' ')FQPY";'; } > $@.tmp
	mv $@.tmp $@

$(BUILD_DIR)/%.o: src/%.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude -MMD -MP -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -f $(TARGET) $(OBJECTS) $(FIT_LIB) $(FIT_OBJECT) $(OBJECTS:.o=.d) $(FIT_OBJECT:.o=.d) $(FIT_SOURCES) $(addsuffix .tmp,$(FIT_SOURCES)) $(PYTHON_BUILD_SETTINGS)
	@rmdir $(BUILD_DIR) 2>/dev/null || true

-include $(OBJECTS:.o=.d) $(FIT_OBJECT:.o=.d)

FORCE:

.PHONY: all clean FORCE
