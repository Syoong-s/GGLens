CXX = mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
PYTHON_CONFIG ?= python3-config

TARGET = esd_bcg_mpi
BUILD_DIR = build
PLOT_OBJECT = $(BUILD_DIR)/esd_plot.o
PLOT_LIB = $(BUILD_DIR)/libesd_plot.a
SOURCES = $(filter-out src/esd_plot.cpp,$(wildcard src/*.cpp))
OBJECTS = $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

PYTHON_PREFIX = $(shell $(PYTHON_CONFIG) --exec-prefix)
PYTHON_INCLUDES = $(shell $(PYTHON_CONFIG) --includes)
PYTHON_LDFLAGS = $(shell $(PYTHON_CONFIG) --ldflags --embed) -Wl,-rpath,$(PYTHON_PREFIX)/lib

all: $(TARGET)

$(TARGET): $(OBJECTS) $(PLOT_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(PLOT_LIB) $(PYTHON_LDFLAGS) $(LDLIBS)

$(PLOT_LIB): $(PLOT_OBJECT)
	$(AR) rcs $@ $<

$(PLOT_OBJECT): src/esd_plot.cpp $(BUILD_DIR)/esd_plot_source.inc Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude -I$(BUILD_DIR) $(PYTHON_INCLUDES) -DFQ_PLOT_PYTHON='"$(PYTHON_PREFIX)/bin/python3"' -MMD -MP -c -o $@ $<

$(BUILD_DIR)/esd_plot_source.inc: script/esd_plot.py Makefile | $(BUILD_DIR)
	{ printf '%s\n' 'static constexpr char plot_source[] = R"FQPLOT('; cat $<; printf '%s\n' ')FQPLOT";'; } > $@.tmp
	mv $@.tmp $@

$(BUILD_DIR)/%.o: src/%.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude -MMD -MP -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -f $(TARGET) $(OBJECTS) $(PLOT_LIB) $(PLOT_OBJECT) $(OBJECTS:.o=.d) $(PLOT_OBJECT:.o=.d) $(BUILD_DIR)/esd_plot_source.inc $(BUILD_DIR)/esd_plot_source.inc.tmp
	@rmdir $(BUILD_DIR) 2>/dev/null || true

-include $(OBJECTS:.o=.d) $(PLOT_OBJECT:.o=.d)

.PHONY: all clean
