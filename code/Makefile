CXX = mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
PYTHON_CONFIG ?= python3-config
PYTHON ?= python3
MPIEXEC ?= mpiexec
MPIEXEC_FLAGS ?=

TARGET = esd_bcg_mpi
BUILD_DIR = build
FIT_OBJECT = $(BUILD_DIR)/esd_fit_backend.o
FIT_LIB = $(BUILD_DIR)/libesd_fit.a
FIT_MODULES = esd_nfw esd_two_halo esd_profiles esd_plot esd_fit
FIT_SOURCES = $(addprefix $(BUILD_DIR)/,$(addsuffix _source.inc,$(FIT_MODULES)))
SOURCES = $(filter-out src/esd_fit_backend.cpp,$(wildcard src/*.cpp))
OBJECTS = $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
CATALOG_TEST = $(BUILD_DIR)/test_catalogs
CATALOG_OBJECTS = $(filter $(BUILD_DIR)/bcg_dr9% $(BUILD_DIR)/esd_config.o $(BUILD_DIR)/esd_lenses.o,$(OBJECTS))

PYTHON_PREFIX = $(shell $(PYTHON_CONFIG) --exec-prefix)
PYTHON_INCLUDES = $(shell $(PYTHON_CONFIG) --includes)
PYTHON_LDFLAGS = $(shell $(PYTHON_CONFIG) --ldflags --embed) -Wl,-rpath,$(PYTHON_PREFIX)/lib

all: $(TARGET)

$(TARGET): $(OBJECTS) $(FIT_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(FIT_LIB) $(PYTHON_LDFLAGS) $(LDLIBS)

$(FIT_LIB): $(FIT_OBJECT)
	$(AR) rcs $@ $<

$(FIT_OBJECT): src/esd_fit_backend.cpp $(FIT_SOURCES) Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude -I$(BUILD_DIR) $(PYTHON_INCLUDES) -DFQ_ESD_PYTHON='"$(PYTHON_PREFIX)/bin/python3"' -MMD -MP -c -o $@ $<

$(BUILD_DIR)/%_source.inc: script/%.py Makefile | $(BUILD_DIR)
	{ printf '%s\n' 'static constexpr char $*_source[] = R"FQPY('; cat $<; printf '%s\n' ')FQPY";'; } > $@.tmp
	mv $@.tmp $@

$(BUILD_DIR)/%.o: src/%.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude -MMD -MP -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $@

$(CATALOG_TEST): tests/catalogs.cpp $(CATALOG_OBJECTS) Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -Iinclude -o $@ $< $(CATALOG_OBJECTS) $(LDLIBS)

test-catalogs: $(CATALOG_TEST)
	$(PYTHON) tests/catalogs.py $(CATALOG_TEST) $(MPIEXEC) $(MPIEXEC_FLAGS)

clean:
	rm -f $(TARGET) $(CATALOG_TEST) $(OBJECTS) $(FIT_LIB) $(FIT_OBJECT) $(OBJECTS:.o=.d) $(FIT_OBJECT:.o=.d) $(FIT_SOURCES) $(addsuffix .tmp,$(FIT_SOURCES))
	@rmdir $(BUILD_DIR) 2>/dev/null || true

-include $(OBJECTS:.o=.d) $(FIT_OBJECT:.o=.d)

.PHONY: all clean test-catalogs
