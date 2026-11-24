include Makefile

PYTHON ?= python3
MPIEXEC ?= mpiexec
MPIEXEC_FLAGS ?=

CATALOG_TEST = $(BUILD_DIR)/test_catalogs
CATALOG_OBJECTS = $(filter $(BUILD_DIR)/bcg_dr9% $(BUILD_DIR)/esd_config.o $(BUILD_DIR)/esd_lenses.o,$(OBJECTS))

$(CATALOG_TEST): tests/catalogs.cpp $(CATALOG_OBJECTS) Makefile GNUmakefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -Iinclude -o $@ $< $(CATALOG_OBJECTS) $(LDLIBS)

test-catalogs: $(CATALOG_TEST)
	$(PYTHON) tests/catalogs.py $(CATALOG_TEST) $(MPIEXEC) $(MPIEXEC_FLAGS)

test-stages: $(TARGET)
	$(PYTHON) tests/stages.py $(TARGET) $(MPIEXEC) $(MPIEXEC_FLAGS)

clean: clean-tests

clean-tests:
	rm -f $(CATALOG_TEST)

.PHONY: test-catalogs test-stages clean-tests
