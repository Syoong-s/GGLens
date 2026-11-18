CXX = mpicxx
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
PYTHON_CONFIG ?= python3-config
TARGET = esd_bcg_mpi
PLOT_LIB = libesd_plot.a
SOURCES = selfcalibration.cpp selfcal_pdf.cpp esd_bcg_mpi.cpp esd_lenses.cpp esd_regions.cpp esd_config.cpp esd_pairs.cpp esd_measure.cpp bcg_shear_catalog.cpp bcg_stream.cpp bcg_pair_store.cpp bcg_quantiles.cpp bcg_pdf_stream.cpp bcg_star_bar.cpp
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = $(wildcard *.hpp)

PYTHON_PREFIX = $(shell $(PYTHON_CONFIG) --exec-prefix)
PYTHON_INCLUDES = $(shell $(PYTHON_CONFIG) --includes)
PYTHON_LDFLAGS = $(shell $(PYTHON_CONFIG) --ldflags --embed) -Wl,-rpath,$(PYTHON_PREFIX)/lib

all: $(TARGET)

$(TARGET): $(OBJECTS) $(PLOT_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(PLOT_LIB) $(PYTHON_LDFLAGS) $(LDLIBS)

$(PLOT_LIB): esd_plot.o
	$(AR) rcs $@ $<

esd_plot.o: esd_plot.cpp esd_plot_source.inc $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PYTHON_INCLUDES) -DFQ_PLOT_PYTHON='"$(PYTHON_PREFIX)/bin/python3"' -c -o $@ $<

esd_plot_source.inc: esd_plot.py Makefile
	{ printf '%s\n' 'static constexpr char plot_source[] = R"FQPLOT('; cat $<; printf '%s\n' ')FQPLOT";'; } > $@.tmp
	mv $@.tmp $@

%.o: %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) $(OBJECTS) $(PLOT_LIB) esd_plot.o esd_plot_source.inc esd_plot_source.inc.tmp

.PHONY: all clean
