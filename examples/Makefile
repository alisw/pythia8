# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2025 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main101
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=$(OBJ_COMMON) -I$(PREFIX_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl
PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)

# Define RIVET options and fix C++ version, rpath, missing HDF5.
ifeq ($(RIVET_USE),true)
  COMMA=,
  RIVET_VERSION=$(shell $(RIVET_BIN)$(RIVET_CONFIG) --version)
  RIVET_LPATH=$(filter -L%,$(shell $(RIVET_BIN)$(RIVET_CONFIG) --ldflags))
  RIVET_FLAGS=$(subst -L,-Wl$(COMMA)-rpath$(COMMA),$(RIVET_LPATH))
  RIVET_FLAGS+= $(shell $(RIVET_BIN)$(RIVET_CONFIG) --cppflags --libs)
  RIVET_CSTD=c++14
  ifeq ("4.0.0","$(word 1, $(sort 4.0.0 $(RIVET_VERSION)))")
    RIVET_CSTD=c++17
    RIVET_LDIR=$(shell $(RIVET_BIN)$(RIVET_CONFIG) --libdir)
    RIVET_HDF5=$(shell nm $(RIVET_LDIR)/libRivet$(LIB_SUFFIX) | grep H5open)
    ifneq ($(strip $(RIVET_HDF5)),)
      RIVET_FLAGS+= -lhdf5
    endif
  endif
  RIVET_OPTS=$(CXX_COMMON:c++11=$(RIVET_CSTD)) $(RIVET_FLAGS)
  RIVET_OPTS+= $(CXX_DTAGS) -DRIVET
endif

# Define additional dependency options.
EVTGEN_OPTS=$(EVTGEN_INCLUDE) $(EVTGEN_LIB) -DEVTGEN_PYTHIA -DEVTGEN_EXTERNAL
HEPMC2_OPTS=$(HEPMC2_INCLUDE) $(HEPMC2_LIB) -DHEPMC2
HEPMC3_OPTS=$(HEPMC3_INCLUDE) $(HEPMC3_LIB) -DHEPMC3
HDF5_OPTS=$(HDF5_INCLUDE) $(HIGHFIVE_INCLUDE) $(MPICH_INCLUDE)
HDF5_OPTS+= $(HDF5_LIB) -DHDF5
FASTJET3_OPTS=$(FASTJET3_INCLUDE) $(FASTJET3_LIB)
ifeq ($(ROOT_USE),true)
  ROOT_OPTS=$(ROOT_LIB) $(shell $(ROOT_CONFIG) --cflags --glibs) -DPY8ROOT
endif

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	$(info Usage: make mainXXX)

# PYTHIA library.
$(PYTHIA):
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

# All programs without external dependencies.
%: $(PYTHIA) %.cc
	$(CXX) $@.cc -o $@ $(CXX_COMMON)

# Plugin libraries.
lib%Lib.so: %Lib.cc
	$(CXX) $< -o $@ -w $(CXX_COMMON) $(CXX_SHARED)\
	 $(CXX_SONAME)$(notdir $@) -Wl,-undefined,dynamic_lookup

# ROOT libraries generated via CLING.
main%Dct.so: main%Dct.cc
	$(CXX) $< -o $@ -w $(CXX_SHARED) $(CXX_COMMON)\
	 $(ROOT_LIB) $(shell $(ROOT_CONFIG) --cflags)
main%Dct.cc: main%Dct.h
ifeq ($(ROOT_USE),true)
	$(ROOT_BIN)rootcling -f $@ -c $^
else
	$(error Error: $@ requires ROOT)
endif

# YODA histogramming. Requires YODA version 2.0.0 or above.
main114: $(PYTHIA) $$@.cc
ifeq ($(YODA_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON:c++11=c++17)\
	 $(YODA_INCLUDE) $(YODA_LIB)
else
	$(error Error: $@ requires YODA)
endif

# HEPMC2 or HEPMC3 (use HEPMC3 if both).
main131 main132 main133 main134 main135: $(PYTHIA) $$@.cc
ifeq ($(HEPMC3_USE),true)
	$(CXX) $@.cc -o $@ $(CXX_COMMON) $(HEPMC3_OPTS)
else ifeq ($(HEPMC2_USE),true)
	$(CXX) $@.cc -o $@ $(CXX_COMMON) $(HEPMC2_OPTS)
else
	$(error Error: $@ requires HEPMC2 or HEPMC3)
endif

# HDF5, HIGHFIVE, and HepMC2 or HepMC3.
main136: $(PYTHIA) $$@.cc
ifeq ($(HDF5_USE)$(HIGHFIVE_USE)$(HEPMC3_USE),truetruetrue)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(HEPMC3_OPTS) $(HDF5_OPTS)
else ifeq ($(HDF5_USE)$(HIGHFIVE_USE)$(HEPMC2_USE),truetruetrue)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(HEPMC2_OPTS) $(HDF5_OPTS)
else
	$(error Error: $@ requires HDF5, HIGHFIVE, and HEPMC2 or HEPMC3)
endif

# General ROOT examples without other external dependencies. 
main141 main143: $(PYTHIA) $$@.cc
ifeq ($(ROOT_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(ROOT_OPTS)
else
	$(error Error: $@ requires ROOT)
endif

# ROOT with FastJet. 
main142: $(PYTHIA) $$@.cc 
ifeq ($(ROOT_USE)$(FASTJET3_USE),truetrue)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(FASTJET3_OPTS) $(ROOT_OPTS)
else	
	$(error Error: $@ requires ROOT and FASTJET3)
endif

# RIVET with optional ROOT and ROOT user library.
main144: $(PYTHIA) $$@.cc $(if $(filter true,$(ROOT_USE)),$$@Dct.so)
	$(CXX) $@.cc -o $@ -w\
	 $(if $(filter true,$(RIVET_USE)),-w $(RIVET_OPTS),$(CXX_COMMON))\
	 $(if $(filter true,$(ROOT_USE)),$@Dct.so $(ROOT_OPTS))

# FASTJET3.
main161 main212 main213 main507: $(PYTHIA) $$@.cc
ifeq ($(FASTJET3_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(FASTJET3_OPTS)
else
	$(error Error: $@ requires FASTJET3)
endif

# Optional HDF5, HEPMC3, and RIVET.
main164: $(PYTHIA) $$@.cc
	$(CXX) $@.cc -o $@ $(if $(filter true,$(HDF5_USE)),$(HDF5_OPTS))\
	 $(if $(filter true,$(HEPMC3_USE)),$(HEPMC3_OPTS))\
	 $(if $(filter true,$(RIVET_USE)),-w $(RIVET_OPTS),$(CXX_COMMON))

# FASTJET3 with recursive tools.
main215: $(PYTHIA) $$@.cc
ifeq ($(FASTJET3_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(FASTJET3_OPTS)\
	 -lfastjettools -lRecursiveTools
else
	$(error Error: $@ requires FASTJET3)
endif

# Optional HEPMC3.
main224: $(PYTHIA) $$@.cc
	$(CXX) $@.cc -o $@ $(CXX_COMMON)\
	 $(if $(filter true,$(HEPMC3_USE)),$(HEPMC3_OPTS))

# MixMax (no dependency, just remove warnings).
main245: $(PYTHIA) $$@.cc
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON)

# EVTGEN.
main364: $(PYTHIA) $$@.cc
ifeq ($(EVTGEN_USE)$(HEPMC2_USE),truetrue)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(EVTGEN_OPTS)
else
	$(error Error: $@ requires EVTGEN)
endif

# Optional RIVET.
main421: $(PYTHIA) $$@.cc
	$(CXX) $@.cc -o $@ $(if $(filter true,$(RIVET_USE)),\
	 -w $(RIVET_OPTS),$(CXX_COMMON))

# Clean.
clean:
	@rm -f main*[0-9]; rm -f mymain*[0-9]; rm -f test*[0-9]; rm -f *.d;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.cc; rm -f *.so;\
	rm -f *.log; rm -f *out[0-9]*; rm -f *.dat; rm -f main*.pyc;\
	rm -f plot*.py; rm -f *plot.py; rm -f fig*.pdf;\
	rm -f *.hepmc; rm -f *.yoda; rm -f *.root; rm -f *.pcm;\
	rm -f weakbosons.lhe; rm -f particles.xml; rm -f settings.xml;\
	rm -rf amcatnlorun; rm -rf helaconiarun; rm -rf madgraphrun;
