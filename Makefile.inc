# PYTHIA configuration file.
# Generated on Tue Jul 14 14:26:37 CEST 2015 with the user supplied options:
# --with-hepmc2=/usr/local/new/
# --enable-shared

# Install directory prefixes.
PREFIX_BIN=/scratch/galette/bierlich/work/pythiaSVN/trunk/bin
PREFIX_INCLUDE=/scratch/galette/bierlich/work/pythiaSVN/trunk/include
PREFIX_LIB=/scratch/galette/bierlich/work/pythiaSVN/trunk/lib
PREFIX_SHARE=/scratch/galette/bierlich/work/pythiaSVN/trunk/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=true
CXX=g++
CXX_COMMON=-O2 -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

# FASTJET3 configuration.
FASTJET3_USE=false
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

# HEPMC2 configuration.
HEPMC2_USE=true
HEPMC2_BIN=/usr/local/new//bin/
HEPMC2_INCLUDE=/usr/local/new//include
HEPMC2_LIB=/usr/local/new//lib

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=
LHAPDF5_PLUGIN=LHAPDF5.h

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=
LHAPDF6_PLUGIN=LHAPDF5.h

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=
PROMC_LIB=

# ROOT configuration.
ROOT_USE=false
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=
BOOST_LIB=
