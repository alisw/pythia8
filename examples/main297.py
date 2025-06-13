# main297.py is a part of the PYTHIA event generator.
# Copyright (C) 2025 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.

# Authors: Philip Ilten <philten@cern.ch>,
#          Artem Havryliuk,
#          Jim Pivarski

# Keywords: python; parallelism

# This example shows how to use the Awkward array interface to
# generating Pythia events.

# To set the path to the Pythia 8 Python interface do either
# (in a shell prompt):
#      export PYTHONPATH=$(PREFIX_LIB):$PYTHONPATH
# or the following which sets the path from within Python.
import sys
cfg = open("Makefile.inc")
lib = "../lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Import the Pythia module.
import pythia8

#==========================================================================

# The following code demonstrates how to use the nextBatch function to
# generate an Awkward array of events using a single Pythia instance.

# Create a Pythia instance.
pythia = pythia8.Pythia()
pythia.readString("Beams:eCM = 8000.")
pythia.readString("HardQCD:all = on")
pythia.readString("PhaseSpace:pTHatMin = 20.")
pythia.init()

# Generate a run of events.
run = pythia.nextBatch(2)

# The following is intended to demonstrate the structure of the
# returned array. Normally, you would not use a for loop like this and
# instead use array slicing, which is shown later.
for event in run:

    # The event object consists of the event info, accessed as
    # event.info, and the particles, accessed as event.prt.

    # Print the ID of the first scattering parton for each event.
    # Most, but not all, of the event info is available this way,
    # where the naming convention follows the methods of the Info
    # class in Pythia.
    print(event.info.id1)

    # Loop over the particles.
    for prt in event.prt:

        # Print the ID of the first scattering parton. The
        # particle information is accessible in a similar fashion to
        # the event info, where all the relevant data members of the
        # Particle class in Pythia are available.
        if (prt.status != -21): continue
        print(prt.id)

        # Print the particle energy. The momentum is accessed via
        # prt.p, and has the relevant Vec4 data members: px, py, pz,
        # and e.
        vec = prt.p
        print(vec.e)
        break

#==========================================================================

# The nextBatch method is also available for PythiaParallel. The
# configuration of PythiaParallel is the same as shown in
# main221.cc. However, an event analysis function is not needed,
# because the result is instead the full event record for each event.

# Create and initialize the parallel instance.
pythia = pythia8.PythiaParallel()
pythia.readString("Beams:eCM = 8000.")
pythia.readString("HardQCD:all = on")
pythia.readString("PhaseSpace:pTHatMin = 20.")
pythia.init()

# Generate events.
events = pythia.nextBatch(1000)

# Generally, we would use array slicing to perform operations on the
# event record. There are some good examples given with the Awkward
# library.
# https://awkward-array.org/doc/main/getting-started/
# jagged-ragged-awkward-arrays.html

# We can use array slicing to quickly access all the charged pions in
# the event.
cpis = events.prt[abs(events.prt.id) == 211]

# We can then plot the energy of these pions. The following requires
# matplotlib to be installed.
#     python -m pip install matplotlib
# To make a histogram of the array, we have to flatten it with the
# awkward.flatten method.
import awkward as ak
try:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.hist(ak.flatten(cpis.p.e))
    fig.savefig("main297_e.pdf")
except:
    pass

# It is possible to also use the scikit-hep vector package to interact
# with the particle momentum as a four-vector. The following requires
# the vector package to be installed.
#     python -m pip install vector
try:
    import vector

    # This tells vector to convert compatible Awkward arrays into
    # appropriate vector objects.
    vector.register_awkward()

    # We can now plot the pseudorapidity of the charged pions. Here,
    # the eta member of the charged pion momenta is now available as
    # p.eta. Similarly, methods like phi, etc. are available.
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.hist(ak.flatten(cpis.p.eta))
    fig.savefig("main297_eta.pdf")
except:
    pass
