# main39.py is a part of the PYTHIA event generator.
# Copyright (C) 2024 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.

# Authors: Philip Ilten <philten@cern.ch>.

# Keywords: particle data; python

# This provides a simple Python class which reads and parses the Pythia 8
# particle database, without requiring the Pythia 8 Python bindings.

#==========================================================================
class ParticleData:
    """
    The 'ParticleData' class stores the necessary information to
    define a particle.
    """

    #----------------------------------------------------------------------
    def __init__(self, id = 0, name = "", m0 = 0, tau0 = float("inf"),
                 spinType = 0, chargeType = 0, colType = 0,
                 mWidth = 0, mMin = 0, mMax = float("inf"), varWidth = False,
                 antiName = ""):
        """
        Initialize the class with the following: 'id' is the particle ID
        number, 'name' the name, 'm0' the mass in GeV, 'tau0' the
        proper lifetime in mm/c, 'spinType' twice the particle spin,
        'chargeType' is three times the electromagnetic chargeType,
        'colType' is the colType type, 'mWidth' is the particle width
        in GeV, 'mMin' is the minimum mass in GeV, 'mMax' is the
        maximum mass in GeV, 'varWidth' is a flag if variable width
        should be used, and 'antiName' is a dummy argument not used.
        """
        self.id = int(id)
        self.name = str(name)
        self.m0 = float(m0)
        self.tau0 = float(tau0)
        self.spinType = int(spinType)
        self.chargeType = int(chargeType)
        self.colType = int(colType)
        self.mWidth = float(mWidth)
        self.mMin = float(mMin)
        self.mMax = float(mMax)
        self.varWidth = varWidth == "on"
        self.anti = antiName if antiName else None
        self.n = tuple([(self.id/pow(10, i)) % 10 for i in range(10)])
        if self.mMax == 0: self.mMax = float("inf")
        
    #----------------------------------------------------------------------
    def isQuark(self): return 0 < abs(self.id) < 10
    def isDiquark(self): return 1000 < abs(self.id) < 10000
    def isBaryon(self):
        id = abs(self.id)
        return not (
            id <= 1000 or 1000000 <= id <= 9000000 or id >= 9900000 or
            self.n[0] == 0 or self.n[1] == 0 or self.n[2] == 0 or
            self.n[3] == 0)
    def isMeson(self):
        id = abs(self.id)
        if id == 130 or id == 310: return True
        return not (
            id <= 100 or 1000000 <= id <= 9000000 or id >= 9900000 or
            self.n[0] == 0 or self.n[1] == 0 or self.n[2] == 0 or
            self.n[3] != 0)
    def isNucleus(self): return abs(self.id) > 1000000000
        
    #----------------------------------------------------------------------
    def __str__(self):
        """
        Return a string to print of this particle data.
        """
        return ("%10s: %s\n"*9 + "%10s: %s") % (
            "id", self.id, "name", self.name, "m0", self.m0,
            "tau0", self.tau0, "spinType", self.spinType, "chargeType",
            self.chargeType, "colType", self.colType, "mWidth", self.mWidth,
            "mMin", self.mMin, "mMax", self.mMax)
    
#==========================================================================
class ParticleDatabase(dict):
    """
    The 'ParticleDatabase' initializes and stores the 'ParticleData' for
    all particles in the 'ParticleData.xml' file from Pythia 8.
    """

    #----------------------------------------------------------------------
    def __init__(self, xmlfile = None):
        """
        Read in the particle data from the XML file 'xmlfile'.
        """
        # Find the XML file, if not provided.
        if xmlfile == None:
            try:
                cfg = open("Makefile.inc")
                for line in cfg:
                    if line.startswith("PREFIX_SHARE="):
                        xmlfile = line.split("=", 1)[-1].strip(); break
                cfg.close()
            except:
                import os
                xmlfile = os.path.dirname(os.path.abspath(__file__))
                xmlfile += "/../share/Pythia8"
            xmlfile += "/xmldoc/ParticleData.xml"
        
        # Instantiate the base class and open the XML file.
        dict.__init__(self)
        xml = open(xmlfile)
        # Loop over the file.
        pstr = ""
        for line in xml:
            line = line.strip()
            if line.startswith("<particle"): pstr = line
            elif pstr:
                pstr += " " + line
                if line.endswith(">"):
                    self.add(pstr)
                    pstr = ""
        xml.close()

    #----------------------------------------------------------------------
    def add(self, pstr):
        """
        Parses the XML for a particle and adds it to the database.
        """
        import shlex
        # Create the dictionary.
        pdct = {pair.split("=", 1)[0]: pair.split("=", 1)[1] for pair
                in shlex.split(pstr[9:-1])}
        # Add the particle.
        pdat = ParticleData(**pdct)
        self[pdat.id] = pdat
        self[pdat.name] = pdat
        # Add the anti-particle if it exists, flip ID and charge.
        if pdat.anti:
            pdct["name"] = pdat.anti
            pdct["id"] = -pdat.id
            pdct["chargeType"] = -pdat.chargeType
            adat = ParticleData(**pdct)
            self[adat.id], self[adat.name] = adat, adat
            pdat.anti, adat.anti = adat, pdat

#==========================================================================
from sys import argv
if __name__== "__main__":
    
    # Parse particle request.
    if len(argv) < 2:
        print("Usage:\n  python main39.py <particle id/name>")
        exit(1)
    try:
        key = int(argv[1])
    except ValueError:
        key = argv[1]

    # Create a ParticleDatabase object and print the electron info.
    pdb = ParticleDatabase()
    try:
        print(pdb[key])
    except KeyError:
        print("Particle id or name not found: " + argv[1])
