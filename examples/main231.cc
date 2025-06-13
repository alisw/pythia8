// main231.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; command line option; command file; analysis

// This is a simple test program.
// It illustrates (a) how to collect the analysis code in a separate class
// and (b) how to provide the .cmnd filename on the command line

// Once you have linked the main program you can run it with a command line
//     ./main231 -c main231.cmnd > main231.log

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/InputParser.h"

using namespace Pythia8;

//==========================================================================

// Put all your own analysis code in the myAnalysis class.

class MyAnalysis {

public:

  // Constructor can be empty.
  MyAnalysis() {}

  // Initialization actions.
  void init();

  // Analysis of each new event.
  void analyze(Event& event);

  // Show final results.
  void finish();

private:

  // Declare variables and objects that span init - analyze - finish.
  int  nEvt;
  Hist brH, yH, etaChg, mult;

};

//--------------------------------------------------------------------------

// The initialization code.

void MyAnalysis::init() {

  // Initialize counter for number of events.
  nEvt = 0;

  // Book histograms.
  brH.book("Higgs branching ratios by flavour", 30, -0.5, 29.5);
  yH.book("Higgs rapidity", 100, -10., 10.);
  etaChg.book("charged pseudorapidity", 100, -10., 10.);
  mult.book( "charged multiplicity", 100, -0.5, 799.5);

}

//--------------------------------------------------------------------------

// The event analysis code.

void MyAnalysis::analyze(Event& event) {

  // Increase counter.
  ++nEvt;

  // Find latest copy of Higgs and plot its rapidity.
  int iH = 0;
  for (int i = 0; i < event.size(); ++i)
    if (event[i].id() == 25) iH = i;
  yH.fill( event[iH].y() );

  // Plot flavour of decay channel.
  int idDau1 = event[ event[iH].daughter1() ].idAbs();
  int idDau2 = event[ event[iH].daughter2() ].idAbs();
  int iChan  = 29;
  if (idDau2 == idDau1 && idDau1 < 25) iChan = idDau1;
  if (min( idDau1, idDau2) == 22 && max( idDau1, idDau2) == 23) iChan = 26;
  brH.fill( iChan);

  // Plot pseudorapidity distribution. Sum up charged multiplicity.
  int nChg = 0;
  for (int i = 0; i < event.size(); ++i)
  if (event[i].isFinal() && event[i].isCharged()) {
    etaChg.fill( event[i].eta() );
    ++nChg;
  }
  mult.fill( nChg );

}

//--------------------------------------------------------------------------

// The finishing code.

void MyAnalysis::finish() {

  // Normalize histograms.
  double binFactor = 5. / nEvt;
  yH     *= binFactor;
  etaChg *= binFactor;

  // Print histograms.
  cout << brH << yH << etaChg << mult;

}

//==========================================================================

// You should not need to touch the main program: its actions are
// determined by the .cmnd file and the rest belongs in MyAnalysis.

int main(int argc, char* argv[]) {

  // Set up command line options.
  InputParser ip("Illustrates how collect analysis code and read commands.",
    {"./main231 -c main231.cmnd"});
  ip.require("c", "Use this user-written command file.", {"-cmnd"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Confirm that external file will be used for settings.
  cout << " PYTHIA settings will be read from file "
       << ip.get<string>("c") << endl;

  // Declare generator. Read in commands from external file.
  Pythia pythia;
  pythia.readFile(ip.get<string>("c"));

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Declare user analysis class. Do initialization part of it.
  MyAnalysis myAnalysis;
  myAnalysis.init();

  // Read in number of event and maximal number of aborts.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");
  bool hasPL = pythia.flag("PartonLevel:all");

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // User Analysis of current event.
    myAnalysis.analyze( (hasPL ? pythia.event : pythia.process) );

  // End of event loop.
  }

  // Final statistics.
  pythia.stat();

  // User finishing.
  myAnalysis.finish();

  // Done.
  return 0;
}
