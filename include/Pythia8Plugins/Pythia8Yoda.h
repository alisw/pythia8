// Pythia8Yoda.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
#include "Pythia8/PythiaStdlib.h"

#ifndef PYTHIA8YODA_H
#define PYTHIA8YODA_H

#include "YODA/Histo.h"
#include "YODA/Profile.h"
#include "YODA/WriterYODA.h"
#include "YODA/ReaderYODA.h"

namespace Pythia8 {

// Simplified interface to booking and writing of YODA histograms.
// Remember to configure Pythia with YODA; --with-yoda.
// YODA version 2 or greater is required.
//
// Usage: 1) Create a Pythia8Yoda object in your main program. 2) Use
// it to book histograms, receive shared pointers to the histograms. 3) Fill
// the histograms while you run, and 4) At the end of the analysis, write out
// the histograms to a file.

// Typedef of shared pointers to analysis objects.
using AnaObjectPtr = shared_ptr<YODA::AnalysisObject>;
using Histo1DPtr   = shared_ptr<YODA::Histo1D>;
using Histo2DPtr   = shared_ptr<YODA::Histo2D>;
using Histo3DPtr   = shared_ptr<YODA::Histo3D>;
using Profile1DPtr = shared_ptr<YODA::Profile1D>;
using Profile2DPtr = shared_ptr<YODA::Profile2D>;

class Pythia8Yoda {

public:

  // The constructor needs a prefix name, which will be put on the histogram
  // path, as well as an output file name.
  Pythia8Yoda(const string& anaNameIn = "/PYTHIA8/", const string& outputIn =
    "pythia") : anaName("/"+anaNameIn+"/"), outName(outputIn),
    finalized(false) { }

  // The destructor will write the histogram file if this has not
  // already been done.
  ~Pythia8Yoda() {
    if (!finalized) write();
  }

  // Write histograms.
  void write() {
    string purge = ".yoda";
    size_t pos = outName.find(purge);
    if (pos != std::string::npos)
      outName.replace(pos, purge.length(), "");
    std::cout << "Writing histograms to: " << outName << ".yoda" << std::endl;
    YODA::WriterYODA::write(outName+".yoda", anaObjects.begin(),
      anaObjects.end());
    finalized = true;
  }

  // Read and return a named YODA object from a file.
  // Input arguments filename.yoda and full object path.
  template<typename T>
  static T read(const string& fName, const string& objPath) {
    // Read the file content.
    vector<YODA::AnalysisObject* > anaObjs;
    YODA::ReaderYODA::read(fName, anaObjs);

    // Find the object and convert if possible.
    T result;
    for (auto ao : anaObjs)
      if (ao->path() == objPath)
        if (T* castObj = dynamic_cast<T*>(ao)) {
          result = *castObj;
          break;
        }

    // Warning if object not found/converted.
    if (!(result.path() == objPath))
      cout << "Warning: " << objPath << " not found/incorrect type," << endl;

    // Clean up
    for (auto* obj : anaObjs) delete obj;

    return result;
  }

  // Book any YODA object and return a shared ptr.
  template<typename T, typename... Args>
  shared_ptr<T> book(Args&&... args) {
    auto h = make_shared<T>(std::forward<Args>(args)...);
    anaObjects.push_back(h);
    return h;
  }

  // Write the most common use cases to give better compiler warnings and
  // easier use.
  Histo1DPtr bookHisto1D(int nBins, double xMin, double xMax,
    const string& title) {
    return book<YODA::Histo1D>(nBins, xMin, xMax, anaName + title, title);
  }

  // Book a 2D histogram.
  Histo2DPtr bookHisto2D(int nBinsX, double xMin, double xMax, int nBinsY,
    double yMin, double yMax, const std::string& title) {
    return book<YODA::Histo2D>(nBinsX, xMin, xMax, nBinsY, yMin, yMax,
      anaName + title, title);
  }

  // Book a 3D histogram.
  Histo3DPtr bookHisto3D(int nBinsX, double xMin, double xMax, int nBinsY,
    double yMin, double yMax, int nBinsZ, double zMin, double zMax,
    const std::string& title) {
    return book<YODA::Histo3D>(nBinsX, xMin, xMax, nBinsY, yMin, yMax,
      nBinsZ, zMin, zMax, anaName + title, title);
  }

  // Book a 1D profile.
  Profile1DPtr bookProfile1D(int nBins, double xMin, double xMax,
    const std::string& title) {
    return book<YODA::Profile1D>(nBins, xMin, xMax, anaName + title, title);
  }

  // Book a 2D profile.
  Profile2DPtr bookProfile2D(int nBinsX, double xMin, double xMax, int nBinsY,
   double yMin, double yMax, const std::string& title) {
    return book<YODA::Profile2D>(nBinsX, xMin, xMax, nBinsY, yMin, yMax,
      anaName + title, title);
    }

private:
  string anaName;
  string outName;
  bool finalized;
  std::vector<AnaObjectPtr> anaObjects;

};

}

#endif
