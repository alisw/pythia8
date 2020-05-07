// VinciaQED.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for Vincia's QED
// shower class and related auxiliary methods. Main author is Rob
// Verheyen.

#include "Pythia8/VinciaQED.h"

namespace Pythia8 {

//==========================================================================

// Class for the "Hungarian" pairing algorithm.

//--------------------------------------------------------------------------

// A single function wrapper for solving assignment problem.

double HungarianAlgorithm::solve(std::vector <std::vector<double> >&
  distMatrix, std::vector<int>& assignment) {

  unsigned int nRows = distMatrix.size();
  unsigned int nCols = distMatrix[0].size();
  double *distMatrixIn = new double[nRows * nCols];
  int *solution = new int[nRows];
  double cost = 0.0;

  // Fill in the distMatrixIn. Mind the index is "i + nRows * j".
  // Here the cost matrix of size MxN is defined as a double precision
  // array of N*M elements. In the solving functions matrices are seen
  // to be saved MATLAB-internally in row-order. (i.e. the matrix [1
  // 2; 3 4] will be stored as a vector [1 3 2 4], NOT [1 2 3 4]).
  for (unsigned int i = 0; i < nRows; i++)
    for (unsigned int j = 0; j < nCols; j++)
      distMatrixIn[i + nRows * j] = distMatrix[i][j];

  // Call solving function.
  optimal(solution, &cost, distMatrixIn, nRows, nCols);
  assignment.clear();
  for (unsigned int r = 0; r < nRows; r++)
    assignment.push_back(solution[r]);
  delete[] distMatrixIn;
  delete[] solution;
  return cost;

}

//--------------------------------------------------------------------------

// Solve optimal solution for assignment.

void HungarianAlgorithm::optimal(int *assignment, double *cost,
  double *distMatrixIn, int nOfRows, int nOfColumns) {

  // Initialization.
  double *distMatrix, *distMatrixTemp, *distMatrixEnd, *columnEnd,
    value, minValue;
  bool *coveredColumns, *coveredRows, *starMatrix, *newStarMatrix,
    *primeMatrix;
  int nOfElements, minDim, row, col;
  *cost = 0;
  for (row = 0; row<nOfRows; row++) assignment[row] = -1;

  // Generate working copy of distance matrix. Check if all matrix
  // elements are positive.
  nOfElements = nOfRows * nOfColumns;
  distMatrix = (double *)malloc(nOfElements * sizeof(double));
  distMatrixEnd = distMatrix + nOfElements;
  for (row = 0; row<nOfElements; row++) {
    value = distMatrixIn[row];
    if (value < 0)
      std::cerr << "HungarianAlgorithm::assigmentoptimal(): All"
                << " matrix elements have to be non-negative." << std::endl;
    distMatrix[row] = value;
  }

  // Memory allocation.
  coveredColumns = (bool *)calloc(nOfColumns, sizeof(bool));
  coveredRows    = (bool *)calloc(nOfRows, sizeof(bool));
  starMatrix     = (bool *)calloc(nOfElements, sizeof(bool));
  primeMatrix    = (bool *)calloc(nOfElements, sizeof(bool));
  // Used in step4.
  newStarMatrix = (bool *)calloc(nOfElements, sizeof(bool));

  // Preliminary steps.
  if (nOfRows <= nOfColumns) {
      minDim = nOfRows;
      for (row = 0; row<nOfRows; row++) {
        // Find the smallest element in the row.
        distMatrixTemp = distMatrix + row;
        minValue = *distMatrixTemp;
        distMatrixTemp += nOfRows;
        while (distMatrixTemp < distMatrixEnd) {
          value = *distMatrixTemp;
          if (value < minValue)
            minValue = value;
          distMatrixTemp += nOfRows;
        }

        // Subtract the smallest element from each element of the row.
        distMatrixTemp = distMatrix + row;
        while (distMatrixTemp < distMatrixEnd) {
          *distMatrixTemp -= minValue;
          distMatrixTemp += nOfRows;
        }
      }

      // Steps 1 and 2a.
      for (row = 0; row<nOfRows; row++)
        for (col = 0; col<nOfColumns; col++)
          if (abs(distMatrix[row + nOfRows*col]) < DBL_EPSILON)
            if (!coveredColumns[col]) {
              starMatrix[row + nOfRows*col] = true;
              coveredColumns[col] = true;
              break;
            }
  } else {
    minDim = nOfColumns;
    for (col = 0; col<nOfColumns; col++) {
        // Find the smallest element in the column.
        distMatrixTemp = distMatrix + nOfRows*col;
        columnEnd = distMatrixTemp + nOfRows;
        minValue = *distMatrixTemp++;
        while (distMatrixTemp < columnEnd) {
          value = *distMatrixTemp++;
          if (value < minValue)
            minValue = value;
        }

        // Subtract the smallest element from each element of the column.
        distMatrixTemp = distMatrix + nOfRows*col;
        while (distMatrixTemp < columnEnd)
          *distMatrixTemp++ -= minValue;
    }

    // Steps 1 and 2a.
    for (col = 0; col<nOfColumns; col++)
      for (row = 0; row<nOfRows; row++)
        if (abs(distMatrix[row + nOfRows*col]) < DBL_EPSILON)
          if (!coveredRows[row]) {
            starMatrix[row + nOfRows*col] = true;
            coveredColumns[col] = true;
            coveredRows[row] = true;
            break;
          }
    for (row = 0; row<nOfRows; row++)
      coveredRows[row] = false;
  }

  // Move to step 2b.
  step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
    coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

  // Compute cost and remove invalid assignments.
  calcCost(assignment, cost, distMatrixIn, nOfRows);

  // Free allocated memory.
  free(distMatrix);
  free(coveredColumns);
  free(coveredRows);
  free(starMatrix);
  free(primeMatrix);
  free(newStarMatrix);
  return;

}

//--------------------------------------------------------------------------

// Build the assignment vector.

void HungarianAlgorithm::vect(int *assignment,
  bool *starMatrix, int nOfRows, int nOfColumns) {
  int row, col;
  for (row = 0; row<nOfRows; row++)
    for (col = 0; col<nOfColumns; col++)
      if (starMatrix[row + nOfRows*col]) {
        assignment[row] = col;
        break;
      }
}

//--------------------------------------------------------------------------

// Calculate the assignment cost.

void HungarianAlgorithm::calcCost(int *assignment, double *cost,
  double *distMatrix, int nOfRows) {
  int row, col;
  for (row = 0; row<nOfRows; row++) {
    col = assignment[row];
    if (col >= 0) *cost += distMatrix[row + nOfRows*col];
  }
}

//--------------------------------------------------------------------------

// Factorized step 2a of the algorithm.

void HungarianAlgorithm::step2a(int *assignment, double *distMatrix,
  bool *starMatrix, bool *newStarMatrix, bool *primeMatrix,
  bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns,
  int minDim) {

  // Cover every column containing a starred zero.
  bool *starMatrixTemp, *columnEnd;
  int col;
  for (col = 0; col<nOfColumns; col++) {
    starMatrixTemp = starMatrix + nOfRows*col;
    columnEnd = starMatrixTemp + nOfRows;
    while (starMatrixTemp < columnEnd) {
      if (*starMatrixTemp++) {
        coveredColumns[col] = true;
        break;
      }
    }
  }

  // Move to step 2b (note, original comment changed by Skands).
  step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
    coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

}

//--------------------------------------------------------------------------

// Factorized step 2b of the algorithm.

void HungarianAlgorithm::step2b(int *assignment, double *distMatrix,
  bool *starMatrix, bool *newStarMatrix, bool *primeMatrix,
  bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns,
  int minDim) {

  // Count covered columns.
  int col, nOfCoveredColumns;
  nOfCoveredColumns = 0;
  for (col = 0; col<nOfColumns; col++)
    if (coveredColumns[col]) nOfCoveredColumns++;

  // Algorithm finished.
  if (nOfCoveredColumns == minDim)
      vect(assignment, starMatrix, nOfRows, nOfColumns);
  // Move to step 3.
  else step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
        coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

}

//--------------------------------------------------------------------------

// Factorized step 3 of the algorithm.

void HungarianAlgorithm::step3(int *assignment, double *distMatrix,
  bool *starMatrix, bool *newStarMatrix, bool *primeMatrix,
  bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns,
  int minDim) {

  bool zerosFound;
  int row, col, starCol;
  zerosFound = true;
  while (zerosFound) {
      zerosFound = false;
      for (col = 0; col<nOfColumns; col++)
        if (!coveredColumns[col])
          for (row = 0; row<nOfRows; row++)
            if ((!coveredRows[row]) && (abs(distMatrix[row + nOfRows*col])
                                        < DBL_EPSILON)) {
              // Prime zero.
              primeMatrix[row + nOfRows*col] = true;

              // Find starred zero in current row.
              for (starCol = 0; starCol<nOfColumns; starCol++)
                if (starMatrix[row + nOfRows*starCol])
                  break;

              // No starred zero found, move to step 4.
              if (starCol == nOfColumns) {
                step4(assignment, distMatrix, starMatrix, newStarMatrix,
                      primeMatrix, coveredColumns, coveredRows, nOfRows,
                      nOfColumns, minDim, row, col);
                return;
              } else {
                coveredRows[row] = true;
                coveredColumns[starCol] = false;
                zerosFound = true;
                break;
              }
            }
  }

  // Move to step 5.
  step5(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
    coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

}

//--------------------------------------------------------------------------

// Factorized step 4 of the algorithm.

void HungarianAlgorithm::step4(int *assignment, double *distMatrix,
  bool *starMatrix, bool *newStarMatrix, bool *primeMatrix,
  bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns,
  int minDim, int row, int col) {

  // Generate temporary copy of starMatrix.
  int n, starRow, starCol, primeRow, primeCol;
  int nOfElements = nOfRows*nOfColumns;
  for (n = 0; n<nOfElements; n++) newStarMatrix[n] = starMatrix[n];
  // Star current zero.
  newStarMatrix[row + nOfRows*col] = true;
  // Find starred zero in current column.
  starCol = col;
  for (starRow = 0; starRow<nOfRows; starRow++)
    if (starMatrix[starRow + nOfRows*starCol])
      break;
  while (starRow < nOfRows) {
      // Unstar the starred zero.
      newStarMatrix[starRow + nOfRows*starCol] = false;
      // Find primed zero in current row.
      primeRow = starRow;
      for (primeCol = 0; primeCol<nOfColumns; primeCol++)
        if (primeMatrix[primeRow + nOfRows*primeCol])
          break;
      // Star the primed zero.
      newStarMatrix[primeRow + nOfRows*primeCol] = true;
      // Find starred zero in current column.
      starCol = primeCol;
      for (starRow = 0; starRow<nOfRows; starRow++)
        if (starMatrix[starRow + nOfRows*starCol])
          break;
  }

  // Use temporary copy as new starMatrix, delete all primes, uncover
  // all rows.
  for (n = 0; n<nOfElements; n++) {
    primeMatrix[n] = false;
    starMatrix[n] = newStarMatrix[n];
  }
  for (n = 0; n<nOfRows; n++) coveredRows[n] = false;

  // Move to step 2a.
  step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
    coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

}

//--------------------------------------------------------------------------

// Factorized step 5 of the algorithm.

void HungarianAlgorithm::step5(int *assignment, double *distMatrix,
  bool *starMatrix, bool *newStarMatrix, bool *primeMatrix,
  bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns,
  int minDim) {

  // Find smallest uncovered element h.
  double h, value;
  int row, col;
  h = DBL_MAX;
  for (row = 0; row<nOfRows; row++)
    if (!coveredRows[row])
      for (col = 0; col<nOfColumns; col++)
        if (!coveredColumns[col]) {
            value = distMatrix[row + nOfRows*col];
            if (value < h)
              h = value;
        }

  // Add h to each covered row.
  for (row = 0; row<nOfRows; row++)
    if (coveredRows[row])
      for (col = 0; col<nOfColumns; col++)
        distMatrix[row + nOfRows*col] += h;

  // Subtract h from each uncovered column.
  for (col = 0; col<nOfColumns; col++)
    if (!coveredColumns[col])
      for (row = 0; row<nOfRows; row++)
        distMatrix[row + nOfRows*col] -= h;

  // Move to step 3.
  step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
    coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

}

//==========================================================================

// Class for QED emissions.

//--------------------------------------------------------------------------

// Initialize the pointers.

void QEDemitElemental::initPtr(Rndm* rndmPtrIn,
  PartonSystems* partonSystemsPtrIn) {
  rndmPtr = rndmPtrIn;
  partonSystemsPtr = partonSystemsPtrIn;
  isInitPtr = true;
}

//--------------------------------------------------------------------------

// Initialize.

void QEDemitElemental::init(Event &event, int xIn, int yIn, double shhIn,
    double verboseIn) {

  if (!isInitPtr) printOut(__METHOD_NAME__, "initPtr not called");
  x = xIn;
  y = yIn;
  shh = shhIn;
  hasTrial = false;
  isII = false;
  isIF = false;
  isFF = false;
  isRF = false;
  isIA = false;
  isDip = false;

  // If an II antenna, make sure x is the positive pz state.
  if (!event[x].isFinal() && !event[y].isFinal() && event[x].pz() < 0)
      swap(x,y);

  // If an IF/RF antenna, make sure x is the initial state.
  if (event[x].isFinal() && !event[y].isFinal()) swap(x,y);

  // If a dipole, make sure x is the emitting object.
  if (event[x].isFinal() && event[y].isFinal())
    if (!event[x].isCharged() || event[y].isCharged()) swap(x,y);

  idx = event[x].id();
  idy = event[y].id();
  mx2 = event[x].m2();
  my2 = event[y].m2();
  ex = event[x].e();
  ey = event[y].e();
  m2Ant = m2(event[x], event[y]);
  sAnt = 2*dot4(event[x], event[y]);
  QQ = - event[x].charge() * event[y].charge();

  // II.
  if (!event[x].isFinal() && !event[y].isFinal()) isII = true;

  // IF/RF.
  if (!event[x].isFinal() && event[y].isFinal()) {
    // QQ is flipped for IF antennae.
    QQ = -QQ;
    // Check if initial state is in a beam.
    int mother1 = event[x].mother1();
    // Check if initial particle is A or B.
    if (mother1 <= 2) {
      isIF = true;
      if (event[x].pz() > 0) isIA = true;
    // Otherwise it's a resonance decay.
    } else isRF = true;
  }

  // FF.
  if (event[x].isFinal() && event[y].isFinal()) isFF = true;
  isInit = true;
  verbose = verboseIn;

}

//--------------------------------------------------------------------------

// Initialize.

void QEDemitElemental::init(Event &event, int xIn, vector<int> iRecoilIn,
    double shhIn, double verboseIn) {

  x = xIn;
  iRecoil = iRecoilIn;
  shh = shhIn;
  hasTrial = false;
  isII = false;
  isIF = false;
  isFF = false;
  isRF = false;
  isIA = false;
  isDip = true;
  idx = event[x].id();
  mx2 = event[x].m2();

  // Compute total recoiler momentum.
  Vec4 pRecoil;
  for (int i = 0; i < (int)iRecoil.size(); i++)
    pRecoil += event[iRecoil[i]].p();
  my2 = pRecoil.m2Calc();
  m2Ant = (pRecoil + event[xIn].p()).m2Calc();
  sAnt = 2*pRecoil*event[xIn].p();
  QQ = 1;
  isInit = true;
  verbose = verboseIn;

}

//--------------------------------------------------------------------------

// Generate a trial point.

double QEDemitElemental::generateTrial(Event &event, double q2Start,
  double q2Low, double alphaIn, double cIn) {

  if (hasTrial) return q2Sav;
  q2Sav = q2Low;
  alpha = alphaIn;
  c = cIn;

  // FF.
  if (isFF || isDip) {
    // Adjust starting scale.
    q2Start = min(q2Start, sAnt/4.);
    if (q2Start < q2Low) return q2Low;

    // Compute phase space constants.
    double lambda = m2Ant*m2Ant + mx2*mx2 + my2*my2
      - 2.*m2Ant*mx2 - 2.*m2Ant*my2 - 2.*mx2*my2;
    // zMin is identical for all instances.
    double zMin = (4*q2Low/sAnt < 1E-8)
      ? q2Low/sAnt
      : 0.5*(1. - sqrt(1. - 4*q2Low/sAnt));

    // Generate scale for eikonal piece.
    if (true) {
      double Iz = (zMin < 1E-8)
        ? -2*log(zMin) - 2*zMin - pow2(zMin)
        : 2*log((1-zMin)/zMin);
      double comFac = 2*M_PI*sqrt(lambda)/alpha/Iz/c/sAnt;
      double q2New  = q2Start*pow(rndmPtr->flat(), comFac);
      if (q2New > q2Sav) {
        q2Sav   = q2New;
        zetaSav = 1/(exp(Iz*(0.5 - rndmPtr->flat())) + 1);
        sxjSav  = sqrt(sAnt*q2Sav*zetaSav/((1-zetaSav)));
        syjSav  = sqrt(sAnt*q2Sav*(1-zetaSav)/(zetaSav));
      }
    }
    // Generate scale for additional W piece on x.
    if (isFF && abs(idx) == 24) {
      double Iz = (zMin < 1E-8) ?
        -log(zMin) - zMin - pow2(zMin)/2. : log((1-zMin)/zMin);
      double comFac = 3.*M_PI*sqrt(lambda)/alpha/Iz/c/sAnt/2.;
      double q2New = q2Start*pow(rndmPtr->flat(), comFac);
      if (q2New > q2Sav) {
        q2Sav    = q2New;
        double r = rndmPtr->flat();
        zetaSav  = (zMin < 1E-8)
          ? 1 - pow(zMin,r)*(1. - (1.-r)*zMin)
          : 1 - pow(zMin,r)*pow(1.-zMin, 1.-r);
        sxjSav   = q2Sav/zetaSav;
        syjSav   = zetaSav*sAnt;
      }
    }
    // Generate scale for additional W piece on y.
    if (isFF && abs(idy) == 24) {
      double Iz = (zMin < 1E-8)
        ? -log(zMin) - zMin - pow2(zMin)/2.
        : log((1-zMin)/zMin);
      double comFac = 3.*M_PI*sqrt(lambda)/alpha/Iz/c/sAnt/2.;
      double q2New  = q2Start*pow(rndmPtr->flat(), comFac);
      if (q2New > q2Sav) {
        q2Sav    = q2New;
        double r = rndmPtr->flat();
        zetaSav  = (zMin < 1E-8)
          ? 1 - pow(zMin,r)*(1. - (1.-r)*zMin)
          : 1 - pow(zMin,r)*pow(1.-zMin, 1.-r);
        sxjSav   = zetaSav*sAnt;
        syjSav   = q2Sav/zetaSav;
      }
    }
  }

  // IF.
  if (isIF) {
    // Compute exmax and sjkMax.
    double exUsed = 0;
    int nSys = partonSystemsPtr->sizeSys();
    for (int i=0; i<nSys; i++) {
      int iEv;
      if (isIA) iEv = partonSystemsPtr->getInA(i);
      else iEv = partonSystemsPtr->getInB(i);
      exUsed += event[iEv].e();
    }
    double exMax = sqrt(shh)/2.0 - (exUsed-ex);
    double sjkMax = sAnt*(exMax-ex)/ex;

    // Adjust starting scale.
    q2Start = min(q2Start, sAnt*(exMax - ex)/ex);
    if (q2Start < q2Low) return q2Low;
    double zMax = sjkMax/(sjkMax + my2);
    double zMin = q2Low/sjkMax;

    // Check if there is any phase space available.
    if (zMin < zMax) {
      // Generate scale for eikonal piece.
      if (true) {
        double Iz     = log(zMax/zMin);
        double Rpdf   = 1.;
        double comFac = M_PI/alpha/Iz/c/Rpdf;
        double q2New  = q2Start*pow(rndmPtr->flat(), comFac);
        if (q2New > q2Sav) {
          q2Sav   = q2New;
          zetaSav = zMin*pow(zMax/zMin, rndmPtr->flat());
          sxjSav  = sAnt*zetaSav + q2Sav;
          syjSav  = q2Sav/zetaSav;
        }
      }

      // Generate scale for additional W piece on y. The veto
      // probability for this antenna piece includes an additional
      // factor which is incorporated by a veto locally.
      if (abs(idy) == 24) {
        double Iz = log((1-zMin)/(1-zMax));
        double Rpdf = 1.;
        double comFac = 3.*M_PI/alpha/Iz/c/Rpdf/2.;
        double q2New = q2Start;
        double zetaNew, sxjNew, syjNew;
        while (true) {
          q2New  *= pow(rndmPtr->flat(), comFac);
          if (q2New < q2Sav) {break;}
          zetaNew = 1. - (1-zMin)*pow((1-zMax)/(1-zMin),rndmPtr->flat());
          sxjNew  = sAnt*zetaNew + q2New;
          syjNew  = q2New/zetaNew;

          // Veto probability.
          double pVeto = sAnt/(sAnt + syjNew);
          if (rndmPtr->flat() < pVeto) {
            q2Sav   = q2New;
            zetaSav = zetaNew;
            sxjSav  = sxjNew;
            syjSav  = syjNew;
            break;
          }
        }
      }
    }
  }

  // II.
  if (isII) {
    // Adjust starting scale.
    q2Start = min(q2Start, pow2(shh-sAnt)/shh/4.);
    if (q2Start < q2Low) {return q2Low;}

    // Generate scale for eikonal piece.
    if (true) {
      double zMin = 0.5*(shh-sAnt -
        sqrt((shh-sAnt)*(shh-sAnt) - (4.*shh*q2Low)))/shh;
      double zMax = 0.5*(shh-sAnt +
        sqrt((shh-sAnt)*(shh-sAnt) - (4.*shh*q2Low)))/shh;
      if (4.*shh*q2Low/pow2((shh-sAnt)) < 1e-8)
        zMin = q2Low/(shh-sAnt);
      double Iz     = log(zMax*(1-zMin)/(1-zMax)/zMin);
      double Rpdf   = 1.;
      double comFac = M_PI/alpha/Iz/c/Rpdf;
      double q2New  = q2Start*pow(rndmPtr->flat(), comFac);
      if (q2New > q2Sav) {
        q2Sav    = q2New;
        double r = rndmPtr->flat();
        double w = pow(zMax/(1-zMax), r) * pow(zMin/(1-zMin), 1.-r);
        zetaSav  = w/(1.+w);
        sxjSav   = (q2Sav + sAnt*zetaSav)/(1.-zetaSav);
        syjSav   = q2Sav/zetaSav;
      }
    }
  }

  // RF.
  if (isRF) {
    // Compute phase space constants.
    double mr2 = abs((event[x].p() - event[y].p()).m2Calc());
    double mx = sqrt(mx2);
    double my = sqrt(my2);
    double mr = sqrt(mr2);
    double lambda = mr2*mr2 + mx2*mx2 + my2*my2
      - 2.*mr2*mx2 - 2.*mr2*my2 - 2.*mx2*my2;
    double sjkMax = pow2(mx - mr) - my2;
    double sajMax = mx2 - pow2(my + mr);
    // Adjust starting scale.
    q2Start = min(q2Start, sajMax*sjkMax/(sAnt + sjkMax));

    // Generate scale for eikonal piece.
    if (true) {
      double zMin   = q2Low/sjkMax;
      double zMax   = sajMax/sAnt;
      if(zMin < zMax){
        double Iz     = log(zMax/zMin);
        double comFac = M_PI*sqrt(lambda)*sAnt/alpha/Iz/c/pow2(sAnt+sjkMax);
        double q2New  = q2Start;
        double zetaNew, sxjNew, syjNew;
        while (true) {
          q2New *= pow(rndmPtr->flat(), comFac);
          if (q2New < q2Sav) {break;}
          zetaNew = zMin*pow(zMax/zMin, rndmPtr->flat());
          sxjNew  = sAnt*zetaNew + q2New;
          syjNew  = q2New/zetaNew;

          // Veto probability.
          double pVeto = pow2(syjNew+sAnt)/pow2(sjkMax+sAnt);
          if (rndmPtr->flat() < pVeto) {
            q2Sav   = q2New;
            zetaSav = zetaNew;
            sxjSav  = sxjNew;
            syjSav  = syjNew;
            break;
          }
        }
      }
    }

    // Generate scale for W in initial state.
    if (abs(idx) == 24) {
      double zMin   = q2Low/(sajMax - q2Low);
      double zMax   = sjkMax/sAnt;
      if(zMin < zMax && zMin > 0){
        double Iz     = pow2(zMax) + (1./3.)*pow3(zMax)
          - pow2(zMin) - (1./3.)*pow3(zMin);
        double comFac = 3.*M_PI*sqrt(lambda)/alpha/Iz/c/sAnt/2.;
        double q2New  = q2Start*pow(rndmPtr->flat(), comFac);

        if (q2New > q2Sav) {
          double a = rndmPtr->flat()*Iz + pow2(zMin) + (1./3.)*pow3(zMin);
          // Solve for zeta using Newton-Raphson.
          int n = 0;
          zetaSav = zMin;
          while(true) {
            n++;
            double f = pow2(zetaSav) + pow3(zetaSav)/3. - a;
            double fPrime  = 2.*zetaSav + pow2(zetaSav);
            double zetaNew = zetaSav - f/fPrime;
            if (zetaNew > zMax) {zetaSav = zMax; continue;}
            if (zetaNew < zMin) {zetaSav = zMin; continue;}
            if (abs(zetaNew - zetaSav) < 1E-8*zetaNew) {
              zetaSav = zetaNew;
              break;
            }
            if (n > 500) {
              printOut(__METHOD_NAME__,
                "RF(W) failed to find zeta with Newton-Raphson");
              break;
            }
            zetaSav = zetaNew;
          }
          q2Sav  = q2New;
          sxjSav = (1.+zetaSav)*q2Sav/zetaSav;
          syjSav = sAnt*zetaSav;
        }
      }
    }

    // Generate scale for W in final state.
    if (abs(idy) == 24) {
      double zMin   = q2Low/sjkMax;
      double zMax   = sajMax/sAnt;
      if (zMin < zMax) {
        double Iz     = log((1-zMin)/(1-zMax));
        double comFac = 3.*M_PI*sqrt(lambda)/alpha/Iz/c/(sAnt+sjkMax)/2.;
        double q2New  = q2Start;
        double zetaNew, sxjNew, syjNew;
        while (true) {
          q2New *= pow(rndmPtr->flat(), comFac);
          if (q2New < q2Sav) {break;}
          zetaNew = 1. - (1-zMin)*pow((1-zMax)/(1-zMin),rndmPtr->flat());
          sxjNew  = sAnt*zetaNew + q2New;
          syjNew  = q2New/zetaNew;

          // Veto probability.
          double pVeto = (syjNew+sAnt)/(sjkMax+sAnt);
          if (rndmPtr->flat() < pVeto) {
            q2Sav = q2New;
            zetaSav = zetaNew;
            sxjSav = sxjNew;
            syjSav = syjNew;
            break;
          }
        }
      }
    }
  }
  phiSav = 2.*M_PI*rndmPtr->flat();
  hasTrial = true;
  return q2Sav;

}

//==========================================================================

// Class for a QED emission system.

//--------------------------------------------------------------------------

// Initialize the pointers.

void QEDemitSystem::initPtr(Info* infoPtrIn, VinciaCommon* vinComPtrIn) {
  infoPtr       = infoPtrIn;
  particleDataPtr  = infoPtr->particleDataPtr;
  partonSystemsPtr = infoPtr->partonSystemsPtr;
  rndmPtr          = infoPtr->rndmPtr;
  settingsPtr      = infoPtr->settingsPtr;
  vinComPtr        = vinComPtrIn;
  isInitPtr = true;
}

//--------------------------------------------------------------------------

// Initialize settings for current run.

void QEDemitSystem::init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    int verboseIn) {

  // Verbose setting.
  if (!isInitPtr)
    printOut(__METHOD_NAME__,"QEDemitSystem:initPtr not called");
  verbose = verboseIn;

  // Set beam pointers.
  beamAPtr = beamAPtrIn;
  beamBPtr = beamBPtrIn;

  // Settings.
  mode = settingsPtr->mode("Vincia:photonEmissionMode");
  useFullWkernel = settingsPtr->flag("Vincia:fullWkernel");
  emitBelowHad = settingsPtr->flag("PartonLevel:Remnants");

  // Constants.
  TINYPDF = pow(10,-10);

  // Initialized.
  isInit = true;

}

//--------------------------------------------------------------------------

// Prepare a QED system.

void QEDemitSystem::prepare(int iSysIn, Event &event, double q2CutIn,
  bool isBelowHadIn, vector<double> evolutionWindowsIn, AlphaEM alIn) {

  if (!isInit) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__+": Not initialised.");
    return;
  }

  // Verbose output.
  if (verbose >= louddebug) printOut(__METHOD_NAME__, "begin --------------");

  // Input.
  iSys = iSysIn;
  shh = infoPtr->s();
  q2Cut = q2CutIn;
  isBelowHad = isBelowHadIn;
  evolutionWindows = evolutionWindowsIn;
  al = alIn;

  // Build internal system.
  buildSystem(event);
  if (verbose >= louddebug) printOut(__METHOD_NAME__, "end --------------");

}

//--------------------------------------------------------------------------

// Trial antenna function.

double QEDemitSystem::aTrial(QEDemitElemental* ele, double sxj, double syj,
  double sxy) {
  int idx = ele->idx;
  int idy = ele->idy;
  double ant = 0;

  // FF.
  if (ele->isFF || ele->isDip) {
    double s = sxj + syj + sxy;
    ant += 4*s/sxj/syj;
    if (ele->isFF && abs(idx) == 24) ant += 8.*s/sxj/(s - syj)/3.;
    if (ele->isFF && abs(idy) == 24) ant += 8.*s/syj/(s - sxj)/3.;
  }

  // IF.
  if (ele->isIF) {
    double s = sxj + sxy - syj;
    ant += 4*pow2(s+syj)/(s*sxj*syj);
    if (abs(idy) == 24) ant += 8.*(s + syj)/syj/(s + syj - sxj)/3.;
  }

  // II.
  if (ele->isII) {
    double s = sxy - sxj - syj;
    ant += 4*sxy*sxy/s/sxj/syj;
  }

  // RF.
  if (ele->isRF) {
    double s = sxj + sxy - syj;
    ant += 4*pow2(s+syj)/s/sxj/syj;
    if (abs(idx) == 24) ant += 8*(2.*syj/s + pow2(syj)/pow2(s))/sxj/3.;
    if (abs(idy) == 24) ant += 8.*(s + syj)/syj/(s + syj - sxj)/3.;
  }
  return ant;

}

//--------------------------------------------------------------------------

// Physical antenna function.

double QEDemitSystem::aPhys(QEDemitElemental* ele, double sxj, double syj,
  double sxy) {
  double mx2 = ele->mx2;
  double my2 = ele->my2;
  int idx = ele->idx;
  int idy = ele->idy;
  double ant = 0;

  // FF.
  if (ele->isFF) {
    double s = sxj + syj + sxy;
    // Eikonal.
    ant += 4.*sxy/sxj/syj - 4.*mx2/sxj/sxj - 4.*my2/syj/syj;

    // Check if x is a W or a fermion.
    if (abs(idx) == 24 && useFullWkernel)
      ant += (4./3.)*(syj/(s - syj) + syj*(s - syj)/s/s)/sxj;
    else
      ant += 2.*syj/sxj/s;

    // Check if y is a W or a fermion.
    if (abs(idy) == 24 && useFullWkernel)
      ant += (4./3.)*(sxj/(s - sxj) + sxj*(s - sxj)/s/s)/syj;
    else
      ant += 2.*sxj/syj/s;
  }

  // FF (dipole).
  if (ele->isDip) {
    double s = sxj + syj + sxy;
    ant += 4.*sxy/sxj/(sxj+syj) - 4.*mx2/sxj/sxj + 2.*syj/sxj/s;
  }

  // IF.
  if (ele->isIF) {
    double s = sxj + sxy - syj;
    // Eikonal + initial state fermion.
    // The initial state is never a W and has no mass.
    ant += 4.*sxy/sxj/syj - 4.*my2/syj/syj + 2.*syj/sxj/s;

    if (abs(idy) == 24 && useFullWkernel)
      ant += (8./3.)*( sxj/(sxy + syj) + sxj/(s + syj)
        - pow2(sxj)/pow2(s + syj) )/syj;
    else
      ant += 2.*sxj/s/syj;
  }

  // II.
  if (ele->isII) {
    double s = sxy - sxj - syj;
    // Eikonal + fermion.
    ant = 4.*sxy/sxj/syj + 2.*(sxj/syj + syj/sxj)/s;
  }

  // RF.
  if (ele->isRF) {
    double s = sxj + sxy - syj;
    // Eikonal.
    ant = 4.*sxy/sxj/syj - 4.*mx2/sxj/sxj - 4.*my2/syj/syj;

    // Check if x is a W or a fermion
    if (abs(idx) == 24 && useFullWkernel)
      ant += (8./3.)*( syj/(s+syj) + syj/s + pow2(syj)/pow2(s) )/sxj;
    else
      ant += 2.*syj/sxj/s;

    // Check if y is a W or a fermion.
    if (abs(idy) == 24 && useFullWkernel)
      ant += (8./3.)*( sxj/(sxy + syj) + sxj/(s + syj)
          - pow2(sxj)/pow2(s + syj) )/syj;
    else
      ant += 2.*sxj/syj/s;
  }
  return ant;

}

//--------------------------------------------------------------------------

// Ratio between PDFs.

double QEDemitSystem::PDFratio(bool isA, double eOld, double eNew, int id,
  double Qt2) {
  double xOld = eOld/(sqrt(shh)/2.0);
  double xNew = eNew/(sqrt(shh)/2.0);
  double newPDF, oldPDF;
  if (isA) {
    newPDF = beamAPtr->xfISR(iSys, id, xNew, Qt2)/xNew;
    oldPDF = beamAPtr->xfISR(iSys, id, xOld, Qt2)/xOld;
    if (abs(newPDF) < TINYPDF) newPDF = TINYPDF;
    if (abs(oldPDF) < TINYPDF) oldPDF = TINYPDF;
  } else {
    newPDF = beamBPtr->xfISR(iSys, id, xNew, Qt2)/xNew;
    oldPDF = beamBPtr->xfISR(iSys, id, xOld, Qt2)/xOld;
    if (abs(newPDF) < TINYPDF) newPDF = TINYPDF;
    if (abs(oldPDF) < TINYPDF) oldPDF = TINYPDF;
  }
  return newPDF/oldPDF;
}

//--------------------------------------------------------------------------

// Set up antenna pairing for incoherent mode.

void QEDemitSystem::buildSystem(Event &event) {

  // Clear previous antennae.
  eleVec.clear();
  eleMat.clear();
  iCoh.clear();

  // Construct hungarian algorithm solver.
  HungarianAlgorithm ha;
  // Below hadronization scale.
  if (isBelowHad && emitBelowHad) {
    map<int, vector<int> > posMap, negMap;
    vector<Vec4> posMoms, negMoms;

    // Find all (final-state) quarks and leptons.
    vector<int> iQuarks, iLeptons;
    int sysSize = partonSystemsPtr->sizeOut(iSys);
    for (int i = 0; i < sysSize; i++) {
      int iEv = partonSystemsPtr->getOut(iSys, i);
      if (event[iEv].col() != 0 && event[iEv].acol()==0 &&
        event[iEv].isFinal()) {
        // For now, ignore quarks that are connected to junctions. In
        // principle, we could add them, and any antijunction dittos.
        bool isJun = false;
        for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
          for (int iLeg = 0; iLeg < 3; ++iLeg) {
            if (event[iEv].col() == event.endColJunction(iJun,iLeg)) {
              isJun = true;
              break;
            }
          }
        }
        if (!isJun) iQuarks.push_back(iEv);
      }
      if (event[iEv].isLepton() && event[iEv].isCharged())
        iLeptons.push_back(iEv);
    }

    // Currently no showering below hadronisation scale if no leptons.
    if (iLeptons.size() == 0) return;

    // Sort all leptons into maps.
    for (int i = 0; i < (int)iLeptons.size(); i++) {
      int iEv = iLeptons[i];
      vector<int> iLeptonVec;
      iLeptonVec.push_back(iEv);
      if (event[iEv].chargeType() == 3) {
        posMoms.push_back(event[iEv].p());
        posMap[posMoms.size()-1] = iLeptonVec;
      }
      if (event[iEv].chargeType() == -3) {
        negMoms.push_back(event[iEv].p());
        negMap[negMoms.size()-1] = iLeptonVec;
      }
    }
    // Find all colour strings.
    for (int i = 0; i < (int)iQuarks.size(); i++) {
      // Get initial quark and add to pseudo particle.
      Vec4 pPseudo;
      int iEv = iQuarks[i];
      vector<int> iPseudoVec;
      iPseudoVec.push_back(iEv);
      pPseudo += event[iEv].p();

      // Find next colour-connected particle.
      do {
        int colTag = event[iEv].col();
        for (int j = 0; j < sysSize; j++) {
          int jEv = partonSystemsPtr->getOut(iSys, j);
          if (event[jEv].acol() == colTag && event[jEv].isFinal()) {
            iEv = jEv;
            break;
          }
        }
        if (iEv == iPseudoVec.back()) {
          infoPtr->errorMsg("Error in "+__METHOD_NAME__
            +": Colour tracing failed.");
          break;
        }
        iPseudoVec.push_back(iEv);
        pPseudo += event[iEv].p();
      }
      while(!event[iEv].isQuark()&&!event[iEv].isDiquark());

      // Get charge of pseudoparticle and sort into maps.
      int chargeTypePseudo = event[iPseudoVec.front()].chargeType()
        + event[iPseudoVec.back()].chargeType();
      // Strings with only quarks are total charge 1 or -1.
      if (chargeTypePseudo == 3) {
        posMoms.push_back(pPseudo);
        posMap[posMoms.size()-1] = iPseudoVec;
      }
      if (chargeTypePseudo == -3) {
        negMoms.push_back(pPseudo);
        negMap[negMoms.size()-1] = iPseudoVec;
      }
      // Strings with a diquark can be charge 2 or -2. Add these
      // twice to list of recoilers.
      if (chargeTypePseudo == 6) {
        posMoms.push_back(pPseudo);
        posMap[posMoms.size()-1] = iPseudoVec;
        posMoms.push_back(pPseudo);
        posMap[posMoms.size()-1] = iPseudoVec;
      }
      if (chargeTypePseudo == -6) {
        negMoms.push_back(pPseudo);
        negMap[negMoms.size()-1] = iPseudoVec;
        negMoms.push_back(pPseudo);
        negMap[negMoms.size()-1] = iPseudoVec;
      }
    }

    // If no leptons and overall hadronic system has charge = 0, do nothing.
    if (posMoms.size() == 0) return;

    // Solve assignment problem.
    vector<vector<double> > weights;
    weights.resize(posMoms.size());
    for (int i=0; i<(int)posMoms.size(); i++) {
      weights[i].resize(negMoms.size());
      for (int j=0; j<(int)negMoms.size(); j++) {
        double w = posMoms[i]*negMoms[j]
          - posMoms[i].mCalc()*negMoms[j].mCalc();
        weights[i][j] = w;
      }
    }
    vector<int> assignment;
    ha.solve(weights, assignment);

    for (int i = 0; i < (int)posMoms.size(); i++) {
      int iPos = i;
      int iNeg = assignment[i];
      // Only keep antennae with at least one lepton.
      if (posMap[iPos].size() == 1 || negMap[iNeg].size() == 1) {
        eleVec.push_back(QEDemitElemental());
        eleVec.back().initPtr(rndmPtr, partonSystemsPtr);
        // If two leptons, add regular antenna.
        if (posMap[iPos].size() == 1 && negMap[iNeg].size() == 1)
          eleVec.back().init(event, posMap[iPos][0], negMap[iNeg][0], shh,
            verbose);
        // If lepton + pseudoparticle, add dipole.
        if (posMap[iPos].size() == 1 && negMap[iNeg].size() != 1)
          eleVec.back().init(event, posMap[iPos][0], negMap[iNeg], shh,
            verbose);
        if (posMap[iPos].size()!=1 && negMap[iNeg].size()==1)
          eleVec.back().init(event, negMap[iNeg][0], posMap[iPos], shh,
            verbose);
      }
    }

  // Above hadronization scale.
  } else if(!isBelowHad) {
    // Collect relevant particles.
    int sysSize = partonSystemsPtr->sizeAll(iSys);
    for (int i = 0; i < sysSize; i++) {
      int iEv = partonSystemsPtr->getAll(iSys, i);
      if (event[iEv].isCharged()) iCoh.push_back(iEv);
    }

    // Catch cases (like hadron->partons decays) where an explicit
    // charged mother may not have been added to the partonSystem as a
    // resonance.
    if (partonSystemsPtr->getInA(iSys) == 0 &&
        partonSystemsPtr->getInB(iSys) == 0 &&
        partonSystemsPtr->getInRes(iSys) == 0) {
      // Guess that the decaying particle is mother of first parton.
      int iRes = event[partonSystemsPtr->getOut(iSys, 0)].mother1();
      if (iRes != 0 && event[iRes].isCharged()) {
        // Check daughter list consistent with whole system.
        int ida1 = event[iRes].daughter1();
        int ida2 = event[iRes].daughter2();
        if (ida2 > ida1) {
          bool isOK = true;
          for (int i=0; i<partonSystemsPtr->sizeOut(iSys); ++i)
            if (partonSystemsPtr->getOut(iSys,i) < ida1
              || partonSystemsPtr->getOut(iSys,i) > ida2) isOK = false;
          if (isOK) {iCoh.push_back(iRes);}
        }
      }
    }

    // First check charge conservation.
    int chargeTypeTot = 0;
    for (int i = 0; i < (int)iCoh.size(); i++) {
      double cType = event[iCoh[i]].chargeType();
      chargeTypeTot += (event[iCoh[i]].isFinal() ? cType : -cType);
    }

    if (chargeTypeTot != 0) {
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": Charge not conserved above hadronization scale");
      if (verbose >= superdebug) {
        printOut(__METHOD_NAME__, "Printing events and systems");
        event.list();
        partonSystemsPtr->list();
      }
    }

    // Pairing algorithm.
    if (mode == 1) {
      vector<vector<int> > posChargeTypes;
      posChargeTypes.resize(3);
      vector<vector<int> > negChargeTypes;
      negChargeTypes.resize(3);

      for (int i = 0; i < (int)iCoh.size(); i++) {
        int iEv = iCoh[i];
        // Separate particles into charge types.
        double Q = event[iEv].charge();
        // Get index in pos/negChargeTypes.
        int n = abs(event[iEv].chargeType()) - 1;
        // Flip charge contribution of initial state.
        if (!event[iEv].isFinal()) {Q = -Q;}
        if (Q > 0)  posChargeTypes[n].push_back(iEv);
        else negChargeTypes[n].push_back(iEv);
      }

      // Clear list of charged particles.
      iCoh.clear();

      // Solve assignment problems.
      for (int i=0; i<3; i++) {
        int posSize = posChargeTypes[i].size();
        int negSize = negChargeTypes[i].size();
        int maxSize = max(posSize,negSize);
        if (maxSize > 0) {
          vector<vector<double> > weights;
          weights.resize(maxSize);
          // Set up matrix of weights.
          for (int x = 0; x < maxSize; x++) {
            weights[x].resize(maxSize);
            for (int y = 0; y < maxSize; y++) {
              // If either index is out of range. Add some random
              // large weight.
              double wIn = (0.9 + 0.2*rndmPtr->flat())*1E300;
              if (x < posSize && y < negSize) {
                int xEv = posChargeTypes[i][x];
                int yEv = negChargeTypes[i][y];
                wIn = event[xEv].p()*event[yEv].p()
                  - event[xEv].m()*event[yEv].m();
              }
              weights[x][y] = wIn;
            }
          }

          // Find solution.
          vector<int> assignment;
          ha.solve(weights, assignment);

          // Add pairings to list of emitElementals.
          // Add unpaired particles to index list for coherent algorithm.
          for (int j = 0; j < maxSize; j++) {
            int x = j;
            int y = assignment[j];
            if (x < posSize && y < negSize) {
              int xEv = posChargeTypes[i][x];
              int yEv = negChargeTypes[i][y];
              eleVec.push_back(QEDemitElemental());
              eleVec.back().initPtr(rndmPtr, partonSystemsPtr);
              eleVec.back().init(event, xEv, yEv, shh, verbose);
            } else if (x < posSize) {
              int xEv = posChargeTypes[i][x];
              iCoh.push_back(xEv);
            } else if (y < negSize) {
              int yEv = negChargeTypes[i][y];
              iCoh.push_back(yEv);
            }
          }
        }
      }
    }

    // Create eleMat.
    eleMat.resize(iCoh.size());
    for (int i = 0; i < (int)iCoh.size(); i++) {
      eleMat[i].resize(i);
      for (int j = 0; j < i; j++) {
        eleMat[i][j].initPtr(rndmPtr, partonSystemsPtr);
        eleMat[i][j].init(event, iCoh[i], iCoh[j], shh, verbose);
      }
    }

    // Compute overestimate constant.
    cMat = 0;
    for (int i = 0; i < (int)eleMat.size(); i++)
      for (int j = 0; j < i; j++) cMat += max(eleMat[i][j].QQ, 0.);
  }

}

//--------------------------------------------------------------------------

// Generate a trial scale.

double QEDemitSystem::generateTrialScale(Event &event, double q2Start) {

  // Check if qTrial is below the cutoff.
  if (q2Start < q2Cut || evolutionWindows.size() == 0) {return 0;}

  // Find lower value from evolution window.
  int iEvol = evolutionWindows.size() - 1;
  while (iEvol >= 1 && q2Start <= evolutionWindows[iEvol]) iEvol--;
  double q2Low = evolutionWindows[iEvol];
  if (q2Low < 0)
    infoPtr->errorMsg("Error in "+__METHOD_NAME__+": Evolution window < 0");
  double q2Trial = 0;

  // Generate a scale.
  double alphaMax = al.alphaEM(q2Start);

  // Pull scales from eleVec.
  for (int i = 0; i < (int)eleVec.size(); i++) {
    double c = eleVec[i].QQ;
    double q2New = eleVec[i].generateTrial(event, q2Start, q2Low, alphaMax, c);
    if (q2New > q2Low && q2New > q2Trial) {
      q2Trial = q2New;
      eleTrial = &eleVec[i];
      trialIsVec = true;
    }
  }

  // Pull scales from eleMat.
  for (int i = 0; i < (int)eleMat.size(); i++) {
    for (int j = 0; j < i; j++) {
      double q2New = eleMat[i][j].generateTrial(event, q2Start, q2Low,
        alphaMax, cMat);
      if (q2New > q2Low && q2New > q2Trial) {
        q2Trial = q2New;
        eleTrial = &eleMat[i][j];
        trialIsVec = false;
      }
    }
  }

  // Check if evolution window was crossed.
  if (q2Trial < q2Low) {
    if (iEvol == 0) {return 0;}
    // Reset all trials.
    for (int i = 0; i < (int)eleVec.size(); i++) eleVec[i].hasTrial = false;
    for (int i=0; i<(int)eleMat.size(); i++)
      for (int j=0; j<i; j++) eleMat[i][j].hasTrial = false;
    return generateTrialScale(event, q2Low);
  }

  // Otherwise return trial scale.
  return q2Trial;

}

//--------------------------------------------------------------------------

// Check the veto.

bool QEDemitSystem::checkVeto(Event &event) {

  // Mark trial as used.
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");
  eleTrial->hasTrial = false;

  // Pre- and post-branching momenta.
  vector<Vec4> pOld, pNew;

  // Global recoil momenta.
  vector<Vec4> pRec;
  vector<int>  iRec;

  // II.
  if (eleTrial->isII) {
    double saj = eleTrial->sxjSav;
    double sbj = eleTrial->syjSav;
    double phi = eleTrial->phiSav;
    double sAB = eleTrial->sAnt;
    double sab = sAB + saj + sbj;

    // Pre-branching momenta.
    pOld.push_back(event[eleTrial->x].p());
    pOld.push_back(event[eleTrial->y].p());

    // Collect the recoiling final state particles.
    int sysSize = partonSystemsPtr->sizeAll(iSys);
    for (int i = 0; i < sysSize; i++) {
      int iEv = partonSystemsPtr->getAll(iSys, i);
      if (iEv < 0 || !event[iEv].isFinal()) continue;
      if (iEv == eleTrial->x || iEv == eleTrial->y) continue;
      pRec.push_back(event[iEv].p());
      iRec.push_back(iEv);
    }

    // Kinematics.
    if (!vinComPtr->map2to3IImassless(pNew, pRec, pOld, sAB,saj,sbj,sab, phi))
      return false;

    // Check if new energies don't exceed hadronic maxima.
    double eaUsed = 0, ebUsed = 0;
    int nSys = partonSystemsPtr->sizeSys();
    for (int i = 0; i < nSys; i++) {
      eaUsed += event[partonSystemsPtr->getInA(i)].e();
      ebUsed += event[partonSystemsPtr->getInB(i)].e();
    }
    if ((eaUsed - pOld[0].e() + pNew[0].e()) > 0.98*sqrt(shh)/2.) return false;
    if ((ebUsed - pOld[1].e() + pNew[2].e()) > 0.98*sqrt(shh)/2.) return false;
  }

  // IF.
  else if (eleTrial->isIF) {
    double saj = eleTrial->sxjSav;
    double sjk = eleTrial->syjSav;
    double phi = eleTrial->phiSav;
    double sAK = eleTrial->sAnt;
    double sak = sAK + sjk - saj;
    double mK2 = eleTrial->my2;

    // Check phase space.
    if (sak < 0 || saj*sjk*sak - saj*saj*mK2 < 0) {return false;}

    // Pre-branching momenta.
    pOld.push_back(event[eleTrial->x].p());
    pOld.push_back(event[eleTrial->y].p());

    // Kinematics. (TODO: check if need for global kinematics map here).
    if (!vinComPtr->map2to3IFlocal(pNew, pOld, sAK, saj, sjk, sak, phi,
        mK2, 0, mK2)) return false;

    // Check if new energy doesn't exceed the hadronic maximum.
    double eaUsed = 0;
    int nSys = partonSystemsPtr->sizeSys();
    for (int i = 0; i < nSys; i++) {
      int iEv;
      if (eleTrial->isIA) iEv = partonSystemsPtr->getInA(i);
      else iEv = partonSystemsPtr->getInB(i);
      eaUsed += event[iEv].e();
    }
    if ((eaUsed - pOld[0].e() + pNew[0].e()) > 0.98*sqrt(shh)/2.) return false;
  }

  // RF.
  else if (eleTrial->isRF) {
    double saj = eleTrial->sxjSav;
    double sjk = eleTrial->syjSav;
    double sAK = eleTrial->sAnt;
    double sak = sAK + sjk - saj;
    double phi = eleTrial->phiSav;
    double mA2 = eleTrial->mx2;
    double mK2 = eleTrial->my2;

    // Check phase space.
    if (sak < 0 || saj*sjk*sak - saj*saj*mK2 - sjk*sjk*mA2 < 0) return false;

    // Pre-branching momenta.
    pOld.push_back(event[eleTrial->x].p());
    pOld.push_back(event[eleTrial->y].p());

    // Collect the recoiling final state particles.
    int sysSize = partonSystemsPtr->sizeAll(iSys);
    for (int i = 0; i < sysSize; i++) {
      int iEv = partonSystemsPtr->getAll(iSys, i);
      if (iEv < 0 || !event[iEv].isFinal()) {continue;}
      if (iEv == eleTrial->x || iEv == eleTrial->y) {continue;}
      pRec.push_back(event[iEv].p());
      iRec.push_back(iEv);
    }

    // Do kinematics.
    vector<double> masses;
    masses.push_back(sqrt(mA2));
    masses.push_back(0.);
    masses.push_back(sqrt(mK2));
    masses.push_back(sqrt(mA2+mK2-sAK));
    vector<double> invariants;
    invariants.push_back(sAK);
    invariants.push_back(saj);
    invariants.push_back(sjk);
    invariants.push_back(sak);
    vector<Vec4> pAfter;
    vector<Vec4> pBefore = pOld;
    pBefore.insert(pBefore.end(), pRec.begin(), pRec.end());
    if (!vinComPtr->map2toNRFmassive(pAfter, pBefore, 0, 1, invariants, phi,
        masses)) return false;
    pNew.push_back(pAfter[0]);
    pNew.push_back(pAfter[1]);
    pNew.push_back(pAfter[2]);

    // Replace momenta with boosted counterpart.
    pRec.clear();
    for (int i = 3; i < (int)pAfter.size(); i++) pRec.push_back(pAfter[i]);

    // Check if nothing got messed up along the way.
    if (pRec.size() != iRec.size()) {
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": inconsistent recoilers in RF kinematics.");
      return false;
    }
  }

  // FF>
  else if (eleTrial->isFF) {
    double sIK = eleTrial->m2Ant - eleTrial->mx2 - eleTrial->my2;
    double sij = eleTrial->sxjSav;
    double sjk = eleTrial->syjSav;
    double sik = sIK - sij - sjk;
    double mi  = sqrt(eleTrial->mx2);
    double mk  = sqrt(eleTrial->my2);
    double phi = eleTrial->phiSav;

    vector<double> invariants;
    invariants.push_back(sIK);
    invariants.push_back(sij);
    invariants.push_back(sjk);

    vector<double> masses;
    masses.push_back(mi);
    masses.push_back(0);
    masses.push_back(mk);

    // Check phase space.
    if (sik < 0) return false;
    if (sij*sjk*sik - pow2(sij)*pow2(mk) - pow2(sjk)*pow2(mi) < 0)
      return false;

    // Pre-branching momenta.
    pOld.push_back(event[eleTrial->x].p());
    pOld.push_back(event[eleTrial->y].p());

    // Kinematics.
    if (!vinComPtr->map2to3FF(pNew, pOld, 3, invariants, phi, masses))
      return false;
  }

  // Dipole.
  else if (eleTrial->isDip) {
    // Construct recoiler momentum.
    Vec4 pk;
    for (int i = 0; i < (int)eleTrial->iRecoil.size(); i++)
      pk += event[eleTrial->iRecoil[i]].p();
    double sIK = eleTrial->m2Ant - eleTrial->mx2 - eleTrial->my2;
    double sij = eleTrial->sxjSav;
    double sjk = eleTrial->syjSav;
    double sik = sIK - sij - sjk;
    double mi  = sqrt(eleTrial->mx2);
    double mk  = pk.mCalc();
    double phi = eleTrial->phiSav;

    vector<double> invariants;
    invariants.push_back(sIK);
    invariants.push_back(sij);
    invariants.push_back(sjk);

    vector<double> masses;
    masses.push_back(mi);
    masses.push_back(0);
    masses.push_back(mk);

    // Check phase space.
    if (sik < 0) {return false;}
    if (sij*sjk*sik - pow2(sij)*pow2(mk) - pow2(sjk)*pow2(mi) < 0)
      return false;

    // Pre-branching momenta.
    pOld.push_back(event[eleTrial->x].p());
    pOld.push_back(pk);

    // Kinematics.
    if (!vinComPtr->map2to3FF(pNew, pOld, 3, invariants, phi, masses))
      return false;
  }
  Vec4 pPhot = pNew[1];
  Vec4 px = pNew[0];
  Vec4 py = pNew[2];
  int x = eleTrial->x;
  int y = eleTrial->y;
  double sxj = eleTrial->sxjSav;
  double syj = eleTrial->syjSav;
  double sxy = px*py*2.;

  // Compute veto probability.
  double pVeto = 1.;

  // Add alpha veto.
  pVeto *= al.alphaEM(eleTrial->q2Sav) / eleTrial->alpha;
  if (pVeto > 1) printOut(__METHOD_NAME__, "Alpha increased");

  // Add antenna veto. Simple veto for eleTrial in eleVec.
  if (trialIsVec) {
    // Note that charge factor is included at generation step.
    double aTrialNow = aTrial(eleTrial, sxj, syj, sxy);
    double aPhysNow = aPhys(eleTrial, sxj, syj, sxy);

    if (aPhysNow/aTrialNow > 1.001) {
      stringstream ss1;
      ss1 << "Incorrect overestimate (eleVec) " << aPhysNow/aTrialNow;
      printOut(__METHOD_NAME__, ss1.str());
      if (verbose > louddebug) {
        stringstream ss2, ss3;
        if (eleTrial->isFF) ss2 << "Antenna is FF";
        if (eleTrial->isIF) ss2 << "Antenna is IF";
        if (eleTrial->isRF) ss2 << "Antenna is RF";
        if (eleTrial->isII) ss2 << "Antenna is II";
        printOut(__METHOD_NAME__, ss2.str());
        printOut(__METHOD_NAME__, ss3.str());
      }
    }
    pVeto *= aPhysNow/aTrialNow;

  // Construct full branching kernel for eleTrial in eleMat. Perform
  // sector check too.
  } else {
    double aSectorNow = aPhys(eleTrial, sxj, syj, sxy);
    double aTrialFull = eleTrial->c*aTrial(eleTrial, sxj, syj, sxy);
    double aPhysFull  = 0;

    // Build map of momenta & invariants with new photon.
    map<int, double> sj;
    map<int, Vec4> p;
    // Loop over the first column in eleMat.
    for (int i = 0; i < (int)iCoh.size(); i++) {
      int iEv = iCoh[i];
      // If the particle is in eleTrial, use shower variables.
      if (iEv == x) {
        p[iEv] = px;
        sj[iEv] = sxj;
      } else if (iEv == y) {
        p[iEv] = py;
        sj[iEv] = syj;
      // Otherwise get the momenta elsewhere
      } else {
        // If global recoil, get them from pRec.
        if (eleTrial->isII) {
          // Find index.
          for (int j = 0; j < (int)iRec.size(); j++) {
            if (iEv == iRec[j]) {
              p[iEv] = pRec[j];
              sj[iEv] = 2.*pRec[j]*pPhot;
              break;
            }
          }
        // Otherwise use momentum from event.
        } else {
          p[iEv] = event[iEv].p();
          sj[iEv] = 2.*event[iEv].p()*pPhot;
        }
      }
    }

    // Then build aPhys.
    for (int v=0; v<(int)eleMat.size(); v++) {
      for (int w=0; w<v; w++) {
        double sxjNow = sj[eleMat[v][w].x];
        double syjNow = sj[eleMat[v][w].y];
        double sxyNow = 2.*p[eleMat[v][w].x]*p[eleMat[v][w].y];
        double aPhysNow = aPhys(&eleMat[v][w], sxjNow, syjNow, sxyNow);

        // Sector veto.
        if (aPhysNow > aSectorNow) return false;

        // Add aPhysNew to aPhys.
        aPhysFull += eleMat[v][w].QQ*aPhysNow;
      }
    }
    // Set aPhys to zeto if below zero.
    if (aPhysFull < 0) {aPhysFull = 0;}

    // Check overestimate.
    if (aPhysFull/aTrialFull > 1) {
      stringstream ss1;
      ss1 << "Incorrect overestimate (eleVec) " << aPhysFull/aTrialFull;
      printOut(__METHOD_NAME__, ss1.str());
    }
    // Add antenna veto.
    pVeto *= aPhysFull/aTrialFull;
  }

  // Add PDF veto.
  if (eleTrial->isIF) {
    pVeto *= PDFratio(eleTrial->isIA, pOld[0].e(), pNew[0].e(),
      eleTrial->idx, eleTrial->q2Sav);
  }
  if (eleTrial->isII) {
    pVeto *= PDFratio(true,  pOld[0].e(), pNew[0].e(),
      eleTrial->idx, eleTrial->q2Sav);
    pVeto *= PDFratio(false, pOld[1].e(), pNew[2].e(),
      eleTrial->idy, eleTrial->q2Sav);
  }

  // Perform veto.
  if (rndmPtr->flat() > pVeto) return false;

  // Accepted, now fix the event. Different procedures for dipole and
  // antennae.

  // If it is a dipole.
  if (eleTrial->isDip) {
    // Set up new particles.
    Particle newPhoton(22, 51, 0, 0, 0, 0, 0, 0, pPhot);
    Particle newPartx = event[x];
    newPartx.p(px);

    // Add to event.
    int xNew = event.append(newPartx);
    int jNew = event.append(newPhoton);

    // Update parton system.
    partonSystemsPtr->replace(iSys, x, xNew);
    partonSystemsPtr->addOut(iSys, jNew);

    // Set old particles to negative.
    event[x].statusNeg();

    // Update mother-daughter structure.
    event[xNew].mothers(x,0);
    event[jNew].mothers(x,0);
    event[x].daughters(xNew, jNew);
    event[xNew].daughters(0,0);
    event[jNew].daughters(0,0);
    event[xNew].statusCode(51);
    event[jNew].statusCode(51);

    // Boost momenta and update.
    for (int i = 0; i < (int)eleTrial->iRecoil.size(); i++) {
      int iDipRec = eleTrial->iRecoil[i];
      Vec4 pDipRec = event[iDipRec].p();
      pDipRec.bstback(pOld[1]);
      pDipRec.bst(pNew[2]);
      // Copy the recoiler.
      int iDipRecNew = event.copy(iDipRec, 52);
      // Change the momentum.
      event[iDipRecNew].p(pDipRec);
      // Update parton system.
      partonSystemsPtr->replace(iSys, iDipRec, iDipRecNew);
    }

  } else if (eleTrial->isRF) {

    // Set up new particles.
    Particle newPhoton(22, 51, 0, 0, 0, 0, 0, 0, pPhot);
    Particle newParty = event[y];
    newParty.p(py);

    // Add branched particles to event.
    int jNew, yNew;
    if (x < y) {
      jNew = event.append(newPhoton);
      yNew = event.append(newParty);
    } else {
      yNew = event.append(newParty);
      jNew = event.append(newPhoton);
    }

    // Update parton system.
    partonSystemsPtr->addOut(iSys, jNew);
    partonSystemsPtr->replace(iSys, y, yNew);

    // Set old particles to negative.
    event[y].statusNeg();
    event[jNew].mothers(y,0);
    event[yNew].mothers(y,0);
    event[jNew].daughters(0,0);
    event[yNew].daughters(0,0);
    event[y].daughters(yNew,jNew);
    event[yNew].statusCode(51);

    // Update event for global recoil.
    vector<pair<int,int> > iRecNew;
    iRecNew.clear();
    iRecNew.resize(0);
    for (int j = 0; j < event.size(); j++)
      if (event[j].isFinal())
        for (int k = 0; k < (int)iRec.size(); k++)
          if (iRec[k] == j) {
            // Copy the recoiler.
            int inew = event.copy(j, 52);
            // Change the momentum.
            event[inew].p(pRec[k]);
            iRecNew.push_back(make_pair(iRec[k], inew));
          }

    // Update parton system.
    for (int k=0; k<(int)iRecNew.size(); k++)
      partonSystemsPtr->replace(iSys, iRecNew[k].first, iRecNew[k].second);

  // If it is an antenna
  } else {
    // Set up new particles.
    Particle newPhoton(22, 51, 0, 0, 0, 0, 0, 0, pPhot);
    Particle newPartx = event[x];
    newPartx.p(px);
    Particle newParty = event[y];
    newParty.p(py);

    // Add branched particles to event.
    int xNew, jNew, yNew;
    if (x < y) {
      xNew = event.append(newPartx);
      jNew = event.append(newPhoton);
      yNew = event.append(newParty);
    } else {
      yNew = event.append(newParty);
      jNew = event.append(newPhoton);
      xNew = event.append(newPartx);
    }

    // Update parton system.
    partonSystemsPtr->replace(iSys, x, xNew);
    partonSystemsPtr->addOut(iSys, jNew);
    partonSystemsPtr->replace(iSys, y, yNew);

    // Set old particles to negative.
    event[x].statusNeg();
    event[y].statusNeg();

    // Update everything.
    if (eleTrial->isII) {
      event[xNew].mothers(event[x].mother1(), event[x].mother2());
      event[jNew].mothers(xNew, yNew);
      event[yNew].mothers(event[y].mother1(), event[y].mother2());
      event[x].mothers(xNew, 0);
      event[y].mothers(yNew, 0);
      event[xNew].daughters(jNew, x);
      event[yNew].daughters(jNew, y);
      event[jNew].daughters(0,0);
      event[xNew].status(-41);
      event[yNew].status(-41);
      event[jNew].status(43);

      // Update beam daughters.
      if (iSys == 0) {
        bool founda = false;
        bool foundb = false;
        for (int i = 0; i < (int)event.size(); i++) {
          if (!founda)
            if (event[i].daughter1() == x) {
              event[i].daughters(xNew, 0);
              founda = true;
            }
          if (!foundb)
            if (event[i].daughter1() == y) {
              event[i].daughters(yNew, 0);
              foundb = true;
            }
          if (founda && foundb) break;
        }
      }

      // Update event for global recoil.
      vector<pair<int,int> > iRecNew;
      iRecNew.clear();
      iRecNew.resize(0);
      for (int j=0; j<event.size(); j++)
        if (event[j].isFinal())
          for (int k=0; k<(int)iRec.size(); k++)
            if (iRec[k] == j) {
              // Copy the recoiler.
              int inew = event.copy(j, 44);
              // Change the momentum.
              event[inew].p(pRec[k]);
              iRecNew.push_back(make_pair(iRec[k], inew));
            }


      // Update parton system.
      for (int k = 0; k < (int)iRecNew.size(); k++)
        partonSystemsPtr->replace(iSys, iRecNew[k].first, iRecNew[k].second);
      partonSystemsPtr->setInA(iSys, xNew);
      partonSystemsPtr->setInB(iSys, yNew);

      // Update beams.
      BeamParticle& beam1 = *beamAPtr;
      BeamParticle& beam2 = *beamBPtr;

      // Check that x is always a with pz>0.
      if (event[xNew].pz() < 0) {
        printOut(__METHOD_NAME__, "Swapped II  antenna");
        beam1 = *beamBPtr;
        beam2 = *beamAPtr;
      }
      beam1[iSys].update(xNew, event[xNew].id(), event[xNew].e()/beam1.e());
      beam2[iSys].update(yNew, event[yNew].id(), event[yNew].e()/beam2.e());
    }

    if (eleTrial->isIF) {
      event[xNew].mothers(event[x].mother1(), event[x].mother2());
      event[jNew].mothers(y,xNew);
      event[yNew].mothers(y,0);
      event[x].mothers(xNew,0);
      event[xNew].daughters(jNew,x);
      event[jNew].daughters(0,0);
      event[yNew].daughters(0,0);
      event[y].daughters(jNew, yNew);
      event[xNew].status(-41);
      event[yNew].status(43);
      event[jNew].status(43);

      // Update beam daughter.
      if (iSys == 0)
        for (int i=0; i<(int)event.size(); i++)
          if (event[i].daughter1() == x) {
            event[i].daughters(xNew, 0);
            break;
          }

      // Update parton system.
      if (eleTrial->isIA) partonSystemsPtr->setInA(iSys, xNew);
      else partonSystemsPtr->setInB(iSys, xNew);

      // Update beams.
      BeamParticle& beam = (eleTrial->isIA ? *beamAPtr : *beamBPtr);
      beam[iSys].update(xNew, event[xNew].id(), event[xNew].e()/beam.e());
    }

    if (eleTrial->isFF) {
      event[xNew].mothers(x,0);
      event[jNew].mothers(x,y);
      event[yNew].mothers(y,0);
      event[x].daughters(xNew, jNew);
      event[y].daughters(yNew, jNew);
      event[xNew].daughters(0,0);
      event[jNew].daughters(0,0);
      event[yNew].daughters(0,0);
      event[xNew].statusCode(51);
      event[jNew].statusCode(51);
      event[yNew].statusCode(51);
    }

    // Update event pointers.
    event.restorePtrs();

    // Fix sHat for parton system.
    double shat = (event[partonSystemsPtr->getInA(iSys)].p() +
      event[partonSystemsPtr->getInB(iSys)].p()).m2Calc();
    partonSystemsPtr->setSHat(iSys, shat);
  }
  if (verbose >= debug) {
    if (verbose >= superdebug) {
      event.list();
      partonSystemsPtr->list();
    }
    printOut(__METHOD_NAME__, "end --------------");
  }
  return true;

}

//--------------------------------------------------------------------------

// Print the QED emit internal system.

void QEDemitSystem::print() {
  cout << "Printing QEDemit internal system" << endl;
  cout << "Pairing elementals" << endl;
  for (int i = 0; i < (int)eleVec.size(); i++) {
    if (eleVec[i].isDip) {
      cout << "Dipole: x = ";
      cout << eleVec[i].x << " Recoilers: (";
      for (int j = 0; j < (int)eleVec[i].iRecoil.size(); j++) {
        cout << eleVec[i].iRecoil[j] << ", ";
        if (j == (int)eleVec[i].iRecoil.size()-1) cout << ")";
        else cout << ", ";
      }
    } else cout << "Antennae: x = " << eleVec[i].x << ", y = " << eleVec[i].y;
    cout << ", QQ = " << eleVec[i].QQ << ", s = " << eleVec[i].sAnt << endl;
  }
  cout << "Coherent elementals" << endl;
  for (int i = 0; i < (int)eleMat.size(); i++)
    for (int j = 0; j < i; j++)
      cout << "x = " << eleMat[i][j].x << ", y = " << eleMat[i][j].y
           << ", QQ = " << eleMat[i][j].QQ << ", s = " << eleMat[i][j].sAnt
           << endl;
}

//==========================================================================

// Class for a QED splitting system.

//--------------------------------------------------------------------------

// Initialize pointers.

void QEDsplitSystem::initPtr(Info* infoPtrIn, VinciaCommon* vinComPtrIn) {
  infoPtr       = infoPtrIn;
  particleDataPtr  = infoPtr->particleDataPtr;
  partonSystemsPtr = infoPtr->partonSystemsPtr;
  rndmPtr          = infoPtr->rndmPtr;
  settingsPtr      = infoPtr->settingsPtr;
  vinComPtr = vinComPtrIn;
  isInitPtr = true;
}

//--------------------------------------------------------------------------

// Initialize.

void QEDsplitSystem::init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  int verboseIn) {
  if (!isInitPtr) printOut(__METHOD_NAME__, "initPtr not called");
  verbose = verboseIn;
  q2Max   = pow2(settingsPtr->parm("Vincia:mMaxGamma"));
  nLepton = settingsPtr->mode("Vincia:nGammaToLepton");
  nQuark  = settingsPtr->mode("Vincia:nGammaToQuark");
  beamAPtr = beamAPtrIn;
  beamBPtr = beamBPtrIn;
  isInit = true;
}

//--------------------------------------------------------------------------

// Prepare list of final-state photons - with recoilers - for splittings.

void QEDsplitSystem::prepare(int iSysIn, Event &event, double q2CutIn,
  bool isBelowHadIn, vector<double> evolutionWindowsIn, AlphaEM alIn) {

  if (!isInit) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__+": Not initialised.");
    return;
  }
  if (verbose >= louddebug) printOut(__METHOD_NAME__, "begin --------------");

  // Input.
  iSys = iSysIn;
  q2Cut = q2CutIn;
  isBelowHad = isBelowHadIn;
  evolutionWindows = evolutionWindowsIn;
  al = alIn;

  // Set up weights for splitting flavours.
  ids.clear();
  idWeights.clear();
  totIdWeight = 0;
  maxIdWeight = 0;

  // Splittings for gamma->lepton+lepton-.
  for (int i = 0; i < nLepton; i++) {
    ids.push_back(11 + 2*i);
    idWeights.push_back(1);
  }
  // Only include gamma->qqbar if above hadronisation scale.
  if (!isBelowHad) {
    for (int i = 1; i <= nQuark; i++) {
      ids.push_back(i);
      idWeights.push_back((i%2==0 ? 4./3. : 1./3.));
    }
  }
  // Total weights.
  for (int i=0; i<(int)ids.size(); i++) {
    totIdWeight += idWeights[i];
    if (idWeights[i] > maxIdWeight) maxIdWeight = idWeights[i];
  }

  // Build internal system.
  buildSystem(event);
  if (verbose >= louddebug) printOut(__METHOD_NAME__, "end --------------");

}

//--------------------------------------------------------------------------

// Build the splitting system.

void QEDsplitSystem::buildSystem(Event &event) {

  // Get rid of saved trial and clear all antennae.
  hasTrial = false;
  eleVec.clear();

  // Build lists of particles.
  vector<int> photList, chSpecList, uchSpecList;
  int sysSize = partonSystemsPtr->sizeAll(iSys);
  for (int i = 0; i < sysSize; i++) {
    int iEv = partonSystemsPtr->getAll(iSys, i);
    if (iEv > 0) {
      // Only involve final state particles.
      if (event[iEv].isFinal()) {
        // Find photons.
        if (event[iEv].id()==22)    photList.push_back(iEv);
        // Find recoilers.
        if (event[iEv].isCharged()) chSpecList.push_back(iEv);
        else                        uchSpecList.push_back(iEv);
      }
    }
  }

  // If no charged and no uncharged spectators, return.
  if (chSpecList.empty() && uchSpecList.empty()) return;

  // Loop over photons.
  for (int i = 0; i < (int)photList.size(); i++) {
    int iPhot = photList[i];
    // If no charged spectators, use uncharged.
    if (chSpecList.empty()) {
      // Check if there is another spectator than the current photon.
      bool otherSpecAvail = false;
      for (int j = 0; j < (int)uchSpecList.size(); j++)
        if (uchSpecList[j] != iPhot) {otherSpecAvail = true; break;}
      // Continue to next photon if no spectator is available.
      if (!otherSpecAvail) continue;

      // Select one at random that's not the photon itself.
      int iSpec;
      while (true) {
        iSpec = uchSpecList[rndmPtr->flat()*uchSpecList.size()];
        if (iSpec != iPhot) break;
      }
      eleVec.push_back(QEDsplitElemental(event, iPhot, iSpec));
      eleVec.back().ariWeight = 1.;

    // Else use charged spectators.
    } else {
      double ariNorm = 0;
      vector<QEDsplitElemental> tempEleVec;
      for (int j = 0; j < (int)chSpecList.size(); j++) {
        int iSpec = chSpecList[j];
        tempEleVec.push_back(QEDsplitElemental(event, iPhot, iSpec));
        ariNorm += 1./tempEleVec.back().m2Ant;
      }
      // Set up Ariadne factors.
      for (int j = 0; j < (int)tempEleVec.size(); j++)
        tempEleVec[j].ariWeight = 1./(tempEleVec[j].m2Ant*ariNorm);
      eleVec.insert(eleVec.end(), tempEleVec.begin(), tempEleVec.end());
    }
  }
}

//--------------------------------------------------------------------------

// Generate a scale for the system.

double QEDsplitSystem::generateTrialScale(Event &event, double q2Start) {

  // Return saved trial.
  if (hasTrial) return q2Trial;

  // Check if there are any photons left.
  if (eleVec.size() == 0) return 0;

  // Starting scale - account for cut on mGammaMax.
  q2Trial = min(q2Max, q2Start);

  // Check if qTrial is below the cutoff.
  if (q2Trial <= q2Cut) return 0;

  // Find lower value from evolution window.
  int iEvol = evolutionWindows.size() - 1;
  while (q2Start <= evolutionWindows[iEvol]) iEvol--;
  double q2Low = evolutionWindows[iEvol];

  // Compute weights.
  vector<double> weightVec;
  double totWeight(0), maxWeight(0);
  for (int i = 0; i < (int)eleVec.size(); i++) {
    double Iz = q2Low > eleVec[i].m2Ant ? 0 : 1. - q2Low/eleVec[i].m2Ant;
    double w = totIdWeight*eleVec[i].ariWeight*Iz*eleVec[i].getKallen();
    weightVec.push_back(w);
    totWeight += w;
    if (w > maxWeight) maxWeight = w;
  }

  // If no antennae are active, don't generate new scale.
  if (totWeight < TINY) q2Trial = 0;

  // Generate scale and do alpha veto.
  else {
    while (true) {
      double alphaMax = al.alphaEM(q2Trial);
      q2Trial *= pow(rndmPtr->flat(), M_PI/totWeight/alphaMax);
      double alphaNew = al.alphaEM(q2Trial);
      if (rndmPtr->flat() < alphaNew/alphaMax) break;
    }
  }

  // Check if evolution window was crossed.
  if (q2Trial < q2Low) {
    if (iEvol == 0) return 0;
    return generateTrialScale(event, q2Low);
  }

  // Select antenna.
  while (true) {
    int i = rndmPtr->flat()*weightVec.size();
    if (rndmPtr->flat() < weightVec[i]/maxWeight) {
      eleTrial = &eleVec[i];
      break;
    }
  }

  // Select splitting ID.
  while (true) {
    int idIndex = rndmPtr->flat()*ids.size();
    idTrial = ids[idIndex];
    if (rndmPtr->flat() < idWeights[idIndex]/maxIdWeight) break;
  }

  // Generate value of zeta and phi.
  zTrial = (1. - q2Low/eleTrial->m2Ant)*rndmPtr->flat();
  phiTrial = rndmPtr->flat()*2*M_PI;
  hasTrial = true;
  return q2Trial;

}

//--------------------------------------------------------------------------

// Check the veto.

bool QEDsplitSystem::checkVeto(Event &event) {

  // Mark trial as used
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");
  hasTrial = false;

  // Set up some shorthands.
  int iPhot = eleTrial->iPhot;
  int iSpec = eleTrial->iSpec;
  double m2Ant = eleTrial->m2Ant;

  // New momenta.
  vector<Vec4> pOld, pNew;
  pOld.push_back(event[iPhot].p());
  pOld.push_back(event[iSpec].p());

  // ij is the new pair, k is the spectator.
  double mFerm = particleDataPtr->m0(idTrial);
  double mSpec = sqrt(eleTrial->m2Spec);
  double sIK = m2Ant - 2*pow2(mFerm) - pow2(mSpec);
  double sij = q2Trial - 2*pow2(mFerm);
  double sjk = zTrial*m2Ant;
  double sik = m2Ant - sij - sjk - 2*pow2(mFerm) - pow2(mSpec);

  // Check phase space.
  if (sik < 0) return false;
  if (sij*sjk*sik - pow2(sij)*pow2(mSpec)
    - (pow2(sjk) + pow2(sik))*pow2(mFerm) < 0) return false;

  // Kernel veto.
  double pVeto = ( (pow2(sik) + pow2(sjk))/m2Ant + 2.*pow2(mFerm)/q2Trial)/2.;
  if (rndmPtr->flat() > pVeto) return false;
  vector<double> invariants;
  invariants.push_back(sIK);
  invariants.push_back(sij);
  invariants.push_back(sjk);
  vector<double> masses;
  masses.push_back(mFerm);
  masses.push_back(mFerm);
  masses.push_back(mSpec);

  // Kinematics.
  if (!vinComPtr->map2to3FF(pNew, pOld, 3, invariants, phiTrial, masses))
    return false;

  // Set up the new fermions. Stochastic colour tag.
  int colTag = idTrial < 10 ? 10*(event.nextColTag()/10 + 1) + 1
                         + rndmPtr->flat()*10 : 0;
  Particle partFermNew(idTrial, 51, iPhot, 0, 0, 0, colTag, 0, pNew[0], mFerm);
  Particle partAFermNew(-idTrial,51,iPhot, 0, 0, 0, 0, colTag, pNew[1], mFerm);
  Particle partSpecNew = event[iSpec];
  partSpecNew.mothers(iSpec, iSpec);
  partSpecNew.p(pNew[2]);
  partSpecNew.statusCode(52);

  // Change the event - add new particles.
  int iFermNew  = event.append(partFermNew);
  int iAFermNew = event.append(partAFermNew);
  int iSpecNew  = event.append(partSpecNew);

  // Adjust old ones.
  event[iPhot].statusNeg();
  event[iPhot].daughters(iFermNew, iAFermNew);
  event[iSpec].statusNeg();
  event[iSpec].daughters(iSpecNew, 0);

  // Update parton system.
  partonSystemsPtr->replace(iSys, iPhot, iFermNew);
  partonSystemsPtr->addOut(iSys, iAFermNew);
  partonSystemsPtr->replace(iSys, iSpec, iSpecNew);
  event.restorePtrs();
  if (verbose >= debug) printOut(__METHOD_NAME__, "end --------------");
  return true;

}

//--------------------------------------------------------------------------

// Print the system.

void QEDsplitSystem::print() {
  cout << "Splitting" << endl;
  for (int i = 0; i < (int)eleVec.size(); i++)
    cout << "(" << eleVec[i].iPhot << " " << eleVec[i].iSpec << ") "
         << "s = " << eleVec[i].m2Ant << " ariFac = " << eleVec[i].ariWeight
         << endl;
}

//==========================================================================

// Class for a QED conversion system.

//--------------------------------------------------------------------------

// Initialize the pointers.

void QEDconvSystem::initPtr(Info* infoPtrIn, VinciaCommon* vinComPtrIn) {
  infoPtr       = infoPtrIn;
  particleDataPtr  = infoPtr->particleDataPtr;
  partonSystemsPtr = infoPtr->partonSystemsPtr;
  rndmPtr          = infoPtr->rndmPtr;
  settingsPtr      = infoPtr->settingsPtr;
  vinComPtr        = vinComPtrIn;
  isInitPtr = true;
}

//--------------------------------------------------------------------------

// Initialize the system.
void QEDconvSystem::init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  int verboseIn) {

  // Verbosity setting.
  if (!isInitPtr) printOut(__METHOD_NAME__, "initPtr not called");
  verbose = verboseIn;

  // Settings, number of incoming flavours to allow conversions to.
  // Could be extended to allow top quarks in future; for now up to b.
  nQuark = 5;
  if (!settingsPtr->flag("Vincia:convertGammaToQuark")) nQuark = 0;

  // Set trial pdf ratios.
  Rhat[1]  = 77;
  Rhat[2]  = 140;
  Rhat[3]  = 30;
  Rhat[4]  = 22;
  Rhat[5]  = 15;
  Rhat[-1] = 63;
  Rhat[-2] = 65;
  Rhat[-3] = 30;
  Rhat[-4] = 30;
  Rhat[-5] = 16;

  // Constants.
  TINYPDF = pow(10, -10);

  // Beam pointers.
  beamAPtr = beamAPtrIn;
  beamBPtr = beamBPtrIn;
  isInit = true;

}

//--------------------------------------------------------------------------

// Prepare for backwards-evolution of photons.

void QEDconvSystem::prepare(int iSysIn, Event &event, double q2CutIn,
  bool isBelowHadIn, vector<double> evolutionWindowsIn, AlphaEM alIn) {

  if (!isInit) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__+": Not initialised.");
    return;
  }
  if (verbose >= louddebug) printOut(__METHOD_NAME__, "begin --------------");

  // Input.
  iSys = iSysIn;
  shh = infoPtr->s();
  isBelowHad = isBelowHadIn;
  q2Cut = q2CutIn;
  evolutionWindows = evolutionWindowsIn;
  al = alIn;

  // Set up weights for conversion flavours.
  ids.clear();
  idWeights.clear();
  totIdWeight = 0;
  maxIdWeight = 0;

  // If switched off, do nothing.
  if (nQuark == 0) return;

  // Only do conversions to quarks if above hadronisation scale.
  if (!isBelowHad)
    for (int i = 1; i <= nQuark; i++) {
      ids.push_back(i);
      ids.push_back(-i);
      idWeights.push_back((i%2==0 ? 4./9. : 1./9.)*Rhat[i]);
      idWeights.push_back((i%2==0 ? 4./9. : 1./9.)*Rhat[-i]);
    }
  // Total weights.
  for (int i = 0; i < (int)idWeights.size(); i++) {
    totIdWeight += idWeights[i];
    if (idWeights[i] > maxIdWeight) maxIdWeight = idWeights[i];
  }

  // Build system.
  buildSystem(event);
  if (verbose >= louddebug) printOut(__METHOD_NAME__, "end --------------");

}

//--------------------------------------------------------------------------

// Build the system.

void QEDconvSystem::buildSystem(Event &event) {

  // Get rid of saved trial.
  hasTrial = false;

  // Get initial states.
  iA = partonSystemsPtr->getInA(iSys);
  iB = partonSystemsPtr->getInB(iSys);
  isAPhot = event[iA].id() == 22;
  isBPhot = event[iB].id() == 22;
  s = (event[iA].p() + event[iB].p()).m2Calc();

}

//--------------------------------------------------------------------------

// Generate a trial scale.

double QEDconvSystem::generateTrialScale(Event &event, double q2Start) {

  // Return if saved trial.
  if (hasTrial) return q2Trial;

  // Check if there are any photons.
  if (!isAPhot && !isBPhot) {return 0;}
  double totWeight = 1.;

  // Select a photon.
  if       (isAPhot && !isBPhot)  {iPhotTrial = iA; iSpecTrial = iB;}
  else  if (isBPhot && !isAPhot)  {iPhotTrial = iB; iSpecTrial = iA;}
  else {
    if (rndmPtr->flat() < 0.5)    {iPhotTrial = iA; iSpecTrial = iB;}
    else                          {iPhotTrial = iB; iSpecTrial = iA;}
    // Two photon antennae -> twice the weight.
    totWeight *= 2.;
  }

  // Starting scale.
  q2Trial = q2Start;

  // Check if qTrial is below the cutoff.
  if (q2Trial <= q2Cut) return 0;

  // Find lower value from evolution window.
  int iEvol = evolutionWindows.size() - 1;
  while(q2Start <= evolutionWindows[iEvol]) {iEvol--;}
  double q2Low = evolutionWindows[iEvol];

  // Iz integral.
  double zPlus = shh/s;
  double zMin = 1 + q2Low/s;
  if (zPlus < zMin) {return 0;}
  double Iz = log(zPlus/zMin);
  totWeight *= totIdWeight*Iz;

  // If no antennae are active, don't generate new scale.
  if (totWeight < TINY) return 0;

  // Generate scale and do alpha veto.
  else
    while(true) {
      double alphaMax = al.alphaEM(q2Trial);
      q2Trial *= pow(rndmPtr->flat(), M_PI/totWeight/alphaMax);
      double alphaNew = al.alphaEM(q2Trial);
      if (rndmPtr->flat() < alphaNew/alphaMax) break;
    }

  // Check if evolution window was crossed.
  if (q2Trial < q2Low) {
    if (iEvol==0) return 0;
    return generateTrialScale(event, q2Low);
  }

  // Select conversion ID.
  while( true) {
    int idIndex = rndmPtr->flat()*ids.size();
    idTrial = ids[idIndex];
    if (rndmPtr->flat() < idWeights[idIndex]/maxIdWeight) break;
  }
  zTrial = zMin*pow(zPlus/zMin, rndmPtr->flat());
  phiTrial = rndmPtr->flat()*2*M_PI;
  hasTrial = true;
  return q2Trial;

}

//--------------------------------------------------------------------------

// Check the veto.

bool QEDconvSystem::checkVeto(Event &event) {

  // Mark trial as used
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");
  hasTrial = false;

  // Conversion mass.
  double mf2 = pow2(particleDataPtr->m0(idTrial));

  // Spectator ID.
  int idSpec = event[iSpecTrial].id();

  // New momenta.
  vector<Vec4> pOld, pNew;
  pOld.push_back(event[iPhotTrial].p());
  pOld.push_back(event[iSpecTrial].p());

  // Note that we treat the initial state as massless and the final
  // state as massive. q2 is defined as saj - 2*mf2, but otherwise we
  // adhere to the proper kinematics.
  double saj = q2Trial + 2*mf2;
  double sbj = (zTrial-1)*s - saj + mf2;
  double sab = s + saj + sbj - mf2;

  // Check kinematic boundaries.
  if (sbj < 0) return false;
  if (sab*saj*sbj - mf2*sab*sab < 0) return false;

  // Check if photon is in beam a or b.
  bool isPhotA = (iPhotTrial == iA) ? true : false;

  // Global recoil momenta.
  vector<Vec4> pRec;
  vector<int>  iRec;
  int sysSize = partonSystemsPtr->sizeAll(iSys);
  for (int i=0; i<sysSize; i++) {
    int iEv = partonSystemsPtr->getAll(iSys, i);
    if (iEv < 0 || !event[iEv].isFinal()) continue;
    pRec.push_back(event[iEv].p());
    iRec.push_back(iEv);
  }

  // Kinematics.
  if (!vinComPtr->map2to3II(pNew, pRec, pOld, s, saj, sbj, sab,
      phiTrial, mf2)) return false;

  // Check if new energies don't exceed hadronic maxima.
  double eaUsed = 0, ebUsed = 0;
  int nSys = partonSystemsPtr->sizeSys();
  for (int i=0; i<nSys; i++) {
    eaUsed += event[partonSystemsPtr->getInA(i)].e();
    ebUsed += event[partonSystemsPtr->getInB(i)].e();
  }
  if (isPhotA) {
    if ((eaUsed - pOld[0].e() + pNew[0].e()) > 0.98*sqrt(shh)/2.) return false;
    if ((ebUsed - pOld[1].e() + pNew[2].e()) > 0.98*sqrt(shh)/2.) return false;
  } else {
    if ((ebUsed - pOld[0].e() + pNew[0].e()) > 0.98*sqrt(shh)/2.) return false;
    if ((eaUsed - pOld[1].e() + pNew[2].e()) > 0.98*sqrt(shh)/2.) return false;
  }

  // Kernel veto probability.
  double pVeto = 0.5*(1. + pow2(sbj)/pow2(sab)
    - 2.*mf2*pow2(s)/pow2(sab)/(saj-2*mf2));

  // Compute pdf ratios.
  double Rpdf = 1.;
  double xPhotOld = pOld[0].e()/(sqrt(shh)/2.);
  double xPhotNew = pNew[0].e()/(sqrt(shh)/2.);
  double xSpecOld = pOld[1].e()/(sqrt(shh)/2.);
  double xSpecNew = pNew[2].e()/(sqrt(shh)/2.);

  if (isPhotA) {
    // Photon pdf.
    double newPDFPhot = beamAPtr->xfISR(iSys, idTrial, xPhotNew, q2Trial);
    double oldPDFPhot = beamAPtr->xfISR(iSys, 22,      xPhotOld, q2Trial);
    if (abs(newPDFPhot) < TINYPDF) newPDFPhot = TINYPDF;
    if (abs(oldPDFPhot) < TINYPDF) oldPDFPhot = TINYPDF;
    Rpdf *= newPDFPhot/oldPDFPhot;

    // Spectator pdf.
    double newPDFSpec = beamBPtr->xfISR(iSys, idSpec, xSpecNew, q2Trial);
    double oldPDFSpec = beamBPtr->xfISR(iSys, idSpec, xSpecOld, q2Trial);
    if (abs(newPDFSpec) < TINYPDF) newPDFSpec = TINYPDF;
    if (abs(oldPDFSpec) < TINYPDF) oldPDFSpec = TINYPDF;
    Rpdf *= newPDFSpec/oldPDFSpec;
  } else {
    // Photon pdf.
    double newPDFPhot = beamBPtr->xfISR(iSys, idTrial, xPhotNew, q2Trial);
    double oldPDFPhot = beamBPtr->xfISR(iSys, 22,      xPhotOld, q2Trial);
    if (abs(newPDFPhot) < TINYPDF) newPDFPhot = TINYPDF;
    if (abs(oldPDFPhot) < TINYPDF) oldPDFPhot = TINYPDF;
    Rpdf *= newPDFPhot/oldPDFPhot;

    // Spectator pdf.
    double newPDFSpec = beamAPtr->xfISR(iSys, idSpec, xSpecNew, q2Trial);
    double oldPDFSpec = beamAPtr->xfISR(iSys, idSpec, xSpecOld, q2Trial);
    if (abs(newPDFSpec) < TINYPDF) newPDFSpec = TINYPDF;
    if (abs(oldPDFSpec) < TINYPDF) oldPDFSpec = TINYPDF;
    Rpdf *= newPDFSpec/oldPDFSpec;
  }

  if (Rpdf > Rhat[idTrial]) {
    stringstream ss;
    ss << "Incorrect pdf overestimate " << "id = " << idTrial << "ratio = "
       << Rpdf/Rhat[idTrial];
    printOut(__METHOD_NAME__, ss.str());
  }

  // Pdf ratio veto probability.
  pVeto *= (Rpdf/Rhat[idTrial]);

  // Do veto.
  if (rndmPtr->flat() > pVeto) return false;

  // Passed veto, set up new particles.
  Particle partSpecNew = event[iSpecTrial];
  partSpecNew.p(pNew[2]);
  partSpecNew.statusCode(41);

  // Stochastic colour tag.
  int colTag = 10*(event.nextColTag()/10 + 1) + 1 + rndmPtr->flat()*10;
  Particle partBeamNew  (idTrial, -41, 0, 0, 0, 0, idTrial > 0 ?
    colTag : 0, idTrial > 0 ? 0 : colTag, pNew[0]);
  Particle partFinalNew (idTrial,  43, 0, 0, 0, 0, idTrial > 0 ?
    colTag : 0, idTrial > 0 ? 0 : colTag, pNew[1], sqrt(mf2));
  int iBeamNew = event.append(partBeamNew);
  int iFinalNew = event.append(partFinalNew);
  int iSpecNew = event.append(partSpecNew);
  event[iPhotTrial].statusNeg();
  event[iSpecTrial].statusNeg();
  event[iBeamNew].mothers(event[iPhotTrial].mother1(),
    event[iPhotTrial].mother2());
  event[iFinalNew].mothers(iBeamNew, 0);
  event[iSpecNew].mothers(event[iSpecTrial].mother1(),
    event[iSpecTrial].mother2());
  event[iPhotTrial].mothers(iBeamNew, 0);
  event[iSpecTrial].mothers(iSpecNew, 0);
  event[iBeamNew].daughters(iFinalNew, iPhotTrial);
  event[iFinalNew].daughters(0,0);
  event[iSpecNew].daughters(iSpecTrial);

  // Change daughters of beams for hard process.
  if (iSys == 0) {
    bool foundPhot = false;
    bool foundSpec = false;
    for (int i = 0; i < (int)event.size(); i++) {
      if (!foundPhot)
        if (event[i].daughter1() == iPhotTrial) {
          event[i].daughters(iBeamNew, 0);
          foundPhot = true;
        }
      if (!foundSpec)
        if (event[i].daughter1() == iSpecTrial) {
          event[i].daughters(iSpecNew, 0);
          foundSpec = true;
        }
      if (foundPhot && foundSpec) break;
    }
  }

  // Update event for global recoil.
  vector<pair<int,int> > iRecNew;
  iRecNew.clear();
  iRecNew.resize(0);
  for (int j = 0; j < event.size(); j++)
    if (event[j].isFinal())
      for (int k = 0; k < (int)iRec.size(); k++)
        if (iRec[k] == j) {
          // Copy the recoiler and change the momentum.
          int inew = event.copy(j,44);
          event[inew].p(pRec[k]);
          iRecNew.push_back(make_pair(iRec[k], inew));
        }

  // Update parton system.
  partonSystemsPtr->replace(iSys, iPhotTrial, iBeamNew);
  partonSystemsPtr->addOut(iSys, iFinalNew);
  partonSystemsPtr->replace(iSys, iSpecTrial, iSpecNew);

  // Recoilers and initial particles.
  for (int k=0; k<(int)iRecNew.size(); k++) {
    partonSystemsPtr->replace(iSys, iRecNew[k].first, iRecNew[k].second);
  }
  if (isPhotA) {
    partonSystemsPtr->setInA(iSys, iBeamNew);
    partonSystemsPtr->setInB(iSys, iSpecNew);
  } else {
    partonSystemsPtr->setInA(iSys, iSpecNew);
    partonSystemsPtr->setInB(iSys, iBeamNew);
  }
  double shat = (event[partonSystemsPtr->getInA(iSys)].p() +
    event[partonSystemsPtr->getInB(iSys)].p()).m2Calc();
  partonSystemsPtr->setSHat(iSys, shat);

  // Update beams.
  BeamParticle& beam1 = *beamAPtr;
  BeamParticle& beam2 = *beamBPtr;
  if (isPhotA) {
    beam1[iSys].update(iBeamNew, event[iBeamNew].id(),
      event[iBeamNew].e()/beam1.e());
    beam2[iSys].update(iSpecNew, event[iSpecNew].id(),
      event[iSpecNew].e()/beam2.e());
  } else {
    beam1[iSys].update(iSpecNew, event[iSpecNew].id(),
      event[iSpecNew].e()/beam1.e());
    beam2[iSys].update(iBeamNew, event[iBeamNew].id(),
      event[iBeamNew].e()/beam2.e());
  }
  if (verbose >= debug) printOut(__METHOD_NAME__, "end --------------");
  return true;

}

//--------------------------------------------------------------------------

// Print.

void QEDconvSystem::print() {
  cout << "Conversion" << endl;
  cout << "s = " << s << endl;
}

//==========================================================================

// Class for performing QED showers.

//--------------------------------------------------------------------------

// Initialize the pointers.

void QEDShower::initPtr(Info* infoPtrIn, VinciaCommon* vinComPtrIn) {
  infoPtr       = infoPtrIn;
  particleDataPtr  = infoPtr->particleDataPtr;
  partonSystemsPtr = infoPtr->partonSystemsPtr;
  rndmPtr          = infoPtr->rndmPtr;
  settingsPtr      = infoPtr->settingsPtr;
  vinComPtr        = vinComPtrIn;
  isInitPtr = true;
}

//--------------------------------------------------------------------------

// Initialize settings for current run.

void QEDShower::init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn) {

  // Initialize alphaEM
  verbose = settingsPtr->mode("Vincia:verbose");
  double alpEM0Vincia = settingsPtr->parm("Vincia:alphaEM0");
  double alpEMmzVincia = settingsPtr->parm("Vincia:alphaEMmz");
  double alpEM0Pythia = settingsPtr->parm("StandardModel:alphaEM0");
  double alpEMmzPythia = settingsPtr->parm("StandardModel:alphaEMmZ");
  int alphaEMorder = settingsPtr->mode("Vincia:alphaEMorder");

  // Change Pythia settings, initialize, then change them back.
  settingsPtr->parm("StandardModel:alphaEM0", alpEM0Vincia);
  settingsPtr->parm("StandardModel:alphaEMmZ", alpEMmzVincia);
  al.init(alphaEMorder, settingsPtr);
  settingsPtr->parm("StandardModel:alphaEM0", alpEM0Pythia);
  settingsPtr->parm("StandardModel:alphaEMmz", alpEMmzPythia);

  // Get settings.
  doQED          = settingsPtr->flag("Vincia:doQED");
  doEmission     = doQED;
  nGammaToLepton = settingsPtr->mode("Vincia:nGammaToLepton");
  nGammaToQuark  = settingsPtr->mode("Vincia:nGammaToQuark") >= 1;
  doConvertGamma = settingsPtr->flag("Vincia:convertGammaToQuark");
  doConvertQuark = settingsPtr->flag("Vincia:convertQuarkToGamma");

  // QED cutoff for coloured particles/hadronisation scale.
  q2minColouredSav = pow2(settingsPtr->parm("Vincia:QminChgQ"));
  q2minSav         = pow2(settingsPtr->parm("Vincia:QminChgL"));

  // Set beam pointers.
  beamAPtr = beamAPtrIn;
  beamBPtr = beamBPtrIn;
  isInitSav = true;

}

//--------------------------------------------------------------------------

// Prepare to shower a system.

void QEDShower::prepare(int iSysIn, Event &event, bool isBelowHad) {

  // Verbose output
  if (!doQED) return;
  if (verbose >= debug) {
    printOut(__METHOD_NAME__, "begin --------------");
    if (verbose >= superdebug) event.list();
  }

  // Above or below hadronisation scale.
  double q2cut = (isBelowHad) ? q2minSav : q2minColouredSav;

  // Initialize windows for the hard systen and
  // the final after beam remnants system.
  if (iSysIn == 0 || iSysIn == -1) {
    // The cutoff scale is the lowest window boundary, then step up to
    // q2Max successively by factor winStep.
    double q2BiggestEver  = infoPtr->s();
    double q2Window       = q2cut;
    double winStep        = 100.0;
    evolutionWindows.clear();
    do {
      evolutionWindows.push_back(q2Window);
      q2Window *= winStep;
    } while(q2Window < q2BiggestEver);
  }

  // Find out if showering a resonance decay.
  bool isResDecay = false;
  if (iSysIn >= 0)
    if (partonSystemsPtr->hasInRes(iSysIn) ||
        (partonSystemsPtr->getInA(iSysIn)==0 &&
         partonSystemsPtr->getInB(iSysIn)==0) ) isResDecay = true;

  // If showering a resonance decay or below hadronization scale clear
  // out all previous systems.
  if ( iSysIn <= 0 || isResDecay || isBelowHad ) {
    iSystems.clear();
    emitSystems.clear();
    splitSystems.clear();
    convSystems.clear();
  }

  // Special case: iSysIn = -1 implies below hadronisation scale.
  // Collect all final-state particles into one new system.  Note:
  // this system will have sum(charge) != 0 if the sum of the beam
  // charges is nonzero.
  if (iSysIn == -1) {
    iSysIn = partonSystemsPtr->addSys();
    // Loop over all particles in event rather than all parton systems
    // since beam remnant partons not part of any partonSystem.
    for (int i = 1; i < event.size(); ++i) {
      if (!event[i].isFinal()) continue;
      partonSystemsPtr->addOut(iSysIn, i);
    }
  }

  // Add new systems.
  iSystems.push_back(iSysIn);
  emitSystems.push_back(QEDemitSystem());
  splitSystems.push_back(QEDsplitSystem());
  convSystems.push_back(QEDconvSystem());

  // Initialze pointers.
  emitSystems.back().initPtr(infoPtr, vinComPtr);
  splitSystems.back().initPtr(infoPtr, vinComPtr);
  convSystems.back().initPtr(infoPtr, vinComPtr);

  // Initialize. Note, these calls read from settings database, should
  // be rewritten to require settings to be read only once.
  emitSystems.back().init(beamAPtr, beamBPtr,verbose);
  splitSystems.back().init(beamAPtr, beamBPtr,verbose);
  convSystems.back().init(beamAPtr, beamBPtr,verbose);

  // Prepare and build QED systems.
  emitSystems.back().prepare(iSysIn, event, q2cut, isBelowHad,
    evolutionWindows, al);
  splitSystems.back().prepare(iSysIn, event, q2cut, isBelowHad,
    evolutionWindows, al);
  convSystems.back().prepare(iSysIn, event, q2cut, isBelowHad,
    evolutionWindows, al);

  // Verbose output.
  if(verbose >= debug) printOut(__METHOD_NAME__, "end --------------");

}

//--------------------------------------------------------------------------

// Update QED shower system(s) each time something has changed in event.

void QEDShower::update(Event &event, int iSys) {

  // Find index of the system.
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");
  for (int i = 0; i < (int)iSystems.size(); i++) {
    if (iSystems[i] == iSys) {
      // Rebuild systems.
      emitSystems[i].buildSystem(event);
      splitSystems[i].buildSystem(event);
      convSystems[i].buildSystem(event);
      break;
    }
  }
  if (verbose >= debug) {
    if (verbose >= superdebug) event.list();
    printOut(__METHOD_NAME__, "end --------------");
  }
}

//--------------------------------------------------------------------------

// Generate a trial scale.

double QEDShower::generateTrialScale(Event &event, double q2Start) {

  // Get a trial from every system.
  q2Trial = 0;
  isTrialEmit = false;
  isTrialSplit = false;
  isTrialConv = false;
  if (!doQED) return 0.0;
  if (verbose >= louddebug) {
    printOut(__METHOD_NAME__, "begin --------------");
    if (verbose >= superdebug)
      cout << " QEDShower::generateTrialScale() q2Start = " << q2Start
           << " doEmit = " << bool2str(doEmission)
           << " nSplitGamToLep = " << num2str(nGammaToLepton)
           << " nSplitGamToQuark = " << num2str(nGammaToQuark)
           << " doConv = " << bool2str(doConvertGamma) << endl;
  }

  // Emissions.
  if (doEmission) {
    for (int i = 0; i < (int)emitSystems.size(); i++) {
      double q2TrialEmitNew =
        emitSystems[i].generateTrialScale(event, q2Start);
      if (q2TrialEmitNew > q2Trial) {
        q2Trial = q2TrialEmitNew;
        iSysTrial = iSystems[i];
        iSysIndexTrial = i;
        isTrialEmit = true;
        isTrialSplit = false;
        isTrialConv = false;
      }
    }
  }

  // Splittings.
  if (nGammaToLepton + nGammaToQuark > 0) {
    for (int i = 0; i < (int)splitSystems.size(); i++) {
      double q2TrialSplitNew =
        splitSystems[i].generateTrialScale(event, q2Start);
      if (q2TrialSplitNew > q2Trial) {
        q2Trial = q2TrialSplitNew;
        iSysTrial = iSystems[i];
        iSysIndexTrial = i;
        isTrialEmit = false;
        isTrialSplit = true;
        isTrialConv = false;
      }
    }
  }

  // Conversions.
  if (doConvertGamma) {
    for (int i = 0; i < (int)convSystems.size(); i++) {
      double q2TrialConvNew =
        convSystems[i].generateTrialScale(event, q2Start);
      if (q2TrialConvNew > q2Trial) {
        q2Trial = q2TrialConvNew;
        iSysTrial = iSystems[i];
        iSysIndexTrial = i;
        isTrialEmit = false;
        isTrialSplit = false;
        isTrialConv = true;
      }
    }
  }
  if (verbose >= louddebug) printOut(__METHOD_NAME__, "end --------------");
  return q2Trial;

}

//--------------------------------------------------------------------------

// Check the veto.

bool QEDShower::checkVeto(Event &event) {
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");
  bool doVeto = false;
  if (isTrialEmit)  doVeto = emitSystems[iSysIndexTrial].checkVeto(event);
  if (isTrialSplit) doVeto = splitSystems[iSysIndexTrial].checkVeto(event);
  if (isTrialConv)  doVeto = convSystems[iSysIndexTrial].checkVeto(event);
  if (verbose >= debug) printOut(__METHOD_NAME__, "end --------------");
  return doVeto;
}

//==========================================================================

} // end namespace Pythia8
