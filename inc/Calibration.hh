#ifndef __CALIBRATION_HH
#define __CALIBRATION_HH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <map>

#include "TSystem.h"
#include "TEnv.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"

#include "Settings.hh"
#ifdef SIMULATION
#include "Gretina.hh"
#include "Miniball.hh"
#include "ZeroDeg.hh"
#include "MINOS.hh"
#include "Tracking.hh"
#else
#include "HiCARI.hh"
#endif

using namespace std;

class Calibration {
public:
  Calibration();
  
//! Constructs the Calibration object based on the given settings file.
/*!
  Constructs the Calibration objects, reading in the various settings.
  @param setting The Settings object from which to read calibration data.
  @param event The event number from which to start counting, used for time-dependent corrections.
 */
  Calibration(Settings*, int event =0);
  ~Calibration();
#ifdef SIMULATION
  //! Read the Miniball coordinates, average first interactions from file
  void ReadMBPositions(const char* filename);
#else
  //! Read the HiCARI coordinates, average first interactions from file
  void ReadHiCARIPositions(const char* filename);
  //! Read the HiCARI calibration parameters
  void ReadHiCARICalibration(const char* filename);
#endif

#ifdef SIMULATION
  //! Construct all calibrated objects.
  void BuildAllCalc(Gretina* inGret, GretinaCalc* outGret, Miniball* inMB, MiniballCalc* outMB, ZeroDeg *zerodeg, MINOS *minos);
  
  //! Build the MiniballCalc object, given a raw Miniball object.
  void BuildMiniballCalc(Miniball* in, MiniballCalc* out);
  //! Build the GretinaCalc object, given a raw Gretina object.
  void BuildGretinaCalc(Gretina* in, GretinaCalc* out);
  //! Build the ZeroDeg object.
  void BuildZeroDeg(ZeroDeg *zerodeg);
  //! Build the MINOS object.
  void BuildMINOS(MINOS *minos);

  void AddBackGretinaCrystal(GretinaCalc* gr);
  void AddBackGretinaCluster(GretinaCalc* gr);
  void AddBackGretinaEverything(GretinaCalc* gr);
  void AddBackMiniballCluster(MiniballCalc* gr);
  void AddBackMiniballEverything(MiniballCalc* gr);
  void ClusterGretina(GretinaCalc* gr, Gretina *in);
  vector<HitCalc*> ExtractAllHits(Gretina* in);
  void AllGretinaHits(GretinaCalc* gr, Gretina *in);

  //! Construct tracked gamma events
  void GammaTrack(GretinaCalc* gr, GretinaEvent* gt);
#else
  //! Build the HiCARICalc object, given a raw HiCARI object.
  void BuildHiCARICalc(HiCARI* in, HiCARICalc* out);
  void AddBackHiCARICluster(HiCARICalc* gr);
  void AddBackHiCARIEverything(HiCARICalc* gr);
#endif

  void PrintCtrs();

  
private:
  void ResetCtrs();

  Settings* fSett;
  int fverbose;

  int fevent;

#ifdef SIMULATION
  //! averaged miniball first interaction positions
  TVector3 fMBpositions[MBCLUST+CLOVERS][CRYST][SEGS];//cluster, crystal, segment
#else
  TRandom* fRand;
  //! HiCARI positions
  TVector3 fHiCARIpositions[12][4][40];//cluster, crystal, segment
  double fCoreGain[12][4];
  double fCoreOffs[12][4];
  double fSegGain[12][4][40];
  double fSegOffs[12][4][40];
#endif

  int fAddBackType;
#ifdef SIMULATION
  Long64_t fgretactr;
  Long64_t fminiballctr;
  Long64_t fzerodegctr;
  Long64_t fminosctr;
  TF1* fMINOSZett;
  
  Tracking* ftracking;
#else
  Long64_t fHiCARIctr;
#endif
};

#endif
