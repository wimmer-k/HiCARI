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
#include "Gretina.hh"
#include "HiCARI.hh"
#include "Trace.hh"

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
  //! Read the HiCARI coordinates, average first interactions from file
  void ReadHiCARIPositions(const char* filename);
  //! Read the HiCARI calibration parameters
  void ReadHiCARICalibration(const char* filename);
  //! Read the HiCARI time offsets
  void ReadHiCARITimeOffset(const char* filename);
  //! Read the matrix file for the position transformation for mode2 data
  void ReadMatrix(const char* filename);
  
  //! Build the GretinaCalc object, given a raw Gretina object.
  void BuildGretinaCalc(Gretina* in, GretinaCalc* out);
  TVector3 TransformCoordinates(int hole, int cry, TVector3 local);
  void CalibrateIPoints(Crystal* cry);
  void AddBackGretinaCrystal(GretinaCalc* gr);
  void AddBackGretinaCluster(GretinaCalc* gr);
  void AddBackGretinaEverything(GretinaCalc* gr);

  void ClusterGretina(GretinaCalc* gr, Gretina *in);
  vector<HitCalc*> ExtractAllHits(Gretina* in);
  void AllGretinaHits(GretinaCalc* gr, Gretina *in);

  //! Construct tracked gamma events
  void GammaTrack(GretinaCalc* gr, GretinaEvent* gt);

  //! Build the HiCARICalc object, given a raw HiCARI object.
  void BuildHiCARICalc(HiCARI* in, HiCARICalc* out);
  void AddBackHiCARICluster(HiCARICalc* gr);
  void AddBackHiCARIEverything(HiCARICalc* gr);

  void PrintCtrs();
  long long int GetBigRIPSCtr(){return fBigRIPSctr;}
  long long int GetBigRIPSHitCtr(){return fBigRIPSHitctr;}
  long long int GetHiCARICtr(){return fHiCARIctr;}
  long long int GetHiCARIHitCtr(){return fHiCARIHitctr;}
  long long int GetHiCARIHitABCtr(){return fHiCARIHitABctr;}
  long long int GetGretinaCtr(){return fGretinactr;}

  long long int GetGretinaHitCtr(){return fGretinaHitctr;}
  long long int GetGretinaHitABCtr(){return fGretinaHitABctr;}
  
private:
  void ResetCtrs();

  Settings* fSett;
  int fverbose;

  int fevent;

  TRandom* fRand;
  //! HiCARI positions
  TVector3 fHiCARIpositions[MAXDETPOS][MAXCRYSTALNO][MAXSEGS];//cluster, crystal, segment
  double fCoreGain[MAXDETPOS][MAXCRYSTALNO];
  double fCoreOffs[MAXDETPOS][MAXCRYSTALNO];
  double fSegGain[MAXDETPOS][MAXCRYSTALNO][MAXSEGS];
  double fSegOffs[MAXDETPOS][MAXCRYSTALNO][MAXSEGS];
  double fCoreTimeOffset[MAXDETPOS][MAXCRYSTALNO];
  float fcrmat[MAXDETPOS][MAXCRYSTALNO][4][4];

  int fAddBackType;
  int fCoincTDiff;
  Long64_t fHiCARIctr;
  Long64_t fBigRIPSctr;
  Long64_t fHiCARIHitctr;
  Long64_t fHiCARIHitABctr;
  Long64_t fBigRIPSHitctr;
  Long64_t fGretinactr;
  Long64_t fGretinaHitctr;
  Long64_t fGretinaHitABctr;
};

#endif
