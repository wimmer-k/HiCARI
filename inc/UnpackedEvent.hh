#ifndef UNPACKED_EVENT_HH__
#define UNPACKED_EVENT_HH__

#include <iostream>
#include <iomanip>

#include "TObject.h"
#include "TTree.h"
#include "TList.h"
#include "TEnv.h"
#include "TRandom.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "Calibration.hh"
#include "Trace.hh"
#include "Gretina.hh"
#include "Miniball.hh"
#include "GammaSim.hh"
#include "ZeroDeg.hh"
#ifdef USEMINOS
#include "MINOS.hh"
#endif
#ifdef USELISA
#include "LISA.hh"
#endif
#include "Settings.hh"
#include "RawHistograms.hh"
#include "CalHistograms.hh"

struct simresolution{
  // resolution = A*sqrt(1 + B*Energy)
  int Detector;
  int Crystal;
  double A;
  double B;
  double C;
};

struct simthreshold{
  // threshold = 0.5*( 1 + tanh( (Energy - E) / dE ) )
  int Detector;
  int Crystal;
  double E;
  double dE;
};

class UnpackedEvent {
public:
  UnpackedEvent(){};
  UnpackedEvent(Settings* settings = NULL);
  ~UnpackedEvent(){
    delete fGretina;
    delete fGretinaCalc;
    delete fZeroDeg;
    delete fMode3Event;
#ifdef USEMINOS
    delete fMINOS;
#endif    
#ifdef USELISA
    delete fLISA;
#endif    
    delete frhist;
    delete fchist;
    delete fMiniball;
    delete fMiniballCalc;
  }
  
  void SetVL(int vl){
    fvl = vl;
  }
  void SetCalibration(Calibration* cal){
    fcal = cal;
  }
  void SetWrite(bool wtree, bool whist, bool wctree, bool wchist, bool wstree){
    fwtree = wtree;
    fwhist = whist;
    fwcaltree = wctree;
    fwcalhist = wchist;
    fwsimtree = wstree;
  }
  void SetRecalibrate(bool recalibrate){
    frecalibrate = recalibrate;
  }
  void SetWrite(bool wtree, bool whist){
    fwtree = wtree;
    fwhist = whist;
  }
  void Init();
  //! Read the resolutions to be applied to the simulated data
  void ReadSimResolution(const char* filename);
  //! Read the energy thresholds to be applied to the simulated data
  void ReadSimThresholds(const char* filename);
  //! Read the specified buffer and make the GretinaG4Sim event.
  int DecodeGammaG4Sim(G4SIM_EGS*, long long int);
  //! Read the specified buffer and make the ZeroDegPhysicsData event.
  int DecodeZeroDegPhysicsData(ZD_PHYSICSDATA*, long long int);
#ifdef USEMINOS
  //! Read the specified buffer and make the MINOSPhysicsData event.
  int DecodeMINOSPhysicsData(MINOS_DATA*, long long int);
#endif
#ifdef USELISA
  //! Read the specified buffer and make the LISAPhysicsData event.
  int DecodeLISAPhysicsData(LISA_DATA*, long long int);
#endif
  //! Passed a gretina Crystal, make a new crystal in the Gretina object.
  int DecodeGretina(Crystal* cryst, long long int gts);
  //! Passed a miniball Crystal, make a new crystal in the Miniball object.
  int DecodeMiniball(MBCrystal* cryst, long long int gts);
  //! Read the specified buffer and make the Mode3Event.
  int DecodeMode3(char* cBuf, int len, long long int ts);

  //! Write the last event to file.
  void WriteLastEvent();
  void PrintHit(struct crys_ips_abcd1234 inbuf);
  void PrintHit(struct crys_ips_abcd5678 inbuf);

  int NrOfEvents(){return fnentries;}
  int NrOfCalEvents(){return fncalentries;}
  int NrOfGRHits(){return fGRhits;}
  int NrOfMBHits(){return fMBhits;}
  int NrOfStrangeHits(){return fstrangehits;}
  //! The number of simulated events
  int NrOfSimEvents(){return fnsimentries;}

  bool SimResolution(Gretina* gr);
  bool SimThresholds(Gretina* gr);
  
  bool SimResolution(Miniball* mb);
  bool SimThresholds(Miniball* mb);
  
  TTree* GetTree(){return ftr;}
  TTree* GetCalTree(){return fcaltr;}
  TTree* GetTraceTree(){return ftr;}
  //! The tree containing the simulated data
  TTree* GetSimTree(){return fsimtr;}

protected:
  //! Create a single trace
  Trace DecodeTrace(unsigned short** wBuf_p, int length, long long int gts);
  //! Performs end of event actions.
  void CloseEvent();
  //! Clears memory of current event.
  void ClearEvent();

  TRandom* fRand;
  TTree *fsimtr;
  TTree *ftr;
  TTree *fcaltr;

  ZeroDeg*      fZeroDeg;
#ifdef USEMINOS
  MINOS*        fMINOS;
#endif
#ifdef USELISA
  LISA*        fLISA;
#endif
  Gretina*      fGretina;
  GretinaCalc*  fGretinaCalc;
  Miniball*     fMiniball;
  MiniballCalc* fMiniballCalc;
  GammaSim*     fGammaSim;
  Mode3Event *fMode3Event;

  bool fhasdata;
  long long int fcurrent_ts;
  long long int ffirst_ts;
  int fvl;
  int fnentries;
  int fncalentries;
  int fnsimentries;
  int fGRhits;
  int fMBhits;
  int fctr;
  bool fwtree;
  bool fwhist;
  bool fwcaltree;
  bool fwcalhist;
  bool fwsimtree;
  int fstrangehits;
  int frecalibrate;
  int fEventTimeDiff;

  Calibration* fcal;

  RawHistograms* frhist;
  CalHistograms* fchist;

  Settings* fSett;

  vector<simresolution> fSimResolutions;
  vector<simthreshold>  fSimThresholds;

};

#endif
