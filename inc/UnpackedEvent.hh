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
#ifdef SIMULATION
#include "Gretina.hh"
#include "Miniball.hh"
#include "GammaSim.hh"
#include "ZeroDeg.hh"
#ifdef USEMINOS
#include "MINOS.hh"
#endif
#else
#include "HiCARI.hh"
#endif
#include "Settings.hh"
#include "RawHistograms.hh"
#include "CalHistograms.hh"

#ifdef SIMULATION
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
#endif

class UnpackedEvent {
public:
  UnpackedEvent(){};
  UnpackedEvent(Settings* settings = NULL);
  ~UnpackedEvent(){
#ifdef SIMULATION
    delete fGretina;
    delete fGretinaCalc;
    delete fZeroDeg;
#ifdef USEMINOS
    delete fMINOS;
#endif    
    delete fMiniball;
    delete fMiniballCalc;
#else
    delete fMode3Event;
    delete fHiCARI;
    delete fHiCARICalc;
#endif
    delete frhist;
    delete fchist;
  }
  
  void SetVL(int vl){
    fvl = vl;
  }
  void SetCalibration(Calibration* cal){
    fcal = cal;
  }
  Calibration* GetCalibration(){
    return fcal;
  }
#ifdef SIMULATION
  void SetWrite(bool wtree, bool whist, bool wctree, bool wchist, bool wstree){
    fwtree = wtree;
    fwhist = whist;
    fwcaltree = wctree;
    fwcalhist = wchist;
    fwsimtree = wstree;
  }
#else
  void SetWrite(bool wtree, bool whist, bool wctree, bool wchist){
    fwtree = wtree;
    fwhist = whist;
    fwcaltree = wctree;
    fwcalhist = wchist;
  }
#endif
  void SetMakeMode2(bool makemode2){
    fmakemode2 = makemode2;
  }
  void SetWrite(bool wtree, bool whist){
    fwtree = wtree;
    fwhist = whist;
  }
  void Init();
#ifdef SIMULATION
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
  //! Passed a gretina Crystal, make a new crystal in the Gretina object.
  int DecodeGretina(Crystal* cryst, long long int gts);
  //! Passed a miniball Crystal, make a new crystal in the Miniball object.
  int DecodeMiniball(MBCrystal* cryst, long long int gts);
#endif

  //! Read the specified buffer and make the Mode3Event.
  int DecodeMode3(char* cBuf, int len, long long int ts);

  //! Write the last event to file.
  void WriteLastEvent();
#ifdef SIMULATION
  void PrintHit(struct crys_ips_abcd1234 inbuf);
  void PrintHit(struct crys_ips_abcd5678 inbuf);
#endif

  int NrOfEvents(){return fnentries;}
  int NrOfCalEvents(){return fncalentries;}
#ifdef SIMULATION
  int NrOfStrangeHits(){return fstrangehits;}
  int NrOfGRHits(){return fGRhits;}
  int NrOfMBHits(){return fMBhits;}
  //! The number of simulated events
  int NrOfSimEvents(){return fnsimentries;}

  bool SimResolution(Gretina* gr);
  bool SimThresholds(Gretina* gr);
  
  bool SimResolution(Miniball* mb);
  bool SimThresholds(Miniball* mb);
#endif  
  TTree* GetTree(){return ftr;}
  TTree* GetCalTree(){return fcaltr;}
#ifdef SIMULATION
  //! The tree containing the simulated data
  TTree* GetSimTree(){return fsimtr;}
#else
  TTree* GetTraceTree(){return ftr;}
  //! Make HiCARI objects from mode3 data
  void MakeMode2();
  void SetMode3(Mode3Event *m3e){fMode3Event = m3e;}
  HiCARI* GetHiCARI(){return fHiCARI;}
  HiCARICalc* GetHiCARICalc(){return fHiCARICalc;}
#endif

protected:
  //! Create a single trace
  Trace DecodeTrace(unsigned short** wBuf_p, int length, long long int gts);
  //! Performs end of event actions.
  void CloseEvent();
  //! Clears memory of current event.
  void ClearEvent();

  TRandom* fRand;
#ifdef SIMULATION
  TTree *fsimtr;
#endif
  TTree *ftr;
  TTree *fcaltr;

#ifdef SIMULATION
  ZeroDeg*      fZeroDeg;
#ifdef USEMINOS
  MINOS*        fMINOS;
#endif
  Gretina*      fGretina;
  GretinaCalc*  fGretinaCalc;
  Miniball*     fMiniball;
  MiniballCalc* fMiniballCalc;
  GammaSim*     fGammaSim;
#else
  Mode3Event *fMode3Event;
  HiCARI *fHiCARI;
  HiCARICalc *fHiCARICalc;
#endif

  bool fhasdata;
  long long int fcurrent_ts;
  long long int ffirst_ts;
  int fvl;
  int fnentries;
  int fncalentries;
#ifdef SIMULATION
  int fnsimentries;
  int fGRhits;
  int fMBhits;
#endif
  int fctr;
  bool fwtree;
  bool fwhist;
  bool fwcaltree;
  bool fwcalhist;
#ifdef SIMULATION
  bool fwsimtree;
  int fstrangehits;
  int frecalibrate;
#endif
  int fEventTimeDiff;

  Calibration* fcal;

  RawHistograms* frhist;
  CalHistograms* fchist;

  Settings* fSett;

#ifdef SIMULATION
  vector<simresolution> fSimResolutions;
  vector<simthreshold>  fSimThresholds;
#endif

  bool fmakemode2;
};

#endif
