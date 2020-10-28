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
#include "Gretina.hh" //now also for data as we take mode2 for P3 and Quad
#include "HiCARI.hh"
#include "Settings.hh"
#include "RawHistograms.hh"
#include "CalHistograms.hh"


class UnpackedEvent {
public:
  UnpackedEvent(){};
  UnpackedEvent(Settings* settings = NULL);
  ~UnpackedEvent(){
    delete fGretina;
    delete fGretinaCalc;
    delete fMode3Event;
    delete fHiCARI;
    delete fHiCARICalc;
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
  void SetWrite(bool wtree, bool whist, bool wctree, bool wchist){
    fwtree = wtree;
    fwhist = whist;
    fwcaltree = wctree;
    fwcalhist = wchist;
  }
  void SetMakeMode2(bool makemode2){
    fmakemode2 = makemode2;
  }
  void SetWrite(bool wtree, bool whist){
    fwtree = wtree;
    fwhist = whist;
  }
  void Init();

  //! Read the specified buffer and make the Mode3Event.
  int DecodeMode3(char* cBuf, int len, long long int ts);
  //! Passed a gretina Crystal, make a new crystal in the Gretina object.
  int DecodeGretina(Crystal* cryst, long long int gts);

  //! Write the last event to file.
  void WriteLastEvent();

  int NrOfEvents(){return fnentries;}
  int NrOfCalEvents(){return fncalentries;}
  TTree* GetTree(){return ftr;}
  TTree* GetCalTree(){return fcaltr;}
  TTree* GetTraceTree(){return ftr;}
  //! Make HiCARI objects from mode3 data
  void MakeMode2();
  void SetMode3(Mode3Event *m3e){fMode3Event = m3e;}
  HiCARI* GetHiCARI(){return fHiCARI;}
  HiCARICalc* GetHiCARICalc(){return fHiCARICalc;}

protected:
  //! Create a single trace
  Trace DecodeTrace(unsigned short** wBuf_p, int length, long long int gts);
  //! Performs end of event actions.
  void CloseEvent();
  //! Clears memory of current event.
  void ClearEvent();

  TRandom* fRand;
  TTree *ftr;
  TTree *fcaltr;

  Gretina*      fGretina;
  GretinaCalc*  fGretinaCalc;
  Mode3Event *fMode3Event;
  HiCARI *fHiCARI;
  HiCARICalc *fHiCARICalc;

  bool fhasdata;
  long long int fcurrent_ts;
  long long int ffirst_ts;
  int fvl;
  int fnentries;
  int fncalentries;
  int fstrangehits;
  int fctr;
  bool fwtree;
  bool fwhist;
  bool fwcaltree;
  bool fwcalhist;
  int fEventTimeDiff;

  Calibration* fcal;

  RawHistograms* frhist;
  CalHistograms* fchist;

  Settings* fSett;

  bool fmakemode2;
};

#endif
