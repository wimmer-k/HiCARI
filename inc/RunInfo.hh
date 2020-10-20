#ifndef __RUNINFO_HH
#define __RUNINFO_HH

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "TObject.h"

#include "Globaldefs.h"
using namespace std;

/*!
  A container for the analysis runinfo
*/
class RunInfo : public TObject {
public:
  //! default constructor
  RunInfo(){
    Clear();
  }
  void Clear();
  void PrintRunInfo();
  void AppendInfo(RunInfo* sec);

  void SetHIRunNumber(int runnr){fHIrunnr = runnr;}
  void SetBRRunNumber(int runnr){fBRrunnr = runnr;}
  void SetM2RunNumber(int runnr){fM2runnr = runnr;}
  void SetHIBytes(long long int bytes){fHIbytes = bytes;}
  void SetHIRawEvents(long long int events){fHIrawevents = events;}
  void SetHICalEvents(long long int events){fHIcalevents = events;}
  void SetBREvents(long long int events){fBRevents = events;}
  void SetMode2Events(long long int events){fMode2events = events;}
  void SetMergedEvents(long long int events){fMergedevents = events;}
  void SetCorrelationRate(double rate){fCorrRate = rate;}
  void SetBigRIPSCtr(long long int ctr){fBigRIPSctr = ctr;}
  void SetBigRIPSHitCtr(long long int ctr){fBigRIPSHitctr = ctr;}
  void SetHiCARICtr(long long int ctr){fHiCARIctr = ctr;}
  void SetHiCARIHitCtr(long long int ctr){fHiCARIHitctr = ctr;}
  void SetMode2HitCtr(long long int ctr){fMode2Hitctr = ctr;}
  

  int GetHIRunNumber(){return fHIrunnr;}
  int GetBRRunNumber(){return fBRrunnr;}
  int GetM2RunNumber(){return fM2runnr;}
  long long int GetHIBytes(){return fHIbytes;}
  long long int GetHIMegaBytes(){return fHIbytes/(1024*1024);}
  long long int GetHIGigaBytes(){return fHIbytes/(1024*1024*1024);}
  long long int GetHIRawEvents(){return fHIrawevents;}
  long long int GetHICalEvents(){return fHIcalevents;}
  long long int GetBREvents(){return fBRevents;}
  long long int GetMode2Events(){return fMode2events;}
  long long int GetMergedEvents(){return fMergedevents;}
  double GetCorrelationRate(){return fCorrRate;}
  
  long long int GetBigRIPSCtr(){return fBigRIPSctr;}
  long long int GetBigRIPSHitCtr(){return fBigRIPSHitctr;}
  long long int GetHiCARICtr(){return fHiCARIctr;}
  long long int GetHiCARIHitCtr(){return fHiCARIHitctr;}
  long long int GetMode2HitCtr(){return fMode2Hitctr;}

protected:
  int fHIrunnr;
  int fBRrunnr;
  int fM2runnr;

  long long int fHIbytes;
  long long int fHIrawevents;
  long long int fHIcalevents;
  long long int fBRevents;
  long long int fMode2events;
  
  long long int fMergedevents;

  long long int fHiCARIctr;
  long long int fBigRIPSctr;
  long long int fHiCARIHitctr;
  long long int fBigRIPSHitctr;

  long long int fMode2Hitctr;
  
  double fCorrRate;
  
  ClassDef(RunInfo, 1)
};

#endif
