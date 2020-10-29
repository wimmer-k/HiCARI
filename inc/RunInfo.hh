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
  void SetM2Bytes(long long int bytes){fM2bytes = bytes;}
  void SetM2RawEvents(long long int events){fM2rawevents = events;}
  void SetM2CalEvents(long long int events){fM2calevents = events;}
  
  void SetMergedEvents(long long int events){fMergedevents = events;}
  void SetCorrelationRate(double rate){fCorrRate = rate;}

  void SetBigRIPSCtr(long long int ctr){fBigRIPSctr = ctr;}
  void SetBigRIPSHitCtr(long long int ctr){fBigRIPSHitctr = ctr;}
  void SetHiCARICtr(long long int ctr){fHiCARIctr = ctr;}
  void SetHiCARIHitCtr(long long int ctr){fHiCARIHitctr = ctr;}
  void SetHiCARIHitABCtr(long long int ctr){fHiCARIHitABctr = ctr;}
  void SetGretinaCtr(long long int ctr){fGretinactr = ctr;}
  void SetGretinaHitCtr(long long int ctr){fGretinaHitctr = ctr;}
  void SetGretinaHitABCtr(long long int ctr){fGretinaHitABctr = ctr;}
  

  int GetHIRunNumber(){return fHIrunnr;}
  int GetBRRunNumber(){return fBRrunnr;}
  int GetM2RunNumber(){return fM2runnr;}
  long long int GetHIBytes(){return fHIbytes;}
  long long int GetHIMegaBytes(){return fHIbytes/(1024*1024);}
  long long int GetHIGigaBytes(){return fHIbytes/(1024*1024*1024);}
  long long int GetHIRawEvents(){return fHIrawevents;}
  long long int GetHICalEvents(){return fHIcalevents;}
  long long int GetM2Bytes(){return fM2bytes;}
  long long int GetM2MegaBytes(){return fM2bytes/(1024*1024);}
  long long int GetM2GigaBytes(){return fM2bytes/(1024*1024*1024);}
  long long int GetM2RawEvents(){return fM2rawevents;}
  long long int GetM2CalEvents(){return fM2calevents;}
  long long int GetBREvents(){return fBRevents;}
  long long int GetMergedEvents(){return fMergedevents;}
  double GetCorrelationRate(){return fCorrRate;}
  
  long long int GetBigRIPSCtr(){return fBigRIPSctr;}
  long long int GetBigRIPSHitCtr(){return fBigRIPSHitctr;}
  long long int GetHiCARICtr(){return fHiCARIctr;}
  long long int GetHiCARIHitCtr(){return fHiCARIHitctr;}
  long long int GetHiCARIHitABCtr(){return fHiCARIHitABctr;}
  long long int GetGretinaCtr(){return fGretinactr;}
  long long int GetGretinaHitCtr(){return fGretinaHitctr;}
  long long int GetGretinaHitABCtr(){return fGretinaHitABctr;}

protected:
  int fHIrunnr;
  int fBRrunnr;
  int fM2runnr;

  long long int fHIbytes;
  long long int fHIrawevents;
  long long int fHIcalevents;
  long long int fBRevents;
  long long int fM2bytes;
  long long int fM2rawevents;
  long long int fM2calevents;
  
  long long int fMergedevents;

  long long int fHiCARIctr;
  long long int fBigRIPSctr;
  long long int fHiCARIHitctr;
  long long int fHiCARIHitABctr;
  long long int fBigRIPSHitctr;

  long long int fGretinactr;
  long long int fGretinaHitctr;
  long long int fGretinaHitABctr;
  
  double fCorrRate;
  
  ClassDef(RunInfo, 1)
};

#endif
