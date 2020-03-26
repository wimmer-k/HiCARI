#ifndef __TRACE_HH
#define __TRACE_HH

#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include "TObject.h"
#include "Tracedefs.h"

using namespace std;

class Settings;

class Trace : public TObject {
public:
  Trace(){
    Clear();
  }
  void Clear(){
    flength = -1;
    fboard = -1;
    fchn = -1;
    fslot = -1;
    fcrystal = -1;
    fhole = -1;
    fLED_ts = -1;
    fen = -1;
    fen_sign = false;
    fpileup = -1;
    fCFD_ts = -1;
    fCFD[0] = -1;
    fCFD[1] = -1;
    fts = -1;
    ftrace.clear();
    fcore = false;
    ftdiff = -1;
  }
  void SetTS(long long int ts){fts = ts;}
  void SetLength(int length){
    flength = length;
    ftrace.resize(length);
  }
  void SetBoard(int boardid){fboard = boardid;}
  void SetChn(int chnid){fchn = chnid;}
  void SetSlot(int slot){fslot = slot;}
  void SetCrystal(int crystal){fcrystal = crystal;}
  void SetHole(int hole){fhole = hole;}
  void SetLED(long long int led){fLED_ts = led;}
  void SetEnergy(int en){fen = en;}
  void SetEnSign(bool sign){fen_sign = sign;}
  void SetPileUp(int pileup){fpileup = pileup;}
  void SetCFD(long long int cfd){fCFD_ts = cfd;}
  void SetCFD(int n, int cfd){fCFD[n] = cfd;}
  void SetTrace(int n, short v){ftrace[n] = v;}
  void SetCore(bool core){fcore = core;}
  void SetTDiff(long long int tdiff){ftdiff = tdiff;}

  long long int GetTS(){return fts;}
  int GetLength(){return flength;}
  int GetBoard(){return fboard;}
  int GetChn(){return fchn;}
  int GetSlot(){return fslot;}
  int GetCrystal(){return fcrystal;}
  int GetHole(){return fhole;}
  int GetID(){return 4*fhole + fcrystal;}
  long long int GetLED(){return fLED_ts;}
  int GetEnergy(){return fen;}
  bool GetEnSign(){return fen_sign;}
  int GetPileUp(){return fpileup;}
  long long int GetCFD(){return fCFD_ts;}
  int GetCFD(int n){return fCFD[n];}
  vector <short> GetTrace(){return ftrace;}
  bool GetCore(){return fcore;}
  long long int GetTDiff(){return ftdiff;}

  void Print(){
    cout << "--------------------- Trace Print ------------------------" << endl;
    cout << "flength  = " << flength  << endl;
    cout << "fboard   = " << fboard   << endl;
    cout << "fchn     = " << fchn     << endl;
    cout << "fslot    = " << fslot    << endl;
    cout << "fcrystal = " << fcrystal << endl;
    cout << "fhole    = " << fhole    << endl;
    cout << "fLED_ts  = " << fLED_ts  << endl;
    cout << "fen      = " << fen      << endl;
    cout << "fen_sign = " << fen_sign << endl;
    cout << "fpileup  = " << fpileup  << endl;
    cout << "fCFD_ts  = " << fCFD_ts  << endl;
    cout << "fCFD[0]  = " << fCFD[0]  << endl;
    cout << "fCFD[1]  = " << fCFD[1]  << endl;
    cout << "fts      = " << fts      << endl;
    cout << "ftdiff   = " << ftdiff   << endl;
  }
  //int GetSegNum(Settings*);
protected:
  //! The number of points in the trace.
  int flength;
  //! The board from which this trace came.
  int fboard;
  //! The channel from which this trace came.
  int fchn;
  int fslot;
  int fcrystal;
  //! The hole of the detector.
  int fhole;
  long long fLED_ts;
  int fen;
  bool fen_sign;
  int fpileup;
  long long fCFD_ts;
  int fCFD[2];
  //if you want to write traces to the root file remove the "//!" in the next line
  vector <short> ftrace;
  long long int fts;
  long long int ftdiff;
  bool fcore;
  ClassDef(Trace, 1);
};


class Mode3Hit : public TObject {
public:
  Mode3Hit(){
    Clear();
  }
  void Clear(){
    ftrace.clear();
    fmult = 0;
    fcore = -1;
    fts = -1;
  }
  void AddTrace(Trace add){
    if(fmult>39){
      cout << "adding segment mult > 39" << endl;
      if (add.GetCore()){
	cout << "    the dropped trace was the core trace" << endl;
      } else {
	cout << "    the dropped trace was not the core trace" << endl;
      }
      return;
    }
    ftrace.push_back(add);
    if(add.GetCore()){
      fcore = fmult;
      fts = add.GetLED();
    }
    fmult++;
  }

  UShort_t GetMult(){return fmult;}
  int GetCoreN(){return fcore;}
  Trace* GetTrace(int n){return &ftrace[n];}
  Trace* GetCoreTrace(){
    if(fcore>-1&&fcore<fmult)
      return &ftrace[fcore];
    return NULL;
  }
  vector<Trace>* GetTrace(){return &ftrace;}
  long long int GetTS(){return fts;}

  //crys_ips_abcd5678 MakeMode2(Settings*);
protected:
  UShort_t fmult;
  long long int fts;
  int fcore;
  vector<Trace> ftrace;

  ClassDef(Mode3Hit, 1);
};


class Mode3Event : public TObject {
public:
  Mode3Event(){
    Clear();
  }
  ~Mode3Event(){
    Clear();
  }
  void Clear(){
    for(vector<Mode3Hit*>::iterator hit = fhit.begin(); hit!=fhit.end(); hit++){
      delete *hit;
    }
    fhit.clear();
    fmult = 0;
    fctr = 0;
  }
  void AddHit(Mode3Hit* add){
    fhit.push_back(add);
    fmult++;
  }
  void AddTrace(Trace add){
    //Look through to see if there are any current hits to which the new trace should be added.
    bool found = false;
    for(vector<Mode3Hit*>::iterator hit_p = fhit.begin(); hit_p!=fhit.end(); hit_p++){
      Mode3Hit* hit = *hit_p;
      if(hit->GetTrace(0)->GetID()==add.GetID() &&
	 abs(hit->GetTrace(0)->GetLED() - add.GetLED()) < MAX_TDIFF_HIT){
	hit->AddTrace(add);
	found = true;
	break;
      }
    }
    //No hit is the same, so we need to make a new hit.
    if (!found){
      Mode3Hit* hit = new Mode3Hit;
      hit->AddTrace(add);
      AddHit(hit);
    }
  }
  void SetCounter(int ctr){fctr = ctr;}
  int GetMult(){return fmult;}
  int GetCounter(){return fctr;}
  Mode3Hit* GetHit(int n){
    if(n<fmult)
      return fhit[n];
    cout << " requesting hit " << n << " mult is " << fmult << endl;
    return NULL;
  }
  vector<Mode3Hit*>* GetHit(){return &fhit;}
protected:
  //! The number of hits in the event.
  int fmult;
  int fctr;
  vector<Mode3Hit*> fhit;

  ClassDef(Mode3Event, 1);
};

#endif

