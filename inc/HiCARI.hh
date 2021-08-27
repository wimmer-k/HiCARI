#ifndef __HICARI_HH
#define __HICARI_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "Settings.hh"
#include "Beam.hh"

using namespace std;

class HiCARIHit : public TObject {
public:
  HiCARIHit();
  HiCARIHit(int clu, int cry,  Short_t nr, double en,  long long int ts, bool tracking);
  ~HiCARIHit();
  void Clear();
  bool InsertCore(int clu, int cry, double en,  long long int ts);
  bool InsertSegment(int clu, int cry, Short_t nr, Float_t en);

  void AddBackHit(HiCARIHit* other);
  void AddSegment(Short_t segnr, Float_t segnen);
  void SetEnergy(Float_t val){fen = val;}
  void SetSegmentEn(int n, Float_t en){
    fsegen[n] = en;
  }

  Short_t GetID(){return fcluster*4 + fcrystal;}
  Short_t GetCluster(){return fcluster;}
  Short_t GetCrystal(){return fcrystal;}
  Float_t GetEnergy(){return fen;}

  Float_t GetSegmentSum();
  int GetMult(){return fmult;}
  long long int GetTS(){return ftimestamp;}

  vector<Float_t> GetSegmentEn(){return fsegen;}
  vector<Short_t> GetSegmentNr(){return fsegnr;}
  
  Float_t GetSegmentEn(int n){
    if(n<fmult)
      return fsegen[n];
    return sqrt(-1);
  }
  Short_t GetSegmentNr(int n){
    if(n<fmult)
      return fsegnr[n];
    return sqrt(-1);
  }

  Short_t GetMaxSegNr(){return fmaxsegnr;}
  Float_t GetMaxSegEn(){return fmaxsegen;}
  
  bool GetTraking(){return ftracking;}
  bool IsTracking(){return ftracking;}
  
  bool IsMiniball(){return (fcluster>-1 && fcluster<6);}
  bool IsSuperClo(){return (fcluster> 5 && fcluster<10);}
  void PrintEvent();

protected:
  //! The cluster number (1-30) of the hit.
  Short_t fcluster;
  //! The crystal number (0-3) of the hit
  Short_t fcrystal;
  //! The energy (keV) of the hit.
  Float_t fen;
  //! The maximum single crystal energy, for addback
  Float_t fmaxsinglecrystal;
  //! The number of segments hit in the crystal.
  Short_t fmult;
  //! The segment energies
  vector<Float_t> fsegen;
  //! The segment numbers
  vector<Short_t> fsegnr;
  //! The maximum segment energy
  Float_t fmaxsegen;
  //! The segment number of the one with the max energy
  Short_t fmaxsegnr;
  //! The GEB header timestamp.
  long long int ftimestamp;
  //! Indicates if it is a tracking detector
  bool ftracking;

  ClassDef(HiCARIHit, 1);
};

class HiCARI : public TObject {
public:
  HiCARI();
  ~HiCARI(){Clear();}
  void Clear();
  void AddHit(HiCARIHit* cry);

  int GetHitPattern(){return fhitpattern;}
  int GetMult(){return fmult;}
  vector<HiCARIHit*> GetHits(){return fhits;}
  HiCARIHit* GetHit(int n){
    if(n<fmult)
      return fhits[n];
    return NULL;
  }
  HiCARIHit* GetHit(int mod, int cry){
    //cout << "trying to get " << mod << " , " << cry << endl;
    for(vector<HiCARIHit*>::iterator his=fhits.begin(); his!=fhits.end(); his++){
      //cout << "finding " << (*his)->GetCluster() << " , " << (*his)->GetCrystal() << endl;
      if((*his)->GetCluster() == mod && (*his)->GetCrystal() == cry)
	return (*his);
    }
    return NULL;
  }
  void PrintEvent();
protected:
  //! An integer whose n-th bit is 1 iff the detector in cluster n fired.
  int fhitpattern;
  //! The crystal multiplicity of the event.
  Short_t fmult;
  vector<HiCARIHit*> fhits;
  ClassDef(HiCARI, 1);
};


class HiCARIHitCalc : public TObject {
public:
  HiCARIHitCalc(){Clear();}
  HiCARIHitCalc(Short_t clu, Short_t cry, Short_t maxseg, Float_t segsum, TVector3 pos, Float_t en, long long int ts);
  HiCARIHitCalc(HiCARIHitCalc* hit);
  ~HiCARIHitCalc(){Clear();}
  void Clear();
  void SetSegments(vector<Short_t> nr, vector<Float_t> en);

  void AddBackHiCARIHitCalc(HiCARIHitCalc* hit);
  void SetDCEnergy(float dcen){fDCen = dcen;}
  void SetPosition(TVector3 in){fposition = in;}
  void AddSegment(Short_t nr, Float_t en);
  
  // ! which cluster
  Short_t GetCluster(){return fcluster;}
  //! The number of the crystal, ranging from 0-3.
  Short_t GetCrystal(){return fcrystal;}
  //! The number of the max segment
  Short_t GetMaxSegment(){return fmaxseg;}
  Float_t GetMaxSegmentEnergy(){
    return *max_element(fsegen.begin(), fsegen.end());
  }
  //! The energy of the hit (keV).
  Float_t GetEnergy(){return fen;}
  //! The sum of the segment energies of the hit (keV).
  Float_t GetSegSum(){return fsegsum;}
  vector<Float_t> GetSegmentEn(){return fsegen;}
  vector<Short_t> GetSegmentNr(){return fsegnr;}
  Float_t GetTime(){return ftime;}

  void CorrectTime(long long int BRTS){
    //cout << "correct " << ftimestamp << "\t" << BRTS << "\t";
    ftime += ftimestamp - BRTS;
    //cout << ftimestamp << endl;
  }
 
  //! The Doppler-corrected energy of the hit.
  Float_t GetDCEnergy(){return fDCen;}
  Float_t GetDCEnergy(float beta, double x=0, double y=0, double z=0){
    TVector3 pos;
    pos.SetXYZ(x,y,z);//in mm
    pos += GetPosition();
    return fen/sqrt(1-beta*beta)*(1-beta*cos(pos.Theta()));
  }
  Float_t GetDCEnergy(float beta, TVector3 gpos){
    return fen/sqrt(1-beta*beta)*(1-beta*cos(gpos.Theta()));
  }
  Float_t GetDCEnergy(float beta, TVector3 gpos, TVector3  ejecdir){
    return fen/sqrt(1-beta*beta)*(1-beta*cos(gpos.Angle(ejecdir)));
  }
  Float_t GetMaxSingleHit(){return fmaxhit;}
  //! The position of the hit in the lab system.
  TVector3 GetPosition(){return fposition;}
  int GetHitsAdded(){return fHitsAdded;}
  unsigned long long int GetTS(){return ftimestamp;}

  void DopplerCorrect(Settings* set){
    fDCen = fen*DopplerCorrectionFactor(GetPosition(),set);
  }
  //! Returns the Doppler-correction factor to correct the energy.
  double DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set);

  void DopplerCorrect(Beam* beam){
    fDCen = fen*DopplerCorrectionFactor(GetPosition(),beam);
  }
  //! Returns the Doppler-correction factor to correct the energy.
  double DopplerCorrectionFactor(TVector3 PosToTarget, Beam* beam);


  
  bool IsMiniball(){return (fcluster>-1 && fcluster<6);}
  bool IsSuperClo(){return (fcluster> 5 && fcluster<10 && fcrystal<4);}
  bool IsTracking(){return (fcluster> 9 && fcluster<12);}
  bool IsBigRIPS(){return (fcluster==9 && fcrystal==9);}
  bool IsHiCARI(){return !(fcluster==9 && fcrystal==9);}

  void Print(){
    cout << "HiCARI: cluster " << fcluster << "\tcrystal " << fcrystal << "\tmaxseg " << fmaxseg << "\tsegsum " << fsegsum << "\ten " << fen << "\tmax hit " << fmaxhit << endl;
    cout << " x " << fposition.X() << " y " << fposition.Y() << " z " << fposition.Z() << endl;
    for(UShort_t s=0;s<fsegnr.size();s++){
      cout << "segment " << fsegnr.at(s) << ", en " << fsegen.at(s) << endl;
    }
    return;
  }


 
protected:
  Short_t fcluster;
  Short_t fcrystal;

  Short_t fmaxseg;
  Float_t fsegsum;
  vector<Float_t> fsegen;
  vector<Short_t> fsegnr;

  Float_t ftime;
  Float_t fen;
  Float_t fDCen;
  TVector3 fposition;
  int fHitsAdded;
  Float_t fmaxhit; // max single hit when addback is used
  long long int ftimestamp;
  
  ClassDef(HiCARIHitCalc, 1)
};

//! The calibrated HiCARI object.
/*!
  Holds all calibrated information of HiCARI.
  Holds three lists, keeping the crystal-by-crystal information,
    the add-backed information
 */
class HiCARICalc : public TObject {
public:
  HiCARICalc(){
    Clear();
  }
  void Clear(){
    ftimestamp = -1;
    fmult = 0;
    fHImult = 0;
    for(vector<HiCARIHitCalc*>::iterator cry=fhits.begin(); cry!=fhits.end(); cry++){
      delete *cry;
    }
    fhits.clear();
    fhadBigRIPS = false;
    ClearAddBack();
  }
  void ClearAddBack(){
    fmult_ab = 0;
    for(vector<HiCARIHitCalc*>::iterator cry=fhits_ab.begin(); cry!=fhits_ab.end(); cry++){
      delete *cry;
    }
    fhits_ab.clear();
  }
  void AddHit(HiCARIHitCalc* cry, bool isBigRIPS){
    fhits.push_back(cry);
    
    ftimestamp = cry->GetTS();
    if(ftimestamp<0){
      cout << "invalid time stamp? " << ftimestamp << " with energy " << cry->GetEnergy() << endl;
      cout << "clu " << cry->GetCluster() << ", cry " << cry->GetCrystal() << endl;
    }

    if(isBigRIPS){
      fhadBigRIPS = true;
      fnBigRIPS = fmult;
    }
    else
      fHImult++;
    fmult++;
  }
  void AddHitAB(HiCARIHitCalc* cry){
    fhits_ab.push_back(cry);
    fmult_ab++;
  }
  bool HadBigRIPS(){return fhadBigRIPS;}
  
  void DopplerCorrect(Settings* set);
  void DopplerCorrect(Beam* beam);
  void CorrectTime(long long int br_TS);

  void Print(){
    cout << " singles mult " <<fmult << endl;
    for(vector<HiCARIHitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      (*hit)->Print();
    }
    if(fmult_ab>0){
      cout << "addback mult " <<fmult_ab <<endl;
      for(vector<HiCARIHitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
	(*hit)->Print();
      }
    }
     cout << "-----------------------------------"<<endl;
  }
  
  int GetMult(){return fmult;}
  vector<HiCARIHitCalc*> GetHits(){return fhits;}
  HiCARIHitCalc* GetBigRIPSHit(){
    if(fhadBigRIPS)
      return fhits[fnBigRIPS];
    return NULL;
  }
  HiCARIHitCalc* GetHit(int n){
    if(n<fmult)
      return fhits[n];
    return NULL;
  }
  vector<HiCARIHitCalc*> GetHitsAB(){return fhits_ab;}
  int GetMultAB(){return fmult_ab;}
  HiCARIHitCalc* GetHitAB(int n){
    if(n<fmult_ab)
      return fhits_ab[n];
    return NULL;
  }
  unsigned long long int GetTS(){return ftimestamp;}

protected:
  long long int ftimestamp;
  Short_t fmult;
  Short_t fHImult;
  vector<HiCARIHitCalc*> fhits;
  Short_t fmult_ab;
  vector<HiCARIHitCalc*> fhits_ab;
  bool fhadBigRIPS;
  Short_t fnBigRIPS;
  ClassDef(HiCARICalc, 1);
};

#endif
