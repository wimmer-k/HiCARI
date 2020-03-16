#ifndef __MINIBALL_HH
#define __MINIBALL_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "ZeroDeg.hh"
#include "MINOS.hh"
#include "Settings.hh"
#include "Miniballdefs.h"

using namespace std;
struct MB_seg{
  int seg_id;
  float e;
};

struct MB_clust{
  int type;
  int crystal_id;
  int num;
  float tot_e;
  MB_seg seg[MBSEGS];
};

class MBCrystal : public TObject {
public:
  MBCrystal();
  MBCrystal(MB_clust inbuf, long long int ts);
  MBCrystal(MBCrystal* old);
  MBCrystal(int clu, int cry,  Short_t nr, double en,  long long int ts);
  ~MBCrystal();
  void Clear();
  bool InsertCore(int clu, int cry, double en,  long long int ts);
  bool InsertSegment(int clu, int cry, Short_t nr, Float_t en);

  void AddBackCrystal(MBCrystal* other);
  void AddSegment(Short_t segnr, Float_t segnen);

  Short_t GetID(){return fcluster*4 + fcrystalid;}
  Short_t GetCluster(){return fcluster;}
  Short_t GetCrystal(){return fcrystalid;}
  Float_t GetEnergy(){return fen;}
  void SetEnergy(Float_t val){fen = val;}
  void SetSegmentEn(int n, Float_t en){
      fsegen[n] = en;
  }

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

  void PrintEvent();

protected:
  //! The cluster number (1-30) of the hit.
  Short_t fcluster;
  //! The crystal number (0-3) of the hit
  Short_t fcrystalid;
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
  
  ClassDef(MBCrystal, 1);
};

class Miniball : public TObject {
public:
  Miniball();
  ~Miniball(){Clear();}
  void Clear();
  void AddHit(MBCrystal* cry);

  int GetHitPattern(){return fhitpattern;}
  int GetMult(){return fmult;}
  vector<MBCrystal*> GetHits(){return fcrystals;}
  MBCrystal* GetHit(int n){
    if(n<fmult)
      return fcrystals[n];
    return NULL;
  }
  MBCrystal* GetHit(int mod, int cry){
    //cout << "trying to get " << mod << " , " << cry << endl;
    for(vector<MBCrystal*>::iterator mbs=fcrystals.begin(); mbs!=fcrystals.end(); mbs++){
      //cout << "finding " << (*mbs)->GetCluster() << " , " << (*mbs)->GetCrystal() << endl;
      if((*mbs)->GetCluster() == mod && (*mbs)->GetCrystal() == cry)
	return (*mbs);
    }
    return NULL;
  }
  void PrintEvent();
protected:
  //! An integer whose n-th bit is 1 iff the detector in cluster n fired.
  int fhitpattern;
  //! The crystal multiplicity of the event.
  Short_t fmult;
  vector<MBCrystal*> fcrystals;
  ClassDef(Miniball, 1);
};


class MBHitCalc : public TObject {
public:
  MBHitCalc(){Clear();}
  MBHitCalc(Short_t clu, Short_t cry, Short_t seg, TVector3 pos, Float_t en, long long int ts);
  MBHitCalc(MBHitCalc* hit);
  ~MBHitCalc(){Clear();}
  void Clear();
  void AddBackMBHitCalc(MBHitCalc* hit);
  void SetDCEnergy(float dcen){fDCen = dcen;}
  void SetPosition(TVector3 in){fposition = in;}

  // ! which cluster
  Short_t GetCluster(){return fcluster;}
  //! The number of the crystal, ranging from 0-3.
  Short_t GetCrystal(){return fcrystal;}
  //! The number of the segment number
  Short_t GetSegment(){return fsegment;}
  //! The energy of the hit (keV).
  Float_t GetEnergy(){return fen;}
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
  long long int GetTS(){return ftimestamp;}

  //void DopplerCorrect(double beta, double z = 0){
  //  fDCen = fen * MBHitCalc::DopplerCorrectionFactor(GetPosition(),beta,z);
  //}
  void DopplerCorrect(Settings* set){
    fDCen = fen*MBHitCalc::DopplerCorrectionFactor(GetPosition(),set);
  }
  //! Returns the Doppler-correction factor to correct the energy.
  static double DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set);

  void DopplerCorrect(Settings* set, ZeroDeg* zerodeg){
    fDCen = fen*MBHitCalc::DopplerCorrectionFactor(GetPosition(),set,zerodeg);
  }
  //! Returns the Doppler-correction factor to correct the energy.
  static double DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set, ZeroDeg* zerodeg);
  
  void DopplerCorrect(Settings* set, ZeroDeg* zerodeg, MINOS* minos){
    fDCen = fen*MBHitCalc::DopplerCorrectionFactor(GetPosition(),set,zerodeg,minos);
  }
  //! Returns the Doppler-correction factor to correct the energy.
  static double DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set, ZeroDeg* zerodeg, MINOS* minos);
  //static double DopplerCorrectionFactor(TVector3 PosToTarget, double beta, double z);

  void Print(){
    cout << "Miniball: cluster " << fcluster << "\tcrystal " << fcrystal << "\tsegment " << fsegment << "\ten " << fen << "\tmax hit " << fmaxhit << endl;
    return;
  }


 
protected:
  Short_t fcluster;
  Short_t fcrystal;
  Short_t fsegment;
  double fen;
  double fDCen;
  TVector3 fposition;
  int fHitsAdded;
  double fmaxhit;
  long long int ftimestamp;
  
  ClassDef(MBHitCalc, 1)
};

//! The calibrated MINIBALL object.
/*!
  Holds all calibrated information of MINIBALL.
  Holds three lists, keeping the crystal-by-crystal information,
    the add-backed information,
    and the interaction point cluster information.
 */
class MiniballCalc : public TObject {
public:
  MiniballCalc(){
    Clear();
  }
  void Clear(){
    fmult = 0;
    for(vector<MBHitCalc*>::iterator cry=fhits.begin(); cry!=fhits.end(); cry++){
      delete *cry;
    }
    fhits.clear();
    ClearAddBack();
  }
  void ClearAddBack(){
    fmult_ab = 0;
    for(vector<MBHitCalc*>::iterator cry=fhits_ab.begin(); cry!=fhits_ab.end(); cry++){
      delete *cry;
    }
    fhits_ab.clear();
  }
  void AddHit(MBHitCalc* cry){
    fhits.push_back(cry);
    fmult++;
  }
  void AddHitAB(MBHitCalc* cry){
    fmult_ab++;
    fhits_ab.push_back(cry);
  }
  void DopplerCorrect(Settings* set);
  void DopplerCorrect(Settings* set, ZeroDeg* zerodeg);
  void DopplerCorrect(Settings* set, ZeroDeg* zerodeg, MINOS* minos);
  void Print(){
    cout << " singles mult " <<fmult << endl;
    for(vector<MBHitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      (*hit)->Print();
    }
    if(fmult_ab>0){
      cout << "addback mult " <<fmult_ab <<endl;
      for(vector<MBHitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
	(*hit)->Print();
      }
    }
     cout << "-----------------------------------"<<endl;
  }
  
  int GetMult(){return fmult;}
  vector<MBHitCalc*> GetHits(){return fhits;}
  MBHitCalc* GetHit(int n){
    if(n<fmult)
      return fhits[n];
    return NULL;
  }
  vector<MBHitCalc*> GetHitsAB(){return fhits_ab;}
  int GetMultAB(){return fmult_ab;}
  MBHitCalc* GetHitAB(int n){
    if(n<fmult_ab)
      return fhits_ab[n];
    return NULL;
  }

protected:
  Short_t fmult;
  vector<MBHitCalc*> fhits;
  Short_t fmult_ab;
  vector<MBHitCalc*> fhits_ab;
  ClassDef(MiniballCalc, 1);
};

#endif
