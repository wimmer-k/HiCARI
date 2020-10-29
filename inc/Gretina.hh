#ifndef __GRETINA_HH
#define __GRETINA_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "Gretinadefs.h"
#include "Settings.hh"

using namespace std;

//! interaction points
struct ip{
  //! coordinates of the interaction point
  float x, y, z;
  //! energy deposition at the interaction point
  float e;
  //! segment ID
  int seg;
  //! energy deposite in that segment
  float seg_ener;
};
//! crystal stuff
struct crys_ips_abcd1234{
  int type;
  int crystal_id;
  int num;
  float tot_e;
  long long int timestamp;
  long long trig_time;
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  int pad;
  ip ips[MAX_INTPTS];
};
//format july 2012
struct crys_ips_abcd5678 {
  int type;          /* defined as abcd5678 */
  int crystal_id;
  int num;           /* # of int pts from decomp, or # of nets on decomp error */
  float tot_e;       /* dnl corrected */
  int core_e[4];     /* 4 raw core energies from FPGA filter (no shift) */
  long long int timestamp;
  long long trig_time;    /* not yet impl */
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  float prestep;    /* avg trace value before step */
  float poststep;   /* avg trace value following step */
  int pad;          /* non-0 on decomp error, value gives error type */
  ip ips[MAX_INTPTS];
};
//format 
struct crys_ips_abcd6789 {
  int type;          /* defined as abcd6789 */
  int crystal_id;
  int num;           /* # of int pts from decomp, or # of nets on decomp error */
  float tot_e;       /* NOT dnl corrected */
  int core_e[4];     /* 4 raw core energies from FPGA filter (no shift) */
  long long int timestamp;
  /* long long trig_time; */   /* not yet impl */
  float tot_e_fixedPickOff_priorEvent1;
  float tot_e_fixedPickOff_priorEvent2;
  float t0;
  /* float cfd; */
  unsigned short int deltaT_priorEvent1;
  unsigned short int deltaT_priorEvent2;
  float chisq;
  float norm_chisq;
  /* float baseline; */
  float tot_e_fixedPickOff_thisEvent;
  float prestep;    /* avg trace value before step */
  float poststep;   /* avg trace value following step */
  int pad;          /* non-0 on decomp error, value gives error type */
  ip ips[MAX_INTPTS];
  // struct {
  //   float x, y, z, e;       /* here e refers to the fraction */
  //   int seg;                /* segment hit */
  //   float seg_ener;         /* energy of hit segment */
  // } intpts[MAX_INTPTS];
};

class IPoint : public TObject {
public:
  IPoint();
  IPoint(Float_t en, Float_t x, Float_t y, Float_t z, Int_t seg, Float_t segen);
  IPoint(IPoint* old);
  ~IPoint(){Clear();}
  void Clear();
  void SetEnergy(Float_t en){fen = en;}
  void SetPosition(Float_t x, Float_t y, Float_t z){fposition.SetXYZ(x,y,z);}
  void SetPosition(TVector3 pos){fposition = pos;}

  Float_t GetEnergy(){return fen;}
  Float_t GetSegEnergy(){return fseg_en;}
  Int_t GetSeg(){return fseg;}
  TVector3 GetPosition(){return fposition;}
  void PrintEvent();

protected:
  //! The energy of the interaction point.
  Float_t fen;
  //! The position relative to the crystal of the interaction point.
  TVector3 fposition;
  //! The segment number (0-35) of the interaction point.
  Int_t fseg;
  //! The total energy deposited in this segment.  (May contain multiple interaction points.)
  Float_t fseg_en;
  ClassDef(IPoint, 1);
};

class Crystal : public TObject {
public:
  Crystal();
  Crystal(crys_ips_abcd1234 inbuf);
  Crystal(crys_ips_abcd5678 inbuf);
  Crystal(crys_ips_abcd6789 inbuf);
  Crystal(Crystal* old);
  ~Crystal();
  void Clear();
  void AddBackCrystal(Crystal* other);
  void AddIP(IPoint *ip);

  void SetTrigTime(long long int trigtime){ftrig_time = trigtime;}
  void SetT0(Float_t t0){ft0 = t0;}
  void SetCFD(Float_t cfd){fcfd = cfd;}
  void SetChiSq(Float_t chisq){fchisq = chisq;}
  void SetNChiSq(Float_t norm_chisq){fnorm_chisq = norm_chisq;}
  void SetBaseline(Float_t baseline){fbaseline = baseline;}
  void SetError(UShort_t error){ferror = error;}
  Short_t GetID(){return fcluster*4 + fcrystalid;}
  Short_t GetCluster(){return fcluster;}
  Short_t GetCrystal(){return fcrystalid;}
  Float_t GetEnergy(){return fen;}
  void SetEnergy(Float_t val){fen = val;}
  Float_t GetIPSum();
  Float_t GetSegmentSum();
  long long int GetTS(){return ftimestamp;}
  long long int GetITS(){return fits;}
  long long int GetTrigTime(){return ftrig_time;}
  Float_t GetT0(){return ft0;}
  Float_t GetCFD(){return fcfd;}
  Float_t GetChiSq(){return fchisq;}
  Float_t GetNChiSq(){return fnorm_chisq;}
  Float_t GetBaseline(){return fbaseline;}
  Short_t GetMult(){return fmult;}
  UShort_t GetError(){return ferror;}
  int GetCoreE(int i){return fcore_e[i];}
  vector<IPoint*> GetIPoints(){return fipoints;}
  IPoint* GetIPoint(int n){return fipoints[n];}
  Short_t GetMaxIPNr(){return fmaxip;}
  IPoint* GetMaxIP(){return fipoints[fmaxip];}
  Float_t GetMaxEn(){return fmaxen;}
  Float_t GetMaxSingleCrystal(){return fMaxSingleCrystal;}
  float GetPreStep(){return fprestep;}
  float GetPostStep(){return fpoststep;}
  void PrintEvent();

protected:
  //! The cluster number (1-30) of the hit.
  Short_t fcluster;
  //! The crystal number (0-3) of the hit
  Short_t fcrystalid;
  //! The energy (keV) of the hit.
  Float_t fen;
  //! No longer used.
  Float_t fMaxSingleCrystal;
  //! The number of interaction points in the crystal.
  Short_t fmult;
  //! Empty as of 2012-08-13
  long long int  ftrig_time;
  //! The CFD time, as calculated by the decomposition.
  Float_t ft0;
  //! Empty as of 2012-08-13
  Float_t fcfd;
  //! chi square of decomposition
  Float_t fchisq;
  //! normalized chi square of decomposition
  Float_t fnorm_chisq;
  //! The baseline for the decomposition.
  Float_t fbaseline;
  vector<IPoint*> fipoints;
  //! The position in fipoints of the IP with the most deposited energy
  Short_t fmaxip;
  //! The energy of the segment with the most deposited energy.
  Float_t fmaxen;
  //! The energy sum of the interaction points.
  Float_t fIPsum;
  //! The 4 core energies.
  /*!
    The 4 core energies.
    According to Mario, for backwards compatibility reasons, the first is always the selected range.
    The remaining 3 are the other ranges, in order.
   */
  int fcore_e[4];
  //! The GEB header timestamp.
  long long int ftimestamp;
  //! The internal timestamp.  (Should always be equal to ftimestamp)
  long long int fits;
  //! A diagnostic added by Mario.
  float fprestep;
  //! A diagnostic added by Mario.
  float fpoststep;
  //! The error code reported by the decomposition.
  UShort_t ferror;

  ClassDef(Crystal, 1);
};

class Gretina : public TObject {
public:
  Gretina();
  ~Gretina(){Clear();}
  void Clear();
  void AddHit(Crystal* cry);

  int GetHitPattern(){return fhitpattern;}
  int GetMult(){return fmult;}
  vector<Crystal*> GetHits(){return fcrystals;}
  Crystal* GetHit(int n){
    if(n<fmult)
      return fcrystals[n];
    return NULL;
  }
  void PrintEvent();
protected:
  //! An integer whose n-th bit is 1 iff the detector in cluster n fired.
  int fhitpattern;
  //! The crystal multiplicity of the event.
  Short_t fmult;
  vector<Crystal*> fcrystals;
  ClassDef(Gretina, 1);
};

//! A class to contain a single calibrated interaction.
/*!
  The HitCalc class is meant to contain a single calibrated interaction.
  This may be the result of multiple addback from multiple sources,
    or may be the result of a single interaction point.
 */
class HitCalc : public TObject {
public:
  //! Required default constructor.
  HitCalc(){
    Clear();
  }
  //! Manually makes a hit.
  HitCalc(Short_t cluster, Short_t crystal, Float_t energy, long long int timestamp, TVector3 pos, Float_t ipsum, Float_t t0, Float_t chisq,int index=-1){
    fcluster = cluster;
    fcrystal = crystal;
    fen = energy;
    fMaxSingleHit = energy;
    fDCen = sqrt(-1);
    ftimestamp = timestamp;
    fposition = pos;
    fIPsum = ipsum;
    ftime = t0;
    fchisq = chisq;
    fHitsAdded = 1;
    fIndex = index;
  }
  //! Creates a hit from the main interaction point of a mode2 crystal.
  HitCalc(Crystal* cry){
    Clear();
    fcluster = cry->GetCluster();
    fcrystal = cry->GetCrystal();
    fen = fMaxSingleHit = cry->GetEnergy();
    fDCen = sqrt(-1);
    ftimestamp = cry->GetTS();
    fIPsum = cry->GetIPSum();
    if (cry->GetMaxIPNr()!=-1){
      fposition = cry->GetMaxIP()->GetPosition();
      fTrueFirst = cry->GetIPoint(0)->GetPosition();
      fUsedTrueFirst = (cry->GetMaxIPNr() == 0);
    }
    ftime = cry->GetT0();
    fchisq = cry->GetChiSq();
    fHitsAdded = 1;
  }
  //! Copies an existing hit.
  HitCalc(HitCalc* hit){
    Clear();
    fcluster = hit->GetCluster();
    fcrystal = hit->GetCrystal();
    fen = hit->GetEnergy();
    fDCen = hit->GetDCEnergy();
    fMaxSingleHit = hit->GetEnergy();
    fIPsum = hit->GetIPSum();
    fposition = hit->GetPosition();
    ftimestamp = hit->GetTS();
    fchisq = hit->GetChiSq();
    ftime = hit->GetTime();
    fHitsAdded = 1;
    fIndex = hit->GetIndex();
  }
  ~HitCalc(){
    Clear();
  }
  void Clear(){
    fcluster = -1;
    fcrystal = -1;
    ftime = fen = fDCen = fMaxSingleHit = fIPsum = sqrt(-1.0);
    ftimestamp = -1;
    fHitsAdded = 0;
    fchisq = 0;
    fposition.SetXYZ(0,0,0);
    fIndex = -1;
    //Simulation parameters
    fTrueFirst.SetXYZ(0,0,0);
    fUsedTrueFirst = false;
  }
  //! Add another HitCalc to the current HitCalc
  /*!
    Add another HitCalc into this one.
    The timestamp and position are taken from the hit with more energy.
   */
  void AddBackHitCalc(HitCalc* other){
    if (other->GetEnergy() > fMaxSingleHit || isnan(fMaxSingleHit)){
      fcluster = other->GetCluster();
      fcrystal = other->GetCrystal();
      ftimestamp = other->GetTS();
      fMaxSingleHit = other->GetEnergy();
      fposition = other->GetPosition();
      fIPsum = other->GetIPSum();
      ftime = other->GetTime();
      fchisq = other->GetChiSq();
    }
    fen += other->GetEnergy();
    fHitsAdded += other->GetHitsAdded();
  }
  void SetDCEnergy(float dcen){fDCen = dcen;}
  void SetPosition(TVector3 in){fposition = in;}
  void SetIPSum(float en){fIPsum = en;}
  void SetChiSq(float chisq){fchisq = chisq;}
  //! An integer identifying the crystal, equal to 4*Cluster+(Position-1)
  Short_t GetID(){return fcluster*4 + fcrystal;}
  //! The number of the cluster, ranging from 1-30.
  Short_t GetCluster(){return fcluster;}
  //! The number of the crystal, ranging from 0-3.
  Short_t GetCrystal(){return fcrystal;}
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
  long long int GetTS(){return ftimestamp;}
  Float_t GetChiSq(){return fchisq;}
  Float_t GetMaxSingleHit(){return fMaxSingleHit;}
  //! The position of the hit in the lab system.
  TVector3 GetPosition(){return fposition;}
  Float_t GetIPSum(){return fIPsum;}
  int GetHitsAdded(){return fHitsAdded;}
  Float_t GetTime(){return ftime;}

  
  void CorrectTime(long long int BRTS){
    ftime += ftimestamp - BRTS;
  }
  //! Apply the Doppler correction using the given settings (for background events without any ZeroDeg information)
  /*!
    Apply the Doppler correction using the given settings.
    Uses the beta and the target position from the settings file.
   */
  void DopplerCorrect(Settings* set);
  //! Returns the Doppler-correction factor to correct the energy.
  static double DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set);
  

  void Print(){
    cout << "cluster " << fcluster << "\tcrystal " << fcrystal << "\ten " << fen << "\tmax hit " << fMaxSingleHit << endl;//"\tipoints " << fipoints.size()<< endl;
    return;
  }

  int GetIndex(){return fIndex;}
  void SetIndex(int i){fIndex = i;}

  //Extra functions for simulated files.
  TVector3 GetTrueFirstPos(){return fTrueFirst;}
  bool GetUsedTrueFirstIP(){return fUsedTrueFirst;}
  double GetThetaFromTrue(){return (fposition.Angle(fTrueFirst));}
  double GetThetaDiff(){return fposition.Theta() - fTrueFirst.Theta();}
  double GetDistanceFromTrue(){return (fposition-fTrueFirst).Mag();}
  double GetPerpDistanceFromTrue(){
    TVector3 diff = fposition - fTrueFirst;
    diff -= diff.Dot(fposition.Unit())*fposition.Unit();
    return diff.Mag();
  }
  double GetDCEnergySimCheat(){return fDCen_simcheat;}
  TVector3 GetTruePos(){return fTrueFirst;}
  TVector3 GetMainPos(){return fposition;}

protected:
  Short_t fcluster;
  Short_t fcrystal;
  Float_t fchisq;
  Float_t fen;
  Float_t fDCen;
  Float_t ftime;
  Float_t fMaxSingleHit; //!
  long long int ftimestamp;
  TVector3 fposition;
  Float_t fIPsum;
  int fHitsAdded;
  int fIndex;
  //Extra variables for simulations
  TVector3 fTrueFirst;
  bool fUsedTrueFirst;
  Float_t fDCen_simcheat;
  ClassDef(HitCalc, 1);
};


//! The calibrated GRETINA object.
/*!
  Holds all calibrated information of GRETINA.
  Holds three lists, keeping the crystal-by-crystal information,
    the add-backed information,
    and the interaction point cluster information.
 */
class GretinaCalc : public TObject {
public:
  GretinaCalc(){
    Clear();
  }
  void Clear(){
    fmult = 0;
    for(vector<HitCalc*>::iterator cry=fhits.begin(); cry!=fhits.end(); cry++){
      delete *cry;
    }
    fhits.clear();
    ClearAddBack();
    ClearCluster();
  }
  void ClearAddBack(){
    fmult_ab = 0;
    for(vector<HitCalc*>::iterator cry=fhits_ab.begin(); cry!=fhits_ab.end(); cry++){
      delete *cry;
    }
    fhits_ab.clear();
  }
  void ClearCluster(){
    for(UShort_t i=0;i<fmult_cl.size();i++){
      for(vector<HitCalc*>::iterator cry=fhits_cl[i].begin(); cry!=fhits_cl[i].end(); cry++){
	delete *cry;
      }
      fhits_cl[i].clear();
    }
    fhits_cl.clear();
    fmult_cl.clear();
    fncluster =0;
  }
  void AddHit(HitCalc* cry){
    fhits.push_back(cry);
    fmult++;
  }
  void AddHitAB(HitCalc* cry){
    fmult_ab++;
    fhits_ab.push_back(cry);
  }
  void AddHitCL(vector<HitCalc*> cry, UShort_t c){
    if(c<fncluster){
      //cout << "invalid " << endl;
      return;
    }
    fncluster++;
    fmult_cl.push_back(cry.size());
    fhits_cl.push_back(cry);
  }
  void AddHitCL(HitCalc* cry, UShort_t c){
    if(c>=fncluster){
      //cout << "invalid " << endl;
      return;
    }
    fmult_cl[c]++;
    fhits_cl[c].push_back(cry);
  }
  void DopplerCorrect(Settings* set);
  void CorrectTime(long long int br_TS);
  void Print(){
    cout << " singles mult " <<fmult << endl;
    for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      (*hit)->Print();
    }
    if(fmult_ab>0){
      cout << "addback mult " <<fmult_ab <<endl;
      for(vector<HitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
	(*hit)->Print();
      }
    }
    if(fncluster>0){
      cout << fncluster << " cluster created" << endl;
      for(UShort_t i=0;i<fncluster;i++){
	cout << i << " mult " << fmult_cl[i] << endl;
	for(vector<HitCalc*>::iterator hit=fhits_cl[i].begin(); hit!=fhits_cl[i].end(); hit++){
	  (*hit)->Print();
	}
      }
    }
    cout << "-----------------------------------"<<endl;
  }
  void PrintCluster(int i){
    cout << i << " mult " << fmult_cl[i] << endl;
    for(vector<HitCalc*>::iterator hit=fhits_cl[i].begin(); hit!=fhits_cl[i].end(); hit++){
      (*hit)->Print();
    }
  }

  int GetMult(){return fmult;}
  vector<HitCalc*> GetHits(){return fhits;}
  HitCalc* GetHit(int n){
    if(n<fmult)
      return fhits[n];
    return NULL;
  }
  unsigned long long int GetTS(){
    if(fmult>0)
      return fhits[0]->GetTS();
    return sqrt(-1);
  }
  vector<HitCalc*> GetHitsAB(){return fhits_ab;}
  int GetMultAB(){return fmult_ab;}
  HitCalc* GetHitAB(int n){
    if(n<fmult_ab)
      return fhits_ab[n];
    return NULL;
  }

  int GetNCluster(){return fncluster;}
  vector<Short_t> GetMultCL(){return fmult_cl;}
  int GetMultCL(int n){return fmult_cl[n];}
  vector<vector<HitCalc*> > GetClusters(){return fhits_cl;}
  vector<HitCalc*> GetClusters(int n){
    if(n<fncluster)
      return fhits_cl[n];
    return vector<HitCalc*>();
  }
  HitCalc* GetHitCL(int n, int m){
    if(n<fncluster && m <fmult_cl.at(n))
      return fhits_cl[n][m];
    return NULL;
  }

  void SetCluster(int n, vector<HitCalc*> hits){
    fhits_cl[n] = hits;
  }

protected:
  Short_t fmult;
  vector<HitCalc*> fhits;
  Short_t fmult_ab;
  vector<HitCalc*> fhits_ab;
  Short_t fncluster;
  vector<Short_t> fmult_cl;
  vector<vector<HitCalc*> > fhits_cl;
  ClassDef(GretinaCalc, 1);
};

//! Allows for comparing of HitCalc objects.
/*!
  Allows for ordering of HitCalc objects by decreasing energy.
  Used as a first guess of an path for tracking.
 */
class HitComparer {
public:
  bool operator() ( HitCalc *lhs, HitCalc *rhs) {
    return (*lhs).GetEnergy() > (*rhs).GetEnergy();
  }
};

//! Class to contain multiple interaction points, used for tracking.
/*!
  Hold an ordered list of HitCalc objects,
    which can be ordered using Tracking::TrackCluster().
  Holds information on the energy of the cluster, and the figure of merit from the tracking.
 */
class GretinaTrack : public TObject {
public:
  GretinaTrack(){
    Clear();
  }
  GretinaTrack(vector<HitCalc*> hits){
    Clear();
    SetHits(hits);
  }
  GretinaTrack(HitCalc* hit){
    Clear();
    AddHit(hit);
  }
  void Clear(){
    fmult = 0;
    fFOM = -2;
    fperm = -1;
    fjumps = 0;
    fesum = -1;
    fDCen = -1;
    for(vector<HitCalc*>::iterator cry=fhits.begin(); cry!=fhits.end(); cry++){
      delete *cry;
    }
    fhits.clear();
  }
  GretinaTrack* Copy(){
    GretinaTrack* output = new GretinaTrack;
    output->fperm = fperm;
    output->fjumps = fjumps;
    output->fesum = fesum;
    output->fDCen = fDCen;
    for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      output->AddHit(new HitCalc(*hit));
    }
    return output;
  }
  void SetHits(vector<HitCalc*> hits){
    fmult = hits.size();
    fhits = hits;
  }
  void AddHit(HitCalc* hit){
    fmult++;
    fhits.push_back(hit);
  }
  void DopplerCorrect(Settings* set){
    fDCen = fesum * HitCalc::DopplerCorrectionFactor(fhits[0]->GetPosition(),set);
  }
  void SetFOM(Double_t fom){ fFOM = fom; }
  void SetPermutation(Int_t perm){ fperm = perm; }
  void SetPermutation(vector<int> perm){ fpermlist = perm; }
  void SetJumps(Double_t jumps){ fjumps = jumps; }
  void SetEsum(Double_t esum){ fesum = esum; }
  void SetDCEnergy(Double_t DCen){ fDCen = DCen; }
  int GetMult(){return fmult;}
  Short_t GetJumps(){return fjumps;}
  Int_t GetPermutation(){return fperm;}
  Double_t GetFOM(){return fFOM;}
  Double_t GetEsum(){return fesum;}
  Double_t GetDCEnergy(){return fDCen;}
  vector<HitCalc*> GetHits(){return fhits;}
  HitCalc* GetHit(int n){return fhits[n];}
  double GetFirstScattering(){
    if (fhits.size()>=2){
      TVector3 pos1 = fhits[0]->GetPosition();
      TVector3 pos2 = fhits[1]->GetPosition();
      return pos1.Angle(pos2-pos1);
    } else {
      return sqrt(-1);
    }
  }
  double GetFirstEnergy(){
    if (fhits.size()>=1){
      return fhits[0]->GetEnergy();
    } else {
      return sqrt(-1);
    }
  }
  double GetFirstEnergyScattering(double initial){
    double en = GetFirstEnergy();
    double costheta = 1 - 511*(1/initial + 1/(initial-en));
    return acos(costheta);
  }
  void WriteIPoints(char* filename){
    FILE* file = fopen(filename,"w");
    for(uint i = 0; i<fhits.size(); i++){
      TVector3 pos =  fhits[i]->GetPosition();
      fprintf(file,
	      "%f\t%f\t%f\t%f\t%d\n",
	      pos.X(),pos.Y(),pos.Z(),
	      fhits[i]->GetEnergy(),
	      fhits[i]->GetIndex());
    }
    fclose(file);
  }
  bool IsTrackedMain(){
    double trackedEn = GetHit(0)->GetEnergy();
    for(int i=1; i<GetMult(); i++){
      if (GetHit(i)->GetEnergy() > trackedEn){
	return false;
      }
    }
    return true;
  }
  void Print(){
    cout << "mult " << fmult << " esum " << fesum << " FOM " << fFOM << endl;
    for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      (*hit)->Print();
    }
  }
protected:
  Short_t fmult;
  Short_t fjumps;
  Int_t fperm;
  Double_t fFOM;
  Double_t fesum;
  Double_t fDCen;
  vector<HitCalc*> fhits;
  vector<int> fpermlist;
  ClassDef(GretinaTrack, 1);
};


class GretinaEvent : public TObject {
public:
  GretinaEvent(){
    Clear();
  }
  void Clear(){
    fmult = 0;
    fmult_ab = 0;
    fmult_raw = 0;
    for(vector<GretinaTrack*>::iterator track=ftracks.begin(); track!=ftracks.end(); track++){
      delete *track;
    }
    ftracks.clear();
    for(vector<HitCalc*>::iterator cry=fhits_ab.begin(); cry!=fhits_ab.end(); cry++){
      delete *cry;
    }
    fhits_ab.clear();
    for(vector<HitCalc*>::iterator cry=fhits_raw.begin(); cry!=fhits_raw.end(); cry++){
      delete *cry;
    }
    fhits_raw.clear();
  }
  void AddTrack(GretinaTrack* track){
    ftracks.push_back(track);
    fmult++;
  }
  void SetHitsRaw(vector<HitCalc*> hits){
    fmult_raw = hits.size();
    fhits_raw = hits;
  }
  void SetHitsAB(vector<HitCalc*> hits){
    fmult_ab = hits.size();
    fhits_ab = hits;
  }

  vector<GretinaTrack*> GetTracks(){return ftracks;}
  GretinaTrack* GetTrack(int n){return ftracks[n];}
  Short_t GetMult(){return fmult;}
  Short_t GetMultAB(){return fmult_ab;}
  Short_t GetMultRaw(){return fmult_raw;}
  void PrintEvent(){
    cout << " mult " << fmult << endl;
    for(vector<GretinaTrack*>::iterator track=ftracks.begin(); track!=ftracks.end(); track++){
      (*track)->Print();
    }
  }
  vector<HitCalc*> GetHitsAB(){return fhits_ab;}
  vector<HitCalc*> GetHitsRaw(){return fhits_raw;}
protected:
  Short_t fmult;
  vector<GretinaTrack*> ftracks;
  Short_t fmult_ab;
  Short_t fmult_raw;
  vector<HitCalc*> fhits_ab;
  vector<HitCalc*> fhits_raw;
  ClassDef(GretinaEvent, 1);
};

#endif
