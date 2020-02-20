#ifndef __GAMMASIM_HH
#define __GAMMASIM_HH
#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"

#include "Simdefs.h"
using namespace std;


class EmittedGamma : public TObject {
public:
  EmittedGamma(){
    Clear();
  }
  EmittedGamma(EG* eg) {
    fenergy = eg->e;
    femissionX = eg->x;
    femissionY = eg->y;
    femissionZ = eg->z;
    femissionPhi = eg->phi;
    femissionTheta = eg->theta;
    femissionBeta = eg->beta;
  }
  EmittedGamma(Float_t e,
	       Float_t x, Float_t y, Float_t z,
	       Float_t ph, Float_t th, Float_t b){
    fenergy = e;
    femissionX = x;
    femissionY = y;
    femissionZ = z;
    femissionPhi = ph;
    femissionTheta = th;
    femissionBeta = b;
  }
  ~EmittedGamma(){
  }
  void Clear(){
    fenergy = 0;
    femissionX = 0;
    femissionY = 0;
    femissionZ = 0;
    femissionPhi = 0;
    femissionTheta = 0;
    femissionBeta = 0;
  }
  Float_t GetEnergy(){ return fenergy; }
  TVector3 GetPos(){
    return TVector3(femissionX, femissionY, femissionZ);
  }
  Float_t GetX(){ return femissionX; }
  Float_t GetY(){ return femissionY; }
  Float_t GetZ(){ return femissionZ; }
  Float_t GetPhi(){ return femissionPhi; }
  Float_t GetTheta(){ return femissionTheta; }
  Float_t GetBeta(){ return femissionBeta; }

  void SetEnergy(Float_t e) { fenergy = e; }
  void SetX(Float_t x) { femissionX = x; }
  void SetY(Float_t y) { femissionY = y; }
  void SetZ(Float_t z) { femissionZ = z; }
  void SetPos(Float_t x, Float_t y, Float_t z) {
    femissionX = x;
    femissionY = y;
    femissionZ = z;
  }
  void SetPhi(Float_t ph) { femissionPhi = ph; }
  void SetTheta(Float_t th) { femissionTheta = th; }
  void SetBeta(Float_t b) { femissionBeta = b; }

protected:
  Float_t fenergy;
  Float_t femissionX;
  Float_t femissionY;
  Float_t femissionZ;
  Float_t femissionPhi;
  Float_t femissionTheta;
  Float_t femissionBeta;

  ClassDef(EmittedGamma, 1);
};

class GammaSim : public TObject {
public:
  GammaSim(){
    Clear();
  }
  ~GammaSim(){
    Clear();
  }
  void Clear(){
    ftimestamp = 0;
    fmult = 0;
    for(vector<EmittedGamma*>::iterator eg=fgammas.begin(); eg!=fgammas.end(); eg++){
      delete *eg;
    }
    fgammas.clear();
  }
  EmittedGamma* GetEmittedGamma(Int_t i){return fgammas[i];}
  Int_t GetMult(){return fmult;}
  long long int GetTimeStamp(){return ftimestamp;}

  void SetTimeStamp(long long int ts){ ftimestamp = ts; }
  void AddEmittedGamma(Float_t e, Float_t x, Float_t y, Float_t z,
		       Float_t ph, Float_t th, Float_t b){
    fmult++;
    fgammas.push_back( new EmittedGamma(e, x, y, z, ph, th, b) );
  }
  void AddEmittedGamma(EG* eg){
    fmult++;
    fgammas.push_back( new EmittedGamma(eg->e,
					eg->x, eg->y, eg->z,
					eg->phi, eg->theta, eg->beta) );
  }

protected:
  long long int ftimestamp;
  Short_t fmult;
  vector<EmittedGamma*> fgammas;
  ClassDef(GammaSim, 1);
};

#endif
