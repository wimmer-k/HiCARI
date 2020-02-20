#ifndef __MINOS_HH
#define __MINOS_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"
#include "MINOSdefs.h"

using namespace std;
//MINOS
class MINOS : public TObject {
public:
  MINOS(){
    Clear();
  }
  ~MINOS(){
    Clear();
  };
  void Clear(){
    fts =0;
    fpos = TVector3(0,0,0);
    fbetare = sqrt(-1.0);
  }

  Float_t GetX(){return fpos.X();}
  Float_t GetY(){return fpos.Y();}
  Float_t GetZ(){return fpos.Z();}
  Float_t GetBetaRE(){return fbetare;}
  Float_t GetPhi(){return fpos.Phi();}
  Float_t GetTheta(){return fpos.Theta();}
  TVector3 GetVertex(){return fpos;}
  long long int GetTimeStamp(){return fts;};
  Float_t GetBeta(){return fbeta;}

  void SetVertex(Float_t x, Float_t y, Float_t z){fpos = TVector3(x,y,z);}
  void SetVertex(TVector3 pos){fpos = pos;}
  void SetBetaRE(Float_t value){fbetare = value;}
  void SetTimeStamp(long long int time){fts = time;}
  void SetBeta(Float_t value){fbeta = value;}

protected:
  //! the reaction vertex
  TVector3 fpos;
  //! The beta at the reaction vertex
  Float_t fbetare;
  //! time stamp
  long long int fts;
  //! beta reconstructed from position
  Float_t fbeta;

  ClassDef(MINOS, 1);
};

#endif
