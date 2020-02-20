#ifndef __LISA_HH
#define __LISA_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"
#include "LISAdefs.h"

using namespace std;
//LISA
class LISA : public TObject {
public:
  LISA(){
    Clear();
  }
  ~LISA(){
    Clear();
  };
  void Clear(){
    fts =0;
    fpos = TVector3(0,0,0);
    fdeltaE.clear();
    fntar = 0;
    frtar = -1;
  }

  Float_t GetX(){return fpos.X();}
  Float_t GetY(){return fpos.Y();}
  Float_t GetZ(){return fpos.Z();}
  Short_t GetNTargets(){return fntar;}
  Short_t GetReaction(){return frtar;}
  TVector3 GetVertex(){return fpos;}
  long long int GetTimeStamp(){return fts;};
  Float_t GetDeltaE(Short_t t){
    if(t>-1 && t<fntar)
      return fdeltaE.at(t);
    return sqrt(-1);
  }
  vector<Float_t> GetDeltaE(){
    return fdeltaE;
  }
  
  void SetVertex(Float_t x, Float_t y, Float_t z){fpos = TVector3(x,y,z);}
  void SetVertex(TVector3 pos){fpos = pos;}
  void SetReaction(Short_t t){frtar = t;}
  void SetTimeStamp(long long int time){fts = time;}
  void AddTarget(Float_t delta){
    fdeltaE.push_back(delta);
    fntar++;
  }

protected:
  //! the reaction vertex
  TVector3 fpos;
  //! the number of targets
  Short_t fntar;
  //! the energy loss in the targets
  vector<Float_t> fdeltaE;
  //! time stamp
  long long int fts;
  //! the target where the reaction happend
  Short_t frtar;
  
  // reconstructed stuff
  ClassDef(LISA, 1);
};

#endif
