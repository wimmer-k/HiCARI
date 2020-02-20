#ifndef __ZERODEG_HH
#define __ZERODEG_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"

using namespace std;
//ZeroDeg
class ZeroDeg : public TObject {
public:
  ZeroDeg(){
    Clear();
  }
  ~ZeroDeg(){
    Clear();
  };
  void Clear(){
    fts =0;
    fata = sqrt(-1.0);
    fbta = sqrt(-1.0);
    fxta = sqrt(-1.0);
    fyta = sqrt(-1.0);
    fbetata = sqrt(-1.0);
    fazita = sqrt(-1.0);
    fscatter = sqrt(-1.0);
    fptot = sqrt(-1.0);
    fppar = sqrt(-1.0);
    fptra = sqrt(-1.0);
    fetot = sqrt(-1.0);
    fekin = sqrt(-1.0);
  }

  Float_t GetATA(){return fata;}
  Float_t GetBTA(){return fbta;}
  Float_t GetXTA(){return fxta;}
  Float_t GetYTA(){return fyta;}
  Float_t GetBetaTA(){return fbetata;}
  Float_t GetPhi(){return fazita;}
  Float_t GetTheta(){return fscatter;}
  Float_t GetAzita(){return fazita;}
  Float_t GetScatter(){return fscatter;}
  Float_t GetPtot(){return fptot;}
  Float_t GetPpar(){return fppar;}
  Float_t GetPtra(){return fptra;}
  Float_t GetEtot(){return fetot;}
  Float_t GetEkin(){return fekin;}
  long long int GetTimeStamp(){return fts;};

  void SetATA(Float_t value){fata = value;}
  void SetBTA(Float_t value){fbta = value;}
  void SetXTA(Float_t value){fxta = value;}
  void SetYTA(Float_t value){fyta = value;}
  void SetBetaTA(Float_t value){fbetata = value;}
  void SetPhi(Float_t value){fazita = value;}
  void SetTheta(Float_t value){fscatter = value;}
  void SetPtot(Float_t value){fptot = value;}
  void SetPpar(Float_t value){fppar = value;}
  void SetPtra(Float_t value){fptra = value;}
  void SetEtot(Float_t value){fetot = value;}
  void SetEkin(Float_t value){fekin = value;}
  void SetTimeStamp(long long int time){fts = time;}

protected:
  //! The x angle at the target. (mrad)
  Float_t fata;
  //! The y angle at the target. (mrad)
  Float_t fbta;
  //! The x position at the target. (mm)
  Float_t fxta;
  //! The y position at the target. (mm)
  Float_t fyta;
  //! The beta after the target
  Float_t fbetata;
  //! The azimuthal angle of the projectile. (deg)
  Float_t fazita;
  //! The scattering angle of the projectile. (deg)
  Float_t fscatter;
  //! The total ejectile momentum
  Float_t fptot;
  //! The ejectile parallel momentum
  Float_t fppar;
  //! The ejectile transversal momentum
  Float_t fptra;
  //! The total ejectile energy
  Float_t fetot;
  //! The ejectile kineticenergy
  Float_t fekin;
  //! time stamp
  long long int fts;

  ClassDef(ZeroDeg, 1);
};

#endif
