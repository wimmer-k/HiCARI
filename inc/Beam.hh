#ifndef __BEAM_HH
#define __BEAM_HH
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "TObject.h"
#include "TVector3.h"
#include "TRotation.h"

using namespace std;


/*!
  Container for the full beam, tof, beta and pid information
*/
class Beam : public TObject {
public:
  //! default constructor
  Beam(){
    Clear();
  };
  //! Clear all information
  void Clear(Option_t *option = ""){
    for(unsigned short j=0;j<6;j++){
      faoq[j] = sqrt(-1.);
      faoqc[j] = sqrt(-1.);
      fzet[j] = sqrt(-1.);
      fzetc[j] = sqrt(-1.);
    }    
    for(unsigned short j=0;j<3;j++){
      ftof[j] = sqrt(-1.);
      fbeta[j] = sqrt(-1.);
    }
    for(unsigned short j=0;j<4;j++){
      fdelta[j] = sqrt(-1.);
    }
    ftargetpos.SetXYZ(0,0,-99999);
    fincdir.SetXYZ(0,0,-99999);
    foutdir.SetXYZ(0,0,-99999);
    fscadir.SetXYZ(0,0,-99999);
    TVector3 X(1,0,0);
    TVector3 Y(0,1,0);
    TVector3 Z(0,0,1);
    //fincrot.RotateAxes(X,Y,Z);
  }
  //! Set the A/Q ratio
  void SetAQ(unsigned short j, double aoq){
    if(j<0 || j>5) return;
    faoq[j] = aoq;
  }
  //! Set the Z number
  void SetZ(unsigned short j, double zet){    
    if(j<0 || j>5) return;
    fzet[j] = zet;
  }
  //! Set both A/Q and Z
  void SetAQZ(unsigned short j, double aoq, double zet){
    if(j<0 || j>5) return;
    faoq[j] = aoq;
    fzet[j] = zet;
  }
  //! Set the time-of-flight
  void SetTOF(unsigned short j, double tof){    
    if(j<0 || j>2) return;
    ftof[j] = tof;
  }
  //! Set the beta
  void SetBeta(unsigned short j, double beta){    
    if(j<0 || j>2) return;
    fbeta[j] = beta;
  }
  //! Set both time-of-flight and beta
  void SetTOFBeta(unsigned short j, double tof, double beta){    
    if(j<0 || j>2) return;
    ftof[j] = tof;
    fbeta[j] = beta;
  }
  //! Set the delta
  void SetDelta(unsigned short j, double delta){    
    if(j<0 || j>3) return;
    fdelta[j] = delta;
  }

  //! Set the target position 
  void SetTargetPosition(TVector3 pos){
    ftargetpos = pos;
  }
  //! Set the target x position 
  void SetTargetXPos(double xpos){
    ftargetpos.SetX(xpos);
  }
  //! Set the target y position 
  void SetTargetYPos(double ypos){
    ftargetpos.SetY(ypos);
  }
  //! Set the direction of the incoming beam
  void SetIncomingDirection(TVector3 dir){
    fincdir = dir;
    fincrot.SetZAxis(fincdir.Unit());
    fincrot = fincrot.Inverse();
  }
  //! Set the direction of the scattered beam
  void SetOutgoingDirection(TVector3 dir){
    foutdir = dir;
    if(fincdir.Z()>-99999){
      fscadir = dir;
      fscadir.Transform(fincrot);
    }
  }

  //! Correct the A/Q ratio based on position
  void CorrectAQ(unsigned short j, double corr){
    if(j<0 || j>5) return;
    faoqc[j] = faoq[j] + corr;
  }

  //! Get the A/Q ratio
  double GetAQ(unsigned short j){
    if(j<0 || j>5) return sqrt(-1.);
    return faoq[j];
  }
  //! Get the corrected A/Q ratio
  double GetCorrAQ(unsigned short j){
    if(j<0 || j>5) return sqrt(-1.);
    return faoqc[j];
  }
  //! Get the Z number
  double GetZ(unsigned short j){
    if(j<0 || j>5) return sqrt(-1.);
    return fzet[j];
  }
  //! Get the time-of-flight
  double GetTOF(unsigned short j){
    if(j<0 || j>2) return sqrt(-1.);
    return ftof[j];
  }
  //! Get beta
  double GetBeta(unsigned short j){
    if(j<0 || j>2) return sqrt(-1.);
    return fbeta[j];
  }
  //! Get Delta
  double GetDelta(unsigned short j){
    if(j<0 || j>3) return sqrt(-1.);
    return fdelta[j];
  }
  //! Get the direction of the incoming beam in lab system
  TVector3 GetIncomingDirection(){
    return fincdir;
  }
  //! Get the direction of the outgoing beam in lab system
  TVector3 GetOutgoingDirection(){
    return foutdir;
  }
  //! Get the direction of the scattering (in beam coordinate system)
  TVector3 GetScatteredDirection(){
    return fscadir;
  }
  //! Get the target position
  TVector3 GetTargetPosition(){
    return ftargetpos;
  }
  //! Get the target X position
  double GetTargetPositionX(){
    return ftargetpos.X();
  }
  //! Get the target Y position
  double GetTargetPositionY(){
    return ftargetpos.Y();
  }

  //! Get the scattering angle phi
  double GetPhi(){
    if(fscadir.Z()>-99999)
      return fscadir.Phi();
    return sqrt(-1.);
  }
  //! Get the scattering angle theta
  double GetTheta(){
    if(fscadir.Z()>-99999)
      return fscadir.Theta();
    return sqrt(-1.);
  }

protected:
  //! A/Q for 3-5, 5-7, 3-7,  8-9, 9-11, 8-11
  double faoq[6];
  //! corrected A/Q
  double faoqc[6];
  //! Z for 3-5, 5-7, 3-7,  8-9, 9-11, 8-11  
  double fzet[6];
  //! corrected Z
  double fzetc[6];

  //! time-of-flight for 3-7, 8-11, 7-8
  double ftof[3];
  //! beta for 3-7, 8-11, 7-8
  double fbeta[3];

  //! delta momentum 3-5, 5-7, 8-9, 9-11
  double fdelta[4];

  //! target position 
  TVector3 ftargetpos;
  //! incoming direction in lab system
  TVector3 fincdir;
  //! outgoing direction in lab system
  TVector3 foutdir;
  //! rotation matrix beam coordinate system
  TRotation fincrot;
  //! ejectile vector in beam coordinate system
  TVector3 fscadir;

  /// \cond CLASSIMP
  ClassDef(Beam,1);
  /// \endcond
};


#endif
