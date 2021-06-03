#ifndef __RECONSTRUCTION_HH
#define __RECONSTRUCTION_HH
#include <iostream>
#include <iomanip>
#include "Calibration.hh"
#include "PPAC.hh"
#include "Beam.hh"
/*!
  A class for reconstruction of the data
*/
class Reconstruction : public Calibration {
public:
  //! default constructor
  Reconstruction(){};
  //! constructor
  Reconstruction(Settings* settings);
  //! dummy destructor
  ~Reconstruction(){
  };
  //! acces the settings
  Settings* GetSettings(){return fSett;}
  /////! sort by energy highest first
  ///vector<DALIHit*> Sort(vector<DALIHit*> dali);
  /////! sort by energy lowest first
  ///vector<DALIHit*> Revert(vector<DALIHit*> dali);

  // ! Align the PPAC3 position
  void AlignPPAC(SinglePPAC* pin0, SinglePPAC* pin1);
  
  // target position and scattering angle
  //! calculate the PPAC position
  TVector3 PPACPosition(SinglePPAC* pina, SinglePPAC* pinb);
  //! calculate the target position
  TVector3 TargetPosition(TVector3 inc, TVector3 ppac);
  
  // gates 
  //!  gate on the F5X position
  bool F5XGate(double f5xpos);
  //!  gate on changing charge state in ZeroDegree
  bool ChargeChange(Beam *beam);
  //!  gate on changing charge state in BigRIPS
  bool ChargeChangeBR(double delta0, double delta1);
  //!  gate on changing charge state in ZeroDegree
  bool ChargeChangeZD(double delta2, double delta3);

  //private:
  
};
#endif
