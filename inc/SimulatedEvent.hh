#ifndef __SIMULATED_EVENT_HH
#define __SIMULATED_EVENT_HH
#include <iostream>
#include <iomanip>

#include "UnpackedEvent.hh"
#include "GammaSim.hh"
#include "ZeroDeg.hh"
#include "Simdefs.h"


/*!
  A class for unpacking simulated data
*/
class SimulatedEvent : public UnpackedEvent {
public:
  SimulatedEvent(Settings* settings = NULL) : UnpackedEvent{settings}{};
  
  void Init();

  //! Read the specified buffer and make the GretinaG4Sim event.                                                                           
  int DecodeGammaG4Sim(G4SIM_EGS*, long long int);
  //! Read the specified buffer and make the ZeroDegPhysicsData event.                                                                     
  int DecodeZeroDegPhysicsData(ZD_PHYSICSDATA*, long long int);
  //! Passed a HiCARI Crystal, make a new crystal in the Miniball object.                                                                
  int DecodeHiCARI(HiCARIHit* hit, long long int gts);

  void SetWrite(bool wtree, bool whist, bool wctree, bool wchist, bool wstree){
    UnpackedEvent::SetWrite(wtree,whist,wctree,wchist);
    fwsimtree = wstree;
  }

  
  //! The tree containing the simulated data                                                                                               
  TTree* GetSimTree(){return fsimtr;}
  int NrOfSimEvents(){return fnsimentries;}

  
private:
  TTree *fsimtr;
  
  GammaSim*     fGammaSim;

  bool fwsimtree;
  int fnsimentries;



};
#endif
