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

  //! Read the specified buffer and make the BeamData event.
  int DecodeBeamData(BEAM_DATA*, long long int);
  
  //! Read the specified buffer and make the ZeroDegPhysicsData event.
  int DecodeZeroDegPhysicsData(ZD_PHYSICSDATA*, long long int);
  
  //! Passed a HiCARI Crystal, make a new crystal in the Miniball object.
  int DecodeHiCARI(HiCARIHit* hit, long long int gts);

  void SetWrite(bool wtree, bool whist, bool wctree, bool wchist, bool wstree){
    UnpackedEvent::SetWrite(wtree,whist,wctree,wchist);
    fwsimtree = wstree;
  }
  void CloseEvent();
  void PrintSimCtrs();
  
  //! The tree containing the simulated data
  TTree* GetSimTree(){return fsimtr;}
  //! The number of simulated events
  int NrOfSimEvents(){return fnsimentries;}
  //! Read with detectors and segments are "broken"
  void ReadSimBroken(const char* filename);
  //! Read the resolutions to be applied to the simulated data
  void ReadSimResolution(const char* filename);
  //! Read the energy thresholds to be applied to the simulated data
  void ReadSimThresholds(const char* filename);
  /*
  //! Set the energy of broken crystals and segments to 0
  void SimBroken(HiCARI* hi);
  //! Apply energy resolution
  void SimResolution(HiCARI* hi);
  //! Apply energy threshold
  void SimThreshold(HiCARI* hi);
  */
  
  //! Apply experimental factors
  void SimSmear(HiCARI* hi);
  //! Apply experimental factors
  void SimSmear(Gretina* gr);
  
private:
  TTree *fsimtr;
  TRandom* fRand;
  
  GammaSim*  fGammaSim;
  ZeroDeg*   fZeroDeg;
  Beam*      fbeam;    
  
  bool fwsimtree;
  int fnsimentries;
  int fnbeamentries;

  // vector<simbroken> fSimBroken;
  // vector<simresolution> fSimResolutions;
  // vector<simthreshold>  fSimThresholds;

  bool fCoreBroken[MAXDETPOS][MAXCRYSTALNO];

  double fCoreThreshE[MAXDETPOS][MAXCRYSTALNO];
  double fCoreThreshdE[MAXDETPOS][MAXCRYSTALNO];

  double fCoreResoA[MAXDETPOS][MAXCRYSTALNO];
  double fCoreResoB[MAXDETPOS][MAXCRYSTALNO];
  double fCoreResoC[MAXDETPOS][MAXCRYSTALNO];

  
  bool fSegBroken[MAXDETPOS][MAXCRYSTALNO][MAXSEGS];

  double fSegResoA[MAXDETPOS][MAXCRYSTALNO][MAXSEGS];
  double fSegResoB[MAXDETPOS][MAXCRYSTALNO][MAXSEGS];
  double fSegResoC[MAXDETPOS][MAXCRYSTALNO][MAXSEGS];



};
#endif
