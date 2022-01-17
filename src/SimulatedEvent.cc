#include "SimulatedEvent.hh"

void SimulatedEvent::Init(){
  UnpackedEvent::Init();

  fGammaSim = new GammaSim;
  if(fwsimtree){
    cout << "setting up simulation tree " << endl;
    fsimtr = new TTree("simtr","Geant4 emitted gamma rays");
    fsimtr->Branch("GammaSim",&fGammaSim, 320000);
    fsimtr->BranchRef();
    if(fvl>1)
      cout << "done setting up simulation tree" << endl;
  }
  fGammaSim->Clear();    
  
}

int SimulatedEvent::DecodeHiCARI(HiCARIHit* hit, long long int gts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << gts << endl;
  if(ffirst_ts<0){
    if(fvl>1)
      cout << "UnpackedEvent: " << "setting first timestamp " << gts << endl;
    ffirst_ts = gts;
  }
  if(fvl>1){
    cout << "UnpackedEvent: " << "-----------------------------next hit: "<< fnentries<< endl;
  }
  int det = hit->GetCluster();
  int cry = hit->GetCrystal();

  if(cry<0 || cry > MBCRYST){
    cout << "UnpackedEvent: " << "invalid MB crystal number " << cry << endl;
    return 12;
  }


  if(fvl>1){
    cout << "SimulatedEvent: mult " << hit->GetMult() <<"\ten " <<  hit->GetEnergy() <<"\tts " <<  hit->GetTS() <<"\tmax seg en " << \
 hit->GetMaxSegEn()<<"\tmax seg nr" <<  hit->GetMaxSegNr()<< endl;
    cout << "SimulatedEvent: " <<"-----------------------------"<< endl;
  }

  //the events which have no good interaction points                                                                                       
  if(hit->GetMaxSegNr()<0){
    fstrangehits++;
  }
  if(fvl>0 && hit->GetMaxSegNr()<0){
    cout << "strange hit " << endl;
    hit->PrintEvent();
  }
  // now check time stamps                                                                                                                 
  long long int deltaEvent = gts - fcurrent_ts;
  if(fcurrent_ts>-1 && deltaEvent < 0 )
    cout << "SimulatedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (hicari) " << gts << " difference " <<\
 deltaEvent<< endl;

  if(fvl>1){
    cout << "SimulatedEvent: " <<fnentries<< "this ts " << gts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "SimulatedEvent: " <<fnentries<< " coincidence difference " << deltaEvent << endl;
  }
  else{
    if(fvl>1)
      cout << "SimulatedEvent: " <<fnentries << " hicari single time difference " << deltaEvent << endl;
    if(fcurrent_ts>-1){
      if(fvl>2)
        cout << "SimulatedEvent: " << "Closing event due to timestamp in HiCARI." << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }

  fHiCARI->AddHit(hit);
  if(fvl>1)
    fHiCARI->PrintEvent();
  //set the current timestamp                                                                                                              
  fcurrent_ts = gts;

  if(fvl>2){
    cout << "SimulatedEvent: " << "Miniball event found with timestamp " << gts << endl;
  }
  return 0;
}
  
int SimulatedEvent::DecodeGammaG4Sim(G4SIM_EGS* g4Sim, long long int ts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << ts << endl;
  if(ffirst_ts<0){
    if(fvl>1)
      cout << "SimulatedEvent: " << "setting first timestamp " << ts << endl;
    ffirst_ts = ts;
  }
  fGammaSim->SetTimeStamp(ts);
  for(Int_t i = 0; i < g4Sim->num; i++){
    fGammaSim->AddEmittedGamma(&g4Sim->gammas[i]);
  }

  if(fvl>2){
    cout << "Simulation: GammaSim: " << endl;
    for(Int_t i = 0; i < fGammaSim->GetMult(); i++)
      cout << "e = " << fGammaSim->GetEmittedGamma(i)->GetEnergy()
           << " (x, y, z) = ("
           << fGammaSim->GetEmittedGamma(i)->GetPos().X() << ", "
           << fGammaSim->GetEmittedGamma(i)->GetPos().Y() << ", "
           << fGammaSim->GetEmittedGamma(i)->GetPos().Z()
           << ")  (phi, theta) = ("
           << fGammaSim->GetEmittedGamma(i)->GetPhi() << ", "
           << fGammaSim->GetEmittedGamma(i)->GetTheta() << ")" << endl;

    cout << "Simulation: GammaSim event found with timestamp " << ts << endl;
   }

  // These events go into a separate tree as singles. No need to pay                                                                       
  // attention to the time window.                                                                                                         

  fsimtr->Fill();
  fnsimentries++;
  fGammaSim->Clear();

  return 0;
}

int SimulatedEvent::DecodeZeroDegPhysicsData(ZD_PHYSICSDATA* zeropd, long long int ts){
  /*
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << ts << endl;
  // now check time stamps                                                                                                                 
  long long int deltaEvent = ts - fcurrent_ts;
  if(fcurrent_ts>-1 && deltaEvent < 0 )
    cout << "SimulatedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (zerodeg) " << ts << " difference " << d\
eltaEvent<< endl;
  if(fvl>1){
    cout << "SimulatedEvent: " <<fnentries<< " this ts " << ts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "SimulatedEvent: " <<fnentries<< " coincidence difference " << deltaEvent << endl;
    if(fhasdata==true){
      if(fvl>2)
        cout << "SimulatedEvent: " << "Closing event due to two ZeroDeg entries." << endl;
      cout << "SimulatedEvent: " << " coincidence with another ZeroDeg! deltaT = " << deltaEvent << " writing and clearing last event" << e\
ndl;
      if(fvl>1&&fhasdata==false)
        cout << "SimulatedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " zerodeg is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
  } else {
    if(fvl>1)
      cout << "SimulatedEvent: " <<fnentries << " zerodeg single event difference " << deltaEvent << endl;
    if(fcurrent_ts>-1&&fhasdata==true){
      if(fvl>2)
        cout << "SimulatedEvent: " << "Closing event due to timestamp in ZeroDeg." << endl;
      if(fvl>1&&fhasdata==false)
        cout << "SimulatedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " zerodeg is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }
  //set the current timestamp                                                                                                              
  fcurrent_ts = ts;

  if(fvl>0)
    cout << "SimulatedEvent: ZD_PHYSICSDATA: ata = " << zeropd->ata
         << " bta = "  << zeropd->bta
         << " xta = "  << zeropd->xta
         << " yta = "  << zeropd->yta
         << " betata = "  << zeropd->betata
         << " ts = "  << ts << endl;

  fZeroDeg->SetTimeStamp(ts);
  fZeroDeg->SetATA(fRand->Gaus(zeropd->ata, fSett->TargetAngleResolution()));
  fZeroDeg->SetBTA(fRand->Gaus(zeropd->bta, fSett->TargetAngleResolution()));
  fZeroDeg->SetXTA(fRand->Gaus(zeropd->xta, fSett->TargetPosResolution()));
  fZeroDeg->SetYTA(fRand->Gaus(zeropd->yta, fSett->TargetPosResolution()));
  fZeroDeg->SetBetaTA(fRand->Gaus(zeropd->betata, fSett->TargetBetaResolution()));
  fhasdata = true;
  if(fvl>2){
    cout << "SimulatedEvent: ZeroDegPhysicsData event found with timestamp " << ts << endl;
  }
  */
  return 0;
}
