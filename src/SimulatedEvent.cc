#include "SimulatedEvent.hh"

void SimulatedEvent::Init(){
  fnbeamentries = 0;
  fnsimentries = 0;
  UnpackedEvent::Init();
  fRand = new TRandom();

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
  ReadSimBroken(fSett->SimBrokenFile());

  fbeam = new Beam;
  fcaltr->Branch("beam",&fbeam,320000);
  fcaltr->BranchRef();

  ReadSimResolution(fSett->SimResolutionFile());
  ReadSimThresholds(fSett->SimThresholdFile());
  return;
}
void SimulatedEvent::CloseEvent(){
  //cout << __PRETTY_FUNCTION__ << endl;
  // SimBroken(fHiCARI);
  // SimThreshold(fHiCARI);
  // SimResolution(fHiCARI);
  SimSmear(fHiCARI);
  
  UnpackedEvent::CloseEvent();
  return;
  //cout << "end " << __PRETTY_FUNCTION__ << endl;
}
void SimulatedEvent::PrintSimCtrs(){
  //cout << "event counters in" << __PRETTY_FUNCTION__ << endl;                                                              
  if(fnbeamentries>0){
    cout << "Beam events \t" << fnbeamentries  << endl;
  }
  if(fnsimentries>0){
    cout << "Emitted events \t" << fnsimentries  << endl;
  }
}
int SimulatedEvent::DecodeHiCARI(HiCARIHit* hit, long long int gts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << gts << endl;
  if(ffirst_ts<0){
    if(fvl>1)
      cout << "SimulatedEvent: " << "setting first timestamp " << gts << endl;
    ffirst_ts = gts;
  }
  if(fvl>1){
    cout << "SimulatedEvent: " << "-----------------------------next hit: "<< fncalentries<< endl;
  }
  int det = hit->GetCluster();
  int cry = hit->GetCrystal();
  if(fvl>2)
    cout << " unpacking det = " << det << ", cry " << cry << endl;

  det = fSett->Clu2Det(det);
  hit->SetCluster(det);
  if(fvl>2)
    cout << " mapping det = " << det << ", cry " << cry << endl;

  
  if(cry<0 || cry > MBCRYST){
    cout << "SimulatedEvent: " << "invalid MB crystal number " << cry << endl;
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
    cout << "SimulatedEvent: " <<fncalentries<< "this ts " << gts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "SimulatedEvent: " <<fncalentries<< " coincidence difference " << deltaEvent << endl;
  }
  else{
    if(fvl>1)
      cout << "SimulatedEvent: " <<fncalentries << " hicari single time difference " << deltaEvent << endl;
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
           << fGammaSim->GetEmittedGamma(i)->GetTheta() << ")"
	   << " beta = " << fGammaSim->GetEmittedGamma(i)->GetBeta() << endl;

    cout << "Simulation: GammaSim event found with timestamp " << ts << endl;
   }

  // These events go into a separate tree as singles. No need to pay                                                                       
  // attention to the time window.                                                                                                         

  fsimtr->Fill();
  fnsimentries++;
  fGammaSim->Clear();

  return 0;
}

int SimulatedEvent::DecodeBeamData(BEAM_DATA* beam, long long int ts){
  if(fvl>1)
    cout << __PRETTY_FUNCTION__  << " time stamp " << ts << endl;
  // now check time stamps                                                                                                                 
  long long int deltaEvent = ts - fcurrent_ts;
  if(fcurrent_ts>-1 && deltaEvent < 0 )
    cout << "SimulatedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (beam) " << ts << " difference " << deltaEvent<< endl;
  if(fvl>1){
    cout << "SimulatedEvent: " <<fncalentries<< " this ts " << ts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "SimulatedEvent: " <<fncalentries<< " coincidence difference " << deltaEvent << endl;
    if(fhasdata==true){
      if(fvl>2)
        cout << "SimulatedEvent: " << "Closing event due to two Beam entries." << endl;
      cout << "SimulatedEvent: " << " coincidence with another Beam! deltaT = " << deltaEvent << " writing and clearing last event" << e\
ndl;
      if(fvl>1&&fhasdata==false)
        cout << "SimulatedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fncalentries << " beam is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
  } else {
    if(fvl>1)
      cout << "SimulatedEvent: " <<fncalentries << " beam single event difference " << deltaEvent << endl;
    if(fcurrent_ts>-1&&fhasdata==true){
      if(fvl>2)
        cout << "SimulatedEvent: " << "Closing event due to timestamp in Beam." << endl;
      if(fvl>1&&fhasdata==false)
        cout << "SimulatedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fncalentries << " beam is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }
  //set the current timestamp                                                                                                              
  fcurrent_ts = ts;

  //cout << "setting beta in " << beam->betain << ", out "  << beam->betata << endl;
  fbeam->SetRIPSBeta(1,fRand->Gaus(beam->betain, fSett->SimTargetBetaResolution()));
  fbeam->SetRIPSBeta(3,fRand->Gaus(beam->betata, fSett->SimTargetBetaResolution()));
  double x = fRand->Gaus(beam->tx,fSett->SimTargetPosResolution()); 
  double y = fRand->Gaus(beam->ty,fSett->SimTargetPosResolution()); 
  //cout << "setting target pos " << x << "\t" << y << "\t" << beam->iz << endl;
  fbeam->SetTargetPosition(TVector3(x,y,beam->tz));

  //cout << "setting incoming dir " << beam->ix << "\t" << beam->iy << "\t" << beam->iz << endl;
  TVector3 v = TVector3(beam->ix,beam->iy,beam->iz);
  fbeam->SetIncomingDirection(v);
  
  v = TVector3(beam->ox,beam->oy,beam->oz);
  fbeam->SetOutgoingDirection(v);
  


  fhasdata = true;
  if(fvl>2){
    cout << "SimulatedEvent: BeamData event found with timestamp " << ts << endl;
  }
  fnbeamentries++;
  return 0;
}

int SimulatedEvent::DecodeZeroDegPhysicsData(ZD_PHYSICSDATA* zeropd, long long int ts){
  if(fvl>0)
    cout << __PRETTY_FUNCTION__  << " time stamp " << ts << endl;
  // now check time stamps                                                                                                                 
  long long int deltaEvent = ts - fcurrent_ts;
  if(fcurrent_ts>-1 && deltaEvent < 0 )
    cout << "SimulatedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (zerodeg) " << ts << " difference " << d\
eltaEvent<< endl;
  if(fvl>1){
    cout << "SimulatedEvent: " <<fncalentries<< " this ts " << ts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "SimulatedEvent: " <<fncalentries<< " coincidence difference " << deltaEvent << endl;
    if(fhasdata==true){
      if(fvl>2)
        cout << "SimulatedEvent: " << "Closing event due to two ZeroDeg entries." << endl;
      cout << "SimulatedEvent: " << " coincidence with another ZeroDeg! deltaT = " << deltaEvent << " writing and clearing last event" << e\
ndl;
      if(fvl>1&&fhasdata==false)
        cout << "SimulatedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fncalentries << " zerodeg is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
  } else {
    if(fvl>1)
      cout << "SimulatedEvent: " <<fncalentries << " zerodeg single event difference " << deltaEvent << endl;
    if(fcurrent_ts>-1&&fhasdata==true){
      if(fvl>2)
        cout << "SimulatedEvent: " << "Closing event due to timestamp in ZeroDeg." << endl;
      if(fvl>1&&fhasdata==false)
        cout << "SimulatedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fncalentries << " zerodeg is empty " << endl;
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
  // fZeroDeg->SetATA(fRand->Gaus(zeropd->ata, fSett->TargetAngleResolution()));
  // fZeroDeg->SetBTA(fRand->Gaus(zeropd->bta, fSett->TargetAngleResolution()));
  // fZeroDeg->SetXTA(fRand->Gaus(zeropd->xta, fSett->TargetPosResolution()));
  // fZeroDeg->SetYTA(fRand->Gaus(zeropd->yta, fSett->TargetPosResolution()));
  // fZeroDeg->SetBetaTA(fRand->Gaus(zeropd->betata, fSett->TargetBetaResolution()));
  fZeroDeg->SetATA(fRand->Gaus(zeropd->ata, 0));
  fZeroDeg->SetBTA(fRand->Gaus(zeropd->bta, 0));
  fZeroDeg->SetXTA(fRand->Gaus(zeropd->xta, 0));
  fZeroDeg->SetYTA(fRand->Gaus(zeropd->yta, 0));
  fZeroDeg->SetBetaTA(fRand->Gaus(zeropd->betata, 0));
  fhasdata = true;
  if(fvl>2){
    cout << "SimulatedEvent: ZeroDegPhysicsData event found with timestamp " << ts << endl;
  }
  return 0;
}


/*!
  Read in which detectors and segments are broken
  Takes the name of the file with the list of broken detectors 
  @param filename: name of the file with the list of broken detectors 
*/
void SimulatedEvent::ReadSimBroken(const char* filename){
  TEnv* env = new TEnv(filename);
  if(fvl>1){
    cout << "reading file " << filename << endl;
    env->Print();
  }
  for(int det=0;det<MBCLUST+CLOVERS+TRACKINGS;det++){
    for(int cr=0;cr<MAXCRYSTALNO;cr++){
      fCoreBroken[det][cr] = env->GetValue(Form("Detector.%d.Crystal.%d.Core",det,cr),0);
      
      // simbroken coreBroken;
      // coreBroken.Detector = det;
      // coreBroken.Crystal  = cr;
      // coreBroken.Segment  = -1;
      // if(env->GetValue(Form("Detector.%d.Crystal.%d.Core",det,cr),0)==true)
      // 	fSimBroken.push_back(coreBroken);
      for(int s=0;s<MAXSEGS;s++){
	fSegBroken[det][cr][s] = env->GetValue(Form("Detector.%d.Crystal.%d.Seg.%d",det,cr,s),0);
 	// simbroken segBroken;
	// segBroken.Detector = det;
	// segBroken.Crystal  = cr;
	// segBroken.Segment  = s;
	// if(env->GetValue(Form("Detector.%d.Crystal.%d.Seg.%d",det,cr,s),0)==true)
	//   fSimBroken.push_back(segBroken);
      }//seg
    }//cr
  }//det
  return;
}

void SimulatedEvent::ReadSimResolution(const char* filename){
  TEnv* env = new TEnv(filename);
  if(fvl>1){
    cout << "reading file " << filename << endl;
    env->Print();
  }
  double globA = env->GetValue("Detector.All.Crystal.All.A",0.0);
  double globB = env->GetValue("Detector.All.Crystal.All.B",0.0);
  double globC = env->GetValue("Detector.All.Crystal.All.C",0.0);
  for(int det=0;det<MBCLUST+CLOVERS+TRACKINGS;det++){
    for(int cr=0;cr<MAXCRYSTALNO;cr++){
      int toRes = env->GetValue(Form("Detector.%d.Crystal.%d.SimResolution",det,cr),0);
      if(toRes==1){
	// simresolution CoreSimres;
	// CoreSimres.Detector = det;
	// CoreSimres.Crystal = cr;
	// CoreSimres.Segment  = -1;
	// CoreSimres.A = env->GetValue(Form("Detector.%d.Crystal.%d.A",det,cr),0.0);
	// CoreSimres.B = env->GetValue(Form("Detector.%d.Crystal.%d.B",det,cr),0.0);
	// CoreSimres.C = env->GetValue(Form("Detector.%d.Crystal.%d.C",det,cr),0.0);
        fCoreResoA[det][cr] = env->GetValue(Form("Detector.%d.Crystal.%d.A",det,cr),0.0);
        fCoreResoB[det][cr] = env->GetValue(Form("Detector.%d.Crystal.%d.B",det,cr),0.0);
        fCoreResoC[det][cr] = env->GetValue(Form("Detector.%d.Crystal.%d.C",det,cr),0.0);
	// if(CoreSimres.A>0)
	//   fSimResolutions.push_back(CoreSimres);
	for(int s=0;s<MAXSEGS;s++){
	  // simresolution SegSimres;
	  // SegSimres.Detector = det;
	  // SegSimres.Crystal = cr;
	  // SegSimres.Segment  = s;
	  // SegSimres.A = env->GetValue(Form("Detector.%d.Crystal.%d.A",det,cr),0.0);
	  // SegSimres.B = env->GetValue(Form("Detector.%d.Crystal.%d.B",det,cr),0.0);
	  // SegSimres.C = env->GetValue(Form("Detector.%d.Crystal.%d.C",det,cr),0.0);
	  fSegResoA[det][cr][s] = env->GetValue(Form("Detector.%d.Crystal.%d.Seg.%d.A",det,cr,s),0.0);
	  fSegResoB[det][cr][s] = env->GetValue(Form("Detector.%d.Crystal.%d.Seg.%d.B",det,cr,s),0.0);
	  fSegResoC[det][cr][s] = env->GetValue(Form("Detector.%d.Crystal.%d.Seg.%d.C",det,cr,s),0.0);
	  //cout << det << "\t" << cr << "\t" << s << "\t" << fSegResoA[det][cr][s] << endl;
	  // if(SegSimres.A>0)
	  //   fSimResolutions.push_back(SegSimres);
	}
      }
      else if (toRes==2){
	// simresolution CoreSimres;
	// CoreSimres.Detector = det;
	// CoreSimres.Crystal = cr;
	// CoreSimres.A = globA;
	// CoreSimres.B = globB;
	// CoreSimres.C = globC;
	// fSimResolutions.push_back(CoreSimres);
        fCoreResoA[det][cr] = globA;
        fCoreResoB[det][cr] = globB;
        fCoreResoC[det][cr] = globC;
	for(int s=0;s<MAXSEGS;s++){
	  // simresolution SegSimres;
	  // SegSimres.Detector = det;
	  // SegSimres.Crystal = cr;
	  // SegSimres.Segment  = s;
	  // SegSimres.A = globA;
	  // SegSimres.B = globB;
	  // SegSimres.C = globC;
	  // fSimResolutions.push_back(SegSimres);
	  fSegResoA[det][cr][s] = globA;
	  fSegResoB[det][cr][s] = globB;
	  fSegResoC[det][cr][s] = globC;
	}
      }
    }
  }
  // cout << "fSimResolutions.size() " << fSimResolutions.size() << endl;
  return;
}
void SimulatedEvent::ReadSimThresholds(const char* filename){
  TEnv* env = new TEnv(filename);
  if(fvl>1){
    cout << "reading file " << filename << endl;
    env->Print();
  }
  double globE = env->GetValue("Detector.All.Crystal.All.E",0.0);
  double globdE = env->GetValue("Detector.All.Crystal.All.dE",0.001);
  for(int det=0;det<MBCLUST+CLOVERS+TRACKINGS;det++){
    for(int cr=0;cr<MAXCRYSTALNO;cr++){
      int toThresh = env->GetValue(Form("Detector.%d.Crystal.%d.SimThreshold",det,cr),0);
      if(toThresh==1){
	// simthreshold newSimthresh;
	// newSimthresh.Detector = det;
	// newSimthresh.Crystal = cr;
	// newSimthresh.E = env->GetValue(Form("Detector.%d.Crystal.%d.E",det,cr),0.0);
	// newSimthresh.dE = env->GetValue(Form("Detector.%d.Crystal.%d.dE",det,cr),0.001);
	// if(fvl>2)
	//   cout << det << "\t" << cr << "\t" << newSimthresh.E << "\t" << newSimthresh.dE << endl;
	// if(newSimthresh.E>0)
	//   fSimThresholds.push_back(newSimthresh);
	fCoreThreshE[det][cr] = env->GetValue(Form("Detector.%d.Crystal.%d.E",det,cr),0.0);
	fCoreThreshdE[det][cr] = env->GetValue(Form("Detector.%d.Crystal.%d.dE",det,cr),0.0001);
      }
      else if(toThresh==2){
	// simthreshold newSimthresh;
	// newSimthresh.Detector = det;
	// newSimthresh.Crystal = cr;
	// newSimthresh.E = globE;
	// newSimthresh.dE = globdE;
	// fSimThresholds.push_back(newSimthresh);
	fCoreThreshE[det][cr] = globE;
	fCoreThreshdE[det][cr] = globdE;
      }
    }
  }
  // cout << "fSimThresholds.size() " << fSimThresholds.size() << endl;
  return;
}


/*! 
  apply the broken file, set energy and positions of broken detectors and segments
  

*/
// void SimulatedEvent::SimBroken(HiCARI* hi){
//   // cout << __PRETTY_FUNCTION__ << endl;
//   for(int i=0; i<hi->GetMult(); i++){
//     HiCARIHit* hit = hi->GetHit(i);
//     //Broken detectors
//     int det = hit->GetCluster();
//     int cry = hit->GetCrystal();
//     for(vector<simbroken>::iterator it = fSimBroken.begin(); it!=fSimBroken.end(); it++){
//       if(det==it->Detector && cry==it->Crystal && -1==it->Segment){
// 	hit->SetEnergy(0);
//       }
//       if(det==it->Detector && cry==it->Crystal){
// 	for(int j=0; j<hit->GetMult(); j++){
// 	  if(hit->GetSegmentNr(j)==it->Segment){
// 	    hit->SetSegmentEn(hit->GetSegmentNr(j),0);
// 	  }
// 	}
//       }
//     }
//   }
//   return;
// }

/*! 
  apply the resolutions
  

*/
// void SimulatedEvent::SimResolution(HiCARI* hi){
//   // cout << __PRETTY_FUNCTION__ << endl;
//   for(int i=0; i<hi->GetMult(); i++){
//     HiCARIHit* hit = hi->GetHit(i);
//     //Broken detectors
//     int det = hit->GetCluster();
//     int cry = hit->GetCrystal();
//     for(vector<simresolution>::iterator it = fSimResolutions.begin(); it!=fSimResolutions.end(); it++){
//       if(det==it->Detector && cry==it->Crystal && -1==it->Segment){
// 	  hit->SetEnergy(fRand->Gaus(hit->GetEnergy(),it->A*sqrt(1.0+hit->GetEnergy()*it->B) + it->C*hit->GetEnergy()));
//       }
//       if(det==it->Detector && cry==it->Crystal){
// 	for(int j=0; j<hit->GetMult(); j++){
// 	  if(hit->GetSegmentNr(j)==it->Segment){
// 	    hit->SetSegmentEn(hit->GetSegmentNr(j),fRand->Gaus(hit->GetSegmentEn(hit->GetSegmentNr(j)),it->A*sqrt(1.0+hit->GetSegmentEn(hit->GetSegmentNr(j))*it->B) + it->C*hit->GetSegmentEn(hit->GetSegmentNr(j))));
// 	  }
// 	}//segment hits
//       }//this crystal
//     }
//   }//hits
//   return;
// }

/*! 
  apply the thresholds
  

*/
// void SimulatedEvent::SimThreshold(HiCARI* hi){
//   //cout << __PRETTY_FUNCTION__ << endl;
//   for(int i=0; i<hi->GetMult(); i++){
//     HiCARIHit* hit = hi->GetHit(i);
//     //Broken detectors
//     int det = hit->GetCluster();
//     int cry = hit->GetCrystal();
//     //cout << "before " << hit->GetEnergy();
//     for(vector<simthreshold>::iterator it = fSimThresholds.begin(); it!=fSimThresholds.end(); it++){
//       //cout << it->Detector << "\t" << it->Crystal << endl;
//       if(det==it->Detector && cry==it->Crystal){
//     	if(fRand->Uniform(0,1) > 0.5*(1.0 + tanh( (hit->GetEnergy() - it->E) / it->dE ) ) ){
//     	  hit->SetEnergy(0);
//     	  for(int j=0; j<hit->GetMult(); j++){
//     	    hit->SetSegmentEn(hit->GetSegmentNr(j),0);
//     	  }
//     	}
//       }//this crystal
//     }
//     //cout << ", after " << hit->GetEnergy() << endl;
//   }//hits
//   return;
// }

/*!
  apply all experimental factors, broken detectors, thresholds, resolution

*/
void SimulatedEvent::SimSmear(HiCARI* hi){
  //cout << __PRETTY_FUNCTION__ << endl;
  for(int i=0; i<hi->GetMult(); i++){
    HiCARIHit* hit = hi->GetHit(i);
    int det = hit->GetCluster();
    int cry = hit->GetCrystal();
    double en = hit->GetEnergy();
    // Broken detectors
    if(fCoreBroken[det][cry]){
      hit->SetEnergy(0);
    }
    // Threshods
    else if(fRand->Uniform(0,1) > 0.5*(1.0 + tanh( (en - fCoreThreshE[det][cry]) / fCoreThreshdE[det][cry] ) ) ){
      hit->SetEnergy(0);
    }
    // Resolution
    else{
      hit->SetEnergy(fRand->Gaus(en,fCoreResoA[det][cry]*sqrt(1.0+en*fCoreResoB[det][cry]) + en*fCoreResoC[det][cry]));
    }
    for(int j=0; j<hit->GetMult(); j++){
      int seg = hit->GetSegmentNr(j);
      double en = hit->GetSegmentEn(j);
      double ensmear = fRand->Gaus(en,fSegResoA[det][cry][seg]*sqrt(1.0+en*fSegResoB[det][cry][seg]) + en*fSegResoC[det\
														    ][cry][seg]);
      hit->SetSegmentEn(j,fRand->Gaus(en,fSegResoA[det][cry][seg]*sqrt(1.0+en*fSegResoB[det][cry][seg]) + en*fSegResoC[det][cry][seg]));
    }
  }//hits
  return;
}
