#include "Calibration.hh"
#ifdef SIMULATION
Double_t inverted(Double_t *x, Double_t *p){
  return ( -p[1] - sqrt(p[1]*p[1] - 4*p[0]*(p[2]-x[0])) ) / (2*p[0]);
}
#endif

using namespace std;

#define uint unsigned int

Calibration::Calibration(){
  ResetCtrs();
#ifdef SIMULATION
  ftracking = new Tracking;
#endif
}

Calibration::Calibration(Settings* setting, int event){
  ResetCtrs();
  fSett = setting;
  fevent =event;


  fverbose = fSett->VLevel();

#ifdef SIMULATION
  fAddBackType = fSett->AddBackType();
  ReadMBPositions(fSett->AveMBPos());
  if(fSett->UseMINOS()){
    fMINOSZett = new TF1("MINOSZett",inverted,fSett->TargetZ() - 200, fSett->TargetZ() + 200, 3);
    for(int i=0;i<3;i++){
      fMINOSZett->SetParameter(i,fSett->MINOSBetaCoefficient(i));
      cout << "fSett->MINOSBetaCoefficient("<<i<<") = " << fSett->MINOSBetaCoefficient(i) << endl;
    }
  }
  ftracking = new Tracking((TrackSettings*)setting);
#else
  fRand = new TRandom();
  ReadGePositions(fSett->AveGePos());
  ReadGeCalibration(fSett->GermaniumCalibrationFile());
#endif
  
}

Calibration::~Calibration(){}

#ifdef SIMULATION
void Calibration::ReadMBPositions(const char* filename){
  TEnv *averagePos = new TEnv(filename);
  for(int m=0;m<MBCLUST;m++){
    for(int c=0;c<MBCRYST;c++){
      for(int s=0;s<MBSEGS;s++){
	// double theta, phi;
	// theta = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Theta",m,c,s),0.0);
	// phi   = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Phi",m,c,s),0.0);
	// fMBpositions[m][c][s].SetMagThetaPhi(1,theta,phi);
	//changed to x,y,z, makes it easier to deal with MINOS and position resolutions
	double x,y,z;
	x = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.X",m,c,s),0.0);
	y = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Y",m,c,s),0.0);
	z = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Z",m,c,s),0.0);
	fMBpositions[m][c][s].SetXYZ(x,y,z);
      }
    }
  }
  for(int m=MBCLUST;m<MBCLUST+CLOVERS;m++){
    for(int c=0;c<CLCRYST;c++){
      for(int s=0;s<CLSEGS;s++){
	// double theta, phi;
	// theta = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Theta",m,c,s),0.0);
	// phi   = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Phi",m,c,s),0.0);
	// fMBpositions[m][c][s].SetMagThetaPhi(1,theta,phi);
	//changed to x,y,z, makes it easier to deal with MINOS and position resolutions
	double x,y,z;
	x = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.X",m,c,s),0.0);
	y = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Y",m,c,s),0.0);
	z = averagePos->GetValue(Form("Miniball.Clu%d.Cry%d.Seg%d.Z",m,c,s),0.0);
	fMBpositions[m][c][s].SetXYZ(x,y,z);
      }
    }
  }
}
#else
void Calibration::ReadGePositions(const char* filename){
  TEnv *averagePos = new TEnv(filename);
  for(int m=0;m<12;m++){
    for(int c=0;c<4;c++){
      for(int s=0;s<40;s++){
	// double theta, phi;
	// theta = averagePos->GetValue(Form("Germanium.Clu%d.Cry%d.Seg%d.Theta",m,c,s),0.0);
	// phi   = averagePos->GetValue(Form("Germanium.Clu%d.Cry%d.Seg%d.Phi",m,c,s),0.0);
	// fGepositions[m][c][s].SetMagThetaPhi(1,theta,phi);
	//changed to x,y,z, makes it easier to deal with MINOS and position resolutions
	double x,y,z;
	x = averagePos->GetValue(Form("Germanium.Clu%d.Cry%d.Seg%d.X",m,c,s),0.0);
	y = averagePos->GetValue(Form("Germanium.Clu%d.Cry%d.Seg%d.Y",m,c,s),0.0);
	z = averagePos->GetValue(Form("Germanium.Clu%d.Cry%d.Seg%d.Z",m,c,s),0.0);
	fGepositions[m][c][s].SetXYZ(x,y,z);
      }
    }
  }
}
void Calibration::ReadGeCalibration(const char* filename){
  //cout << "filename " << filename << endl;
  TEnv *calF = new TEnv(filename);
  for(int m=0;m<12;m++){
    for(int c=0;c<4;c++){
      fGeCoreGain[m][c] = calF->GetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",m,c),0.0);
      fGeCoreOffs[m][c] = calF->GetValue(Form("Core.Clu.%02d.Cry.%02d.Offset",m,c),0.0);      
      for(int s=0;s<40;s++){
	fGeSegGain[m][c][s] = calF->GetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Gain",m,c,s),0.0);
	fGeSegOffs[m][c][s] = calF->GetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Offset",m,c,s),0.0);
      }
      //cout << fGeCoreGain[m][c] << "\t" << fGeCoreOffs[m][c] << endl;
    }
  }
}

#endif

#ifdef SIMULATION
/*!
  First three parameters are the three raw data structures, Miniball, Gretina
  Last three parameters are the three calibrated data structures to write to.
  Internally, calls BuildMiniballCalc(), BuildGretinaCalc()
  Additionally, performs calibrations that depend on lisaple systems, such as the doppler correction.
  @param inGret The raw Gretina object as input
  @param outGret A pointer to the GretinaCalc object to be built.
  @param zerodeg pointer to the zerodeg data
  @param minos pointer to the MINOS data
*/
void Calibration::BuildAllCalc(Gretina* inGret, GretinaCalc* outGret, Miniball* inMB, MiniballCalc* outMB, ZeroDeg* zerodeg, MINOS* minos){
  if(fverbose>2)
    cout <<__PRETTY_FUNCTION__ << endl;
  //Determine which of the components are present in this event.
  bool hasgret = inGret->GetMult()>0;
  bool hasmini = inMB->GetMult()>0;
  bool haszero = zerodeg->GetBetaTA()>0;
  bool hasmino = minos->GetBetaRE()>0;

  outGret->Clear();

  if(haszero){
    if(fverbose>2)
      cout << "has zerodeg " << endl;
    BuildZeroDeg(zerodeg);
    fzerodegctr++;
  }
  if(hasmino && fSett->UseMINOS()){
    if(fverbose>2)
      cout << "has minos " << endl;
    BuildMINOS(minos);
    fminosctr++;
  }
  if(hasgret){
    if(fverbose>2)
      cout << "has tracking " << endl;
    BuildGretinaCalc(inGret,outGret);
    fgretactr++;
    if(hasmino && haszero && fSett->UseMINOS())
      outGret->DopplerCorrect(fSett,zerodeg,minos);
    else if(haszero)
      outGret->DopplerCorrect(fSett,zerodeg);
    else
      outGret->DopplerCorrect(fSett);
  }
  if(hasmini){
    if(fverbose>2)
      cout << "has miniball " << endl;
    BuildMiniballCalc(inMB, outMB);
    fminiballctr++;
    if(hasmino && haszero && fSett->UseMINOS())
      outMB->DopplerCorrect(fSett,zerodeg,minos);
    else if(haszero)
      outMB->DopplerCorrect(fSett,zerodeg);
    else
      outMB->DopplerCorrect(fSett);
  }
  fevent++;
}

void Calibration::BuildMiniballCalc(Miniball* in, MiniballCalc* out){
  //cout <<__PRETTY_FUNCTION__ << endl;
  out->Clear();
  if(fverbose)
    in->PrintEvent();
  if(in->GetMult()==0){
    return;
  }
  vector<MBCrystal*> cr= in->GetHits();
  for(vector<MBCrystal*>::iterator hit = cr.begin(); hit!=cr.end(); hit++){
    Short_t clu = (*hit)->GetCluster();
    Short_t cry = (*hit)->GetCrystal();
    Short_t seg = (*hit)->GetMaxSegNr();
    long long int ts = (*hit)->GetTS();
    out->AddHit(new MBHitCalc(clu,cry,seg,fMBpositions[clu-FIRSTMB][cry][seg],(*hit)->GetEnergy(),ts));

  }
  if(fverbose)
    out->Print();
  //Perform the addback, which fills the add-backed vector.
  if (fAddBackType == 0){
    // do nothing
  } else if (fAddBackType == 1){
    AddBackMiniballCluster(out);
  } else if (fAddBackType == 2 || fAddBackType == 3){
    AddBackMiniballEverything(out);
  } else {
    cout << "unknown addback type: " << fAddBackType << endl;
  }
}

void Calibration::BuildGretinaCalc(Gretina* in, GretinaCalc* out){

  //For no add-back, we don't need to do anything beyond copying relevant parameters.
  out->Clear();
  if (in->GetMult()==0){
    return;
  }
  for(int i=0; i<in->GetMult(); i++){
    if(in->GetHit(i)->GetError()>0)
      continue;
    out->AddHit(new HitCalc(in->GetHit(i)));
  }
  //Perform the addback, which fills the add-backed vector.
  if (fAddBackType == 0){
    AddBackGretinaCrystal(out);
  } else if (fAddBackType == 1){
    AddBackGretinaCluster(out);
  } else if (fAddBackType == 2){
    AddBackGretinaEverything(out);
  } else if (fAddBackType == 3){
    ClusterGretina(out,in);
  } else {
    cout << "unknown addback type: " << fAddBackType << endl;
  }
  if(fSett->StoreAllIPoints())
    AllGretinaHits(out,in);

}

vector<HitCalc*> Calibration::ExtractAllHits(Gretina* in){
  vector<HitCalc*> output;
  for(int i=0; i<in->GetMult(); i++){
    Crystal* cry = in->GetHit(i);
    if(cry->GetEnergy()<fSett->OverflowThreshold()){
      Short_t cluster = cry->GetCluster();
      Short_t crystal = cry->GetCrystal();
      long long int ts = cry->GetTS();
      Float_t t0 = cry->GetT0();
      Float_t chisq =  cry->GetChiSq();
      for(UShort_t j=0;j<cry->GetIPoints().size();j++){
	Float_t en = cry->GetIPoints()[j]->GetEnergy();
	TVector3 pos = cry->GetIPoints()[j]->GetPosition();
	output.push_back(new HitCalc(crystal,crystal,en,ts,pos,t0,chisq));
      }
    }
  }
  return output;
}

void Calibration::AddBackGretinaCrystal(GretinaCalc* gr){
  vector<HitCalc*> hits= gr->GetHits();
  for(vector<HitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    gr->AddHitAB(new HitCalc(*iter));
  }
}

void Calibration::AddBackGretinaCluster(GretinaCalc* gr){
  //All hits within a cluster are summed.
  vector<HitCalc*> hits= gr->GetHits();
  for(vector<HitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    HitCalc* hit = *iter;
    bool addbacked = false;
    for (int j=0; j<gr->GetMultAB(); j++){
      if (gr->GetHitAB(j)->GetCluster() == hit->GetCluster()){
	gr->GetHitAB(j)->AddBackHitCalc(hit);
	addbacked = true;
	break;
      }
    }
    if (!addbacked){
      gr->AddHitAB(new HitCalc(*hit));
    }
  }
}

void Calibration::AddBackGretinaEverything(GretinaCalc* gr){
  //All hits within Gretina are summed.
  vector<HitCalc*> hits = gr->GetHits();
  if(hits.size() < 1)
    return;
  for(vector<HitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    HitCalc* hit = *iter;
    if(iter==hits.begin()){
      gr->AddHitAB(new HitCalc(*hit));
    } else {
      gr->GetHitAB(0)->AddBackHitCalc(hit);
    }
  }

}

void Calibration::AddBackMiniballCluster(MiniballCalc* gr){
  //All hits within a cluster 
  vector<MBHitCalc*> hits= gr->GetHits();
  for(vector<MBHitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    MBHitCalc* hit = *iter;
    bool addbacked = false;
    for (int j=0; j<gr->GetMultAB(); j++){
      if (gr->GetHitAB(j)->GetCluster() == hit->GetCluster()){
	gr->GetHitAB(j)->AddBackMBHitCalc(hit);
	addbacked = true;
	break;
      }
    }
    if (!addbacked){
      gr->AddHitAB(new MBHitCalc(*hit));
    }
  }
}

void Calibration::AddBackMiniballEverything(MiniballCalc* gr){
  //All hits within Gretina are summed.
  vector<MBHitCalc*> hits = gr->GetHits();
  if(hits.size() < 1)
    return;
  for(vector<MBHitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    MBHitCalc* hit = *iter;
    if(iter==hits.begin()){
      gr->AddHitAB(new MBHitCalc(*hit));
    } else {
      gr->GetHitAB(0)->AddBackMBHitCalc(hit);
    }
  }

}

void Calibration::AllGretinaHits(GretinaCalc* gr, Gretina* in){
  //monstercluster for Ragnar
  vector<HitCalc*> cluster = ExtractAllHits(in);
  if (cluster.size() > 0 )
    gr->AddHitCL(cluster,0);
}
void Calibration::GammaTrack(GretinaCalc* gr, GretinaEvent* gt){
  ftracking->SetGretina(gr);
  ftracking->SortInClusters();    
  gt = ftracking->GetEvent();
}

void Calibration::ClusterGretina(GretinaCalc* gr, Gretina* in){
  vector<HitCalc*> unusedHits = ExtractAllHits(in);
  int c =0; //clustercounter
  //All hits within a cone should be summed.
  //cone adapts for each hit.
  while (unusedHits.size() > 0){
     //First hit is automatically good.
    vector<HitCalc*> curCluster;
    curCluster.push_back(unusedHits.back());

    //also add-back clusters
    HitCalc* curHit = new HitCalc(*unusedHits.back());
    gr->AddHitAB(curHit);

    unusedHits.pop_back();

    int added =0;
    do{
      added =0;
      int incone = -1;
      bool found = false;
      for(UShort_t k=0;k<curCluster.size();k++){
	for(UShort_t j=0;j<unusedHits.size();j++){
	  if(curCluster[k]->GetPosition().Angle(unusedHits[j]->GetPosition())<fSett->ClusterAngle()*TMath::Pi()/180.){
	    incone = j;
	    found = true;
	    break;
	  }
	}
	if(found)
	  break;
      }
      if(found&&incone>-1){
	added++;
	curCluster.push_back(unusedHits[incone]);
	curHit->AddBackHitCalc(unusedHits[incone]);
	unusedHits.erase(unusedHits.begin() + incone);
      }
      if(unusedHits.size()==0)
	break;
    }while(added>0);
    gr->AddHitCL(curCluster,c);
    c++;
  }
}
void Calibration::BuildZeroDeg(ZeroDeg* zerodeg){
  //azita azimutal angle at target = phi
  double xsin, ysin;
  xsin = sin(zerodeg->GetATA()/1000);
  ysin = sin(zerodeg->GetBTA()/1000);
  zerodeg->SetPhi(atan2(ysin,xsin));//in rad


  //scatter polar angle at traget = theta
  zerodeg->SetTheta(asin(sqrt(xsin*xsin + ysin*ysin)));//in rad
  double beta =  zerodeg->GetBetaTA();
  double gamma = 1./sqrt(1-beta*beta);
    
  zerodeg->SetPtot(gamma*beta*AMU*fSett->EjectileMass());
  zerodeg->SetPpar(zerodeg->GetPtot()*cos(zerodeg->GetTheta()));
  zerodeg->SetPtra(zerodeg->GetPtot()*sin(zerodeg->GetTheta()));
  zerodeg->SetEtot(gamma*AMU*fSett->EjectileMass());
  zerodeg->SetEkin((gamma-1)*AMU*fSett->EjectileMass());
}
void Calibration::BuildMINOS(MINOS* minos){
  double beta = fMINOSZett->Eval(minos->GetZ());
  minos->SetBeta(beta);
}
#else
void Calibration::BuildGermaniumCalc(Germanium* in, GermaniumCalc* out){
  //cout <<__PRETTY_FUNCTION__ << endl;
  out->Clear();
  if(fverbose)
    in->PrintEvent();
  if(in->GetMult()==0){
    return;
  }
  vector<GeCrystal*> cr= in->GetHits();
  for(vector<GeCrystal*>::iterator hit = cr.begin(); hit!=cr.end(); hit++){
    Short_t clu = (*hit)->GetCluster();
    Short_t cry = (*hit)->GetCrystal();
    Short_t seg = (*hit)->GetMaxSegNr();
    long long int ts = (*hit)->GetTS();
    double en = (*hit)->GetEnergy() + fRand->Uniform(0,1);
    en = en*fGeCoreGain[clu][cry] + fGeCoreOffs[clu][cry];
    out->AddHit(new GeHitCalc(clu,cry,seg,fGepositions[clu][cry][seg],en,ts));

  }
  if(fverbose)
    out->Print();
  //Perform the addback, which fills the add-backed vector.
  if (fAddBackType == 0){
    // do nothing
  } else if (fAddBackType == 1){
    AddBackGermaniumCluster(out);
  } else if (fAddBackType == 2 || fAddBackType == 3){
    AddBackGermaniumEverything(out);
  } else {
    cout << "unknown addback type: " << fAddBackType << endl;
  }
}
void Calibration::AddBackGermaniumCluster(GermaniumCalc* gr){
  //All hits within a cluster 
  vector<GeHitCalc*> hits= gr->GetHits();
  for(vector<GeHitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    GeHitCalc* hit = *iter;
    bool addbacked = false;
    for (int j=0; j<gr->GetMultAB(); j++){
      if (gr->GetHitAB(j)->GetCluster() == hit->GetCluster()){
	gr->GetHitAB(j)->AddBackGeHitCalc(hit);
	addbacked = true;
	break;
      }
    }
    if (!addbacked){
      gr->AddHitAB(new GeHitCalc(*hit));
    }
  }
}

void Calibration::AddBackGermaniumEverything(GermaniumCalc* gr){
  //All hits within Gretina are summed.
  vector<GeHitCalc*> hits = gr->GetHits();
  if(hits.size() < 1)
    return;
  for(vector<GeHitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    GeHitCalc* hit = *iter;
    if(iter==hits.begin()){
      gr->AddHitAB(new GeHitCalc(*hit));
    } else {
      gr->GetHitAB(0)->AddBackGeHitCalc(hit);
    }
  }

}
#endif
void Calibration::ResetCtrs(){
#ifdef SIMULATION
  fgretactr = 0;
  fminiballctr = 0;
  fzerodegctr = 0;
  fminosctr = 0;
#else
  fgectr = 0;
#endif
}

void Calibration::PrintCtrs(){
  cout << "event counters in" << __PRETTY_FUNCTION__ << endl;
  cout << "fevent\t" << fevent << endl;
#ifdef SIMULATION
  cout << "fgretactr  \t" << fgretactr  << endl;
  cout << "fminiballctr  \t" << fminiballctr  << endl;
  cout << "fzerodegctr  \t" << fzerodegctr  << endl;
  cout << "fminosctr  \t" << fminosctr  << endl;
#else
  cout << "fgectr  \t" << fgectr  << endl;
#endif
}
