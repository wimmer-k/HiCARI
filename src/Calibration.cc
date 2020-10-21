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
  fAddBackType = fSett->AddBackType();
  fCoincTDiff = fSett->CoincTimeDiff();

#ifdef SIMULATION
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
  ReadHiCARIPositions(fSett->HiCARIPos());
  ReadHiCARICalibration(fSett->HiCARICalibrationFile());
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
void Calibration::ReadHiCARIPositions(const char* filename){
  TEnv *positions = new TEnv(filename);
  for(int m=0;m<12;m++){
    for(int c=0;c<4;c++){
      for(int s=0;s<40;s++){
	// double theta, phi;
	// theta = positions->GetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Theta",m,c,s),0.0);
	// phi   = positions->GetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Phi",m,c,s),0.0);
	// fHiCARIpositions[m][c][s].SetMagThetaPhi(1,theta,phi);
	//changed to x,y,z, makes it easier to deal with MINOS and position resolutions
	double x,y,z;
	x = positions->GetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.X",m,c,s),0.0);
	y = positions->GetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Y",m,c,s),0.0);
	z = positions->GetValue(Form("HiCARI.Clu%d.Cry%d.Seg%d.Z",m,c,s),0.0);
	fHiCARIpositions[m][c][s].SetXYZ(x,y,z);
      }
    }
  }
}
void Calibration::ReadHiCARICalibration(const char* filename){
  //cout << "filename " << filename << endl;
  TEnv *calF = new TEnv(filename);
  for(int m=0;m<12;m++){
    for(int c=0;c<4;c++){
      fCoreGain[m][c] = calF->GetValue(Form("Core.Clu.%02d.Cry.%02d.Gain",m,c),0.0);
      fCoreOffs[m][c] = calF->GetValue(Form("Core.Clu.%02d.Cry.%02d.Offset",m,c),0.0);      
      for(int s=0;s<40;s++){
	fSegGain[m][c][s] = calF->GetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Gain",m,c,s),0.0);
	fSegOffs[m][c][s] = calF->GetValue(Form("Clu.%02d.Cry.%02d.Seg.%02d.Offset",m,c,s),0.0);
      }
      //cout << fCoreGain[m][c] << "\t" << fCoreOffs[m][c] << endl;
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
#endif

void Calibration::BuildGretinaCalc(Gretina* in, GretinaCalc* out){

  //For no add-back, we don't need to do anything beyond copying relevant parameters.
  out->Clear();
  fevent++;
  if (in->GetMult()==0){
    return;
  }
  for(int i=0; i<in->GetMult(); i++){
    if(in->GetHit(i)->GetError()>0)
      continue;
    out->AddHit(new HitCalc(in->GetHit(i)));
    fGretinaHitctr++;
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

  fGretinactr++;
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
    fGretinaHitABctr++;
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
      fGretinaHitABctr++;
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
      fGretinaHitABctr++;
    } else {
      gr->GetHitAB(0)->AddBackHitCalc(hit);
    }
  }

}

#ifdef SIMULATION

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

void Calibration::GammaTrack(GretinaCalc* gr, GretinaEvent* gt){
  ftracking->SetGretina(gr);
  ftracking->SortInClusters();    
  gt = ftracking->GetEvent();
}

#endif

void Calibration::AllGretinaHits(GretinaCalc* gr, Gretina* in){
  //monstercluster for Ragnar
  vector<HitCalc*> cluster = ExtractAllHits(in);
  if (cluster.size() > 0 )
    gr->AddHitCL(cluster,0);
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

#ifdef SIMULATION

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
void Calibration::BuildHiCARICalc(HiCARI* in, HiCARICalc* out){
  //cout <<__PRETTY_FUNCTION__ << endl;
  out->Clear();
  if(fverbose)
    in->PrintEvent();
  if(in->GetMult()==0){
    return;
  }
  vector<HiCARIHit*> hits= in->GetHits();
  for(vector<HiCARIHit*>::iterator hit = hits.begin(); hit!=hits.end(); hit++){
    Short_t clu = (*hit)->GetCluster();
    Short_t cry = (*hit)->GetCrystal();
    long long int ts = (*hit)->GetTS();
    Float_t en = (*hit)->GetEnergy() + fRand->Uniform(0,1);
    en = en*fCoreGain[clu][cry] + fCoreOffs[clu][cry];
    //Short_t seg = (*hit)->GetMaxSegNr();
    if(fverbose)
      cout << "calibrating hit clu = " << clu << ", cry " << cry << ", en " << en << ", ts " << ts << endl;
    //here calibrate segments and find largest one
    Float_t maxen = 0;
    Float_t sumen = 0;
    Short_t maxnr = -1;
    vector<Short_t> nrs;
    vector<Float_t> ens;
    for(UShort_t s=0;s<(*hit)->GetSegmentNr().size();s++){
      Short_t segnr = (*hit)->GetSegmentNr().at(s);
      Float_t segen = (*hit)->GetSegmentEn().at(s) + fRand->Uniform(0,1);
      segen = segen*fSegGain[clu][cry][segnr] + fSegOffs[clu][cry][segnr];
      if(fverbose)
	cout << "segment " << segnr << " energy " << (*hit)->GetSegmentEn().at(s) << " cal " << segen<< endl;
      if(segen>maxen){
	maxen = segen;
	maxnr = segnr;
      }
      if(segen>0){
	sumen += segen;
	ens.push_back(segen);
	nrs.push_back(segnr);
      }
    }//uncalibrated segs
    HiCARIHitCalc* newHit = new HiCARIHitCalc(clu,cry,maxnr,sumen,fHiCARIpositions[clu][cry][maxnr],en,ts);
    newHit->SetSegments(nrs,ens);
    bool isBigRIPS = false;
    if(clu==fSett->BigRIPSCluster() && cry==fSett->BigRIPSCrystal()){
      fBigRIPSHitctr++;
      isBigRIPS = true;
    }
    else if(en>0 || sumen>0)
      fHiCARIHitctr++;
    out->AddHit(newHit,isBigRIPS);
  }
  
  if(fverbose)
    out->Print();
  //Perform the addback, which fills the add-backed vector.
  if (fAddBackType == 0){
    // do nothing
  } else if (fAddBackType == 1){
    AddBackHiCARICluster(out);
  } else if (fAddBackType == 2 || fAddBackType == 3){
    AddBackHiCARIEverything(out);
  } else {
    cout << "unknown addback type: " << fAddBackType << endl;
  }
  if(out->GetMult()>0)
    fHiCARIctr++;
  if(out->HadBigRIPS())
    fBigRIPSctr++;
  if(fverbose>2 && out->HadBigRIPS())
    cout << out->GetMult() << "\t" << (out->GetHits()).size() << endl;
  fevent++;
}

void Calibration::AddBackHiCARICluster(HiCARICalc* gr){
  //All hits within a cluster 
  vector<HiCARIHitCalc*> hits= gr->GetHits();
  for(vector<HiCARIHitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    HiCARIHitCalc* hit = *iter;
    bool addbacked = false;
    for (int j=0; j<gr->GetMultAB(); j++){
      if (gr->GetHitAB(j)->GetCluster() == hit->GetCluster() && (fabs(gr->GetHitAB(j)->GetTS()-hit->GetTS()) < fCoincTDiff || fCoincTDiff<0)){
	gr->GetHitAB(j)->AddBackHiCARIHitCalc(hit);
	addbacked = true;
	break;
      }
    }
    if (!addbacked){
      gr->AddHitAB(new HiCARIHitCalc(*hit));
    }
  }
}

void Calibration::AddBackHiCARIEverything(HiCARICalc* gr){
  //All hits within HiCARI are summed.
  vector<HiCARIHitCalc*> hits = gr->GetHits();
  if(hits.size() < 1)
    return;
  for(vector<HiCARIHitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    HiCARIHitCalc* hit = *iter;
    if(iter==hits.begin()){
      gr->AddHitAB(new HiCARIHitCalc(*hit));
    }
    else if(fabs(gr->GetHitAB(0)->GetTS()-hit->GetTS()) < fCoincTDiff || fCoincTDiff<0){
      gr->GetHitAB(0)->AddBackHiCARIHitCalc(hit);
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
  fHiCARIctr = 0;
  fBigRIPSctr = 0;
  fHiCARIHitctr = 0;
  fBigRIPSHitctr = 0;
  fGretinactr = 0;
  fGretinaHitctr = 0;
  fGretinaHitABctr = 0;
#endif
}

void Calibration::PrintCtrs(){
  //cout << "event counters in" << __PRETTY_FUNCTION__ << endl;
  cout << endl << "event counters: " << endl;
  cout << "events   \t" << fevent << endl;
#ifdef SIMULATION
  cout << "fgretactr  \t" << fgretactr  << endl;
  cout << "fminiballctr  \t" << fminiballctr  << endl;
  cout << "fzerodegctr  \t" << fzerodegctr  << endl;
  cout << "fminosctr  \t" << fminosctr  << endl;
#else
  if(fHiCARIctr>0){
    cout << "HiCARI events \t" << fHiCARIctr  << endl;
    cout << "HiCARI hits   \t" << fHiCARIHitctr  << endl;
  }
  if(fBigRIPSctr>0){
    cout << "BigRIPS events\t" << fBigRIPSctr  << endl;
    cout << "BigRIPS hits  \t" << fBigRIPSHitctr  << endl;
  }
  if(fGretinactr>0){
    cout << "Gretina events \t" << fGretinactr  << endl;
    cout << "Gretina hits   \t" << fGretinaHitctr  << endl;
    cout << "Gretina AB hits\t" << fGretinaHitABctr  << endl;
  }
#endif
}
