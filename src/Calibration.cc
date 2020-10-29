#include "Calibration.hh"

using namespace std;

#define uint unsigned int

Calibration::Calibration(){
  ResetCtrs();
}

Calibration::Calibration(Settings* setting, int event){
  ResetCtrs();
  fSett = setting;
  fevent = event;


  fverbose = fSett->VLevel();
  fAddBackType = fSett->AddBackType();
  fCoincTDiff = fSett->CoincTimeDiff();

  fRand = new TRandom();
  ReadHiCARIPositions(fSett->HiCARIPos());
  ReadHiCARICalibration(fSett->HiCARICalibrationFile());
  ReadMatrix(fSett->MatrixFile());
  
}

Calibration::~Calibration(){}

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
void Calibration::ReadMatrix(const char* filename){
  if(fSett->VLevel()>0)
    cout <<__PRETTY_FUNCTION__  << " " << filename << endl;
  
  ifstream infile;
  infile.open(filename);
  if(!infile.is_open()){
    cout << "no matrix found" << endl;
    return;
  }
  int hole,cry;
  while(!infile.eof()){
    //int hole,cry;
    infile >> hole >> cry;
    infile.ignore(100,'\n');
    for(int i=0;i<4;i++){
      for(int j=0;j<4;j++){
	infile >> fcrmat[hole][cry][i][j];
      }
      infile.ignore(100,'\n');
    }
    if(infile.eof())
      break;
  }
  if(fverbose>2){
    for(int hole=0;hole<MAXDETPOS;hole++){
      for(int cry=0;cry<MAXCRYSTALNO;cry++){
	if(fcrmat[hole][cry][3][3] > 0){
	  for(int i=0;i<4;i++){
	    for(int j=0;j<4;j++){
	      cout << fcrmat[hole][cry][i][j] << "\t";
	    }
	    cout << endl;
	  }
	}
      }//crystals
    }//holes
  }//verbose
}


void Calibration::BuildGretinaCalc(Gretina* in, GretinaCalc* out){

  //For no add-back, we don't need to do anything beyond copying relevant parameters, and calibrate the interaction points
  out->Clear();
  fevent++;
  if(in->GetMult()==0){
    return;
  }
  for(int i=0; i<in->GetMult(); i++){
    if(in->GetHit(i)->GetError()>0)
      continue;
    CalibrateIPoints(in->GetHit(i));
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
  
  out->DopplerCorrect(fSett);
  fGretinactr++;
}

TVector3 Calibration::TransformCoordinates(int hole, int cry, TVector3 local){
  /* Need to convert from mm to cm for this to actually work properly. NOT for Aoi-san's file!*/
  double x = local.X();///10;
  double y = local.Y();///10;
  double z = local.Z();///10;
  double xt = fcrmat[hole][cry][0][0] * x + fcrmat[hole][cry][0][1] * y + fcrmat[hole][cry][0][2] * z + fcrmat[hole][cry][0][3];
  double yt = fcrmat[hole][cry][1][0] * x + fcrmat[hole][cry][1][1] * y + fcrmat[hole][cry][1][2] * z + fcrmat[hole][cry][1][3];
  double zt = fcrmat[hole][cry][2][0] * x + fcrmat[hole][cry][2][1] * y + fcrmat[hole][cry][2][2] * z + fcrmat[hole][cry][2][3];
  // xt*=10.;
  // yt*=10.;
  // zt*=10.; //in mm
  return TVector3(xt,yt,zt) - fSett->TargetPos();
}

void Calibration::CalibrateIPoints(Crystal* cry){
  double sum =0;
  for(int j=0; j<cry->GetMult(); j++){
    IPoint* ipoint = cry->GetIPoint(j);
    ipoint->SetPosition(TransformCoordinates(cry->GetCluster(),
					     cry->GetCrystal(),
					     ipoint->GetPosition()));
    sum+=ipoint->GetEnergy();
  }
  double core = cry->GetEnergy();
  for(int j=0; j<cry->GetMult(); j++){
    IPoint* ipoint = cry->GetIPoint(j);
    double en=ipoint->GetEnergy();
    en*=core/sum;
    ipoint->SetEnergy(en);
  }
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
	Float_t ipsum = cry->GetIPSum();
	output.push_back(new HitCalc(crystal,crystal,en,ts,pos,ipsum,t0,chisq));
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
    
    
    HiCARIHitCalc* newHit = new HiCARIHitCalc();
    bool isBigRIPS = false;
    if(clu==fSett->BigRIPSCluster() && cry==fSett->BigRIPSCrystal()){
      fBigRIPSHitctr++;
      isBigRIPS = true;
      newHit = new HiCARIHitCalc(clu,cry,maxnr,sumen,TVector3(0,0,0),en,ts);
      out->AddHit(newHit,isBigRIPS);
    }
    else if(en>0 || sumen>0){
      newHit = new HiCARIHitCalc(clu,cry,maxnr,sumen,fHiCARIpositions[clu][cry][maxnr],en,ts);
      if(fverbose)
	cout << "x = "<< fHiCARIpositions[clu][cry][maxnr].X() << ", y = "<< fHiCARIpositions[clu][cry][maxnr].Y()<< ", z = "<< fHiCARIpositions[clu][cry][maxnr].Z() << endl;
      newHit->SetSegments(nrs,ens);

      if(fSett->ExcludeTracking() && newHit->IsTracking())
	continue;
      
      fHiCARIHitctr++;
      out->AddHit(newHit,isBigRIPS);
      if(fverbose)
	cout << "hit x = "<< newHit->GetPosition().X() << ", y = "<< newHit->GetPosition().Y()<< ", z = "<< newHit->GetPosition().Z() << endl;
    }
  }//hit iterator
  
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

  out->DopplerCorrect(fSett);
  
  fevent++;
}

void Calibration::AddBackHiCARICluster(HiCARICalc* gr){
  //All hits within a cluster 
  vector<HiCARIHitCalc*> hits= gr->GetHits();
  for(vector<HiCARIHitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    HiCARIHitCalc* hit = *iter;
    if(hit->IsBigRIPS())
      continue;
	      
    bool addbacked = false;
    for(int j=0; j<gr->GetMultAB(); j++){
      //if (gr->GetHitAB(j)->GetCluster() == hit->GetCluster() && (fabs(gr->GetHitAB(j)->GetTS()-hit->GetTS()) < fCoincTDiff || fCoincTDiff<0)){
      if(gr->GetHitAB(j)->GetCluster() == hit->GetCluster()){
	gr->GetHitAB(j)->AddBackHiCARIHitCalc(hit);
	addbacked = true;
	break;
      }
    }
    if(!addbacked){
      gr->AddHitAB(new HiCARIHitCalc(*hit));
      fHiCARIHitABctr++;
    }
  }
  //gr->Print();
}

void Calibration::AddBackHiCARIEverything(HiCARICalc* gr){
  //All hits within HiCARI are summed.
  vector<HiCARIHitCalc*> hits = gr->GetHits();
  if(hits.size() < 1)
    return;
  for(vector<HiCARIHitCalc*>::iterator iter = hits.begin(); iter!=hits.end(); iter++){
    HiCARIHitCalc* hit = *iter;
    if(hit->IsBigRIPS())
      continue;
    if(iter==hits.begin()){
      gr->AddHitAB(new HiCARIHitCalc(*hit));
    }
    else if(fabs(gr->GetHitAB(0)->GetTS()-hit->GetTS()) < fCoincTDiff || fCoincTDiff<0){
      gr->GetHitAB(0)->AddBackHiCARIHitCalc(hit);
    }
  }

}
void Calibration::ResetCtrs(){
  fHiCARIctr = 0;
  fBigRIPSctr = 0;
  fHiCARIHitctr = 0;
  fHiCARIHitABctr = 0;
  fBigRIPSHitctr = 0;
  fGretinactr = 0;
  fGretinaHitctr = 0;
  fGretinaHitABctr = 0;
}

void Calibration::PrintCtrs(){
  //cout << "event counters in" << __PRETTY_FUNCTION__ << endl;
  cout << endl << "event counters: " << endl;
  cout << "events   \t" << fevent << endl;
  if(fHiCARIctr>0){
    cout << "HiCARI events \t" << fHiCARIctr  << endl;
    cout << "HiCARI hits   \t" << fHiCARIHitctr  << endl;
    cout << "HiCARI AB hits\t" << fHiCARIHitABctr  << endl;
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
}
