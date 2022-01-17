#include "HiCARI.hh"

HiCARIHit::HiCARIHit(){
  Clear();
}

HiCARIHit::HiCARIHit(int clu, int cry, Short_t nr, double en,  long long int ts, bool tracking){
  Clear();
  fcluster = clu;
  fcrystal = cry;
  ftracking = tracking;
  if( (!tracking&&nr==9) || (tracking&&nr==39)){
    fen = en;
    fmaxsinglecrystal = fen;
    ftimestamp = ts;
  }
  else{
    if(tracking && (nr%10)==9)
      return;
    AddSegment(nr,en);
    fen=0;
    ftimestamp = ts;
  }
}

#ifdef WITHSIM
HiCARIHit::HiCARIHit(sim_clust inbuf, long long int ts){
  Clear();
  if (inbuf.type!=(int)HICARI_TAG){
    cout << "Error, unexpected mode 2 hicari data format:" << hex << inbuf.type << dec
         << " in " << __PRETTY_FUNCTION__ << endl;
  }
  fcluster = inbuf.crystal_id / 4;
  fcrystal = inbuf.crystal_id % 4;
  fen = inbuf.tot_e;
  fmaxsinglecrystal = fen;
  //i am here                                                                                                                              
  for(int i=0; i<inbuf.num; i++){
    AddSegment(inbuf.seg[i].seg_id, inbuf.seg[i].e);
  }
  ftimestamp = ts;
}
#endif

bool HiCARIHit::InsertCore(int clu, int cry, double en,  long long int ts){
  if(fcluster != clu || fcrystal != cry){
    cout << "wrong cluster or crystal id! " << fcluster <<"!=" << clu <<"||"<< fcrystal <<"!="<< cry << endl;
    return false;
  }
  fen = en;
  fmaxsinglecrystal = fen;
  ftimestamp = ts;
  return true;
}

bool HiCARIHit::InsertSegment(int clu, int cry, Short_t nr, Float_t en){
  if(fcluster != clu || fcrystal != cry){
    cout << "wrong cluster or crystal id! " << fcluster <<"!=" << clu <<"||"<< fcrystal <<"!="<< cry << endl;
    return false;
  }
  AddSegment(nr,en);
  return true;
}


HiCARIHit::~HiCARIHit(){
  Clear();
}

void HiCARIHit::Clear(){
  fcluster = -1;
  fcrystal = -1;
  fen = sqrt(-1.0);
  fmult = 0;
  fmaxsinglecrystal = sqrt(-1);
  fmaxsegen = sqrt(-1);
  fmaxsegnr = -1;
  fsegen.clear();
  fsegnr.clear();
  ftimestamp = -1;
}

//Returns the sum of the segment energies
//Should be approximately equal to fen.
Float_t HiCARIHit::GetSegmentSum(){
  Float_t sum = 0;
  for(vector<Float_t>::iterator iter=fsegen.begin(); iter!=fsegen.end(); iter++){
    sum += (*iter);
  }
  return sum;
}


void HiCARIHit::AddBackHit(HiCARIHit* other){
  if(other->GetEnergy() > fmaxsinglecrystal){
    fcluster = other->GetCluster();
    fcrystal = other->GetCrystal();
    fmaxsinglecrystal =other->GetEnergy();
    ftimestamp = other->GetTS();
  }
  fen += other->GetEnergy();
  for(int i=0; i<other->GetMult(); i++){
    AddSegment(other->GetSegmentNr(i),other->GetSegmentEn(i));
  }
}

void HiCARIHit::AddSegment(Short_t nr, Float_t en){
  if(isnan(fmaxsegen) || en>fmaxsegen){
    fmaxsegen = en;
    fmaxsegnr = nr;
  }
  fsegen.push_back(en);
  fsegnr.push_back(nr);
  fmult++;
}

void HiCARIHit::PrintEvent(){
  cout << __PRETTY_FUNCTION__ << endl;
  cout << "cluster " << fcluster << ", crystal " << fcrystal << ", ts " << ftimestamp<< endl;
  cout << "energy " << fen << " keV, nr of segment hits " << fmult << endl;
  for(UShort_t i=0;i<fsegen.size();i++)
    cout << fsegnr[i] << ": " << fsegen[i] << " keV" << endl;
  cout << "max " << fmaxsegnr << " maxen " << fmaxsegen << endl;
}


HiCARI::HiCARI(){
  Clear();
}

void HiCARI::Clear(){
  fhitpattern = 0;
  fmult = 0;
  for(vector<HiCARIHit*>::iterator cry=fhits.begin(); cry!=fhits.end(); cry++){
    delete *cry;
  }
  fhits.clear();
}

void HiCARI::AddHit(HiCARIHit* cry){
  fhits.push_back(cry);
  fmult++;
  fhitpattern |= 1 << (int)cry->GetCluster();
}

void HiCARI::PrintEvent(){
  cout << __PRETTY_FUNCTION__ << endl;
  cout << "mult " << fmult << endl;
  cout << "hitpattern " << fhitpattern << endl;
  for(UShort_t i=0;i<fhits.size();i++)
    fhits[i]->PrintEvent();
}

HiCARIHitCalc::HiCARIHitCalc(Short_t clu, Short_t cry, Short_t maxseg, Float_t segsum, TVector3 pos, Float_t en, long long int ts){
  Clear();
  fcluster = clu;
  fcrystal = cry;
  fmaxseg = maxseg;
  fsegsum = segsum;
  fen = en;
  fDCen = sqrt(-1.0);
  fposition = pos;
  fHitsAdded = 1;
  fmaxhit = en;
  ftimestamp = ts;
  ftime = 0;
}

void HiCARIHitCalc::Clear(){
  fcluster = -1;
  fcrystal = -1;
  fmaxseg = -1;
  fsegsum = sqrt(-1);
  fen = sqrt(-1);
  fDCen = sqrt(-1.0);
  fposition.SetXYZ(0,0,0);
  fHitsAdded = 0;
  fmaxhit = sqrt(-1);
  ftimestamp = -1;
  fsegnr.clear();
  fsegen.clear();
  ftime = sqrt(-1);
}

HiCARIHitCalc::HiCARIHitCalc(HiCARIHitCalc* hit){
  Clear();
  fcluster = hit->GetCluster();
  fcrystal = hit->GetCrystal();
  fmaxseg = hit->GetMaxSegment();
  fsegsum = hit->GetSegSum();
  fen = hit->GetEnergy();
  fDCen = hit->GetDCEnergy();
  fposition = hit->GetPosition();
  ftimestamp = hit->GetTS();
  fHitsAdded = 1;
  fmaxhit = fen;
  fsegnr = hit->GetSegmentNr();
  fsegen = hit->GetSegmentEn();
  ftime = hit->GetTime();
}

void HiCARIHitCalc::SetSegments(vector<Short_t>nr, vector<Float_t> en){
  fsegnr = nr;
  fsegen = en;
}
void HiCARIHitCalc::AddSegment(Short_t nr, Float_t en){
  // int savemaxseg = fmaxseg;
  if(fmaxseg<0 || en>GetMaxSegmentEnergy())
    fmaxseg = nr;
  fsegnr.push_back(nr);
  fsegen.push_back(en);
  fsegsum += en;

  // if(fmaxseg!=savemaxseg)
  //   Print();
  
}

void HiCARIHitCalc::AddBackHiCARIHitCalc(HiCARIHitCalc* hit){
  if(!isnan(fmaxhit) && hit->GetEnergy() > fmaxhit){
    fcluster = hit->GetCluster();
    fcrystal = hit->GetCrystal();
    fmaxseg = hit->GetMaxSegment();
    fposition = hit->GetPosition();
    fmaxhit = hit->GetEnergy();
    ftimestamp = hit->GetTS();
    fsegnr = hit->GetSegmentNr();
    fsegen = hit->GetSegmentEn();
    ftime = hit->GetTime();
  }
  fen += hit->GetEnergy();
  fsegsum += hit->GetSegSum();
  fHitsAdded += hit->GetHitsAdded();
 
}


// fixed beta
double HiCARIHitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set){

  PosToTarget.SetX(PosToTarget.X() - set->TargetX());
  PosToTarget.SetY(PosToTarget.Y() - set->TargetY());
  PosToTarget.SetZ(PosToTarget.Z() - set->TargetZ());
  double CosDop = cos(PosToTarget.Theta());

  double beta;
  beta = set->TargetBeta();
  double gamma = 1/sqrt(1.0 - beta*beta);
  return gamma*(1-beta*CosDop);

}
void HiCARICalc::DopplerCorrect(Settings* set){
  for(vector<HiCARIHitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    if((*hit)->IsHiCARI() )
      (*hit)->DopplerCorrect(set);
  }
  for(vector<HiCARIHitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    if((*hit)->IsHiCARI() )
      (*hit)->DopplerCorrect(set);
  }
}

// event by event beta and position
// this is the function that is actually used right now
double HiCARIHitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, Beam* beam){
  TVector3 tp = beam->GetTargetPosition();
  
  PosToTarget.SetX(PosToTarget.X() - tp.X());
  PosToTarget.SetY(PosToTarget.Y() - tp.Y());
  PosToTarget.SetZ(PosToTarget.Z() - tp.Z());
  //double CosDop = cos(PosToTarget.Theta());
  double CosDop = cos(PosToTarget.Angle(beam->GetOutgoingDirection()));

  double beta;
  beta = beam->GetDopplerBeta();
  double gamma = 1/sqrt(1.0 - beta*beta);
  return gamma*(1-beta*CosDop);

}
void HiCARICalc::DopplerCorrect(Beam* beam){
  for(vector<HiCARIHitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    if((*hit)->IsHiCARI() )
      (*hit)->DopplerCorrect(beam);
  }
  for(vector<HiCARIHitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    if((*hit)->IsHiCARI() )
      (*hit)->DopplerCorrect(beam);
  }
}
  
void HiCARICalc::CorrectTime(long long int br_TS){
  for(vector<HiCARIHitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    if((*hit)->IsHiCARI() )
      (*hit)->CorrectTime(br_TS);
  }
  for(vector<HiCARIHitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    if((*hit)->IsHiCARI() )
      (*hit)->CorrectTime(br_TS);
  }
}
