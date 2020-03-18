#include "Germanium.hh"

GeCrystal::GeCrystal(){
  Clear();
}

GeCrystal::GeCrystal(int clu, int cry, Short_t nr, double en,  long long int ts, bool tracking){
  Clear();
  fcluster = clu;
  fcrystalid = cry;
  ftracking = tracking;
  if( (!tracking&&nr==9) || (tracking&&nr==39)){
    fen = en;
    fmaxsinglecrystal = fen;
    ftimestamp = ts;
  }
  else{
    AddSegment(nr,en);
  }
}

bool GeCrystal::InsertCore(int clu, int cry, double en,  long long int ts){
  if(fcluster != clu || fcrystalid != cry){
    cout << "wrong cluster or crystal id! " << fcluster <<"!=" << clu <<"||"<< fcrystalid <<"!="<< cry << endl;
    return false;
  }
  fen = en;
  fmaxsinglecrystal = fen;
  ftimestamp = ts;
  return true;
}

bool GeCrystal::InsertSegment(int clu, int cry, Short_t nr, Float_t en){
  if(fcluster != clu || fcrystalid != cry){
    cout << "wrong cluster or crystal id! " << fcluster <<"!=" << clu <<"||"<< fcrystalid <<"!="<< cry << endl;
    return false;
  }
  AddSegment(nr,en);
  return true;
}


GeCrystal::~GeCrystal(){
  Clear();
}

void GeCrystal::Clear(){
  fcluster = -1;
  fcrystalid = -1;
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
Float_t GeCrystal::GetSegmentSum(){
  Float_t sum = 0;
  for(vector<Float_t>::iterator iter=fsegen.begin(); iter!=fsegen.end(); iter++){
    sum += (*iter);
  }
  return sum;
}


void GeCrystal::AddBackCrystal(GeCrystal* other){
  if(other->GetEnergy() > fmaxsinglecrystal){
    fcluster = other->GetCluster();
    fcrystalid = other->GetCrystal();
    fmaxsinglecrystal =other->GetEnergy();
    ftimestamp = other->GetTS();
  }
  fen += other->GetEnergy();
  for(int i=0; i<other->GetMult(); i++){
    AddSegment(other->GetSegmentNr(i),other->GetSegmentEn(i));
  }
}

void GeCrystal::AddSegment(Short_t nr, Float_t en){
  if(isnan(fmaxsegen) || en>fmaxsegen){
    fmaxsegen = en;
    fmaxsegnr = nr;
  }
  fsegen.push_back(en);
  fsegnr.push_back(nr);
  fmult++;
}

void GeCrystal::PrintEvent(){
  cout << __PRETTY_FUNCTION__ << endl;
  cout << "cluster " << fcluster << ", crystal " << fcrystalid << ", ts " << ftimestamp<< endl;
  cout << "energy " << fen << " keV, nr of segment hits " << fmult << endl;
  for(UShort_t i=0;i<fsegen.size();i++)
    cout << fsegnr[i] << ": " << fsegen[i] << " keV" << endl;
  cout << "max " << fmaxsegnr << " maxen " << fmaxsegen << endl;
}


Germanium::Germanium(){
  Clear();
}

void Germanium::Clear(){
  fhitpattern = 0;
  fmult = 0;
  for(vector<GeCrystal*>::iterator cry=fcrystals.begin(); cry!=fcrystals.end(); cry++){
    delete *cry;
  }
  fcrystals.clear();
}

void Germanium::AddHit(GeCrystal* cry){
  fcrystals.push_back(cry);
  fmult++;
  fhitpattern |= 1 << (int)cry->GetCluster();
}

void Germanium::PrintEvent(){
  cout << __PRETTY_FUNCTION__ << endl;
  cout << "mult " << fmult << endl;
  cout << "hitpattern " << fhitpattern << endl;
  for(UShort_t i=0;i<fcrystals.size();i++)
    fcrystals[i]->PrintEvent();
}

GeHitCalc::GeHitCalc(Short_t clu, Short_t cry, Short_t seg, TVector3 pos, Float_t en, long long int ts){
  Clear();
  fcluster = clu;
  fcrystal = cry;
  fsegment = seg;
  fen = en;
  fDCen = sqrt(-1.0);
  fposition = pos;
  fHitsAdded = 1;
  fmaxhit = en;
  ftimestamp = ts;
}

void GeHitCalc::Clear(){
  fcluster = -1;
  fcrystal = -1;
  fsegment = -1;
  fen = sqrt(-1);
  fDCen = sqrt(-1.0);
  fposition.SetXYZ(0,0,0);
  fHitsAdded = 0;
  fmaxhit = sqrt(-1);
  ftimestamp = -1;
}

GeHitCalc::GeHitCalc(GeHitCalc* hit){
  Clear();
  fcluster = hit->GetCluster();
  fcrystal = hit->GetCrystal();
  fsegment = hit->GetSegment();
  fen = hit->GetEnergy();
  fDCen = hit->GetDCEnergy();
  fposition = hit->GetPosition();
  ftimestamp = hit->GetTS();
  fHitsAdded = 1;
  fmaxhit = fen;
}

void GeHitCalc::AddBackGeHitCalc(GeHitCalc* hit){
  if(!isnan(fmaxhit) && hit->GetEnergy() > fmaxhit){
    fcluster = hit->GetCluster();
    fcrystal = hit->GetCrystal();
    fsegment = hit->GetSegment();
    fposition = hit->GetPosition();
    fmaxhit = hit->GetEnergy();
    ftimestamp = hit->GetTS();
  }
  fen += hit->GetEnergy();
  fHitsAdded += hit->GetHitsAdded();
 
}

double GeHitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set){

  PosToTarget.SetX(PosToTarget.X() - set->TargetX());
  PosToTarget.SetY(PosToTarget.Y() - set->TargetY());
  PosToTarget.SetZ(PosToTarget.Z() - set->TargetZ());
  double CosDop = cos(PosToTarget.Theta());

  double beta;
  beta = set->TargetBeta();
  double gamma = 1/sqrt(1.0 - beta*beta);
  return gamma*(1-beta*CosDop);

}

void GermaniumCalc::DopplerCorrect(Settings* set){
  for(vector<GeHitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    (*hit)->DopplerCorrect(set);
  }
  for(vector<GeHitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    (*hit)->DopplerCorrect(set);
  }
}
/*  
double GeHitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set, ZeroDeg* zerodeg){

  PosToTarget.SetX(PosToTarget.X() - zerodeg->GetXTA() - set->TargetX());
  PosToTarget.SetY(PosToTarget.Y() - zerodeg->GetYTA() - set->TargetY());
  PosToTarget.SetZ(PosToTarget.Z() - set->TargetZ());
  
  TVector3 BeamDir;
  BeamDir.SetMagThetaPhi(1,zerodeg->GetTheta(),zerodeg->GetPhi());

  double CosDop = cos(PosToTarget.Angle(BeamDir));

  double beta;
  beta = set->TargetBeta() * ( 1 + (zerodeg->GetBetaTA() - set->AverageAfterBeta())/set->AverageAfterBeta());
  double gamma = 1/sqrt(1.0 - beta*beta);
  return gamma*(1-beta*CosDop);

}

void GermaniumCalc::DopplerCorrect(Settings* set, ZeroDeg* zerodeg){
  for(vector<GeHitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    (*hit)->DopplerCorrect(set,zerodeg);
  }
  for(vector<GeHitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    (*hit)->DopplerCorrect(set,zerodeg);
  }
}
  
  
double GeHitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set, ZeroDeg* zerodeg, MINOS* minos){

  PosToTarget = PosToTarget - minos->GetVertex();

  TVector3 BeamDir;
  BeamDir.SetMagThetaPhi(1,zerodeg->GetTheta(),zerodeg->GetPhi());

  double CosDop = cos(PosToTarget.Angle(BeamDir));

  double beta = minos->GetBeta();
  //cout << beta << "\t" << minos->GetBetaRE() << endl;;
  double gamma = 1/sqrt(1.0 - beta*beta);
  return gamma*(1-beta*CosDop);

}

void GermaniumCalc::DopplerCorrect(Settings* set, ZeroDeg* zerodeg, MINOS* minos){
  for(vector<GeHitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    (*hit)->DopplerCorrect(set,zerodeg,minos);
  }
  for(vector<GeHitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    (*hit)->DopplerCorrect(set,zerodeg,minos);
  }
}
*/
