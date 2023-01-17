#include "Gretina.hh"

#include <map>

IPoint::IPoint(Float_t en, Float_t x, Float_t y, Float_t z, Int_t seg,
	       Float_t segen){
  Clear();
  fen = en;
  fseg = seg;
  fseg_en = segen;
  fposition.SetXYZ(x,y,z);
}
IPoint::IPoint(IPoint* old){
  Clear();
  fen = old->GetEnergy();
  fposition = old->GetPosition();
  fseg = old->GetSeg();
  fseg_en = old->GetSegEnergy();
}
void IPoint::Clear(){
  fen = sqrt(-1.0);
  fposition.SetXYZ(0.,0.,0.);
  fseg_en = sqrt(-1.0);
  fseg = -1;
}
void IPoint::PrintEvent(){
  cout << "energy " << fen
       << " x " << fposition.X() << " y " << fposition.Y() << " z " << fposition.Z()
       << " theta " << fposition.Theta() << " phi " << fposition.Phi() << endl;
}


Crystal::Crystal(crys_ips_abcd1234 inbuf){
  Clear();
  if (inbuf.type!=(int)0xabcd1234){
    cout << "Error, unexpected mode 2 data format in " << __PRETTY_FUNCTION__ << endl;
  }

  fcluster = inbuf.crystal_id / 4;
  fcrystalid = inbuf.crystal_id % 4;
  fen = inbuf.tot_e;
  fMaxSingleCrystal = fen;
  ftimestamp = inbuf.timestamp;
  fits = inbuf.timestamp;
  ftrig_time = inbuf.trig_time;
  ft0 = (inbuf.t0==0) ? sqrt(-1.0) : inbuf.t0;

  if (inbuf.pad>0){
    SetError(inbuf.pad);
  }

  for (int i=0; i<inbuf.num; i++){
    AddIP(new IPoint(inbuf.ips[i].e,
		     inbuf.ips[i].x,inbuf.ips[i].y,inbuf.ips[i].z,
		     inbuf.ips[i].seg,inbuf.ips[i].seg_ener));
  }
}

Crystal::Crystal(crys_ips_abcd5678 inbuf){
  Clear();
  if (inbuf.type!=(int)0xabcd5678){
    cout << "Error, unexpected mode 2 data format:" << hex << inbuf.type << dec
	 << " in " << __PRETTY_FUNCTION__ << endl;
  }
  fcluster = inbuf.crystal_id / 4;
  fcrystalid = inbuf.crystal_id % 4;
  fen = inbuf.tot_e;
  fMaxSingleCrystal = fen;
  ftimestamp = inbuf.timestamp;
  fits = inbuf.timestamp;
  ftrig_time = inbuf.trig_time;
  ft0 = (inbuf.t0==0) ? sqrt(-1.0) : inbuf.t0;
  for(int i=0; i<4; i++){
    fcore_e[i] = inbuf.core_e[i];
  }
  fprestep = inbuf.prestep;
  fpoststep = inbuf.poststep;

  if (inbuf.pad>0){
    SetError(inbuf.pad);
  }

  for (int i=0; i<inbuf.num; i++){
    AddIP(new IPoint(inbuf.ips[i].e,
		     inbuf.ips[i].x,inbuf.ips[i].y,inbuf.ips[i].z,
		     inbuf.ips[i].seg,inbuf.ips[i].seg_ener));
  }

}

Crystal::Crystal(crys_ips_abcd6789 inbuf){
  Clear();
  if (inbuf.type!=(int)0xabcd6789){
    cout << "Error, unexpected mode 2 data format:" << hex << inbuf.type << dec
	 << " in " << __PRETTY_FUNCTION__ << endl;
  }
  fcluster = inbuf.crystal_id / 4;
  fcrystalid = inbuf.crystal_id % 4;
  fen = inbuf.tot_e;
  fMaxSingleCrystal = fen;
  ftimestamp = inbuf.timestamp;
  fits = inbuf.timestamp;
  //ftrig_time = inbuf.trig_time;
  ft0 = (inbuf.t0==0) ? sqrt(-1.0) : inbuf.t0;
  for(int i=0; i<4; i++){
    fcore_e[i] = inbuf.core_e[i];
  }
  fprestep = inbuf.prestep;
  fpoststep = inbuf.poststep;

  if (inbuf.pad>0){
    SetError(inbuf.pad);
  }

  for (int i=0; i<inbuf.num; i++){
    AddIP(new IPoint(inbuf.ips[i].e,
		     inbuf.ips[i].x,inbuf.ips[i].y,inbuf.ips[i].z,
		     inbuf.ips[i].seg,inbuf.ips[i].seg_ener));
  }
  //cout << inbuf.tot_e << "\t" << inbuf.tot_e_fixedPickOff_thisEvent << endl;
}

Crystal::Crystal(Crystal* old){
  Clear();
  fcluster = old->GetCluster();
  fcrystalid = old->GetCrystal();
  fen = old->GetEnergy();
  fMaxSingleCrystal = old->GetMaxSingleCrystal();
  ftimestamp = old->GetTS();
  fits = old->GetITS();
  ftrig_time = old->GetTrigTime();
  fcfd = old->GetCFD();
  fbaseline = old->GetBaseline();
  for(int i=0; i<4; i++){
    fcore_e[i] = old->GetCoreE(i);
  }
  for(int i=0; i<old->GetMult(); i++){
    AddIP(new IPoint(old->GetIPoint(i)));
  }
}

void Crystal::Clear(){
  fcluster = -1;
  fcrystalid = -1;
  fen = sqrt(-1.0);
  fMaxSingleCrystal = fen;
  ftimestamp = -1;
  fits = -1;
  ftrig_time = -1;
  fcfd = sqrt(-1.0);
  fbaseline = sqrt(-1.0);
  fmult = 0;
  ferror = 0;
  for(vector<IPoint*>::iterator ip = fipoints.begin(); ip!=fipoints.end(); ip++){
    delete *ip;
  }
  for(int i=0; i<4; i++){
    fcore_e[i] = -1;
  }
  fipoints.clear();
  fmaxip = -1;
  fmaxen = 0.0;
}

//Returns the sum of the segment energies
//Should be approximately equal to fen.
//If two interaction points are in the same segment,
// the segment energy will be the full energy for each.
//Therefore, keeping track of which segments have been recorded previously.
Float_t Crystal::GetSegmentSum(){
  Float_t sum = 0;
  std::map<int,bool> visited;
  for (vector<IPoint*>::iterator ipIter=fipoints.begin(); ipIter!=fipoints.end(); ipIter++){
    IPoint* ip = *ipIter;
    if (visited.count(ip->GetSeg())==0){
      sum += ip->GetSegEnergy();
      visited[ip->GetSeg()] = true;
    }
  }
  return sum;
}

Float_t Crystal::GetIPSum(){
  Float_t sum = 0;
  for (vector<IPoint*>::iterator ipIter=fipoints.begin(); ipIter!=fipoints.end(); ipIter++){
    sum += (*ipIter)->GetEnergy();
  }
  return sum;
}

void Crystal::AddBackCrystal(Crystal* other){
  if (other->GetEnergy() > fMaxSingleCrystal){
    fcluster = other->GetCluster();
    fcrystalid = other->GetCrystal();
    ftimestamp = other->GetTS();
    fits = other->GetITS();
    ftrig_time = other->GetTrigTime();
    fcfd = other->GetCFD();
    fbaseline = other->GetBaseline();
    ferror = other->GetError();
  }
  fen += other->GetEnergy();
  for(int i=0; i<other->GetMult(); i++){
    AddIP(new IPoint(other->GetIPoint(i)));
  }
}

void Crystal::AddIP(IPoint *ip){
  if(ip->GetEnergy()>fmaxen){
    fmaxen = ip->GetEnergy();
    fmaxip = fmult;
  }
  fipoints.push_back(ip);
  fmult++;
}

void Crystal::PrintEvent(){
  cout << "Cluster " << fcluster << " crystal " << fcrystalid << endl;
  cout << "energy " << fen << " ts " << ftimestamp << endl;
  cout << "t0 " << ft0 << " cfd " << fcfd << endl;
  cout << "nr of interaction points " << fmult << endl;
  for(UShort_t i=0;i<fipoints.size();i++)
    fipoints[i]->PrintEvent();
  cout << "max " << fmaxip << " maxen " << fmaxen << endl;
}


void Gretina::Clear(){
  fhitpattern = 0;
  fmult = 0;
  for(vector<Crystal*>::iterator cry=fcrystals.begin(); cry!=fcrystals.end(); cry++){
    delete *cry;
  }
  fcrystals.clear();
}

void Gretina::AddHit(Crystal* cry){
  fcrystals.push_back(cry);
  fmult++;
  fhitpattern |= 1 << (int)cry->GetCluster();
}

void Gretina::PrintEvent(){
  cout << "mult " << fmult << endl;
  cout << "hitpattern " << fhitpattern << endl;
  for(UShort_t i=0;i<fcrystals.size();i++)
    fcrystals[i]->PrintEvent();
}


void HitCalc::DopplerCorrect(Settings* set){
  fDCen = fen*HitCalc::DopplerCorrectionFactor(GetPosition(),set);
  fDCen_simcheat = fen*HitCalc::DopplerCorrectionFactor(fTrueFirst,set);
}

double HitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set){

  PosToTarget.SetX(PosToTarget.X() - set->TargetX());
  PosToTarget.SetY(PosToTarget.Y() - set->TargetY());
  PosToTarget.SetZ(PosToTarget.Z() - set->TargetZ());

  double CosDop = cos(PosToTarget.Theta());

  double beta;
  beta = set->TargetBeta();
  double gamma = 1/sqrt(1.0 - beta*beta);
  return gamma*(1-beta*CosDop);

}

// event by event beta and position
// this is the function that is actually used right now
double HitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, Beam* beam){
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
void GretinaCalc::DopplerCorrect(Beam* beam){
  for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    (*hit)->DopplerCorrect(beam);
  }
  for(vector<HitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    (*hit)->DopplerCorrect(beam);
  }
}

void GretinaCalc::DopplerCorrect(Settings* set){
  for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    (*hit)->DopplerCorrect(set);
  }
  for(vector<HitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    (*hit)->DopplerCorrect(set);
  }
}

void GretinaCalc::CorrectTime(long long int br_TS){
  for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    (*hit)->CorrectTime(br_TS);
  }
  for(vector<HitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    (*hit)->CorrectTime(br_TS);
  }
}
