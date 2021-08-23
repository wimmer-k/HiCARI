#include "Reconstruction.hh"
using namespace std;
/*!
  Constructor, reads in the settings for the reconstruction
  \param settings the settings file
*/
Reconstruction::Reconstruction(Settings* setting){
  //get positions for gamma
  fSett = setting;
  ReadHiCARIPositions(fSett->HiCARIPos());
  ReadMatrix(fSett->MatrixFile());
  
}

/*!
  Align the F8 PPAC3 with respect to the others, calibration parameters from the empty target runs
  \param pin0 uncalibrated PPAC 0 position
  \param pin1 uncalibrated PPAC 1 position
*/
void Reconstruction::AlignPPAC(SinglePPAC* pin0, SinglePPAC* pin1){
  double x = pin0->GetX() + fSett->PPAC3PositionX0();
  double y = pin0->GetY() + fSett->PPAC3PositionY0();
  pin0->SetXY(x,y);
  x = pin1->GetX() + fSett->PPAC3PositionX1();
  y = pin1->GetY() + fSett->PPAC3PositionY1();
  pin1->SetXY(x,y);
}

/*!
  Calculate the ppac position as the average of A and B PPACs
  \param inc incoming beam
  \param ppac position from which to extrapolate  
  \return vector to the target position
*/
TVector3 Reconstruction::PPACPosition(SinglePPAC* pina, SinglePPAC* pinb){
  double x = sqrt(-1.);
  double y = sqrt(-1.);
  double z = sqrt(-1.);
  if(pina->Fired() && pinb->Fired()){
    x = (pina->GetX()+pinb->GetX())/2;
    y = (pina->GetY()+pinb->GetY())/2;
    z = (pina->GetZ()+pinb->GetZ())/2;
  }
  else if(pina->Fired() && !pinb->Fired()){
    x = pina->GetX();
    y = pina->GetY();
    z = pina->GetZ();
  }
  else if(!pina->Fired() && pinb->Fired()){
    x = pinb->GetX();
    y = pinb->GetY();
    z = pinb->GetZ();
  }
  return TVector3(x,y,z);
}

/*!
  Calculate the ppac position as the average of A and B PPACs
  \param inc incoming beam
  \param ppac position from which to extrapolate  
  \return vector to the target position, with respect to the nominal focus
*/
TVector3 Reconstruction::TargetPosition(TVector3 inc, TVector3 ppac){
  double a = inc.X()/inc.Z();
  double b = inc.Y()/inc.Z();
  
  double x = ppac.X() + a * (fSett->TargetPosition()-ppac.Z());
  double y = ppac.Y() + b * (fSett->TargetPosition()-ppac.Z());
  
  return TVector3(x,y,fSett->TargetPosition());
}

/*!
  Calculate the event-by-event beta
  \param beam
  \return value of beta 
*/
double Reconstruction::EventBeta(Beam* beam){
#ifdef BRHOBETA
  cout << "here" << endl;
  return fSett->TargetBeta() * ( 1 + (beam->GetRIPSBeta(3) - fSett->AverageAfterBeta())/fSett->AverageAfterBeta());
#else
  return fSett->TargetBeta() * ( 1 + (beam->GetBeta(1) - fSett->AverageAfterBeta())/fSett->AverageAfterBeta());
#endif
}

/*!
  Gate on the F5X position
  \return is inside the gate
*/
bool Reconstruction::F5XGate(double f5xpos){
  if(f5xpos < fSett->F5XGate(0) || f5xpos > fSett->F5XGate(1))
      return false;
  return true;
}

/*!
  Gate on changing charge states in BigRIPS and ZeroDegree
  \return changes in charge state
*/
bool Reconstruction::ChargeChange(Beam* beam){
  double ddelta = beam->GetDelta(1) - beam->GetDelta(0);
  bool br = false;
  if(fSett->DeltaGate(0)<-998)//no gate set
    br = false;
  if(ddelta < fSett->DeltaGate(0) || ddelta > fSett->DeltaGate(1))
    br =  true;

  ddelta = beam->GetDelta(3) - beam->GetDelta(2);
  bool zd = false;
  if(fSett->DeltaGate(2)<-998)//no gate set
    zd = false;
  if(ddelta < fSett->DeltaGate(2) || ddelta > fSett->DeltaGate(3))
    zd =  true;

  if(br || zd)
    return true;
  return false;
}

/*!
  Gate on changing charge states in BigRIPS
  \return changes in charge state
*/
bool Reconstruction::ChargeChangeBR(double delta0, double delta1){
  if(fSett->DeltaGate(0)<-998)//no gate set
    return false;
  if((delta0-delta1) < fSett->DeltaGate(0) || (delta0-delta1) > fSett->DeltaGate(1))
      return true;
  return false;
}

/*!
  Gate on changing charge states in ZeroDegree
  \param delta2 deltaBrho in F8-9
  \param delta3 deltaBrho in F9-11
  \return changes in charge state
*/
bool Reconstruction::ChargeChangeZD(double delta2, double delta3){
  if(fSett->DeltaGate(2)<-998)//no gate set
    return false;
  if((delta2-delta3) < fSett->DeltaGate(2) || (delta2-delta3) > fSett->DeltaGate(3))
      return true;
  return false;
}

/*!
  Re-set the gamma positions
  \param HiCARI event
*/
void Reconstruction::SetGammaPositions(HiCARICalc* hi){
  //new HiCARI positions
  vector<HiCARIHitCalc*> hits = hi->GetHits();
  for(vector<HiCARIHitCalc*>::iterator hit=hits.begin(); hit!=hits.end(); hit++){
    if((*hit)->IsBigRIPS()){
      continue;
    }
    int cl = (*hit)->GetCluster();
    int cr = (*hit)->GetCrystal();
    int se = (*hit)->GetMaxSegment();
    //TVector3 pos = (*hit)->GetPosition();
    //TVector3 newpos = GammaPosition(cl,cr,se);
    //cout << setw(7) << setprecision(5) << pos.Theta() << "\t" << newpos.Theta() << "\t" << pos.Theta()- newpos.Theta() << endl;
    (*hit)->SetPosition(GammaPosition(cl,cr,se));
  }


  return;
}



///
////*!
///  Sort the hits by their energy, first the one with the lowest energy
///  \param hits a vector with the unsorted hits
///  \return a vector with the sorted hits
///*/
///vector<DALIHit*> Reconstruction::Revert(vector<DALIHit*> hits){
///  sort(hits.begin(), hits.end(), HitComparer());
///  reverse(hits.begin(), hits.end());
///  return hits;
///}
///
////*!
///  Sort the hits by their energy, first the one with the highest energy
///  \param hits a vector with the unsorted hits
///  \return a vector with the sorted hits
///*/
///vector<DALIHit*> Reconstruction::Sort(vector<DALIHit*> hits){
///  sort(hits.begin(), hits.end(), HitComparer());
///  return hits;
///}
///
