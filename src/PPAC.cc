#include "PPAC.hh"
/*!
  Calculate the ppac position as the average of A and B PPACs
  \param pina A PPAC
  \param pinb B PPAC
  \return vector to the position
*/
TVector3 PPAC::PPACPosition(SinglePPAC* pina, SinglePPAC* pinb){
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

////*!
///  Align the F8 PPAC3 with respect to the others, calibration parameters from the empty target runs
///  \param pin0 uncalibrated PPAC 0 position
///  \param pin1 uncalibrated PPAC 1 position
///*/
///void PACC::AlignPPAC(SinglePPAC* pin0, SinglePPAC* pin1){
///  double x = pin0->GetX() + fset->PPAC3PositionX0();
///  double y = pin0->GetY() + fset->PPAC3PositionY0();
///  pin0->SetXY(x,y);
///  x = pin1->GetX() + fset->PPAC3PositionX1();
///  y = pin1->GetY() + fset->PPAC3PositionY1();
///  pin1->SetXY(x,y);
///}
