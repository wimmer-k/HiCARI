#ifndef __FOCALPLANE_HH
#define __FOCALPLANE_HH
#include <iostream>
#include <vector>
#include <math.h>

#include "TObject.h"
#include "FocalPlanedefs.h"
using namespace std;
extern int fpID[6];
int fpNr(int id);
int firstPPAC(int id);

/*!
  Container for the MUSIC information
*/
class MUSIC : public TObject {
public:
  //! default constructor
  MUSIC(){
    Clear();
  };
  //! Clear the music information
  void Clear(Option_t *option = ""){
    fnumhits = -1.;
    fenergy = sqrt(-1.);
    fsquaresum = sqrt(-1.);
    fADC.clear();
    fGMADC.clear();
    fch.clear();
  }
  //! Set the number of hits
  void SetNHits(int nhits){fnumhits = nhits;}
  //! Set the energy
  void SetEnergy(double energy, double squaresum){
    fenergy = energy;
    fsquaresum = squaresum;
  }
  //! Set the raw ADC values and channels
  void SetRawADC(double adc, int ch){
    fADC.push_back(adc);
    fch.push_back(ch);
  }
  //! Set the gain matched ADC values
  void SetGainMatchADC(double adc){
    fGMADC.push_back(adc);
  }
  
  //! Get the number of hits
  double GetNHits(){return fnumhits;}
  //! Get the charge
  double GetEnergy(){return fenergy;}
  //! Get the channel numbers
  vector<int> GetChan(){return fch;}
  //! Get the raw ADC values
  vector<int> GetADC(){return fADC;}
  //! Get the raw ADC value of a specific channel
  int GetADC(int ch){
    for(unsigned short c=0; c<fch.size(); c++){
      if(fch.at(c)==ch)
	return fADC.at(c);
    }
    return sqrt(-1);
  }
  //! Get the gain matched ADC values
  vector<double> GetGainMatchADC(){return fGMADC;}
  //! Get the raw ADC value of a specific channel
  double GetGainMatchADC(int ch){
    for(unsigned short c=0; c<fch.size(); c++){
      if(fch.at(c)==ch)
	return fGMADC.at(c);
    }
    return sqrt(-1);
  }
protected:
  //! the number of hits
  int fnumhits;
  //! the energy deposit in the ion chamber (average)
  double fenergy;
  //! alternative way for energy calculation
  double fsquaresum;
  //! raw ADC values
  vector<int> fADC;
  //! gain matched ADC values
  vector<double> fGMADC;
  //! ADC channels
  vector<int> fch;
  

  /// \cond CLASSIMP
  ClassDef(MUSIC,1);
  /// \endcond
};

/*!
  Container for the plastic information
*/
class Plastic : public TObject {
public:
  //! default constructor
  Plastic(){
    Clear();
  };
  //! Clear the plastic information
  void Clear(Option_t *option = ""){
    ftime = sqrt(-1.);
    fcharge = sqrt(-1.);
    ftimeL = sqrt(-1.);
    ftimeR = sqrt(-1.);
    fchargeL = sqrt(-1.);
    fchargeR = sqrt(-1.);
  }
  //! Set the time
  void SetTime(Float_t timeL, Float_t timeR){
    ftimeL = timeL;
    ftimeR = timeR;
    if( !isnan(timeL) && !isnan(timeR))
      ftime = (timeL + timeR)/2;
  }
  //! Set the charge
  void SetCharge(double chargeL, double chargeR){
    fchargeL = chargeL;
    fchargeR = chargeR;
    fcharge = sqrt(chargeL * chargeR);
  }
  //! Get the time
  double GetTime(){return ftime;}
  //! Get the charge
  double GetCharge(){return fcharge;}
  //! Get the time left
  double GetTimeL(){return ftimeL;}
  //! Get the charge left
  double GetChargeL(){return fchargeL;}
  //! Get the time right
  double GetTimeR(){return ftimeR;}
  //! Get the charge right
  double GetChargeR(){return fchargeR;}
  
protected:
  //! timing of the plastic (TL + TR)
  Float_t ftime;
  //! charge deposit in the plastic (sqrt(QL*QR))
  Float_t fcharge;
  //! timing of left PMT
  Float_t ftimeL;
  //! timing of right PMT
  Float_t ftimeR;
  //! charge of left PMT
  Float_t fchargeL;
  //! charge of right PMT
  Float_t fchargeR;

  /// \cond CLASSIMP
  ClassDef(Plastic,1);
  /// \endcond
};

/*!
  Container for the Track information
*/
class Track : public TObject {
public:
  //! default constructor
  Track(){
    Clear();
  };
  //! Clear the track information
  void Clear(Option_t *option = ""){
    fx = sqrt(-1.);
    fy = sqrt(-1.);
    fa = sqrt(-1.);
    fb = sqrt(-1.);
  }
  //! Set the track information
  void Set(double x, double y, double a, double b){
    fx = x;
    fy = y;
    fa = a;
    fb = b;
  }
  //! Set x
  void SetX(double x){fx = x;}
  //! Set y
  void SetY(double y){fy = y;}
  //! Set a
  void SetA(double a){fa = a;}
  //! Set b
  void SetB(double b){fb = b;}

  //! Get x
  double GetX(){return fx;}
  //! Get y
  double GetY(){return fy;}
  //! Get a
  double GetA(){return fa;}
  //! Get b
  double GetB(){return fb;}

protected:
  //! horizontal position x
  double fx;
  //! vertical position y
  double fy;
  //! horizontal angle a
  double fa;
  //! vertical angle b
  double fb;

  /// \cond CLASSIMP
  ClassDef(Track,1);
  /// \endcond
};


/*!
  Container for the focal plane information of an event
*/
class FocalPlane : public TObject {
public:
  //! default constructor
  FocalPlane(){
    Clear();
  };
  //! Clear the track information
  void Clear(Option_t *option = ""){
    ftrack.Clear();
    fplastic.Clear();
    fmusic.Clear();
  }
  //! set the track
  void SetTrack(Track track){ftrack = track;}
  //! return the track
  Track* GetTrack(){return &ftrack;}
  //! set the plastic
  void SetPlastic(Plastic plastic){fplastic = plastic;}
  //! return the plastic
  Plastic* GetPlastic(){return &fplastic;}
  //! set the music
  void SetMUSIC(MUSIC music){fmusic = music;}
  //! return the music
  MUSIC* GetMUSIC(){return &fmusic;}
  
protected:
  //! Track of that focal plane
  Track ftrack;
  //! Plastic for that focal plane
  Plastic fplastic;
  //! Ionchamber for that focal plane
  MUSIC fmusic;

  /// \cond CLASSIMP
  ClassDef(FocalPlane,1);
  /// \endcond
};

#endif
