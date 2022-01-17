#ifndef __PPAC_HH
#define __PPAC_HH
#include <iostream>
#include <vector>
#include <math.h>

#include "TObject.h"
#include "TVector3.h"
#include "PPACdefs.h"
using namespace std;

/*!
  Container for the information of a single PPAC
*/
class SinglePPAC : public TObject {
public:
  //! default constructor
  SinglePPAC(){
    Clear();
  };
  //! default destructor
  ~SinglePPAC(){
    Clear();
  }
  //! Clear the ppac information
  void Clear(Option_t *option = ""){
    fID = -1;
    fx = sqrt(-1.);
    fy = sqrt(-1.);
    fz = sqrt(-1.);
    ftsumx = sqrt(-1.);
    ftsumy = sqrt(-1.);
  }
  //! Set the ID
  void SetID(short id){
    if(id>0 && id <NPPACS+1) // id runs from 1
      fID = id-1; //fID runs from 0
  }
  //! Set the x position
  void SetX(double x){fx = x;}
  //! Set the y position
  void SetY(double y){fy = y;}
  //! Set the z position
  void SetZ(double z){fz = z;}
  //! Set the timing sum x
  void SetTsumX(double tsumx){ftsumx = tsumx;}
  //! Set the timing sum y
  void SetTsumY(double tsumy){ftsumy = tsumy;}
  //! Set the position
  void SetXY(double x, double y){
    fx = x;
    fy = y;
  }
  //! Set the position
  void SetXYZ(double x, double y, double z){
    fx = x;
    fy = y;
    fz = z;
  }
  //! Set the timing sum
  void SetTsumXY(double tsumx, double tsumy){
    ftsumx = tsumx;
    ftsumy = tsumy;
  }
  //! Set everything (x,y only for backwards compatability)
  void Set(short id, double x, double y, double tsumx, double tsumy){
    fID = id;
    fx = x;
    fy = y;
    ftsumx = tsumx;
    ftsumy = tsumy;
  }
  //! Set everything
  void Set(short id, double x, double y, double z, double tsumx, double tsumy){
    fID = id;
    fx = x;
    fy = y;
    fz = z;
    ftsumx = tsumx;
    ftsumy = tsumy;
  }
  //! Set everything
  void Set(short id, double x, double y, double xz, double yz, double tsumx, double tsumy){
    fID = id;
    fx = x;
    fy = y;
    fxz = xz;
    fyz = yz;
    fz = (xz+yz)/2;
    ftsumx = tsumx;
    ftsumy = tsumy;
  }

  //! Get the ID
  short GetID(){return fID;}
  //! Get the x position
  double GetX(){return fx;}
  //! Get the y position
  double GetY(){return fy;}
  //! Get the z position
  double GetZ(){return fz;}
  //! Get the z position
  double GetXZ(){return fxz;}
  //! Get the z position
  double GetYZ(){return fyz;}
  //! Get the timing sum x
  double GetTsumX(){return ftsumx;}
  //! Get the timing sum y
  double GetTsumY(){return ftsumy;}
  
  //! Has the PPAC fired
  bool Fired(){
    if( (!isnan(fx) && fx>-500) && (!isnan(fy) && fy>-500) )
      return true;
    return false;
  }
  //! Has the PPAC fired in X
  bool FiredX(){
    if(!isnan(fx) && fx>-500)
      return true;
    return false;
  }
  //! Has the PPAC fired in Y
  bool FiredY(){
    if(!isnan(fy) && fy>-500)
      return true;
    return false;
  }
protected:
  //! ID
  Short_t fID;
  //! x position
  Float_t fx;
  //! y position
  Float_t fy;
  //! z position
  Float_t fz;
  //! xz position
  Float_t fxz;
  //! yz position
  Float_t fyz;
  //! timing sum x
  Float_t ftsumx;
  //! timing sum y
  Float_t ftsumy;

  /// \cond CLASSIMP
  ClassDef(SinglePPAC,1);
  /// \endcond
};

/*!
  Container for the PPAC information
*/
class PPAC : public TObject {
public:
  //! default constructor
  PPAC(){
    Clear();
  };
  //! default destructor
  ~PPAC(){
    Clear();
  }
  //! Clear all ppac information
  void Clear(Option_t *option = ""){
    fnppacs = 0;
    for(vector<SinglePPAC*>::iterator ppac=fppacs.begin(); ppac!=fppacs.end(); ppac++){
      delete *ppac;
    }
    fppacs.clear();
  }
  //! Add a ppac
  void AddPPAC(SinglePPAC* ppac){
    fppacs.push_back(ppac);
    fnppacs++;
  }
  //! Returns the number of hit ppacs
  unsigned short GetN(){return fnppacs;}
  //! Returns the whole vector of ppacs
  vector<SinglePPAC*> GetPPACS(){return fppacs;}
  //! Returns the ppac number n in the vector
  SinglePPAC* GetPPAC(unsigned short n){return fppacs.at(n);}
  //! Returns the ppac with ID n
  SinglePPAC* GetPPACID(unsigned short n){
    for(vector<SinglePPAC*>::iterator ppac=fppacs.begin(); ppac!=fppacs.end(); ppac++){
      if((*ppac)->GetID()==n)
	return (*ppac);
    }
    return new SinglePPAC;
  }

  //! reconstruct the position for a double PPAC
  TVector3 PPACPosition(SinglePPAC* pina, SinglePPAC* pinb);
  ////! align F8 PPAC
  // needs settings
  //void AlignPPAC(SinglePPAC* pin0, SinglePPAC* pin1);


  
protected:
  //! number of ppacs hit
  unsigned short fnppacs;
  //! vector of ppacs
  vector<SinglePPAC*> fppacs;

  /// \cond CLASSIMP
  ClassDef(PPAC,1);
  /// \endcond
};
#endif
