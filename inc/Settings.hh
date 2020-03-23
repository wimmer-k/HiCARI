#ifndef __SETTINGS_HH
#define __SETTINGS_HH

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "TObject.h"
#include "TEnv.h"
#include "TVector3.h"

#ifdef SIMULATION
#include "Simdefs.h"
#endif
#include "Globaldefs.h"
using namespace std;

/*!
  A container for the analysis settings
*/
class Settings : public TObject {
public:
  //! default constructor
  Settings(){};
  //! constructor
  Settings(const char* filename);
  //! constructor
  Settings(vector<char*> files);
  //! dummy destructor
  ~Settings(){
  };
  
  void ReadSettings(TEnv* set);
  void PrintSettings();


  const string InputFile(){return fInputFiles[0];}

  void SetVLevel(int vl){fVerboseLevel = vl;}
  int VLevel(){return fVerboseLevel;}
  int EventTimeDiff(){return fEventTimeDiff;}
  Float_t TargetX(){return ftargetX;}
  Float_t TargetY(){return ftargetY;}
  Float_t TargetZ(){return ftargetZ;}
  TVector3 TargetPos(){return TVector3(ftargetX,ftargetY,ftargetZ);}
  Float_t TargetBeta(){return ftargetBeta;}
  Float_t AverageAfterBeta(){return fAveAfterBeta;}
#ifdef SIMULATION
  const char* SimResolutionFile(){return fResFile.c_str();}
  const char* SimThresholdFile(){return fThreshFile.c_str();}
  double SimGretinaPositionResolution(){return fGretPosRes;}
  double EjectileMass(){return fEjectileMass;}

  double TargetAngleResolution(){return fTargetAngleRes;}
  double TargetPosResolution(){return fTargetPosRes;}
  double TargetBetaResolution(){return fTargetBetaRes;}

  double MINOSXYResolution(){return fMINOSXYRes;}
  double MINOSZResolution(){return fMINOSZRes;}

  int UseMINOS(){return fUseMINOS;}
  double MINOSBetaCoefficient(int c){return fMINOSBetaCoeff[c];}
  
  const char* MatrixFile(){return fMatrixFile.c_str();}
  const char* NeighborFile(){return fNeighborFile.c_str();}

  int AddBackType(){return fAddBackType;}

  double ClusterAngle(){return fClusterAngle;}
  int StoreAllIPoints(){return fStoreAllIPoints;}
  void SetAddBackType(int val){fAddBackType = val;}
  double OverflowThreshold(){return fOverflowThreshold;}

  const char* AveMBPos(){return fAveMBPos.c_str();}

  int Clu2Det(int clu);
  int Det2Clu(int det);

  void SetTracking(bool tracking){fTracking = tracking;}
  bool GetTracking(){return fTracking;}
#else
  const char* HiCARIPos(){return fHiCARIPos.c_str();}
  const char* HiCARIMappingTable(){return fHiCARImapping.c_str();}
  void ReadHiCARIMappingTable();
  void PrintHiCARIMappingTable();
  int HiCARIModule(int hole, int cry, int slot){return fHiCARImap[hole][cry][slot]/10;}
  int HiCARICrystal(int hole, int cry, int slot){return fHiCARImap[hole][cry][slot]%10;}
  int RawThresh(){return fRawThresh;}
  const char* HiCARICalibrationFile(){return fHiCARIcalfile.c_str();}
  int BaselineLength(){return fBaselineLength;}
  bool TracePlots(){return fTracePlots;}

  //! Get the beta for the Doppler correction
  double Beta(){return fbeta;}
  //! Get the BigRIPS PPAC xml file
  char *PPACFile(){return (char*)fPPACfile.c_str();}
  //! Get the BigRIPS PPAC Default xml file
  char *PPACDefFile(){return (char*)fPPACdefaultfile.c_str();}
  //! Get the BigRIPS plastic xml file
  char *PlasticFile(){return (char*)fplasticfile.c_str();}
  //! Get the BigRIPS IC xml file
  char *ICFile(){return (char*)fICfile.c_str();}
  //! Get the BigRIPS Focalplane xml file
  char *FocalFile(){return (char*)ffocalfile.c_str();}
  //! Get the BigRIPS/ZeroDegree matrix file
  char *MatrixFile(int i){return (char*)fmatrixfile[i].c_str();}
  //! Get the time of flight offsets for the A/Q
  double TimeOffset(int b){return ftoffset[b];}

  //! Get the alignment shift (X0) for the PPAC3 at F8 after the target
  double PPAC3PositionX0(){return fppac3align[0];}
  //! Get the alignment shift (Y0) for the PPAC3 at F8 after the target
  double PPAC3PositionY0(){return fppac3align[1];}
  //! Get the alignment shift (X1) for the PPAC3 at F8 after the target
  double PPAC3PositionX1(){return fppac3align[2];}
  //! Get the alignment shift (Y1) for the PPAC3 at F8 after the target
  double PPAC3PositionY1(){return fppac3align[3];}
  //! Get the target position with respect to nominal focus
  double TargetPosition(){return ftargetposition;}
  //! Get the gate on the F5X position
  double F5XGate(int i){return ff5xgate[i];}
  //! Get the gate on the change in delta for charge changes
  double DeltaGate(int i){return fdeltagate[i];}

#endif

protected:
  int fEventTimeDiff;
  vector<string> fInputFiles;

  int fVerboseLevel;

  Float_t ftargetX;
  Float_t ftargetY;
  Float_t ftargetZ;
  Float_t ftargetBeta;
  Float_t fAveAfterBeta;

#ifdef SIMULATION
  double fTargetAngleRes;
  double fTargetPosRes;
  double fTargetBetaRes;

  int fUseMINOS;
  double fMINOSXYRes;
  double fMINOSZRes;
  double fMINOSBetaCoeff[3];
  
  string fResFile;
  string fThreshFile;
  double fGretPosRes;
  double fEjectileMass;

  string fAveMBPos;
  //gretina
  string fMatrixFile;
  string fNeighborFile;

  map<int,int> fdet2clu;
  map<int,int> fclu2det;

  int fAddBackType;
  double fClusterAngle;
  int fStoreAllIPoints;
  double fOverflowThreshold;
  bool fTracking;

#else
  string fHiCARIPos;
  string fHiCARImapping;
  map<int, map<int, map<int, int> > > fHiCARImap;
  int fRawThresh;
  string fHiCARIcalfile;
  int fBaselineLength;
  bool fTracePlots;

  //! BigRIPS PPAC xml file
  string fPPACfile;
  //! BigRIPS PPAC default xml file
  string fPPACdefaultfile;
  //! BigRIPS Plastic xml file
  string fplasticfile;
  //! BigRIPS IC xml file
  string fICfile;
  //! BigRIPS Focalplane xml file
  string ffocalfile;
  //! BigRIPS/ZeroDegree matrix files
  string fmatrixfile[4];
  //! time offsets for A/Q calibration
  double ftoffset[6];
  //! averge beta for Doppler correction
  double fbeta;
  //! alignment of PPAc at F8
  double fppac3align[4];
  //! target position with respect to nominal focus
  double ftargetposition;
  //! gate on the F5X position
  double ff5xgate[2];
  //! gate on the delta change
  double fdeltagate[4];
#endif 

  ClassDef(Settings, 1)
};

#endif
