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
  int AddBackType(){return fAddBackType;}
  int CoincTimeDiff(){return fCoincTimeDiff;}
  int StoreAllIPoints(){return fStoreAllIPoints;}

  double ClusterAngle(){return fClusterAngle;}
  void SetAddBackType(int val){fAddBackType = val;}
  double OverflowThreshold(){return fOverflowThreshold;}

  bool IgnoreTrace(){return fIgnoreTrace;}
  const char* HiCARIPos(){return fHiCARIPos.c_str();}
  const char* HiCARIMappingTable(){return fHiCARImapping.c_str();}
  void ReadHiCARIMappingTable();
  void PrintHiCARIMappingTable();
  int HiCARICluster(int hole, int cry, int slot){return fHiCARImap[hole][cry][slot]/10;}
  int HiCARICrystal(int hole, int cry, int slot){return fHiCARImap[hole][cry][slot]%10;}
  int RawThresh(){return fRawThresh;}
  const char* HiCARICalibrationFile(){return fHiCARIcalfile.c_str();}
  int BaselineLength(){return fBaselineLength;}
  bool TracePlots(){return fTracePlots;}
  bool Mode3Histos(){return fMode3Histos;}

  bool ExcludeTracking(){return fExcludeTracking;}
  const char* MatrixFile(){return fMatrixfile.c_str();}
  
  
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

  //! Set the event correlation mode
  void SetCorrelationMode(int mode){fcorrelationMode = mode;}
  //! Get the event correlation mode
  int CorrelationMode(){return fcorrelationMode;}
  //! Get the cluster used for correlation checking
  int CorrelationCluster(){return fcorrelationCluster;}
  //! Get the crystal used for correlation checking
  int CorrelationCrystal(){return fcorrelationCrystal;}

  //! Get the ADC channel for checkADC
  int CorrelationChannel(){return fcorrelationChannel;}
  
  //! Get the cluster used for BigRIPS time checking
  int BigRIPSCluster(){return fBigRIPSCluster;}
  //! Get the crystal used for BigRIPS time checking
  int BigRIPSCrystal(){return fBigRIPSCrystal;}
  //! Get the channel used for BigRIPS time checking
  int BigRIPSChannel(){return fBigRIPSChannel;}
  //! Set the detail level for the BigRIPSTree
  void SetBigRIPSDetail(int det){fBigRIPSDetail = det;}
  //! Get the detail level for the BigRIPSTree
  int BigRIPSDetail(){return fBigRIPSDetail;}

  double GetAoQCorrection(int sp, int fp, int tr){return faoq_corr[sp][fp][tr];}
  double GetBRAoQCorrection_F3X(){return faoq_corr[0][0][0];}
  double GetBRAoQCorrection_F3A(){return faoq_corr[0][0][1];}
  double GetBRAoQCorrection_F5X(){return faoq_corr[0][1][0];}
  double GetBRAoQCorrection_F5A(){return faoq_corr[0][1][1];}
  double GetBRAoQCorrection_F7X(){return faoq_corr[0][2][0];}
  double GetBRAoQCorrection_F7A(){return faoq_corr[0][2][1];}
  double GetZDAoQCorrection_F8X(){return faoq_corr[1][0][0];}
  double GetZDAoQCorrection_F8A(){return faoq_corr[1][0][1];}
  double GetZDAoQCorrection_F9X(){return faoq_corr[1][1][0];}
  double GetZDAoQCorrection_F9A(){return faoq_corr[1][1][1];}
  double GetZDAoQCorrection_F11X(){return faoq_corr[1][2][0];}
  double GetZDAoQCorrection_F11A(){return faoq_corr[1][2][1];}

  
protected:
  int fEventTimeDiff;
  vector<string> fInputFiles;

  int fVerboseLevel;

  Float_t ftargetX;
  Float_t ftargetY;
  Float_t ftargetZ;
  Float_t ftargetBeta;
  Float_t fAveAfterBeta;
  int fAddBackType;
  int fCoincTimeDiff;
  int fStoreAllIPoints;
  double fOverflowThreshold;
  string fMatrixFile;
  string fNeighborFile;

  map<int,int> fdet2clu;
  map<int,int> fclu2det;

  double fClusterAngle;
  bool fTracking;


  bool fIgnoreTrace;
  string fHiCARIPos;
  string fHiCARImapping;
  map<int, map<int, map<int, int> > > fHiCARImap;
  int fRawThresh;
  string fHiCARIcalfile;
  int fBaselineLength;
  bool fTracePlots;
  bool fMode3Histos;
  bool fExcludeTracking;
  string fMatrixfile;
  
  int fcorrelationMode;
  int fcorrelationCluster;
  int fcorrelationCrystal;
  int fcorrelationChannel;

  int fBigRIPSCluster;
  int fBigRIPSCrystal;
  int fBigRIPSChannel;
  int fBigRIPSDetail;
  
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


  double faoq_corr[2][3][2];  // BR/ZD, focal plane, x,angle

  ClassDef(Settings, 1)
};

#endif
