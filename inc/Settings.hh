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

#include "Simdefs.h"

using namespace std;

class Settings : public TObject {
public:
  Settings(){};//default ctor
  Settings(const char* filename);
  Settings(vector<char*> files);
  ~Settings(){};
  
  void ReadSettings(TEnv* set);
  void PrintSettings();


  const string InputFile(){return fInputFiles[0];}

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
  const char* AveGePos(){return fAveGePos.c_str();}
  const char* GermaniumMappingTable(){return fGemapping.c_str();}
  void ReadGermaniumMappingTable();
  void PrintGermaniumMappingTable();
  int GermaniumModule(int hole, int cry, int slot){return fGemap[hole][cry][slot]/10;}
  int GermaniumCrystal(int hole, int cry, int slot){return fGemap[hole][cry][slot]%10;}
  int RawThresh(){return fRawThresh;}
  const char* GermaniumCalibrationFile(){return fGecalfile.c_str();}
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
  string fAveGePos;
  string fGemapping;
  map<int, map<int, map<int, int> > > fGemap;
  int fRawThresh;
  string fGecalfile;
#endif 

  ClassDef(Settings, 1)
};

#endif
