#include "Settings.hh"
using namespace std;
Settings::Settings(const char* filename){
  fInputFiles.push_back(filename);
  TEnv set(filename);
  ReadSettings(&set);
  if(fVerboseLevel>1){
    PrintSettings();
    PrintMiniballMappingTable();
  }
}

Settings::Settings(vector<char*> files){
  TEnv set;
  //Reverse order, because duplicate entries are ignored, instead of overwriting.
  for(vector<char*>::reverse_iterator it = files.rbegin(); it!=files.rend(); it++){
    set.ReadFile(*it,kEnvLocal);
    fInputFiles.push_back(*it);
  }
  ReadSettings(&set);
  if(fVerboseLevel>1){
    PrintSettings();
    PrintMiniballMappingTable();
  }

}

void Settings::ReadSettings(TEnv* set){
  char* defaultfile = (char*)"~/analysis/settings/nocal.dat";

  fVerboseLevel = set->GetValue("VerboseLevel",0);
  fEventTimeDiff = set->GetValue("EventTimeDiff", 500);

  fResFile = set->GetValue("Sim.Resolution.File",defaultfile);
  fThreshFile = set->GetValue("Sim.Threshold.File",defaultfile);
  fGretPosRes = set->GetValue("Sim.Gretina.Position.Resolution",0.0);
  fEjectileMass = set->GetValue("Ejectile.Mass",0.0);

  ftargetX = set->GetValue("Target.X",0.0);
  ftargetY = set->GetValue("Target.Y",0.0);
  ftargetZ = set->GetValue("Target.Z",0.0);
  ftargetBeta = set->GetValue("Target.Beta",0.0);
  fAveAfterBeta = set->GetValue("Average.Beta.After",0.0);

  fTargetAngleRes = set->GetValue("Target.Angle.Resolution",0.0);
  fTargetPosRes   = set->GetValue("Target.Pos.Resolution",0.0);
  fTargetBetaRes  = set->GetValue("Target.Beta.Resolution",0.0);

  fMINOSXYRes = set->GetValue("MINOS.XY.Resolution",0.0);
  fMINOSZRes = set->GetValue("MINOS.Z.Resolution",0.0);


  fUseMINOS = set->GetValue("Use.MINOS",0);
  for(int i=0;i<3;i++){
    fMINOSBetaCoeff[i] = set->GetValue(Form("MINOS.Beta.Coefficient.%d",i),0.0);
  }
  
  
  fAveMBPos = set->GetValue("Average.MBPositions",defaultfile);
  
  fAddBackType = set->GetValue("AddBackType",0);
  fClusterAngle = set->GetValue("ClusterAngle",20);
  fStoreAllIPoints = set->GetValue("StoreAllIPoints",0);

  fOverflowThreshold = set->GetValue("OverflowThreshold",6000);

  fMatrixFile = set->GetValue("Gretina.Matrix.File",defaultfile);
  fNeighborFile = set->GetValue("Gretina.Neighbor.File",defaultfile);

  fdet2clu.clear();
  fclu2det.clear();
  
  for(int det=0; det<20; det++){
    int clu = set->GetValue(Form("Detector.%d",det),-1);
    if (clu!=-1){
      fdet2clu[det] = clu;
      fclu2det[clu] = det;
    }
  }

  fTracking = set->GetValue("DoTracking",false);
  fMBmapping = set->GetValue("Miniball.Mapping.Table",defaultfile);
  ReadMiniballMappingTable();
}

int Settings::Det2Clu(int det){
  if (fdet2clu.count(det)==1){
    return fdet2clu[det];
  } else {
    return -1;
  }
}

int Settings::Clu2Det(int clu){
  if (fclu2det.count(clu)==1){
    return fclu2det[clu];
  } else {
    return -1;
  }
}

void Settings::PrintSettings(){
  cout << __PRETTY_FUNCTION__ << endl;
  cout << "Sim.Resolution.File\t" << fResFile << endl;
  cout << "Sim.Threshold.File\t" << fThreshFile << endl;
  cout << "Sim.Gretina.Position.Resolution\t" << fGretPosRes << endl;
  cout << "Ejectile.Mass\t" << fEjectileMass << endl;


  cout << "Target.X\t" << ftargetX << endl;
  cout << "Target.Y\t" << ftargetY << endl;
  cout << "Target.Z\t" << ftargetZ << endl;
  cout << "Target.Beta\t" << ftargetBeta << endl;
  cout << "Average.Beta.After\t" << fAveAfterBeta << endl;


  cout << "MINOS.XY.Resolution\t" << fMINOSXYRes << endl;
  cout << "MINOS.Z.Resolution\t" << fMINOSZRes << endl;

  cout << "Use.MINOS\t" << fUseMINOS << endl;
  for(int i=0;i<3;i++)
    cout << Form("MINOS.Beta.Coefficient.%d\t",i) << fMINOSBetaCoeff[i] << endl;  

  cout << "Target.Angle.Resolution\t" << fTargetAngleRes << endl;
  cout << "Target.Pos.Resolution\t"   << fTargetPosRes << endl;
  cout << "Target.Beta.Resolution\t"   << fTargetBetaRes << endl;

  cout << "Average.MBPositions\t" << fAveMBPos << endl;

  cout << "Gretina.Matrix.File\t" << fMatrixFile << endl;
  cout << "Gretina.Neighbor.File\t" << fNeighborFile << endl;

  cout << "AddBackType\t" << fAddBackType << endl;
  cout << "ClusterAngle\t" << fClusterAngle << endl;
  cout << "StoreAllIPoints\t" << fStoreAllIPoints << endl;
  cout << "OverflowThreshold\t"<< fOverflowThreshold << endl;

  cout << "DoTracking\t"<< fTracking << endl;
  cout << "Miniball.Mapping.Table\t"<< fMBmapping << endl;
  
}


void Settings::ReadMiniballMappingTable(){
  TEnv* mapenv = new TEnv(fMBmapping.c_str());
  for(int h=0;h<20;h++){
    for(int c=0;c<4;c++){
      for(int s=0;s<4;s++){
	fMBmap[h][c][s] = mapenv->GetValue(Form("Hole.%d.Crystal.%d.Slot.%d",h,c,s),-1);
      }
    }
  }
}

void Settings::PrintMiniballMappingTable(){
  // for accessing first map 
  map<int, map<int, map<int, int> > >::iterator ftr;   
  // for accessing second map 
  map<int, map<int, int> >::iterator str; 
  // for accessing inner map 
  map<int, int>::iterator ptr;

  char* abc = "ABC";
  for(ftr = fMBmap.begin(); ftr != fMBmap.end(); ftr++) {   
    for(str = ftr->second.begin(); str != ftr->second.end(); str++) { 
      for(ptr = str->second.begin(); ptr != str->second.end(); ptr++) {
	if(ptr->second>-1)
	  cout << "Hole " << ftr->first 
	       << ", Crystal " << str->first 
	       << ", Slot " << ptr->first 
	       << " is MB " << ptr->second/10 << abc[ptr->second%10] << ", or " << fMBmap[ftr->first][str->first][ptr->first] << ", or " << MiniballModule(ftr->first,str->first,ptr->first)<< " and " << MiniballCrystal(ftr->first,str->first,ptr->first) << endl; 
	  
      } 
    } 
  }
}
