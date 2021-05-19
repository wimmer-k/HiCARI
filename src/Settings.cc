#include "Settings.hh"
using namespace std;
Settings::Settings(const char* filename){
  fInputFiles.push_back(filename);
  TEnv set(filename);
  ReadSettings(&set);
  if(fVerboseLevel>1){
    PrintSettings();
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
  }
}

void Settings::ReadSettings(TEnv* set){
  char* defaultfile = (char*)"~/analysis/settings/nocal.dat";

  fVerboseLevel = set->GetValue("VerboseLevel",0);
  fEventTimeDiff = set->GetValue("EventTimeDiff", 500);

  ftargetX = set->GetValue("Target.X",0.0);
  ftargetY = set->GetValue("Target.Y",0.0);
  ftargetZ = set->GetValue("Target.Z",0.0);
  ftargetBeta = set->GetValue("Target.Beta",0.0);
  fAveAfterBeta = set->GetValue("Average.Beta.After",0.0);
  fAddBackType = set->GetValue("AddBackType",0);
  fCoincTimeDiff = set->GetValue("CoincTimeDiff",-1);

  fIgnoreTrace = set->GetValue("HiCARI.IgnoreTrace",0);
  fHiCARIPos = set->GetValue("HiCARI.Positions",defaultfile);
  fHiCARImapping = set->GetValue("HiCARI.Mapping.Table",defaultfile);
  fRawThresh = set->GetValue("HiCARI.Raw.Thresh",0);
  ReadHiCARIMappingTable();
  fHiCARIcalfile = set->GetValue("HiCARI.Calibration.File",defaultfile);
  fBaselineLength = set->GetValue("BaseLine.Length",60);
  fTracePlots = set->GetValue("Trace.Plots",0);
  fMode3Histos = set->GetValue("Mode3.Histos",0);
  fExcludeTracking = set->GetValue("Merge.ExcludeTracking",0);
  fMatrixfile = set->GetValue("Matrix.File",defaultfile);

  fPPACfile = set->GetValue("BigRIPS.PPAC.File","/home/gamma20/exp/db/db/BigRIPSPPAC.xml");
  fPPACdefaultfile = set->GetValue("BigRIPS.PPAC.Def.File","/home/gamma20/exp/db/db/BigRIPSPPAC.xml");
  fplasticfile = set->GetValue("BigRIPS.Plastic.File","/home/gamma20/exp/db/db/BigRIPSPlastic.xml");
  fICfile = set->GetValue("BigRIPS.IC.File","/home/gamma20/exp/db/db/BigRIPSIC.xml");
  ffocalfile = set->GetValue("BigRIPS.Focal.File","/home/gamma20/exp/db/db/FocalPlane.xml");
  fmatrixfile[0] = set->GetValue("Matrix.35.File","/home/gamma20/exp/db/matrix/mat1.mat");
  fmatrixfile[1] = set->GetValue("Matrix.57.File","/home/gamma20/exp/db/matrix/mat2.mat");
  fmatrixfile[2] = set->GetValue("Matrix.89.File","/home/gamma20/exp/db/matrix/F8F9_LargeAccAchr.mat");
  fmatrixfile[3] = set->GetValue("Matrix.911.File","/home/gamma20/exp/db/matrix/F9F11_LargeAccAchr.mat");
  for(int i=0;i<6;i++)
    ftoffset[i] = set->GetValue(Form("TOF.Offset.%d",i),300.0);

  fbeta = set->GetValue("AverageBeta",0.5);

  fppac3align[0] = set->GetValue("PPAC3.Align.X0",0.0);
  fppac3align[1] = set->GetValue("PPAC3.Align.Y0",0.0);
  fppac3align[2] = set->GetValue("PPAC3.Align.X1",0.0);
  fppac3align[3] = set->GetValue("PPAC3.Align.Y1",0.0);
  ftargetposition = set->GetValue("Target.Position",129.5);
  ff5xgate[0] = set->GetValue("F5X.Gate.Low", -200.);
  ff5xgate[1] = set->GetValue("F5X.Gate.High", 200.);
  fdeltagate[2] = set->GetValue("Delta.Gate.Low", -999.);
  fdeltagate[3] = set->GetValue("Delta.Gate.High", 999.);
  fdeltagate[0] = set->GetValue("Delta.Gate.BR.Low", -999.);
  fdeltagate[1] = set->GetValue("Delta.Gate.BR.High", 999.);

  fcorrelationMode = set->GetValue("Correlation.Mode",0);
  fcorrelationCluster = set->GetValue("Correlation.Cluster",3);
  fcorrelationCrystal = set->GetValue("Correlation.Crystal",0);
  fcorrelationChannel = set->GetValue("Correlation.Channel",0);
  
  fBigRIPSCluster = set->GetValue("BigRIPS.Cluster",9);
  fBigRIPSCrystal = set->GetValue("BigRIPS.Crystal",9);
  fBigRIPSChannel = set->GetValue("BigRIPS.Channel",9);
  fBigRIPSDetail = set->GetValue("BigRIPS.Detail",-1);

  fTimeCorFile = set->GetValue("Time.Corrections.File","settings.dat");
  fEvtNrFile = set->GetValue("Event.Number.File","settings.dat");

  faoq_corr[0][0][0] = set->GetValue("BigRIPS.AoQCorr_F3X",0.0);
  faoq_corr[0][0][1] = set->GetValue("BigRIPS.AoQCorr_F3A",0.0);
  faoq_corr[0][1][0] = set->GetValue("BigRIPS.AoQCorr_F5X",0.0);
  faoq_corr[0][1][1] = set->GetValue("BigRIPS.AoQCorr_F5A",0.0);
  faoq_corr[0][2][0] = set->GetValue("BigRIPS.AoQCorr_F7X",0.0);
  faoq_corr[0][2][1] = set->GetValue("BigRIPS.AoQCorr_F7A",0.0);
  
  faoq_corr[1][0][0] = set->GetValue("ZeroDeg.AoQCorr_F8X",0.0);
  faoq_corr[1][0][1] = set->GetValue("ZeroDeg.AoQCorr_F8A",0.0);
  faoq_corr[1][1][0] = set->GetValue("ZeroDeg.AoQCorr_F9X",0.0);
  faoq_corr[1][1][1] = set->GetValue("ZeroDeg.AoQCorr_F9A",0.0);
  faoq_corr[1][2][0] = set->GetValue("ZeroDeg.AoQCorr_F11X",0.0);
  faoq_corr[1][2][1] = set->GetValue("ZeroDeg.AoQCorr_F11A",0.0);

}

/*
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
*/

void Settings::PrintSettings(){
  //cout << __PRETTY_FUNCTION__ << endl;
  cout << "Target.X\t" << ftargetX << endl;
  cout << "Target.Y\t" << ftargetY << endl;
  cout << "Target.Z\t" << ftargetZ << endl;
  cout << "Target.Beta\t" << ftargetBeta << endl;
  cout << "Average.Beta.After\t" << fAveAfterBeta << endl;

  cout << "AddBackType\t" << fAddBackType << endl;
  cout << "CoincTimeDiff\t" << fCoincTimeDiff << endl;
  cout << "ClusterAngle\t" << fClusterAngle << endl;
  cout << "StoreAllIPoints\t" << fStoreAllIPoints << endl;
  cout << "OverflowThreshold\t"<< fOverflowThreshold << endl;

  cout << "DoTracking\t"<< fTracking << endl;

  cout << "HiCARI.IgnoreTrace\t" << fIgnoreTrace << endl;
  cout << "HiCARI.Positions\t" << fHiCARIPos << endl;
  cout << "HiCARI.Mapping.Table\t"<< fHiCARImapping << endl;
  PrintHiCARIMappingTable();
  cout << "HiCARI.Raw.Thresh\t"<< fRawThresh << endl;
  cout << "HiCARI.Calibration.File\t"<< fHiCARIcalfile << endl;
  cout << "BaseLine.Length\t"<< fBaselineLength << endl;
  cout << "Trace.Plots\t"<< fTracePlots << endl;
  cout << "Mode3.Histos\t"<< fMode3Histos << endl;
  cout << "ExcludeTracking\t"<< fExcludeTracking << endl;
  cout << "MAtrixFile\t"<< fMatrixfile << endl;

  for(int i=0;i<6;i++)
    cout << Form("TOF offset.%d\t",i) << ftoffset[i] << endl;

  cout << "BigRIPS.PPAC.File\t"	<< fPPACfile << endl;
  cout << "BigRIPS.PPAC.Def.File\t"	<< fPPACdefaultfile << endl;
  cout << "BigRIPS.Plastic.File\t" << fplasticfile << endl;
  cout << "BigRIPS.IC.File\t" << fICfile << endl;
  cout << "BigRIPS.Focal.File\t" << ffocalfile	<< endl;
  cout << "Matrix.35.File\t" << fmatrixfile[0] << endl;
  cout << "Matrix.57.File\t" << fmatrixfile[1] << endl;
  cout << "Matrix.89.File\t" << fmatrixfile[2] << endl;
  cout << "Matrix.911.File\t" << fmatrixfile[3] << endl;

  cout << "beta\t" << fbeta << endl;  
  cout << "align PPAC 3 x0 = " << fppac3align[0] << " , y0 = " << fppac3align[1] << endl;
  cout << "align PPAC 3 x1 = " << fppac3align[2] << " , y1 = " << fppac3align[3] << endl;
  cout << "target position\t" << ftargetposition << endl;
  cout << "gate on F5X position\t" <<ff5xgate[0] << " to " << ff5xgate[1] << endl;
  cout << "gate on delta for charge changes BR\t" <<fdeltagate[0] << " to " << fdeltagate[1] << endl;
  cout << "gate on delta for charge changes ZD\t" <<fdeltagate[2] << " to " << fdeltagate[3] << endl;

  cout << "Correlation.Mode\t" << fcorrelationMode << endl;  
  cout << "Correlation.Cluster\t" << fcorrelationCluster << endl;  
  cout << "Correlation.Crystal\t" << fcorrelationCrystal << endl;  
  cout << "Correlation.Channel\t" << fcorrelationChannel << endl;  
 
  cout << "BigRIPS.Cluster\t" << fBigRIPSCluster << endl;  
  cout << "BigRIPS.Crystal\t" << fBigRIPSCrystal << endl;  
  cout << "BigRIPS.Channel\t" << fBigRIPSChannel << endl;  
  cout << "BigRIPS.Tree.Detail\t" << fBigRIPSDetail << endl;  

  cout << "Time.Corrections.File\t" << fTimeCorFile << endl; 
  cout << "Event.Number.File\t" << fEvtNrFile << endl;

  cout << "BigRIPS.AoQCorr_F3X\t" << faoq_corr[0][0][0] << endl;  
  cout << "BigRIPS.AoQCorr_F3A\t" << faoq_corr[0][0][1] << endl;  
  cout << "BigRIPS.AoQCorr_F5X\t" << faoq_corr[0][1][0] << endl;  
  cout << "BigRIPS.AoQCorr_F5A\t" << faoq_corr[0][1][1] << endl;  
  cout << "BigRIPS.AoQCorr_F7X\t" << faoq_corr[0][2][0] << endl;  
  cout << "BigRIPS.AoQCorr_F7A\t" << faoq_corr[0][2][1] << endl;  
  cout << "ZeroDeg.AoQCorr_F8X\t" << faoq_corr[1][0][0] << endl;  
  cout << "ZeroDeg.AoQCorr_F8A\t" << faoq_corr[1][0][1] << endl;  
  cout << "ZeroDeg.AoQCorr_F9X\t" << faoq_corr[1][1][0] << endl;  
  cout << "ZeroDeg.AoQCorr_F9A\t" << faoq_corr[1][1][1] << endl;  
  cout << "ZeroDeg.AoQCorr_F11X\t" << faoq_corr[1][2][0] << endl;  
  cout << "ZeroDeg.AoQCorr_F11A\t" << faoq_corr[1][2][1] << endl;  


}
void Settings::ReadHiCARIMappingTable(){
  TEnv* mapenv = new TEnv(fHiCARImapping.c_str());
  for(int h=0;h<20;h++){
    for(int c=0;c<4;c++){
      for(int s=0;s<4;s++){
	fHiCARImap[h][c][s] = mapenv->GetValue(Form("Hole.%02d.Crystal.%d.Slot.%d",h,c,s),-1);
      }
    }
  }
}

void Settings::PrintHiCARIMappingTable(){
  // for accessing first map 
  map<int, map<int, map<int, int> > >::iterator ftr;   
  // for accessing second map 
  map<int, map<int, int> >::iterator str; 
  // for accessing inner map 
  map<int, int>::iterator ptr;

  char* abc = "ABCD";
  for(ftr = fHiCARImap.begin(); ftr != fHiCARImap.end(); ftr++) {   
    for(str = ftr->second.begin(); str != ftr->second.end(); str++) { 
      for(ptr = str->second.begin(); ptr != str->second.end(); ptr++) {
	if(ptr->second>-1)
	  cout << "Hole " << ftr->first 
	       << ", Crystal " << str->first 
	       << ", Slot " << ptr->first 
	       << " is Ge " << ptr->second/10 << abc[ptr->second%10] << ", or " << fHiCARImap[ftr->first][str->first][ptr->first] << ", or " << HiCARICluster(ftr->first,str->first,ptr->first)<< " and " << HiCARICrystal(ftr->first,str->first,ptr->first) << endl; 
	  
      } 
    } 
  }
}
