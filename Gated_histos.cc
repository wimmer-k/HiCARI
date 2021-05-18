#include <iostream>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <signal.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TCutG.h"
#include "TStopwatch.h"

#include "CommandLineInterface.hh"
#include "RunInfo.hh"
#include "Beam.hh"
#include "HiCARI.hh"
#include "Gretina.hh"
#include "FocalPlane.hh"

#include "Globaldefs.h"
using namespace TMath;
using namespace std;

bool signal_received = false;
void signalhandler(int sig);
double get_time();
int main(int argc, char* argv[]){
  double time_start = get_time();  
  TStopwatch timer;
  timer.Start();
  signal(SIGINT,signalhandler);
  vector<char*> InputFiles;char* BigRIPSFile = NULL;
  char* OutputFile = NULL;
  char* SettingFile = NULL;
  int nmax =0;
  int vl =0;
  double targetz = 0;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-v", "verbose", &vl);
  interface->Add("-t", "target z position", &targetz);
  interface->CheckFlags(argc, argv);

  if(InputFiles.size() == 0 || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }
  cout<<"input file(s):"<<endl;
  for(unsigned int i=0; i<InputFiles.size(); i++){
    cout<<InputFiles[i]<<endl;
  }
  cout<<"output file: "<<OutputFile<< endl;

  int incuts = 0;
  int outcuts = 0;

  TCutG* TimeCut[2] = {NULL,NULL};
  
  vector<TCutG*> InCut;
  vector<TCutG*> OutCut;
  vector<int> InCut_sel;
  vector<double> beta;

  int useCorrected = 0;
  int usePlastic = 0;
  int BR_AoQ = 2;
  int ZD_AoQ = 5;
  int nbins = 4000;
  double erange = 4000;
  double zrange[2] = {10,30};
  double aoqrange[2] = {2.2,2.8};
  


  if(SettingFile != NULL){
    TEnv* set = new TEnv(SettingFile);
    char* cfilename = (char*)set->GetValue("CutFile","file");
    incuts = set->GetValue("NInCuts",0);
    outcuts = set->GetValue("NOutCuts",0);

    useCorrected = set->GetValue("UseAoQCorr",0);
    BR_AoQ = set->GetValue("UseBRAoQ",2);
    ZD_AoQ = set->GetValue("UseZDAoQ",5);

    nbins = set->GetValue("Energy.Bins",4000);
    erange = set->GetValue("Energy.Range",4000);
    zrange[0] = set->GetValue("Z.Range.Min",10.);
    zrange[1] = set->GetValue("Z.Range.Max",30.);
    aoqrange[0] = set->GetValue("AoQ.Range.Min",2.2);
    aoqrange[1] = set->GetValue("AoQ.Range.Max",2.8);

    cout<<"Using A/q ["<<BR_AoQ<<"] (BigRIPS) and ["<<ZD_AoQ<<"] (ZD) " << endl;
    
    if(useCorrected)
      cout << "using corrected A/q values for PID gates" << endl;
    else
      cout << "using raw A/q values for PID gates" << endl;
    usePlastic = set->GetValue("UsePlastic",0);
    if(usePlastic)
      cout << "using F7 Plastic for BR PID gates" << endl;
    else
      cout << "using F7 IC for BR PID gates" << endl;
    
    TFile* cFile = new TFile(cfilename);

    char* tcutname = (char*)set->GetValue("TimeCut.HiCARI","name");
    if(cFile->GetListOfKeys()->Contains(tcutname)){
      cout << "time cut found "<< tcutname << endl;
      TimeCut[0] = (TCutG*)cFile->Get(tcutname);
    }
    tcutname = (char*)set->GetValue("TimeCut.Mode2","name");
    if(cFile->GetListOfKeys()->Contains(tcutname)){
      cout << "time cut found "<< tcutname << endl;
      TimeCut[1] = (TCutG*)cFile->Get(tcutname);
    }
    for(int i=0;i<incuts;i++){
      char* iname = (char*)set->GetValue(Form("InCut.%d",i),"file");
      if(cFile->GetListOfKeys()->Contains(iname)){
	cout << "incoming cut found "<< iname << endl;
	TCutG* icg = (TCutG*)cFile->Get(iname);
	InCut.push_back(icg);
      }
      else{
	cout << "Incoming cut " << iname << " does not exist! Error in settings file" << endl;
	return 77;
      }
    }
    for(int i=0;i<outcuts;i++){
      char* oname = (char*)set->GetValue(Form("OutCut.%d",i),"file");
      if(cFile->GetListOfKeys()->Contains(oname)){
	cout << "outgoing cut found "<< oname << endl;
	TCutG* ocg = (TCutG*)cFile->Get(oname);
	OutCut.push_back(ocg);
	double b = set->GetValue(Form("Beta.%d",i),0.57);
	beta.push_back(b);
	int s = set->GetValue(Form("InCutSel.%d",i),0);
	if(s<(int)InCut.size()){
	  InCut_sel.push_back(s);
	}
	else{
	  cout << "Selected incoming cut Nr. " << s << " does not exist! Error in settings file" << endl;
	  return 77;
	}
      }
      else{
	cout << "Outgoing cut " << oname << " does not exist! Error in settings file" << endl;
	return 77;
      }
    }
    
    cFile->Close();
  }//settings present
  
  TChain* tr;
  tr = new TChain("tr");
  for(unsigned int i=0; i<InputFiles.size(); i++){
    tr->Add(InputFiles[i]);
  }
  if(tr == NULL){
    cout << "could not find tree build in file " << endl;
    for(unsigned int i=0; i<InputFiles.size(); i++){
      cout<<InputFiles[i]<<endl;
    }
    return 3;
  }
  HiCARICalc* hi = new HiCARICalc;
  tr->SetBranchAddress("hicari",&hi);
  GretinaCalc* gr = new GretinaCalc;
  tr->SetBranchAddress("mode2",&gr);

  int trigbit = -1;
  tr->SetBranchAddress("trigbit",&trigbit);
  Beam* bz = new Beam;
  tr->SetBranchAddress("beam",&bz);
  FocalPlane* fp[NFPLANES];
  for(unsigned short f=0;f<NFPLANES;f++){
    fp[f] = new FocalPlane;
  }
  for(unsigned short f=0;f<NFPLANES;f++){
    tr->SetBranchAddress(Form("fp%d",fpID[f]),&fp[f]);
  }
  Double_t nentries = tr->GetEntries();


  TList *hlist = new TList();
  TH1F* trigger = new TH1F("trigger","trigger",10,0,10);hlist->Add(trigger);
  TH2F* bigrips = new TH2F("bigrips","bigrips",1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);hlist->Add(bigrips);
  TH2F* bigrips_TB1 = new TH2F("bigrips_TB1","bigrips_TB1",1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);hlist->Add(bigrips_TB1);
  TH2F* zerodeg = new TH2F("zerodeg","zerodeg",1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);hlist->Add(zerodeg);
  TH2F* bigrips_Pl = new TH2F("bigrips_Pl","bigrips_Pl",1000,aoqrange[0],aoqrange[1],2000,0,2000);hlist->Add(bigrips_Pl);
  TH1F* h_egamdc = new TH1F("h_egamdc","h_egamdc",nbins,0,erange);hlist->Add(h_egamdc);
  TH1F* g_egamdc = new TH1F("g_egamdc","g_egamdc",nbins,0,erange);hlist->Add(g_egamdc);
  TH1F* h_egamABdc = new TH1F("h_egamABdc","h_egamABdc",nbins,0,erange);hlist->Add(h_egamABdc);
  TH1F* g_egamABdc = new TH1F("g_egamABdc","g_egamABdc",nbins,0,erange);hlist->Add(g_egamABdc);
  TH2F* h_egamtgam = new TH2F("h_egamtgam","h_egamtgam",1000,-500,500,nbins,0,erange);hlist->Add(h_egamtgam);
  TH2F* g_egamtgam = new TH2F("g_egamtgam","g_egamtgam",1000,-500,500,nbins,0,erange);hlist->Add(g_egamtgam);
  TH2F* h_egamtgamAB = new TH2F("h_egamtgamAB","h_egamtgamAB",1000,-500,500,nbins,0,erange);hlist->Add(h_egamtgamAB);
  TH2F* g_egamtgamAB = new TH2F("g_egamtgamAB","g_egamtgamAB",1000,-500,500,nbins,0,erange);hlist->Add(g_egamtgamAB);
  TH2F* h_egamtgamdc = new TH2F("h_egamtgamdc","h_egamtgamdc",1000,-500,500,nbins,0,erange);hlist->Add(h_egamtgamdc);
  TH2F* g_egamtgamdc = new TH2F("g_egamtgamdc","g_egamtgamdc",1000,-500,500,nbins,0,erange);hlist->Add(g_egamtgamdc);
  TH2F* h_egamtgamABdc = new TH2F("h_egamtgamABdc","h_egamtgamABdc",1000,-500,500,nbins,0,erange);hlist->Add(h_egamtgamABdc);
  TH2F* g_egamtgamABdc = new TH2F("g_egamtgamABdc","g_egamtgamABdc",1000,-500,500,nbins,0,erange);hlist->Add(g_egamtgamABdc);

  TH2F* target = new TH2F("target","target",1000,-50,50,1000,-50,50);hlist->Add(target);
  TH2F* target_hit = new TH2F("target_hit","target_hit",1000,-50,50,1000,-50,50);hlist->Add(target_hit);
  TH2F* target_frame = new TH2F("target_frame","target_frame",1000,-50,50,1000,-50,50);hlist->Add(target_frame);
  
  
  vector<TH2F*> zerodeg_c;
  vector<TH1F*> h_egamdc_c;
  vector<TH1F*> g_egamdc_c;
  vector<TH1F*> h_egamABdc_c;
  vector<TH1F*> g_egamABdc_c;
  
  vector<TH2F*> h_egam_summary_c;
  vector<TH2F*> h_tgam_summary_c;
  vector<TH2F*> h_egamdc_summary_c;
  vector<TH2F*> h_egamAB_summary_c;
  vector<TH2F*> h_tgamAB_summary_c;
  vector<TH2F*> h_egamABdc_summary_c;

  vector<TH2F*> g_egamdc_summary_c;
  vector<TH2F*> g_egamABdc_summary_c;
  
  vector<TH2F*> h_egam_theta_c;
  vector<TH2F*> h_egamdc_theta_c;
  vector<TH2F*> h_egamAB_theta_c;
  vector<TH2F*> h_egamABdc_theta_c;

  vector<TH2F*> g_egam_theta_c;
  vector<TH2F*> g_egamdc_theta_c;
  vector<TH2F*> g_egamAB_theta_c;
  vector<TH2F*> g_egamABdc_theta_c;

  vector<TH2F*> h_egamdc_hmult_c;
  vector<TH2F*> h_egamABdc_hmultAB_c;
  vector<TH2F*> h_egamdc_tmult_c;
  vector<TH2F*> h_egamABdc_tmultAB_c;

  vector<TH2F*> g_egamdc_gmult_c;
  vector<TH2F*> g_egamABdc_gmultAB_c;
  vector<TH2F*> g_egamdc_tmult_c;
  vector<TH2F*> g_egamABdc_tmultAB_c;
  
  vector<TH2F*> g_egamdc_depth_c;
  vector<TH2F*> g_egamABdc_depth_c;
  
  vector<TH2F*> h_egamtgamdc_c;
  vector<TH2F*> g_egamtgamdc_c;
  vector<TH2F*> h_egamtgamABdc_c;
  vector<TH2F*> g_egamtgamABdc_c;
  vector<TH2F*> h_egamegamdc_c;
  vector<TH2F*> g_egamegamdc_c;
  vector<TH2F*> gh_egamegamdc_c;
  vector<TH2F*> h_egamegamABdc_c;
  vector<TH2F*> g_egamegamABdc_c;
  vector<TH2F*> gh_egamegamABdc_c;
  for(int i=0;i<incuts;i++){
    TH2F* h = new TH2F(Form("zerodeg_%s",InCut[i]->GetName()),
		       Form("zerodeg_%s",InCut[i]->GetName()),
		       1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);
    zerodeg_c.push_back(h);
    hlist->Add(zerodeg_c.back());
  }
  
  for(int o=0;o<outcuts;o++){
    TH1F* h;
    
    h = new TH1F(Form("h_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 Form("h_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 nbins,0,erange);
    h_egamdc_c.push_back(h);
    hlist->Add(h_egamdc_c.back());
    h = new TH1F(Form("g_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 Form("g_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 nbins,0,erange);
    g_egamdc_c.push_back(h);
    hlist->Add(g_egamdc_c.back());
    h = new TH1F(Form("h_egamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 Form("h_egamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 nbins,0,erange);
    h_egamABdc_c.push_back(h);
    hlist->Add(h_egamABdc_c.back());
    h = new TH1F(Form("g_egamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 Form("g_egamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 nbins,0,erange);
    g_egamABdc_c.push_back(h);
    hlist->Add(g_egamABdc_c.back());


    TH2F* h2;
    //summary plots
    h2 = new TH2F(Form("h_egam_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egam_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,nbins,0,erange);
    h_egam_summary_c.push_back(h2);
    hlist->Add(h_egam_summary_c.back());
    
    h2 = new TH2F(Form("h_tgam_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_tgam_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,1000,-500,500);
    h_tgam_summary_c.push_back(h2);
    hlist->Add(h_tgam_summary_c.back());
    
    h2 = new TH2F(Form("h_egamdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,nbins,0,erange);
    h_egamdc_summary_c.push_back(h2);
    hlist->Add(h_egamdc_summary_c.back());

    h2 = new TH2F(Form("h_egamAB_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamAB_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,nbins,0,erange);
    h_egamAB_summary_c.push_back(h2);
    hlist->Add(h_egamAB_summary_c.back());
    
    h2 = new TH2F(Form("h_tgamAB_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_tgamAB_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,1000,-500,500);
    h_tgamAB_summary_c.push_back(h2);
    hlist->Add(h_tgamAB_summary_c.back());
    
    h2 = new TH2F(Form("h_egamABdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamABdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,nbins,0,erange);
    h_egamABdc_summary_c.push_back(h2);
    hlist->Add(h_egamABdc_summary_c.back());


    h2 = new TH2F(Form("g_egamdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  64,0,64,nbins,0,erange);
    g_egamdc_summary_c.push_back(h2);
    hlist->Add(g_egamdc_summary_c.back());

    h2 = new TH2F(Form("g_egamABdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamABdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  64,0,64,nbins,0,erange);
    g_egamABdc_summary_c.push_back(h2);
    hlist->Add(g_egamABdc_summary_c.back());


    //multiplicity
    h2 = new TH2F(Form("h_egamdc_hmult_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamdc_hmult_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  40,0,40,nbins,0,erange);
    h_egamdc_hmult_c.push_back(h2);
    hlist->Add(h_egamdc_hmult_c.back());
    
    h2 = new TH2F(Form("h_egamABdc_hmultAB_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamABdc_hmultAB_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  40,0,40,nbins,0,erange);
    h_egamABdc_hmultAB_c.push_back(h2);
    hlist->Add(h_egamABdc_hmultAB_c.back());
    
    h2 = new TH2F(Form("h_egamdc_tmult_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamdc_tmult_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  40,0,40,nbins,0,erange);
    h_egamdc_tmult_c.push_back(h2);
    hlist->Add(h_egamdc_tmult_c.back());
    
    h2 = new TH2F(Form("h_egamABdc_tmultAB_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamABdc_tmultAB_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  40,0,40,nbins,0,erange);
    h_egamABdc_tmultAB_c.push_back(h2);
    hlist->Add(h_egamABdc_tmultAB_c.back());
    
    h2 = new TH2F(Form("g_egamdc_gmult_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamdc_gmult_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  40,0,40,nbins,0,erange);
    g_egamdc_gmult_c.push_back(h2);
    hlist->Add(g_egamdc_gmult_c.back());
    
    h2 = new TH2F(Form("g_egamABdc_gmultAB_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamABdc_gmultAB_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  40,0,40,nbins,0,erange);
    g_egamABdc_gmultAB_c.push_back(h2);
    hlist->Add(g_egamABdc_gmultAB_c.back());
    
    h2 = new TH2F(Form("g_egamdc_tmult_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamdc_tmult_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  40,0,40,nbins,0,erange);
    g_egamdc_tmult_c.push_back(h2);
    hlist->Add(g_egamdc_tmult_c.back());
    
    h2 = new TH2F(Form("g_egamABdc_tmultAB_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamABdc_tmultAB_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  40,0,40,nbins,0,erange);
    g_egamABdc_tmultAB_c.push_back(h2);
    hlist->Add(g_egamABdc_tmultAB_c.back());
    
    //depth in tracking detectors
    h2 = new TH2F(Form("g_egamdc_depth_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamdc_depth_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  300,0,300,nbins,0,erange);
    g_egamdc_depth_c.push_back(h2);
    hlist->Add(g_egamdc_depth_c.back());
    
    h2 = new TH2F(Form("g_egamABdc_depth_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamABdc_depth_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  300,0,300,nbins,0,erange);
    g_egamABdc_depth_c.push_back(h2);
    hlist->Add(g_egamABdc_depth_c.back());
    


    //verus theta 
    h2 = new TH2F(Form("h_egam_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egam_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,nbins,0,erange);
    h_egam_theta_c.push_back(h2);
    hlist->Add(h_egam_theta_c.back());
    
    h2 = new TH2F(Form("h_egamdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,nbins,0,erange);
    h_egamdc_theta_c.push_back(h2);
    hlist->Add(h_egamdc_theta_c.back());

    h2 = new TH2F(Form("h_egamAB_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamAB_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,nbins,0,erange);
    h_egamAB_theta_c.push_back(h2);
    hlist->Add(h_egamAB_theta_c.back());
    
    h2 = new TH2F(Form("h_egamABdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamABdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,nbins,0,erange);
    h_egamABdc_theta_c.push_back(h2);
    hlist->Add(h_egamABdc_theta_c.back());


    h2 = new TH2F(Form("g_egam_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egam_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,nbins,0,erange);
    g_egam_theta_c.push_back(h2);
    hlist->Add(g_egam_theta_c.back());
    
    h2 = new TH2F(Form("g_egamdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,nbins,0,erange);
    g_egamdc_theta_c.push_back(h2);
    hlist->Add(g_egamdc_theta_c.back());

    h2 = new TH2F(Form("g_egamAB_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamAB_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,nbins,0,erange);
    g_egamAB_theta_c.push_back(h2);
    hlist->Add(g_egamAB_theta_c.back());
    
    h2 = new TH2F(Form("g_egamABdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamABdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,nbins,0,erange);
    g_egamABdc_theta_c.push_back(h2);
    hlist->Add(g_egamABdc_theta_c.back());



    
    //e gamma versus time diff to BR
    h2 = new TH2F(Form("h_egamtgamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamtgamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,-500,500,nbins,0,erange);
    h_egamtgamdc_c.push_back(h2);
    hlist->Add(h_egamtgamdc_c.back());

    h2 = new TH2F(Form("g_egamtgamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamtgamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,-500,500,nbins,0,erange);
    g_egamtgamdc_c.push_back(h2);
    hlist->Add(g_egamtgamdc_c.back());

    h2 = new TH2F(Form("h_egamtgamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamtgamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,-500,500,nbins,0,erange);
    h_egamtgamABdc_c.push_back(h2);
    hlist->Add(h_egamtgamABdc_c.back());

    h2 = new TH2F(Form("g_egamtgamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamtgamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,-500,500,nbins,0,erange);
    g_egamtgamABdc_c.push_back(h2);
    hlist->Add(g_egamtgamABdc_c.back());


    //e gamma versus e gamma
    h2 = new TH2F(Form("h_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  nbins/4,0,erange,nbins/4,0,erange);
    h_egamegamdc_c.push_back(h2);
    hlist->Add(h_egamegamdc_c.back());

    h2 = new TH2F(Form("g_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  nbins/4,0,erange,nbins/4,0,erange);
    g_egamegamdc_c.push_back(h2);
    hlist->Add(g_egamegamdc_c.back());

    h2 = new TH2F(Form("gh_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("gh_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  nbins/4,0,erange,nbins/4,0,erange);
    gh_egamegamdc_c.push_back(h2);
    hlist->Add(gh_egamegamdc_c.back());

    h2 = new TH2F(Form("h_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  nbins/4,0,erange,nbins/4,0,erange);
    h_egamegamABdc_c.push_back(h2);
    hlist->Add(h_egamegamABdc_c.back());

    h2 = new TH2F(Form("g_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  nbins/4,0,erange,nbins/4,0,erange);
    g_egamegamABdc_c.push_back(h2);
    hlist->Add(g_egamegamABdc_c.back());

    h2 = new TH2F(Form("gh_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("gh_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  nbins/4,0,erange,nbins/4,0,erange);
    gh_egamegamABdc_c.push_back(h2);
    hlist->Add(gh_egamegamABdc_c.back());


    
  }
  
  Int_t nbytes = 0;
  Int_t status;
  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    hi->Clear();
    gr->Clear();
    bz->Clear();
    for(unsigned short f=0;f<NFPLANES;f++){
      fp[f]->Clear();
    }
    if(vl>2)
      cout << "getting entry " << i << endl;
    status = tr->GetEvent(i);
    if(vl>2)
      cout << "status " << status << endl;
    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;

    trigger->Fill(trigbit);
    if(!useCorrected){
      bigrips->Fill(bz->GetAQ(BR_AoQ),bz->GetZ(BR_AoQ));
      zerodeg->Fill(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ));
      bigrips_Pl->Fill(bz->GetAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge());
      if((trigbit&1)==1)
	bigrips_TB1->Fill(bz->GetAQ(BR_AoQ),bz->GetZ(BR_AoQ));
      
    }
    else{
      bigrips->Fill(bz->GetCorrAQ(BR_AoQ),bz->GetZ(BR_AoQ));
      zerodeg->Fill(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ));
      bigrips_Pl->Fill(bz->GetCorrAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge());
      if((trigbit&1)==1)
	bigrips_TB1->Fill(bz->GetAQ(BR_AoQ),bz->GetZ(BR_AoQ));
    }
    bool ingood = false;
    bool outgood = false;
    for(int in=0;in<incuts;in++){
      if(!useCorrected && !usePlastic && InCut[in]->IsInside(bz->GetAQ(BR_AoQ),bz->GetZ(BR_AoQ))){
	zerodeg_c[in]->Fill(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ));
	ingood = true;
      }
      else if(useCorrected && !usePlastic && InCut[in]->IsInside(bz->GetCorrAQ(BR_AoQ),bz->GetZ(BR_AoQ))){
	zerodeg_c[in]->Fill(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ));
	ingood = true;
      }
      else if(!useCorrected && usePlastic && InCut[in]->IsInside(bz->GetAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge())){
	zerodeg_c[in]->Fill(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ));
	ingood = true;
      }
      else if(useCorrected && usePlastic && InCut[in]->IsInside(bz->GetCorrAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge())){
	zerodeg_c[in]->Fill(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ));
	ingood = true;
      }
    }
    for(int o=0;o<outcuts;o++){
      if(!useCorrected && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ)))
	outgood = true;
      else if(useCorrected && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ)))
	outgood = true;	
    }
    if(ingood && outgood){
      target->Fill(fp[fpNr(8)]->GetTrack()->GetX(),fp[fpNr(8)]->GetTrack()->GetY());
      if(bz->GetBeta(0) - bz->GetBeta(1) < 0.114)
	target_frame->Fill(fp[fpNr(8)]->GetTrack()->GetX(),fp[fpNr(8)]->GetTrack()->GetY());
      else
	target_hit->Fill(fp[fpNr(8)]->GetTrack()->GetX(),fp[fpNr(8)]->GetTrack()->GetY());
    }

    
    int hmult = hi->GetMult() + gr->GetMult();
    int tmult = hi->GetMult() + gr->GetMult();
    if(hi->HadBigRIPS()){
      hmult-=1;
      tmult-=1;
    }
      
    for(int h=0; h<hi->GetMult(); h++){
      HiCARIHitCalc* hit = hi->GetHit(h);
      if(hit->IsBigRIPS()){
	continue;
      }
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	h_egamtgam->Fill(hit->GetTime(),hit->GetEnergy());
	h_egamtgamdc->Fill(hit->GetTime(),hit->GetDCEnergy());
	if(TimeCut[0] && TimeCut[0]->IsInside(hit->GetTime(),hit->GetDCEnergy()))
	  h_egamdc->Fill(hit->GetDCEnergy());
	else if(TimeCut[0]==NULL)
	  h_egamdc->Fill(hit->GetDCEnergy());
	for(int o=0;o<outcuts;o++){	
	  if( (!useCorrected && !usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(BR_AoQ),bz->GetZ(BR_AoQ)) && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ||
	      (useCorrected && !usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(BR_AoQ),bz->GetZ(BR_AoQ)) && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) || 
	      (!useCorrected && usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge()) && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ||
	      (useCorrected && usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge()) && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ){
	    double edc = hit->GetDCEnergy(beta[o],0,0,targetz);
	    h_tgam_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetTime());
	    h_egamtgamdc_c[o]->Fill(hit->GetTime(),edc);

	    // gate on timing
	    if(TimeCut[0]==NULL || TimeCut[0]->IsInside(hit->GetTime(),edc)){
	      h_egamdc_c[o]->Fill(edc);
	      h_egam_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetEnergy());
	      h_egamdc_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), edc);
	      h_egam_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetEnergy());
	      h_egamdc_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, edc);
	      h_egamdc_hmult_c[o]->Fill(hmult, edc);
	      h_egamdc_tmult_c[o]->Fill(tmult, edc);

	      for(int l=h+1; l<hi->GetMult(); l++){
		HiCARIHitCalc* hit2 = hi->GetHit(l);
		if(hit2->IsBigRIPS()){
		  continue;
		}
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
		if(TimeCut[0] && !TimeCut[0]->IsInside(hit2->GetTime(),edc2))
		  continue;
		if(edc>edc2)
		  h_egamegamdc_c[o]->Fill(edc,edc2);
		else
		  h_egamegamdc_c[o]->Fill(edc2,edc);
	      }//second hit
	      // gretina hits
	      for(int l=0; l<gr->GetMult(); l++){
		HitCalc* hit2 = gr->GetHit(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,0);
		if(TimeCut[1] && !TimeCut[1]->IsInside(hit2->GetTime(),edc2))
		  continue;
		gh_egamegamdc_c[o]->Fill(edc,edc2);
	      }//gretina hits


	      
	    }
	  }// in cut ZeroDeg
	}// loop over cuts
      }// energy and position good
    }// hits
    for(int g=0; g<gr->GetMult(); g++){
      HitCalc* hit = gr->GetHit(g);
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	g_egamtgam->Fill(hit->GetTime(),hit->GetEnergy());
	g_egamtgamdc->Fill(hit->GetTime(),hit->GetDCEnergy());
	if(TimeCut[1] && TimeCut[1]->IsInside(hit->GetTime(),hit->GetDCEnergy()))
	  g_egamdc->Fill(hit->GetDCEnergy());
	else if(TimeCut[1]==NULL)
	  g_egamdc->Fill(hit->GetDCEnergy());
	for(int o=0;o<outcuts;o++){	
	  if( (!useCorrected && !usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(BR_AoQ),bz->GetZ(BR_AoQ)) && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ||
	      (useCorrected && !usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(BR_AoQ),bz->GetZ(BR_AoQ)) && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) || 
	      (!useCorrected && usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge()) && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ||
	      (useCorrected && usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge()) && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ){
	    double edc = hit->GetDCEnergy(beta[o],0,0,0);
	    g_egamtgamdc_c[o]->Fill(hit->GetTime(),edc);
	    
	    // gate on timing
	    if(TimeCut[1]==NULL || TimeCut[1]->IsInside(hit->GetTime(),edc)){
	      g_egamdc_c[o]->Fill(edc);
	      g_egamdc_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), edc);
	      g_egam_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetEnergy());
	      g_egamdc_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, edc);
	      g_egamdc_gmult_c[o]->Fill(gr->GetMult(), edc);
	      g_egamdc_tmult_c[o]->Fill(tmult, edc);
	      g_egamdc_depth_c[o]->Fill(hit->GetPosition().Mag(), edc);
	      for(int l=g+1; l<gr->GetMult(); l++){
		HitCalc* hit2 = gr->GetHit(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,0);
		if(TimeCut[1] && !TimeCut[1]->IsInside(hit2->GetTime(),edc2))
		  continue;
		if(edc>edc2)
		  g_egamegamdc_c[o]->Fill(edc,edc2);
		else
		  g_egamegamdc_c[o]->Fill(edc2,edc);
	      }//second hit

	    }
	  }// in cut ZeroDeg
	}// loop over cuts
      }// energy and position good
    }// hits
    int tmultAB = hi->GetMultAB() + gr->GetMultAB();
    for(int h=0; h<hi->GetMultAB(); h++){
      HiCARIHitCalc* hit = hi->GetHitAB(h);
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	h_egamtgamAB->Fill(hit->GetTime(),hit->GetEnergy());
	h_egamtgamABdc->Fill(hit->GetTime(),hit->GetDCEnergy());
	if(TimeCut[0] && TimeCut[0]->IsInside(hit->GetTime(),hit->GetDCEnergy()))
	  h_egamABdc->Fill(hit->GetDCEnergy());
	else if(TimeCut[0]==NULL)
	  h_egamABdc->Fill(hit->GetDCEnergy());
	for(int o=0;o<outcuts;o++){	
	  if( (!useCorrected && !usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(BR_AoQ),bz->GetZ(BR_AoQ)) && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ||
	      (useCorrected && !usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(BR_AoQ),bz->GetZ(BR_AoQ)) && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) || 
	      (!useCorrected && usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge()) && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ||
	      (useCorrected && usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge()) && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ){
	    double edc = hit->GetDCEnergy(beta[o],0,0,targetz);
	    h_tgamAB_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetTime());
	    h_egamtgamABdc_c[o]->Fill(hit->GetTime(),edc);
	    // gate on timing
	    if(TimeCut[0]==NULL || TimeCut[0]->IsInside(hit->GetTime(),edc)){
	      h_egamABdc_c[o]->Fill(edc);
	      h_egamAB_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetEnergy());
	      h_egamABdc_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), edc);
	      h_egamAB_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetEnergy());
	      h_egamABdc_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, edc);
	      h_egamABdc_hmultAB_c[o]->Fill(hi->GetMultAB(), edc);
	      h_egamABdc_tmultAB_c[o]->Fill(tmultAB, edc);

	      for(int l=h+1; l<hi->GetMultAB(); l++){
		HiCARIHitCalc* hit2 = hi->GetHitAB(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
		if(TimeCut[0] && !TimeCut[0]->IsInside(hit2->GetTime(),edc2))
		  continue;
		if(edc>edc2)
		  h_egamegamABdc_c[o]->Fill(edc,edc2);
		else
		  h_egamegamABdc_c[o]->Fill(edc2,edc);
	      }//second hit
	      // gretina hits
	      for(int l=0; l<gr->GetMultAB(); l++){
		HitCalc* hit2 = gr->GetHitAB(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,0);
		if(TimeCut[1] && !TimeCut[1]->IsInside(hit2->GetTime(),edc2))
		  continue;
		gh_egamegamABdc_c[o]->Fill(edc,edc2);
	      }//gretina hits
	    }
	  }// in cut ZeroDeg
	}// loop over cuts
      }// energy and position good
    }// hits
    for(int g=0; g<gr->GetMultAB(); g++){
      HitCalc* hit = gr->GetHitAB(g);
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	g_egamtgamAB->Fill(hit->GetTime(),hit->GetEnergy());
	g_egamtgamABdc->Fill(hit->GetTime(),hit->GetDCEnergy());
	if(TimeCut[1] && TimeCut[1]->IsInside(hit->GetTime(),hit->GetDCEnergy()))
	  g_egamABdc->Fill(hit->GetDCEnergy());
	else if(TimeCut[1]==NULL)
	  g_egamABdc->Fill(hit->GetDCEnergy());
	for(int o=0;o<outcuts;o++){	
	  if( (!useCorrected && !usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(BR_AoQ),bz->GetZ(BR_AoQ)) && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ||
	      (useCorrected && !usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(BR_AoQ),bz->GetZ(BR_AoQ)) && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) || 
	      (!useCorrected && usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge()) && OutCut[o]->IsInside(bz->GetAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ||
	      (useCorrected && usePlastic && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(BR_AoQ),fp[fpNr(7)]->GetPlastic()->GetCharge()) && OutCut[o]->IsInside(bz->GetCorrAQ(ZD_AoQ),bz->GetZ(ZD_AoQ))) ){
	    double edc = hit->GetDCEnergy(beta[o],0,0,0);
	    g_egamtgamABdc_c[o]->Fill(hit->GetTime(),edc);
	    // gate on timing
	    if(TimeCut[1]==NULL || TimeCut[1]->IsInside(hit->GetTime(),edc)){
	      g_egamABdc_c[o]->Fill(edc);
	      g_egamABdc_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), edc);
	      g_egamAB_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetEnergy());
	      g_egamABdc_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, edc);
	      g_egamABdc_gmultAB_c[o]->Fill(gr->GetMultAB(), edc);
	      g_egamABdc_tmultAB_c[o]->Fill(tmultAB, edc);
	      g_egamABdc_depth_c[o]->Fill(hit->GetPosition().Mag(), edc);
	      for(int l=g+1; l<gr->GetMultAB(); l++){
		HitCalc* hit2 = gr->GetHitAB(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,0);
		if(TimeCut[1] && !TimeCut[1]->IsInside(hit2->GetTime(),edc2))
		  continue;
		if(edc>edc2)
		  g_egamegamABdc_c[o]->Fill(edc,edc2);
		else
		  g_egamegamABdc_c[o]->Fill(edc2,edc);
	      }//second hit
	    }
	  }// in cut ZeroDeg
	}// loop over cuts
      }// energy and position good
    }// hits
    
    if(i%10000 == 0){
      double time_end = get_time();
      cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*i)/nentries<<" % done\t"<<(Float_t)i/(time_end - time_start)<<" events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i<<"s to go \r"<<flush;
    }

  }
  cout << endl;


  cout << "writing histograms to file" << endl;
  cout << endl;
  cout << "creating outputfile " << OutputFile << endl;
  TFile* ofile = new TFile(OutputFile,"recreate");
  ofile->cd();
  TH1F* h1;
  TH2F* h2;
  TIter next(hlist);
  while( (h1 = (TH1F*)next()) ){
    if(h1->GetEntries()>0)
      h1->Write("",TObject::kOverwrite);
  }
  while( (h2 = (TH2F*)next()) ){
    if(h2->GetEntries()>0)
      h2->Write("",TObject::kOverwrite);
  }
  ofile->Close();
  delete tr;

  double time_end = get_time();
  cout << "Run time " << time_end - time_start << " s." << endl;

  return 0;
}
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec/1000000;
    return d;
}
