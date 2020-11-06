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

  if(SettingFile != NULL){
    TEnv* set = new TEnv(SettingFile);
    char* cfilename = (char*)set->GetValue("CutFile","file");
    incuts = set->GetValue("NInCuts",0);
    outcuts = set->GetValue("NOutCuts",0);
    
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
  Double_t nentries = tr->GetEntries();


  TList *hlist = new TList();
  TH1F* trigger = new TH1F("trigger","trigger",10,0,10);hlist->Add(trigger);
  TH2F* bigrips = new TH2F("bigrips","bigrips",1000,2.2,2.8,1000,20,40);hlist->Add(bigrips);
  TH2F* zerodeg = new TH2F("zerodeg","zerodeg",1000,2.2,2.8,1000,20,40);hlist->Add(zerodeg);
  TH1F* h_egamdc = new TH1F("h_egamdc","h_egamdc",4000,0,4000);hlist->Add(h_egamdc);
  TH1F* g_egamdc = new TH1F("g_egamdc","g_egamdc",4000,0,4000);hlist->Add(g_egamdc);
  TH1F* h_egamABdc = new TH1F("h_egamABdc","h_egamABdc",4000,0,4000);hlist->Add(h_egamABdc);
  TH1F* g_egamABdc = new TH1F("g_egamABdc","g_egamABdc",4000,0,4000);hlist->Add(g_egamABdc);
  TH2F* h_egamtgam = new TH2F("h_egamtgam","h_egamtgam",1000,-500,500,4000,0,4000);hlist->Add(h_egamtgam);
  TH2F* g_egamtgam = new TH2F("g_egamtgam","g_egamtgam",1000,-500,500,4000,0,4000);hlist->Add(g_egamtgam);
  TH2F* h_egamtgamAB = new TH2F("h_egamtgamAB","h_egamtgamAB",1000,-500,500,4000,0,4000);hlist->Add(h_egamtgamAB);
  TH2F* g_egamtgamAB = new TH2F("g_egamtgamAB","g_egamtgamAB",1000,-500,500,4000,0,4000);hlist->Add(g_egamtgamAB);
  TH2F* h_egamtgamdc = new TH2F("h_egamtgamdc","h_egamtgamdc",1000,-500,500,4000,0,4000);hlist->Add(h_egamtgamdc);
  TH2F* g_egamtgamdc = new TH2F("g_egamtgamdc","g_egamtgamdc",1000,-500,500,4000,0,4000);hlist->Add(g_egamtgamdc);
  TH2F* h_egamtgamABdc = new TH2F("h_egamtgamABdc","h_egamtgamABdc",1000,-500,500,4000,0,4000);hlist->Add(h_egamtgamABdc);
  TH2F* g_egamtgamABdc = new TH2F("g_egamtgamABdc","g_egamtgamABdc",1000,-500,500,4000,0,4000);hlist->Add(g_egamtgamABdc);

  
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
  
  vector<TH2F*> h_egam_theta_c;
  vector<TH2F*> h_egamdc_theta_c;
  vector<TH2F*> h_egamAB_theta_c;
  vector<TH2F*> h_egamABdc_theta_c;
  
  vector<TH2F*> h_egamtgamdc_c;
  vector<TH2F*> g_egamtgamdc_c;
  vector<TH2F*> h_egamtgamABdc_c;
  vector<TH2F*> g_egamtgamABdc_c;
  vector<TH2F*> h_egamegamdc_c;
  vector<TH2F*> g_egamegamdc_c;
  vector<TH2F*> h_egamegamABdc_c;
  vector<TH2F*> g_egamegamABdc_c;
  for(int i=0;i<incuts;i++){
    TH2F* h = new TH2F(Form("zerodeg_%s",InCut[i]->GetName()),
		       Form("zerodeg_%s",InCut[i]->GetName()),
		       1000,2.2,2.8,1000,20,40);
    zerodeg_c.push_back(h);
    hlist->Add(zerodeg_c.back());
  }
  
  for(int o=0;o<outcuts;o++){
    TH1F* h;
    
    h = new TH1F(Form("h_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 Form("h_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 4000,0,4000);
    h_egamdc_c.push_back(h);
    hlist->Add(h_egamdc_c.back());
    h = new TH1F(Form("g_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 Form("g_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 4000,0,4000);
    g_egamdc_c.push_back(h);
    hlist->Add(g_egamdc_c.back());
    h = new TH1F(Form("h_egamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 Form("h_egamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 4000,0,4000);
    h_egamABdc_c.push_back(h);
    hlist->Add(h_egamABdc_c.back());
    h = new TH1F(Form("g_egamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 Form("g_egamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		 4000,0,4000);
    g_egamABdc_c.push_back(h);
    hlist->Add(g_egamABdc_c.back());


    TH2F* h2;
    //summary plots
    h2 = new TH2F(Form("h_egam_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egam_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,4000,0,4000);
    h_egam_summary_c.push_back(h2);
    hlist->Add(h_egam_summary_c.back());
    
    h2 = new TH2F(Form("h_tgam_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_tgam_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,1000,-500,500);
    h_tgam_summary_c.push_back(h2);
    hlist->Add(h_tgam_summary_c.back());
    
    h2 = new TH2F(Form("h_egamdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,4000,0,4000);
    h_egamdc_summary_c.push_back(h2);
    hlist->Add(h_egamdc_summary_c.back());

    h2 = new TH2F(Form("h_egamAB_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamAB_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,4000,0,4000);
    h_egamAB_summary_c.push_back(h2);
    hlist->Add(h_egamAB_summary_c.back());
    
    h2 = new TH2F(Form("h_tgamAB_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_tgamAB_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,1000,-500,500);
    h_tgamAB_summary_c.push_back(h2);
    hlist->Add(h_tgamAB_summary_c.back());
    
    h2 = new TH2F(Form("h_egamABdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamABdc_summary_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  48,0,48,4000,0,4000);
    h_egamABdc_summary_c.push_back(h2);
    hlist->Add(h_egamABdc_summary_c.back());


    //verus theta 
    h2 = new TH2F(Form("h_egam_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egam_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,4000,0,4000);
    h_egam_theta_c.push_back(h2);
    hlist->Add(h_egam_theta_c.back());
    
    h2 = new TH2F(Form("h_egamdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,4000,0,4000);
    h_egamdc_theta_c.push_back(h2);
    hlist->Add(h_egamdc_theta_c.back());

    h2 = new TH2F(Form("h_egamAB_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamAB_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,4000,0,4000);
    h_egamAB_theta_c.push_back(h2);
    hlist->Add(h_egamAB_theta_c.back());
    
    h2 = new TH2F(Form("h_egamABdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamABdc_theta_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  180,0,180,4000,0,4000);
    h_egamABdc_theta_c.push_back(h2);
    hlist->Add(h_egamABdc_theta_c.back());



    
    //e gamma versus time diff to BR
    h2 = new TH2F(Form("h_egamtgamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamtgamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,-500,500,4000,0,4000);
    h_egamtgamdc_c.push_back(h2);
    hlist->Add(h_egamtgamdc_c.back());

    h2 = new TH2F(Form("g_egamtgamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamtgamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,-500,500,4000,0,4000);
    g_egamtgamdc_c.push_back(h2);
    hlist->Add(g_egamtgamdc_c.back());

    h2 = new TH2F(Form("h_egamtgamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamtgamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,-500,500,4000,0,4000);
    h_egamtgamABdc_c.push_back(h2);
    hlist->Add(h_egamtgamABdc_c.back());

    h2 = new TH2F(Form("g_egamtgamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamtgamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,-500,500,4000,0,4000);
    g_egamtgamABdc_c.push_back(h2);
    hlist->Add(g_egamtgamABdc_c.back());


    //e gamma versus e gamma
    h2 = new TH2F(Form("h_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,0,4000,1000,0,4000);
    h_egamegamdc_c.push_back(h2);
    hlist->Add(h_egamegamdc_c.back());

    h2 = new TH2F(Form("g_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamegamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,0,4000,1000,0,4000);
    g_egamegamdc_c.push_back(h2);
    hlist->Add(g_egamegamdc_c.back());

    h2 = new TH2F(Form("h_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("h_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,0,4000,1000,0,4000);
    h_egamegamABdc_c.push_back(h2);
    hlist->Add(h_egamegamABdc_c.back());

    h2 = new TH2F(Form("g_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  Form("g_egamegamABdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		  1000,0,4000,1000,0,4000);
   g_egamegamABdc_c.push_back(h2);
    hlist->Add(g_egamegamABdc_c.back());


    
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
    bigrips->Fill(bz->GetAQ(2),bz->GetZ(2));
    zerodeg->Fill(bz->GetAQ(5),bz->GetZ(5));
    for(int in=0;in<incuts;in++){
      if(InCut[in]->IsInside(bz->GetAQ(2),bz->GetZ(2))){
	zerodeg_c[in]->Fill(bz->GetAQ(5),bz->GetZ(5));
      }
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
	  if(InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))){
	    double edc = hit->GetDCEnergy(beta[o],0,0,targetz);
	    h_tgam_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetTime());
	    h_egamtgamdc_c[o]->Fill(hit->GetTime(),edc);

	    // gate on timing
	    if(TimeCut[0] && TimeCut[0]->IsInside(hit->GetTime(),edc)){
	      h_egamdc_c[o]->Fill(edc);
	      h_egam_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetEnergy());
	      h_egamdc_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), edc);
	      h_egam_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetEnergy());
	      h_egamdc_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, edc);

	      for(int l=h+1; l<hi->GetMult(); l++){
		HiCARIHitCalc* hit2 = hi->GetHit(l);
		if(hit2->IsBigRIPS()){
		  continue;
		}
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
		if(edc>edc2)
		  h_egamegamdc_c[o]->Fill(edc,edc2);
		else
		  h_egamegamdc_c[o]->Fill(edc2,edc);
	      }//second hit
	    }
	    // if not time cut defined take all
	    else if(TimeCut[0]==NULL){
	      h_egamdc_c[o]->Fill(edc);
	      h_egam_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetEnergy());
	      h_egamdc_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), edc);
	      h_egam_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetEnergy());
	      h_egamdc_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, edc);
	      
	      for(int l=h+1; l<hi->GetMult(); l++){
		HiCARIHitCalc* hit2 = hi->GetHit(l);
		if(hit2->IsBigRIPS()){
		  continue;
		}
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
		if(edc>edc2)
		  h_egamegamdc_c[o]->Fill(edc,edc2);
		else
		  h_egamegamdc_c[o]->Fill(edc2,edc);
	      }//second hit
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
	  if(InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))){
	    double edc = hit->GetDCEnergy(beta[o],0,0,targetz);
	    g_egamtgamdc_c[o]->Fill(hit->GetTime(),edc);
	    
	    // gate on timing
	    if(TimeCut[1] && TimeCut[1]->IsInside(hit->GetTime(),edc)){
	      g_egamdc_c[o]->Fill(edc);
	      for(int l=g+1; l<gr->GetMult(); l++){
		HitCalc* hit2 = gr->GetHit(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
		if(edc>edc2)
		  g_egamegamdc_c[o]->Fill(edc,edc2);
		else
		  g_egamegamdc_c[o]->Fill(edc2,edc);
	      }//second hit

	    }
	    // if not time cut defined take all
	    else if(TimeCut[1]==NULL){
	      g_egamdc_c[o]->Fill(edc);
	      for(int l=g+1; l<gr->GetMult(); l++){
		HitCalc* hit2 = gr->GetHit(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
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
	  if(InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))){
	    double edc = hit->GetDCEnergy(beta[o],0,0,targetz);
	    h_tgamAB_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetTime());
	    h_egamtgamABdc_c[o]->Fill(hit->GetTime(),edc);
	    // gate on timing
	    if(TimeCut[0] && TimeCut[0]->IsInside(hit->GetTime(),edc)){
	      h_egamABdc_c[o]->Fill(edc);
	      h_egamAB_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetEnergy());
	      h_egamABdc_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), edc);
	      h_egamAB_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetEnergy());
	      h_egamABdc_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, edc);

	      for(int l=h+1; l<hi->GetMultAB(); l++){
		HiCARIHitCalc* hit2 = hi->GetHitAB(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
		if(edc>edc2)
		  h_egamegamABdc_c[o]->Fill(edc,edc2);
		else
		  h_egamegamABdc_c[o]->Fill(edc2,edc);
	      }//second hit
	    }
	    // if not time cut defined take all
	    else if(TimeCut[0]==NULL){
	      h_egamABdc_c[o]->Fill(edc);
	      h_egamAB_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), hit->GetEnergy());
	      h_egamABdc_summary_c[o]->Fill(hit->GetCluster()*4+hit->GetCrystal(), edc);
	      h_egamAB_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetEnergy());
	      h_egamABdc_theta_c[o]->Fill(hit->GetPosition().Theta()*180./3.1415, edc);
	      for(int l=h+1; l<hi->GetMultAB(); l++){
		HiCARIHitCalc* hit2 = hi->GetHitAB(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
		if(edc>edc2)
		  h_egamegamABdc_c[o]->Fill(edc,edc2);
		else
		  h_egamegamABdc_c[o]->Fill(edc2,edc);
	      }//second hit
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
	  if(InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))){
	    double edc = hit->GetDCEnergy(beta[o],0,0,targetz);
	    g_egamtgamABdc_c[o]->Fill(hit->GetTime(),edc);
	    // gate on timing
	    if(TimeCut[1] && TimeCut[1]->IsInside(hit->GetTime(),edc)){
	      g_egamABdc_c[o]->Fill(edc);
	      for(int l=g+1; l<gr->GetMultAB(); l++){
		HitCalc* hit2 = gr->GetHitAB(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
		if(edc>edc2)
		  g_egamegamABdc_c[o]->Fill(edc,edc2);
		else
		  g_egamegamABdc_c[o]->Fill(edc2,edc);
	      }//second hit
	    }
	    // if not time cut defined take all
	    else if(TimeCut[1]==NULL){
	      g_egamABdc_c[o]->Fill(edc);
	      for(int l=g+1; l<gr->GetMultAB(); l++){
		HitCalc* hit2 = gr->GetHitAB(l);
		double edc2 = hit2->GetDCEnergy(beta[o],0,0,targetz);
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
