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
  char* SettingFile;
  int nmax =0;
  int vl =0;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-v", "verbose", &vl);
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
  if(SettingFile == NULL){
    cout << "Please give a settings file to use" << endl;
    return 5;
  }


  int incuts = 0;
  int outcuts = 0;
  
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
    for(int i=0;i<incuts;i++){
      char* iname = (char*)set->GetValue(Form("InCut.%d",i),"file");
      TCutG* icg = (TCutG*)cFile->Get(iname);
      cout << "incoming cut found "<< iname << endl;
      InCut.push_back(icg);
    }
    for(int i=0;i<outcuts;i++){
      char* oname = (char*)set->GetValue(Form("OutCut.%d",i),"file");
      TCutG* ocg = (TCutG*)cFile->Get(oname);
      cout << "outgoing cut found "<< oname << endl;
      OutCut.push_back(ocg);
      double b = set->GetValue(Form("Beta.%d",i),0.57);
      beta.push_back(b);
      int s = set->GetValue(Form("InCutSel.%d",i),0);
      InCut_sel.push_back(s);
      
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

  
  vector<TH2F*> zerodeg_c;
  vector<TH1F*> h_egamdc_c;
  vector<TH1F*> g_egamdc_c;
  vector<TH1F*> h_egamABdc_c;
  vector<TH1F*> g_egamABdc_c;
  for(int i=0;i<incuts;i++){
    TH2F* h = new TH2F(Form("zerodeg_%s",InCut[i]->GetName()),
		       Form("zerodeg_%s",InCut[i]->GetName()),
		       1000,2.2,2.8,1000,20,40);
    zerodeg_c.push_back(h);
    hlist->Add(zerodeg_c.back());
  }
  
  for(int o=0;o<outcuts;o++){
    TH1F* h = new TH1F(Form("h_egamdc_%s_%s",InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
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
    for(int i=0;i<incuts;i++){
      if(InCut[i]->IsInside(bz->GetAQ(2),bz->GetZ(2))){
	zerodeg_c[i]->Fill(bz->GetAQ(5),bz->GetZ(5));
      }
    }
    for(int h=0; h<hi->GetMult(); h++){
      HiCARIHitCalc* hit = hi->GetHit(h);
      if(hit->IsBigRIPS()){
	continue;
      }
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	h_egamdc->Fill(hit->GetDCEnergy());
	for(int o=0;o<outcuts;o++){	
	  if(InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))){
	    h_egamdc_c[o]->Fill(hit->GetDCEnergy(beta[o]));
	  }
	}
      }
    }
    for(int g=0; g<gr->GetMult(); g++){
      HitCalc* hit = gr->GetHit(g);
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	g_egamdc->Fill(hit->GetDCEnergy());
	for(int o=0;o<outcuts;o++){	
	  if(InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))){
	    g_egamdc_c[o]->Fill(hit->GetDCEnergy(beta[o]));
	  }
	}
      }
    }
    for(int h=0; h<hi->GetMultAB(); h++){
      HiCARIHitCalc* hit = hi->GetHit(h);
      if(hit->IsBigRIPS()){
	continue;
      }
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	h_egamABdc->Fill(hit->GetDCEnergy());
	for(int o=0;o<outcuts;o++){	
	  if(InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))){
	    h_egamABdc_c[o]->Fill(hit->GetDCEnergy(beta[o]));
	  }
	}
      }
    }
    for(int g=0; g<gr->GetMultAB(); g++){
      HitCalc* hit = gr->GetHit(g);
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	g_egamABdc->Fill(hit->GetDCEnergy());
	for(int o=0;o<outcuts;o++){	
	  if(InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))){
	    g_egamABdc_c[o]->Fill(hit->GetDCEnergy(beta[o]));
	  }
	}
      }
    }
    
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
