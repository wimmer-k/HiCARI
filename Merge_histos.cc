#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>
#include <signal.h>


#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCutG.h"
#include "TEnv.h"
#include "TKey.h"

#include "CommandLineInterface.hh"
#include "MergeHistograms.hh"
using namespace TMath;
using namespace std;

bool signal_received = false;
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
int main(int argc, char* argv[]){
  signal(SIGINT,signalhandler);

  double time_start = get_time();
  vector<char*> InputFiles;
  char* OutputFile = NULL;
  //vector<char*> SettingFile;
  int nmax =0;
  int vl =0;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
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
  int checkADC;
  tr->SetBranchAddress("checkADC",&checkADC);
  Beam* beam = new Beam;
  tr->SetBranchAddress("beam",&beam);
  FocalPlane* fp[NFPLANES];
  for(unsigned short f=0;f<NFPLANES;f++){
    fp[f] = new FocalPlane;
    tr->SetBranchAddress(Form("fp%d",fpID[f]),&fp[f]);
  }
  PPAC* ppacs = new PPAC;
  tr->SetBranchAddress("ppacs",&ppacs);
  unsigned long long int brTS;
  tr->SetBranchAddress("brTS",&brTS);
  int brentry;
  tr->SetBranchAddress("brentry",&brentry);

  HiCARICalc* hi = new HiCARICalc;
  tr->SetBranchAddress("hicari",&hi);
  unsigned long long int hiTS;
  tr->SetBranchAddress("hiTS",&hiTS);
  int hientry;
  tr->SetBranchAddress("hientry",&hientry);
  

  GretinaCalc* m2 = new GretinaCalc;
  tr->SetBranchAddress("mode2",&m2);
  unsigned long long int m2TS;
  tr->SetBranchAddress("m2TS",&m2TS);
  int m2entry;
  tr->SetBranchAddress("m2entry",&m2entry);

  Double_t nentries = tr->GetEntries();

  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");

  if(outfile->IsZombie()){
    return 4;
  }

  Int_t nbytes = 0;
  Int_t status;

  // if(SettingFile.size() == 0){
  //   cout << "Please give a settings file to use" << endl;
  //   return 5;
  // }
  // MergeHistograms* hists = new MergeHistograms(new Settings(SettingFile),nentries);
  TFile *firstinfile = new TFile(InputFiles[0]);
  Settings* set = (Settings*)firstinfile->Get("settings");
  MergeHistograms* hists = new MergeHistograms(set,nentries);

  cout << nentries << " entries in tree " << endl;

  if(nmax>0)
    nentries = nmax;

  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    //trigbit = 0;
    checkADC = -1;
    brTS = 0;
    brentry = 0;
    for(int f=0;f<NFPLANES;f++){
      fp[f]->Clear();
    }
    ppacs->Clear();
    beam->Clear();

    hi->Clear();
    hiTS = 0;
    hientry = 0;
    
    m2->Clear();
    m2TS = 0;
    m2entry = 0;
    

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

    hists->FillHistograms(checkADC,hi,m2,brTS,hiTS,m2TS);

    if(i%10000 == 0){
      double time_end = get_time();
      cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*i)/nentries<<" % done\t"<<(Float_t)i/(time_end - time_start)<<" events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i<<"s to go \r"<<flush;
    }

  }
  cout << endl;


  cout << "writing histograms to file" << endl;
  outfile->cd();
  hists->Write(-1,-1);
  outfile->Close();
  delete tr;

  double time_end = get_time();
  cout << "Run time " << time_end - time_start << " s." << endl;

  return 0;
}
