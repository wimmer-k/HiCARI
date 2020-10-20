#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>
#include <signal.h>


#include "TFile.h"
#include "TTree.h"

#include "CommandLineInterface.hh"
#include "HiCARI.hh"
#include "Calibration.hh"
#include "CalHistograms.hh"
#include "RunInfo.hh"

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
  char* InputFile = NULL;
  char* OutputFile = NULL;
  vector<char*> SettingFile;
  int nmax = 0;
  int vl = 0;
  int whist = 0;
  int wtree = 1;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfile", &InputFile);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-h", "write histos", &whist);
  interface->Add("-t", "write tree", &wtree);
  interface->Add("-v", "verbose", &vl);
  interface->CheckFlags(argc, argv);

  if(InputFile == NULL || OutputFile == NULL){
    cerr<<"You have to provide at least the input file and the output file!"<<endl;
    exit(1);
  }
  cout<<"input file: " << InputFile <<endl;
  cout<<"output file: "<<OutputFile<< endl;

  TFile* infile = new TFile(InputFile);
  TTree* tr = (TTree*) infile->Get("build");

  if(tr == NULL){
    cout << "could not find tree build in file " <<InputFile<<endl;
    return 3;
  }
  RunInfo *info = (RunInfo*)infile->Get("info");
  cout << "calibrating data from run " << GREEN << info->GetHIRunNumber() << DEFCOLOR << endl;

  HiCARI* hi = new HiCARI;
  tr->SetBranchAddress("hicari",&hi);

  Double_t nentries = tr->GetEntries();
  Int_t ncalentries = 0;
  cout << nentries << " entries in tree " << endl;

  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");
  cout << "setting up calibrated tree " << endl;
  HiCARICalc* hicalc = new HiCARICalc;
  TTree* caltr = new TTree("calib","calibrated and built events");
  caltr->Branch("hicaricalc",&hicalc, 320000);
  caltr->BranchRef();

  if(outfile->IsZombie()){
    return 4;
  }

  unsigned long long nbytes = 0;
  Int_t status;

  if(SettingFile.size() == 0){
    cout << "Please give a settings file to use" << endl;
    return 5;
  }
  //Read the settings
  Settings* set = new Settings(SettingFile);

  Calibration *cal = new Calibration(set);
  CalHistograms* chist = new CalHistograms(set);

  if(nmax>0)
    nentries = nmax;

  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    hi->Clear();
    hicalc->Clear();

    if(vl>2)
      cout << "getting entry " << i << endl;
    status = tr->GetEvent(i);
    if(vl>2)
      cout << "status " << status << endl;
    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<tr->GetName()<<" in file "<<InputFile<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<tr->GetName()<<" in file "<<InputFile<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;

    cal->BuildHiCARICalc(hi,hicalc);
    if(hicalc->GetMult()>0||hicalc->HadBigRIPS()){
      if(wtree)
	caltr->Fill();
      ncalentries++;
    }
    if(whist){
      chist->FillHistograms(hicalc);
    }

    if(i%10000 == 0){
      double time_end = get_time();
      cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*i)/nentries<<" % done\t"<<(Float_t)i/(time_end - time_start)<<" events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i<<"s to go\r"<<flush;
    }

  }
  cout << endl;
  cout << "Total of " << BLUE << nentries << DEFCOLOR << " entries ("<< BLUE <<nbytes/(1024*1024) << DEFCOLOR << " MB unzipped) read." << endl;
  cout << "Total of " << BLUE << ncalentries << DEFCOLOR << " calibrated events ("<< BLUE << caltr->GetZipBytes()/(1024*1024)<< DEFCOLOR << " MB) written."  << endl;
  if(wtree)
    caltr->Write();
  info->SetHICalEvents(ncalentries);
  info->SetBigRIPSCtr(cal->GetBigRIPSCtr());
  info->SetBigRIPSHitCtr(cal->GetBigRIPSHitCtr());
  info->SetHiCARICtr(cal->GetHiCARICtr());
  info->SetHiCARIHitCtr(cal->GetHiCARIHitCtr());
  if(whist){
    cout << "writing histograms to file" << endl;
    chist->Write();
  }
  info->Write("info");
  outfile->Close();
  delete tr;

  double time_end = get_time();
  cout << "Run time " << time_end - time_start << " s." << endl;

  return 0;
}
