#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/time.h>
#include <signal.h>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "CommandLineInterface.hh"
#include "BuildEvents.hh"
#include "RunInfo.hh"

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
  char* BigRIPSFile = NULL;
  char* OutputFile = NULL;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-b", "BigRIPS tree input file", &BigRIPSFile);
  interface->Add("-o", "output file", &OutputFile);
  interface->CheckFlags(argc, argv);

  if(BigRIPSFile == NULL || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }

  TTree* trbigrips = NULL;
  TFile* inbigrips = NULL;

  /// temp KW proc info
  ProcInfo_t pinfo;

  RunInfo *info = NULL;
  if(BigRIPSFile == NULL){
    cout << "No BigRIPS input file given " << endl;
  }
  else{
    cout<<"BigRIPS file: "<<BigRIPSFile<<endl;
    inbigrips = new TFile(BigRIPSFile);
    trbigrips = (TTree*) inbigrips->Get("BRZDtr");
    if(trbigrips == NULL){
      cout << "could not find BigRIPS tree in file " << inbigrips->GetName() << endl;
    }
    //trbigrips->SetCacheSize(1);
    //trbigrips->SetMaxVirtualSize(100);
    info = (RunInfo*) inbigrips->Get("info");
  }
  //local copy for intermediate storing
  unsigned long long int localBRts = 0;
  int localcheckADC = -1;
  int localtrigbit = -1;
  Beam* localbeam = new Beam;
  FocalPlane* localfp[NFPLANES];
  for(unsigned short f=0;f<NFPLANES;f++){
    localfp[f] = new FocalPlane;
  }

  PPAC *localppacs = new PPAC;
  
  trbigrips->SetBranchAddress("timestamp",&localBRts);
  trbigrips->SetBranchAddress("checkADC",&localcheckADC);
  trbigrips->SetBranchAddress("trigbit",&localtrigbit);
  trbigrips->SetBranchAddress("beam",&localbeam);
  trbigrips->SetBranchAddress("ppacs",&localppacs);
  for(unsigned short f=0;f<NFPLANES;f++){
    trbigrips->SetBranchAddress(Form("fp%d",fpID[f]),&localfp[f]);
  }




  unsigned long long int BRts = 0;
  int checkADC = -1;
  int trigbit = -1;
  Beam* beam = new Beam;
  FocalPlane *fp[NFPLANES];
  for(unsigned short f=0;f<NFPLANES;f++){
    fp[f] = new FocalPlane;
  }
  PPAC* ppacs = new PPAC;

  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");
    
  if(outfile->IsZombie()){
    return 4;
  }
  outfile->cd();
  TTree* mtr = new TTree("tr","merged tree");
  //mtr->SetCacheSize(1);
  //mtr->SetMaxVirtualSize(100);
  mtr->Branch("beam",&beam,320000);
  for(unsigned short f=0;f<NFPLANES;f++){
    mtr->Branch(Form("fp%d",fpID[f]),&fp[f],320000);
  }
  mtr->Branch("ppacs",&ppacs,320000);
  mtr->Branch("checkADC",&checkADC,320000);
  mtr->Branch("trigbit",&trigbit,320000);
  mtr->Branch("brTS",&BRts,320000);
  mtr->BranchRef();
  //mtr->SetAutoFlush( 0 );
  mtr->SetAutoFlush(0);
  mtr->SetCacheSize(0);
  mtr->StopCacheLearningPhase();


  
  
  int nentries = trbigrips->GetEntries();
  cout << nentries << " entries in BigRIPS tree" << endl;

  double time_last = get_time();


  int ctr=0;
  int fnbytes = 0;
  for(int i=0;i<nentries;i++){

    localcheckADC = -1;
    localtrigbit = -1;
    localbeam->Clear();
    for(unsigned short f=0;f<NFPLANES;f++){
      localfp[f]->Clear();
    }
    localppacs->Clear();
    localBRts = 0;
    


    Int_t status = trbigrips->GetEvent(i);
    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<trbigrips->GetName()<<endl;
      return false;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<trbigrips->GetName()<<" in file doesn't exist"<<endl;
      return false;
    }
    fnbytes += status;


    BRts = localBRts;
    checkADC = localcheckADC;
    trigbit = localtrigbit;
    //beam = (Beam*)localbeam->Clone();
    for(unsigned short f=0;f<NFPLANES;f++){
    	fp[f] = (FocalPlane*)localfp[f]->Clone(Form("fp%d",f));
    }
    //ppacs = (PPAC*)localppacs->Clone();
    
    mtr->Fill();
    
    //beam->Clear();
    //delete beam;
    
    if(i%10000 == 0){
      gSystem->GetProcInfo(&pinfo);
      cout << i << " memory usage " << pinfo.fMemResident << " " << pinfo.fMemVirtual << endl;
    }
    if(signal_received){
      break;
    }
  }
  mtr->Write("",TObject::kOverwrite);
  outfile->Close();

  double time_end = get_time();
  cout << "Program Run time: " << time_end - time_start << " s." << endl;
  timer.Stop();
  cout << "CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;
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
