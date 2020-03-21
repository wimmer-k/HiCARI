#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/time.h>
#include <signal.h>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "CommandLineInterface.hh"
#include "BuildEvents.hh"

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
  char* BigRIPSFile;
  char* HiCARIFile;
  char* OutputFile = NULL;
  char* SetFile = NULL;
  int LastEvent = -1;
  int vl = 0;

  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-b", "BigRIPS tree input file", &BigRIPSFile);
  interface->Add("-h", "HiCARI tree input file", &HiCARIFile);
  interface->Add("-o", "output file", &OutputFile);
  interface->Add("-s", "settings file", &SetFile);
  interface->Add("-le", "last event to be read", &LastEvent);
  interface->Add("-v", "verbose", &vl);
  interface->CheckFlags(argc, argv);

  if((BigRIPSFile == NULL && HiCARIFile == NULL) || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }

  TTree* trbigrips = NULL;
  TTree* trhicari = NULL;
  TFile* inbigrips = NULL;
  TFile* inhicari = NULL;

  if(BigRIPSFile == NULL){
    cout << "No BigRIPS input file given " << endl;
  }
  else{
    inbigrips = new TFile(BigRIPSFile);
    trbigrips = (TTree*) inbigrips->Get("BRZDtr");
    if(trbigrips == NULL){
      cout << "could not find BigRIPS tree in file " << inbigrips->GetName() << endl;
    }
  }
  if(HiCARIFile == NULL){
    cout << "No HiCARI input file given " << endl;
  }
  else{
    inhicari = new TFile(HiCARIFile);
    trhicari = (TTree*) inhicari->Get("calib");
    if(trhicari == NULL){
      cout << "could not find HiCARI tree in file " << inhicari->GetName() << endl;
    }
  }

  cout<<"BigRIPS file: "<<BigRIPSFile<<endl;
  cout<<"HiCARI file: "<<HiCARIFile<<endl;
  cout<<"output file: "<<OutputFile<< endl;
  if(SetFile == NULL)
    cerr<<"No settings file! Using standard values"<<endl;
  else
    cout<<"settings file:"<<SetFile<<endl;

  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");
    
  if(outfile->IsZombie()){
    return 4;
  }
  outfile->cd();
  Settings* set = new Settings(SetFile);

  BuildEvents* evts = new BuildEvents();
  evts->SetVerbose(vl);
  evts->SetWindow(set->EventTimeDiff());
  evts->SetCoincMode(0); // for p-gamma coinc  
  evts->Init(trbigrips,trhicari);
  evts->SetLastEvent(LastEvent);

  double time_last = get_time();
  evts->ReadEach();
  int ctr=0;
  int total = evts->GetNEvents();
  
  while(evts->Merge()){
    if(ctr%10000 == 0){
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*ctr)/total<<" % done\t" << 
	(Float_t)ctr/(time_end - time_start) << " events/s (average) " <<
	10000./(time_end - time_last) << " events/s (current) " <<
	(total-ctr)*(time_end - time_start)/(Float_t)ctr << "s to go \r" << flush;
      time_last = time_end;
    }
    if(signal_received){
      break;
    }
    ctr++;
  }
  evts->CloseEvent();
  evts->GetTree()->Write("",TObject::kOverwrite);
  outfile->Close();
  if(inbigrips!=NULL)
    inbigrips->Close();
  if(inhicari!=NULL)
    inhicari->Close();
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
