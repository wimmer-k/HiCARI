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
  char* BigRIPSFile;
  char* HiCARIFile;
  char* OutputFile = NULL;
  char* SetFile = NULL;
  int LastEvent = -1;
  int vl = 0;
  int mode = -1;

  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-b", "BigRIPS tree input file", &BigRIPSFile);
  interface->Add("-h", "HiCARI tree input file", &HiCARIFile);
  interface->Add("-o", "output file", &OutputFile);
  interface->Add("-s", "settings file", &SetFile);
  interface->Add("-le", "last event to be read", &LastEvent);
  interface->Add("-v", "verbose", &vl);
  interface->Add("-m", "coincidence mode", &mode);
  interface->CheckFlags(argc, argv);

  if((BigRIPSFile == NULL && HiCARIFile == NULL) || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }

  TTree* trbigrips = NULL;
  TTree* trhicari = NULL;
  TFile* inbigrips = NULL;
  TFile* inhicari = NULL;
  int BRrunnr = -1;
  int HIrunnr = -1;

  RunInfo *info = NULL;
  if(BigRIPSFile == NULL){
    cout << "No BigRIPS input file given " << endl;
  }
  else{
    inbigrips = new TFile(BigRIPSFile);
    trbigrips = (TTree*) inbigrips->Get("BRZDtr");
    if(trbigrips == NULL){
      cout << "could not find BigRIPS tree in file " << inbigrips->GetName() << endl;
    }
    info = (RunInfo*) inbigrips->Get("info");
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
    if(info==NULL)
      info = (RunInfo*) inhicari->Get("info");
    else
      info->AppendInfo((RunInfo*) inhicari->Get("info"));
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
  set->Write("settings");
  BuildEvents* evts = new BuildEvents(set);
  evts->SetVerbose(vl);
  if(mode>-1)
    evts->SetCorrelationMode(mode);
  evts->Init(trbigrips,trhicari);
  evts->SetLastEvent(LastEvent);

  double time_last = get_time();
  evts->ReadEach();
  int ctr=0;
  int total = evts->GetNEvents();

  while(evts->Merge()){
    if(ctr%10000 == 0){
      double time_end = get_time();
      double r = evts->GetHistos()->GetCorrRate();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*ctr)/total<<" % done\t" << 
	"correlation rate "; 
      if(r>90)
	cout << GREEN;                                                              
      else if(r>70)
	cout << ORANGE;                                                             
      else
	cout << RED;                                                                
      cout << r << DEFCOLOR << " %\t";                                        

      cout <<(Float_t)ctr/(time_end - time_start) << " events/s (average) " <<
	10000./(time_end - time_last) << " events/s (current)\t" <<
	(total-ctr)*(time_end - time_start)/(Float_t)ctr << "s to go \r" << flush;
      time_last = time_end;
    }
    if(signal_received){
      break;
    }
    ctr++;
  }
  info->SetMergedEvents(ctr);
  evts->CloseEvent();
  evts->GetTree()->Write("",TObject::kOverwrite);
  evts->GetHistos()->Write();
  double r = evts->GetHistos()->GetCorrRate();
  info->SetCorrelationRate(r);
  cout << endl << "Total correlation rate from checkADC ";
  if(r>90)
    cout << GREEN;
  else if(r>70)
    cout << ORANGE;
  else
    cout << RED;
  cout << r << DEFCOLOR << " %" << endl;
  info->Write("info");
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
