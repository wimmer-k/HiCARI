#include <iostream>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <signal.h>

#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TStopwatch.h"

#include "CommandLineInterface.hh"
#include "RunInfo.hh"
#include "Beam.hh"
#include "FocalPlane.hh"

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
  int LastEvent = -1;
  int Verbose = 0;
  char* InputFile = NULL;
  char* OutFile = NULL;
  char* SettingFile = NULL;
  int br = 2;
  int zd = 5;
  char* tname = "tr";
  
  // Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "input file", &InputFile);
  interface->Add("-o", "output file", &OutFile);    
  interface->Add("-s", "settings file", &SettingFile);
  interface->Add("-le", "last event to be read", &LastEvent);  
  interface->Add("-v", "verbose level", &Verbose);
  interface->Add("-t", "name of the tree, default \"tr\"", &tname);
  interface->CheckFlags(argc, argv);

  // Complain about missing mandatory arguments
  if(InputFile == NULL){
    cout << "No input file given " << endl;
    return 1;
  }
  if(OutFile == NULL){
    cout << "No output ROOT file given " << endl;
    return 2;
  }
  TFile* infile = new TFile(InputFile);
  TTree* tr = (TTree*)infile->Get(tname);

  if(tr == NULL){
    cout << "could not find tree tr in file " << infile->GetName() << endl;
    return 4;
  }

  //check what the contents of the input file were
  Settings* inset = (Settings*)infile->Get("settings");
  
  HiCARICalc* hi = new HiCARICalc;
  tr->SetBranchAddress("hicari",&hi);
  unsigned long long int hits;
  tr->SetBranchAddress("hiTS",&hits);
  unsigned int hientry;
  tr->SetBranchAddress("hientry",&hientry);
 
  GretinaCalc* m2 = new GretinaCalc;
  tr->SetBranchAddress("mode2",&m2);
  unsigned long long int m2ts;
  tr->SetBranchAddress("m2TS",&m2ts);
  unsigned int m2entry;
  tr->SetBranchAddress("m2entry",&m2entry);

  int trigbit = -1;
  tr->SetBranchAddress("trigbit",&trigbit);
  unsigned long long int brts;
  tr->SetBranchAddress("brTS",&brts);
  unsigned int brentry;
  tr->SetBranchAddress("brentry",&brentry);
   
  
  Beam* beam = new Beam;
  tr->SetBranchAddress("beam",&beam);
  FocalPlane* fp[NFPLANES];
  for(unsigned short f=0;f<NFPLANES;f++){
    fp[f] = new FocalPlane;
  }

  if(inset->BigRIPSDetail()>0){
    for(unsigned short f=0;f<NFPLANES;f++){
      tr->SetBranchAddress(Form("fp%d",fpID[f]),&fp[f]);
    }
  }
  PPAC* ppacs = new PPAC;
  if(inset->BigRIPSDetail()>1)
    tr->SetBranchAddress("ppacs",&ppacs);
  
  int incuts = 0;
  int outcuts = 0;
  
  vector<TCutG*> InCut;
  vector<TCutG*> OutCut;
  vector<int> InCut_sel;
  vector<TTree*> InTree;
  vector<TTree*> OutTree;
 
  int useCorrected = 0;
  int BR_AoQ = 2;
  int ZD_AoQ = 5;
  
  if(SettingFile != NULL){
    TEnv* set = new TEnv(SettingFile);
    char* cfilename = (char*)set->GetValue("CutFile","file");
    incuts = set->GetValue("NInCuts",0);
    outcuts = set->GetValue("NOutCuts",0);

    useCorrected = set->GetValue("UseAoQCorr",0);
    BR_AoQ = set->GetValue("UseBRAoQ",2);
    ZD_AoQ = set->GetValue("UseZDAoQ",5);

    cout<<"Using A/q ["<<BR_AoQ<<"] (BigRIPS) and ["<<ZD_AoQ<<"] (ZD) " << endl;
    
    if(useCorrected)
      cout << "using corrected A/q values for PID gates" << endl;
    else
      cout << "using raw A/q values for PID gates" << endl;
    
    TFile* cFile = new TFile(cfilename);

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

  TFile* ofile = new TFile(OutFile,"recreate");
  ofile->cd();
  InTree.resize(InCut.size());
  for(UShort_t in=0;in<InCut.size();in++){ // loop over incoming cuts
    InTree[in] = new TTree(Form("tr_%s",InCut[in]->GetName()),Form("Data Tree with cut on %s",InCut[in]->GetName()));
    InTree[in]->Branch("hicari",&hi,320000);
    InTree[in]->Branch("hiTS",&hits,320000);
    InTree[in]->Branch("hientry",&hientry,320000);
       
    InTree[in]->Branch("mode2",&m2,320000);
    InTree[in]->Branch("m2TS",&m2ts,320000);
    InTree[in]->Branch("m2entry",&m2entry,320000);
       
    InTree[in]->Branch("trigbit",&trigbit,"trigbit/I");
    InTree[in]->Branch("brTS",&brts,320000);
    InTree[in]->Branch("brentry",&brentry,320000);
      
    InTree[in]->Branch("beam",&beam,320000);

    if(inset->BigRIPSDetail()>0){
      for(unsigned short f=0;f<NFPLANES;f++)
	      InTree[in]->Branch(Form("fp%d",fpID[f]),&fp[f],320000);
    }
    if(inset->BigRIPSDetail()>1)
      InTree[in]->Branch("ppacs",&ppacs,320000);

  }
  
  OutTree.resize(OutCut.size());
  for(UShort_t out=0;out<OutCut.size();out++){ // loop over outgoing cuts
    
    OutTree[out] = new TTree(Form("tr_%s_%s",InCut[InCut_sel[out]]->GetName(),OutCut[out]->GetName()),Form("Data Tree with cut on %s and %s",InCut[InCut_sel[out]]->GetName(),OutCut[out]->GetName()));
    OutTree[out]->Branch("hicari",&hi,320000);
    OutTree[out]->Branch("hiTS",&hits,320000);
    OutTree[out]->Branch("hientry",&hientry,320000);
       
    OutTree[out]->Branch("mode2",&m2,320000);
    OutTree[out]->Branch("m2TS",&m2ts,320000);
    OutTree[out]->Branch("m2entry",&m2entry,320000);
       
    OutTree[out]->Branch("trigbit",&trigbit,"trigbit/I");
    OutTree[out]->Branch("brTS",&brts,320000);
    OutTree[out]->Branch("brentry",&brentry,320000);
      
    OutTree[out]->Branch("beam",&beam,320000);

    if(inset->BigRIPSDetail()>0){
      for(unsigned short f=0;f<NFPLANES;f++)
	      OutTree[out]->Branch(Form("fp%d",fpID[f]),&fp[f],320000);
    }
    if(inset->BigRIPSDetail()>1)
      OutTree[out]->Branch("ppacs",&ppacs,320000);
  }
  
    
  
  Double_t nentries = tr->GetEntries();
  cout << nentries << " entries in tree" << endl;
  if(LastEvent>0)
    nentries = LastEvent;
  
  Int_t nbytes = 0;
  Int_t status;
  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    hi->Clear();
    hits = 0;
    hientry = -1;
    m2->Clear();
    m2ts = 0;
    m2entry = -1;
    
    trigbit = 0;
    brts = 0;
    brentry = -1;
    for(int f=0;f<NFPLANES;f++){
      fp[f]->Clear();
    }
    beam->Clear();
    if(inset->BigRIPSDetail()>0){
      for(int f=0;f<NFPLANES;f++){
	fp[f]->Clear();
      }
    }
    if(inset->BigRIPSDetail()>1)
      ppacs->Clear();


    if(Verbose>2)
      cout << "getting entry " << i << endl;
    status = tr->GetEvent(i);
    if(Verbose>2)
      cout << "status " << status << endl;
    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<tr->GetName()<<" in file "<<infile->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<tr->GetName()<<" in file "<<infile->GetName()<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;
    //start analysis
    for(UShort_t in=0;in<InCut.size();in++){ // loop over incoming cuts
      if((!useCorrected && InCut[in]->IsInside(beam->GetAQ(BR_AoQ),beam->GetZ(BR_AoQ))) ||
	 (useCorrected && InCut[in]->IsInside(beam->GetCorrAQ(BR_AoQ),beam->GetZ(BR_AoQ)))){
	InTree[in]->Fill();
      }
    }
    for(UShort_t out=0;out<OutCut.size();out++){ // loop over outcomoutg cuts
      if((!useCorrected && InCut[InCut_sel[out]]->IsInside(beam->GetAQ(BR_AoQ),beam->GetZ(BR_AoQ)) && OutCut[out]->IsInside(beam->GetAQ(ZD_AoQ),beam->GetZ(ZD_AoQ))) ||
	 (useCorrected && InCut[InCut_sel[out]]->IsInside(beam->GetCorrAQ(BR_AoQ),beam->GetZ(BR_AoQ)) && OutCut[out]->IsInside(beam->GetCorrAQ(ZD_AoQ),beam->GetZ(ZD_AoQ)))){
	OutTree[out]->Fill();
      }
    }
    

    if(i%10000 == 0){
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/nentries <<
	" % done\t" << (Float_t)i/(time_end - time_start) << " events/s " << 
	(nentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
    if(i%10000000 == 0){
      for(UShort_t in=0;in<InCut.size();in++) // loop over incoming cuts
	InTree[in]->AutoSave();
      for(UShort_t out=0;out<OutCut.size();out++) // loop over outcomoutg cuts
	OutTree[out]->AutoSave();
    }
  }
  cout << endl;
  cout << "creating outputfile " << endl;
  Long64_t filesize =0;
  for(UShort_t in=0;in<InCut.size();in++){ // loop over incoming cuts
    InTree[in]->Write("",TObject::kOverwrite);
    filesize += InTree[in]->GetZipBytes();
  }
  for(UShort_t out=0;out<OutCut.size();out++){ // loop over outcomoutg cuts
    OutTree[out]->Write("",TObject::kOverwrite);
    filesize += OutTree[out]->GetZipBytes();
  }
  cout<<"Size of input tree  "<<setw(7)<<tr->GetZipBytes()/(1024*1024)<<" MB"<<endl
      <<"size of splitted trees "<<setw(7)<<filesize/(1024*1024)<<" MB"<<endl
      <<"=> size of splitted tree(s) is "<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*filesize)/tr->GetZipBytes()<<" % of size of input tree"<<endl;
  ofile->Close();
  infile->Close();
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
