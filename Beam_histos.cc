#include <iostream>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <signal.h>

#include "TFile.h"
#include "TEnv.h"
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
#include "FocalPlane.hh"
#include "PPAC.hh"

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
  char* tname = "BRZDtr";
  
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-v", "verbose", &vl);
  interface->Add("-t", "name of the tree, default \"tr\"", &tname);
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

  vector<TCutG*> InCut;
  vector<TCutG*> OutCut;
  vector<int> InCut_sel;
  vector<double> beta;

  int useCorrected = 0;
  int BR_AoQ = 2;
  int ZD_AoQ = 5;
  double zrange[2] = {10,30};
  double aoqrange[2] = {2.2,2.8};
  double timerange = 50e6;
  
  
  if(SettingFile != NULL){
    TEnv* set = new TEnv(SettingFile);
    char* cfilename = (char*)set->GetValue("CutFile","file");
    incuts = set->GetValue("NInCuts",0);
    outcuts = set->GetValue("NOutCuts",0);
    
    useCorrected = set->GetValue("UseAoQCorr",0);
    if(useCorrected)
      cout << "using corrected A/q values for PID gates" << endl;
    else
      cout << "using raw A/q values for PID gates" << endl;

    BR_AoQ = set->GetValue("UseBRAoQ",2);
    ZD_AoQ = set->GetValue("UseZDAoQ",5);

    zrange[0] = set->GetValue("Z.Range.Min",10.);
    zrange[1] = set->GetValue("Z.Range.Max",30.);
    aoqrange[0] = set->GetValue("AoQ.Range.Min",2.2);
    aoqrange[1] = set->GetValue("AoQ.Range.Max",2.8);

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
  tr = new TChain(tname);
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
  PPAC* ppac = new PPAC;
  tr->SetBranchAddress("ppacs",&ppac);
  Double_t nentries = tr->GetEntries();


  TList *hlist = new TList();
  TH1F* trigger = new TH1F("trigger","trigger",10,0,10);hlist->Add(trigger);
  TH2F* bigrips = new TH2F("bigrips","bigrips",1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);hlist->Add(bigrips);
  TH2F* zerodeg = new TH2F("zerodeg","zerodeg",1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);hlist->Add(zerodeg);
  TH2F* bigripsC = new TH2F("bigripsC","bigripsC",1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);hlist->Add(bigripsC);
  TH2F* zerodegC = new TH2F("zerodegC","zerodegC",1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);hlist->Add(zerodegC);

  TH2F* RawF7IC = new TH2F("RawF7IC","RawF7IC",6,0,6,4096,0,4096);hlist->Add(RawF7IC);
  TH2F* RawF11IC = new TH2F("RawF11IC","RawF11IC",6,0,6,4096,0,4096);hlist->Add(RawF11IC);
  TH2F* MatchF7IC = new TH2F("MatchF7IC","MatchF7IC",6,0,6,4096,0,4096);hlist->Add(MatchF7IC);
  TH2F* MatchF11IC = new TH2F("MatchF11IC","MatchF11IC",6,0,6,4096,0,4096);hlist->Add(MatchF11IC);


  TH2F* bigrips_Z_time = new TH2F("bigrips_Z_time","bigrips_Z_time",
				  timerange/10000,0,timerange,1000,zrange[0],zrange[1]);hlist->Add(bigrips_Z_time);
  TH2F* bigrips_AoQ_time = new TH2F("bigrips_AoQ_time","bigrips_AoQ_time",
				    timerange/10000,0,timerange,1000,aoqrange[0],aoqrange[1]);hlist->Add(bigrips_AoQ_time);
  TH2F* zerodeg_Z_time = new TH2F("zerodeg_Z_time","zerodeg_Z_time",
				  timerange/10000,0,timerange,1000,zrange[0],zrange[1]);hlist->Add(zerodeg_Z_time);
  TH2F* zerodeg_AoQ_time = new TH2F("zerodeg_AoQ_time","zerodeg_AoQ_time",
				    timerange/10000,0,timerange,1000,aoqrange[0],aoqrange[1]);hlist->Add(zerodeg_AoQ_time);
  
  
  // ppacs
  TH2F* tsumx_id = new TH2F("tsumx_id","tsumx_id",NPPACS,0,NPPACS,2500,0,250);hlist->Add(tsumx_id);
  TH2F* tsumy_id = new TH2F("tsumy_id","tsumy_id",NPPACS,0,NPPACS,2500,0,250);hlist->Add(tsumy_id);
  TH1F* tsumx[NPPACS];
  TH1F* tsumy[NPPACS];
  for(unsigned short p=0;p<NPPACS;p++){
    tsumx[p] = new TH1F(Form("tsumx_%d",p),Form("tsumx_%d",p),1000,-200,800);hlist->Add(tsumx[p]);
    tsumy[p] = new TH1F(Form("tsumy_%d",p),Form("tsumy_%d",p),1000,-200,800);hlist->Add(tsumy[p]);
  }
  
  // check F8 and beam directions
  TH1F* f8ppacX[6];
  TH1F* f8ppacY[6];
  TH2F* f8ppacXY[6];
  for(int p=0;p<6;p++){
    f8ppacX[p] = new TH1F(Form("f8ppacX_%d",p),Form("f8ppacX_%d",p),200,-100,100);hlist->Add(f8ppacX[p]);
    f8ppacY[p] = new TH1F(Form("f8ppacY_%d",p),Form("f8ppacY_%d",p),200,-100,100);hlist->Add(f8ppacY[p]);
    f8ppacXY[p] = new TH2F(Form("f8ppacXY_%d",p),Form("f8ppacXY_%d",p),200,-100,100,200,-100,100);hlist->Add(f8ppacXY[p]);
  }
  TH2F* compareX[2];
  TH2F* compareY[2];
  TH1F* compare1dX[2];
  TH1F* compare1dY[2];
  for(int p=0;p<2;p++){
    compareX[p] = new TH2F(Form("compareX_%d",p),Form("compareX_%d",p),200,-100,100,200,-100,100);hlist->Add(compareX[p]);
    compareY[p] = new TH2F(Form("compareY_%d",p),Form("compareY_%d",p),200,-100,100,200,-100,100);hlist->Add(compareY[p]);
    compare1dX[p] = new TH1F(Form("compare1dX_%d",p),Form("compare1dX_%d",p),1000,-100,100);hlist->Add(compare1dX[p]);
    compare1dY[p] = new TH1F(Form("compare1dY_%d",p),Form("compare1dY_%d",p),1000,-100,100);hlist->Add(compare1dY[p]);
  }
  TH2F* ppacZpos = new TH2F("ppacZpos","ppacZpos",40,0,40,3000,-1500,1500);hlist->Add(ppacZpos);

  //background inspection
  TH2F* PL3PPAC3   = new TH2F("PL3PPAC3","PL3PPAC3",    200,  0, 5, 200, -50, 50);hlist->Add(PL3PPAC3);
  TH2F* PL7PPAC7   = new TH2F("PL7PPAC7","PL7PPAC7",    200,  0, 5, 200, -50, 50);hlist->Add(PL7PPAC7);
  TH2F* PL8PPAC8   = new TH2F("PL8PPAC8","PL8PPAC8",    200,-10,10, 200,-100,100);hlist->Add(PL8PPAC8);
  TH2F* PL11PPAC11 = new TH2F("PL11PPAC11","PL11PPAC11",200,  0,15, 200, -50, 50);hlist->Add(PL11PPAC11);
  TH2F* PL3PL3     = new TH2F("PL3PL3","PL3PL3",        200,  0, 5, 200,-2,2);hlist->Add(PL3PL3);
  TH2F* PL7PL7     = new TH2F("PL7PL7","PL7PL7",        200,  0, 5, 200,-2,2);hlist->Add(PL7PL7);
  TH2F* PL8PL8     = new TH2F("PL8PL8","PL8PL8",        200,-10,10, 200,-2,2);hlist->Add(PL8PL8);
  TH2F* PL11PL11   = new TH2F("PL11PL11","PL11PL11",    200,  0,15, 200,-5,2);hlist->Add(PL11PL11);

    
  TH1F* deltadiff[2];
  for(unsigned short b=0;b<2;b++){
    deltadiff[b] = new TH1F(Form("deltadiff_%d",b),Form("deltadiff_%d",b),1000,-10,10);hlist->Add(deltadiff[b]);
  }
  TH1F* delta[4];
  for(unsigned short b=0;b<4;b++){
    delta[b] = new TH1F(Form("delta_%d",b),Form("delta_%d",b),1000,-10,10);hlist->Add(delta[b]);
  }

  // PID stuff
  vector<TH2F*> zerodeg_b;
  vector<vector<TH2F*> > position_b_z;
  vector<vector<TH2F*> > position_b;
  position_b.resize(NFPLANES);
  position_b_z.resize(NFPLANES);


  vector<TH2F*> MatchF7IC_b;
  vector<TH2F*> MatchF11IC_z;
  TH2F* h;
  for(int i=0;i<incuts;i++){
    h = new TH2F(Form("zerodeg_%s",InCut[i]->GetName()),
		 Form("zerodeg_%s",InCut[i]->GetName()),
		 1000,aoqrange[0],aoqrange[1],1000,zrange[0],zrange[1]);
    zerodeg_b.push_back(h);
    hlist->Add(zerodeg_b.back());
    
    h = new TH2F(Form("MatchF7IC_%s",InCut[i]->GetName()),
		 Form("MatchF7IC_%s",InCut[i]->GetName()),
		 6,0,6,4096,0,4096);
    MatchF7IC_b.push_back(h);
    hlist->Add(MatchF7IC_b.back());
    
    
    for(unsigned short f=0;f<NFPLANES;f++){
      h = new TH2F(Form("F%dXY_%s",fpID[f],InCut[i]->GetName()),
		   Form("F%dXY_%s",fpID[f],InCut[i]->GetName()),
		   500,-120,120,500,-120,120);
      
      position_b[f].push_back(h);
      hlist->Add(position_b[f].back());
    }
  }

  for(int o=0;o<outcuts;o++){
    h = new TH2F(Form("MatchF11IC_%s",OutCut[o]->GetName()),
		 Form("MatchF11IC_%s",OutCut[o]->GetName()),
		 6,0,6,4096,0,4096);
    MatchF11IC_z.push_back(h);
    hlist->Add(MatchF11IC_z.back());
    
    
    for(unsigned short f=0;f<NFPLANES;f++){
      h = new TH2F(Form("F%dXY_%s_%s",fpID[f],InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		   Form("F%dXY_%s_%s",fpID[f],InCut[InCut_sel[o]]->GetName(),OutCut[o]->GetName()),
		   500,-120,120,500,-120,120);
      
      position_b_z[f].push_back(h);
      hlist->Add(position_b_z[f].back());
      
    }
    
  }
  
  Int_t nbytes = 0;
  Int_t status;
  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    bz->Clear();
    for(unsigned short f=0;f<NFPLANES;f++){
      fp[f]->Clear();
    }
    ppac->Clear();
    
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
    bigripsC->Fill(bz->GetCorrAQ(2),bz->GetZ(2));
    zerodegC->Fill(bz->GetCorrAQ(5),bz->GetZ(5));

    bigrips_Z_time->Fill(i,bz->GetZ(2));
    bigrips_AoQ_time->Fill(i,bz->GetAQ(2));
    zerodeg_Z_time->Fill(i,bz->GetZ(5));
    zerodeg_AoQ_time->Fill(i,bz->GetAQ(2));
  

    // plastics
    Plastic *pl3 = fp[fpNr(3)]->GetPlastic();
    PL3PPAC3->Fill(pl3->GetTimeL() - pl3->GetTimeR(), fp[fpNr(3)]->GetTrack()->GetX());
    PL3PL3->Fill(pl3->GetTimeL() - pl3->GetTimeR(), log(pl3->GetChargeL()/pl3->GetChargeR()));
    Plastic *pl7 = fp[fpNr(7)]->GetPlastic();
    PL7PPAC7->Fill(pl7->GetTimeL() - pl7->GetTimeR(), fp[fpNr(7)]->GetTrack()->GetX());
    PL7PL7->Fill(pl7->GetTimeL() - pl7->GetTimeR(), log(pl7->GetChargeL()/pl7->GetChargeR()));
    Plastic *pl8 = fp[fpNr(8)]->GetPlastic();
    PL8PPAC8->Fill(pl8->GetTimeL() - pl8->GetTimeR(), fp[fpNr(8)]->GetTrack()->GetX());
    PL8PL8->Fill(pl8->GetTimeL() - pl8->GetTimeR(), log(pl8->GetChargeL()/pl8->GetChargeR()));
    Plastic *pl11 = fp[fpNr(11)]->GetPlastic();
    PL11PPAC11->Fill(pl11->GetTimeL() - pl11->GetTimeR(), fp[fpNr(11)]->GetTrack()->GetX());
    PL11PL11->Fill(pl11->GetTimeL() - pl11->GetTimeR(), log(pl11->GetChargeL()/pl11->GetChargeR()));
     
    // ppacs
    for(unsigned short p=0;p<ppac->GetN();p++){
      SinglePPAC *sp = ppac->GetPPAC(p);
      tsumx_id->Fill(sp->GetID(),sp->GetTsumX());
      tsumy_id->Fill(sp->GetID(),sp->GetTsumY());
      if(sp->GetID()>-1 && sp->GetID()<NPPACS){
	tsumx[sp->GetID()]->Fill(sp->GetTsumX());
	tsumy[sp->GetID()]->Fill(sp->GetTsumY());
      }
    }
    TVector3 ppacpos[3];
    ppacpos[0] = ppac->PPACPosition(ppac->GetPPACID(19),ppac->GetPPACID(20));
    ppacpos[1] = ppac->PPACPosition(ppac->GetPPACID(21),ppac->GetPPACID(22));
    ppacpos[2] = ppac->PPACPosition(ppac->GetPPACID(35),ppac->GetPPACID(36));

    bz->SetIncomingDirection(ppacpos[1]-ppacpos[0]);
    TVector3 inc = bz->GetIncomingDirection();
    double a = inc.X()/inc.Z();
    double b = inc.Y()/inc.Z();


    double x = ppacpos[1].X() + a * (ppac->GetPPACID(35)->GetZ()-ppacpos[1].Z());
    double y = ppacpos[1].Y() + b * (ppac->GetPPACID(35)->GetZ()-ppacpos[1].Z());
    compareX[0]->Fill(ppac->GetPPACID(35)->GetX(),x);
    compareY[0]->Fill(ppac->GetPPACID(35)->GetY(),y);
    compare1dX[0]->Fill(ppac->GetPPACID(35)->GetX()-x);
    compare1dY[0]->Fill(ppac->GetPPACID(35)->GetY()-y);

    x = ppacpos[1].X() + a * (ppac->GetPPACID(36)->GetZ()-ppacpos[1].Z());
    y = ppacpos[1].Y() + b * (ppac->GetPPACID(36)->GetZ()-ppacpos[1].Z());
    compareX[1]->Fill(ppac->GetPPACID(36)->GetX(),x);
    compareY[1]->Fill(ppac->GetPPACID(36)->GetY(),y);
    compare1dX[1]->Fill(ppac->GetPPACID(36)->GetX()-x);
    compare1dY[1]->Fill(ppac->GetPPACID(36)->GetY()-y);

    // PPACs
    for(unsigned short p=0;p<ppac->GetN();p++){
      SinglePPAC *sp = ppac->GetPPAC(p);
      ppacZpos->Fill(sp->GetID(),sp->GetZ());
      if(sp->GetID()>=19 && sp->GetID()<=22){
	f8ppacX[sp->GetID()-19]->Fill(sp->GetX());
	f8ppacY[sp->GetID()-19]->Fill(sp->GetY());
	if(sp->Fired())
	  f8ppacXY[sp->GetID()-19]->Fill(sp->GetX(),sp->GetY());
      }
      if(sp->GetID()>=35 && sp->GetID()<=36){
	f8ppacX[sp->GetID()-35+4]->Fill(sp->GetX());
	f8ppacY[sp->GetID()-35+4]->Fill(sp->GetY());
	if(sp->Fired())
	  f8ppacXY[sp->GetID()-35+4]->Fill(sp->GetX(),sp->GetY());
      }
    
    }
    
    // musics
    MUSIC* f7ic = fp[fpNr(7)]->GetMUSIC();
    for(ushort i=0; i<f7ic->GetChan().size();i++){
      RawF7IC->Fill(f7ic->GetChan().at(i), f7ic->GetADC().at(i));
      MatchF7IC->Fill(f7ic->GetChan().at(i), f7ic->GetGainMatchADC().at(i));
    }
    MUSIC* f11ic = fp[fpNr(11)]->GetMUSIC();
    for(ushort i=0; i<f11ic->GetChan().size();i++){
      RawF11IC->Fill(f11ic->GetChan().at(i), f11ic->GetADC().at(i));
      MatchF11IC->Fill(f11ic->GetChan().at(i), f11ic->GetGainMatchADC().at(i));
    }

    // charge state change
    deltadiff[0]->Fill(bz->GetDelta(1) - bz->GetDelta(0));
    deltadiff[1]->Fill(bz->GetDelta(3) - bz->GetDelta(2));
    for(unsigned short b=0;b<4;b++){
      delta[b]->Fill(bz->GetDelta(b));
    }

    
    // gated
    for(int in=0;in<incuts;in++){
      if( (!useCorrected && InCut[in]->IsInside(bz->GetAQ(2),bz->GetZ(2))) ||
	  (useCorrected && InCut[in]->IsInside(bz->GetCorrAQ(2),bz->GetZ(2))) ){

	zerodeg_b[in]->Fill(bz->GetAQ(5),bz->GetZ(5));

	for(ushort i=0; i<f7ic->GetChan().size();i++){
	  MatchF7IC_b[in]->Fill(f7ic->GetChan().at(i), f7ic->GetGainMatchADC().at(i));
	}
	
	for(unsigned short f=0;f<NFPLANES;f++){
	  position_b[f][in]->Fill(fp[f]->GetTrack()->GetX(), fp[f]->GetTrack()->GetY());
	}
      }//isinside
    }//incuts
    for(int o=0;o<outcuts;o++){
      if( (!useCorrected && InCut[InCut_sel[o]]->IsInside(bz->GetAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetAQ(5),bz->GetZ(5))) ||
    	  (useCorrected && InCut[InCut_sel[o]]->IsInside(bz->GetCorrAQ(2),bz->GetZ(2)) && OutCut[o]->IsInside(bz->GetCorrAQ(5),bz->GetZ(5)))){

	for(ushort i=0; i<f11ic->GetChan().size();i++){
	  MatchF11IC_z[o]->Fill(f11ic->GetChan().at(i), f11ic->GetGainMatchADC().at(i));
	}
	
    	for(unsigned short f=0;f<NFPLANES;f++){
	  position_b_z[f][o]->Fill(fp[f]->GetTrack()->GetX(), fp[f]->GetTrack()->GetY());
	}
      }//isinside
    }//outcuts
    
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
