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
#include "Reconstruction.hh"
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
double deg2rad = TMath::Pi()/180.;
double rad2deg = 180./TMath::Pi();
int main(int argc, char* argv[]){
  double time_start = get_time();  
  TStopwatch timer;
  timer.Start();
  signal(SIGINT,signalhandler);
  vector<char*> InputFiles;
  char* OutputFile = NULL;
  char* SettingFile = NULL;
  int nmax =0;
  int vl =0;
  char* tname = "tr";
  int writeTree = 0;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-v", "verbose", &vl);
  interface->Add("-t", "name of the tree, default \"tr\"", &tname);
  interface->Add("-w", "write also a tree", &writeTree);
  interface->CheckFlags(argc, argv);

  if(InputFiles.size() == 0 || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }
  if(SettingFile == NULL){
    cout << "No settings file given " << endl;
    return 2;
  }
  cout<<"input file(s):"<<endl;
  for(unsigned int i=0; i<InputFiles.size(); i++){
    cout<<InputFiles[i]<<endl;
  }
  cout<<"output file: "<<OutputFile<< endl;
  
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
  //PPAC* ppac = new PPAC;
  //tr->SetBranchAddress("ppacs",&ppac);

  Settings* set = new Settings(SettingFile);
  Reconstruction* rec = new Reconstruction(set);
  cout << "creating outputfile " << OutputFile << endl;
  TFile* ofile = new TFile(OutputFile,"recreate");
  ofile->cd();
  set->Write();
  

  TTree* rtr = new TTree("rtr","Reconstructed events");
  rtr->Branch("hicari",&hi,320000);
  rtr->Branch("mode2",&gr,320000);
  rtr->Branch("trigbit",&trigbit,"trigbit/I");
  rtr->Branch("beam",&bz,320000);
  

  TList *hlist = new TList();

  //histograms
  TH1F* trigger = new TH1F("trigger","trigger",10,0,10);hlist->Add(trigger);
  TH2F* bigrips = new TH2F("bigrips","bigrips",1000,1.8,2.3,1000,20,40);hlist->Add(bigrips);
  TH2F* zerodeg = new TH2F("zerodeg","zerodeg",1000,1.8,2.3,1000,20,40);hlist->Add(zerodeg);
  TH2F* bigrips_tr[10];
  TH2F* zerodeg_tr[10];
  TH1F* f5X_tr[10];
  for(int i=0;i<10;i++){
    bigrips_tr[i] = new TH2F(Form("bigrips_tr%d",i),Form("bigrips_tr%d",i),1000,1.8,2.3,1000,20,40);hlist->Add(bigrips_tr[i]);
    zerodeg_tr[i] = new TH2F(Form("zerodeg_tr%d",i),Form("zerodeg_tr%d",i),1000,1.8,2.3,1000,20,40);hlist->Add(zerodeg_tr[i]);
    f5X_tr[i] = new TH1F(Form("f5X_tr%d",i),Form("f5X_tr%d",i),3000,-150,150);hlist->Add(f5X_tr[i]);
  }

  TH2F* thetaphi = new TH2F("thetaphi","thetaphi",800,-4,4,1500,0,150);hlist->Add(thetaphi);
  TH2F* thetaphideg = new TH2F("thetaphideg","thetaphideg",800,-180,180,1500,0,5);hlist->Add(thetaphideg);
  TH2F* thetaphi_tr[10];
  TH2F* thetaphideg_tr[10];
  for(int i=0;i<10;i++){
    thetaphi_tr[i] = new TH2F(Form("thetaphi_tr%d",i),Form("thetaphi_tr%d",i),800,-4,4,1500,0,150);hlist->Add(thetaphi_tr[i]);
    thetaphideg_tr[i] = new TH2F(Form("thetaphideg_tr%d",i),Form("thetaphideg_tr%d",i),800,-180,180,1500,0,5);hlist->Add(thetaphideg_tr[i]);
  }
  TH2F* F8xy[3];
  for(int i=0;i<3;i++){
    F8xy[i] = new TH2F(Form("F8_%dxy",i),Form("F8_%dxy",i),100,-50,50,100,-50,50);hlist->Add(F8xy[i]);
  }  
  TH2F* targetxy = new TH2F("targetxy","targetxy",100,-50,50,100,-50,50);hlist->Add(targetxy);
  TH2F* targetxz = new TH2F("targetxz","targetxz",100,-50,50,100,-50,50);hlist->Add(targetxz);
  TH2F* targetyz = new TH2F("targetyz","targetyz",100,-50,50,100,-50,50);hlist->Add(targetyz);
  
  TH2F* h_egam_tgam = new TH2F("h_egam_tgam","h_egam_tgam",1000,-500,500,4000,0,4000);hlist->Add(h_egam_tgam);
  TH2F* h_egamdc_tgam = new TH2F("h_egamdc_tgam","h_egamdc_tgam",1000,-500,500,4000,0,4000);hlist->Add(h_egamdc_tgam);
  TH1F* h_egamdc = new TH1F("h_egamdc","h_egamdc",4000,0,4000);hlist->Add(h_egamdc);
  TH2F* h_egamdc_theta = new TH2F("h_egamdc_theta","h_egamdc_theta",180,0,180,4000,0,4000);hlist->Add(h_egamdc_theta);
  TH2F* h_egam_beta = new TH2F("h_egam_beta","h_egam_beta",1000,0.55,0.65,4000,0,4000);hlist->Add(h_egam_beta);
  TH2F* h_egamdc_beta = new TH2F("h_egamdc_beta","h_egamdc_beta",1000,0.55,0.65,4000,0,4000);hlist->Add(h_egamdc_beta);
  TH2F* h_egamdc_summary = new TH2F("h_egamdc_summary","h_egamdc_summary",60,0,60,4000,0,4000);hlist->Add(h_egamdc_summary);
  TH2F* h_tgam_summary = new TH2F("h_tgam_summary","h_tgam_summary",60,0,60,2000,-1000,1000);hlist->Add(h_tgam_summary);
  TH2F* h_tgam_summary_HE = new TH2F("h_tgam_summary_HE","h_tgam_summary_HE",60,0,60,2000,-1000,1000);hlist->Add(h_tgam_summary_HE);
  TH2F* h_tgam_trigbit = new TH2F("h_tgam_trigbit","h_tgam_trigbit",10,0,10,2000,-1000,1000);hlist->Add(h_tgam_trigbit);
  TH2F* h_egamdc_egamdc = new TH2F("h_egamdc_egamdc","h_egamdc_egamdc",1000,0,4000,1000,0,4000);hlist->Add(h_egamdc_egamdc);
  
  TH2F* h_egamdc_segsummary = new TH2F("h_egamdc_segsummary","h_egamdc_segsummary",400,0,400,4000,0,4000);hlist->Add(h_egamdc_segsummary);
  
  
  Double_t nentries = tr->GetEntries();
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

    //gate on F5X position 
    if(!isnan(fp[fpNr(5)]->GetTrack()->GetX()) && !rec->F5XGate(fp[fpNr(5)]->GetTrack()->GetX()))
      continue;

    //gate on charge changes in BigRIPS and Zerodeg
    if(trigbit>1 && rec->ChargeChange(bz))
      continue;

    
    //beam direction, scattering angle
    //align
    rec->AlignPPAC(bz->GetF8PPAC3A(), bz->GetF8PPAC3B());
    TVector3 ppacpos[3];
    for(int j=0;j<3;j++){
      ppacpos[j] = rec->PPACPosition(bz->GetF8PPAC(j,0),bz->GetF8PPAC(j,1));
      bz->SetF8Position(j,ppacpos[j]);
      F8xy[j]->Fill(ppacpos[j].X(), ppacpos[j].Y());
    }
    bz->SetIncomingDirection(ppacpos[1]-ppacpos[0]);
    TVector3 inc = bz->GetIncomingDirection();

    // target position with respect to the nominal focal point
    TVector3 targ = rec->TargetPosition(inc,ppacpos[1]);
    targ.SetZ(set->TargetZ());
    bz->SetTargetPosition(targ);
    TVector3 out, sca;
    if(trigbit>1){
      bz->SetOutgoingDirection(ppacpos[2]-targ);
      out = bz->GetOutgoingDirection();
      sca = bz->GetScatteredDirection();
    }

    //beta
    //bz->SetDopplerBeta(set->TargetBeta());
    //bz->SetDopplerBeta(bz->GetRIPSBeta(3));
    bz->SetDopplerBeta(rec->EventBeta(bz));

    
    // //new HiCARI positions
    // for(int h=0; h<hi->GetMult(); h++){
    //   HiCARIHitCalc* hit = hi->GetHit(h);
    //   if(hit->IsBigRIPS()){
    // 	continue;
    //   }
    //   TVector3 pos = hit->GetPosition();
    //   int cl = hit->GetCluster();
    //   int cr = hit->GetCrystal();
    //   int se = hit->GetMaxSegment();
    //   TVector3 newpos = rec->GammaPosition(cl,cr,se);
    //   cout << setw(7) << setprecision(5) << pos.Theta()*rad2deg << "\t" << newpos.Theta()*rad2deg << "\t" << pos.Theta()*rad2deg - newpos.Theta()*rad2deg << endl;
    // }
    
    
    rec->SetGammaPositions(hi);

    hi->DopplerCorrect(bz);

    // histos
    // BEAM
    trigger->Fill(trigbit);
    bigrips->Fill(bz->GetCorrAQ(2),bz->GetZ(2));
    zerodeg->Fill(bz->GetCorrAQ(5),bz->GetZ(5));
    if(trigbit>-1 && trigbit<10){
      bigrips_tr[trigbit]->Fill(bz->GetCorrAQ(2),bz->GetZ(2));
      zerodeg_tr[trigbit]->Fill(bz->GetCorrAQ(5),bz->GetZ(5));
      f5X_tr[trigbit]->Fill(fp[fpNr(5)]->GetTrack()->GetX());
      thetaphi_tr[trigbit]->Fill(bz->GetPhi(),bz->GetTheta()*1000);
      thetaphideg_tr[trigbit]->Fill(bz->GetPhi()*rad2deg,bz->GetTheta()*rad2deg);
    }
    if(trigbit>1){
      //cout << bz->GetPhi() <<"\t" << sca.Phi() - inc.Phi() << endl;
      thetaphi->Fill(bz->GetPhi(),bz->GetTheta()*1000);
      thetaphideg->Fill(bz->GetPhi()*rad2deg,bz->GetTheta()*rad2deg);
    }

    //tp position with respect to HiCARI center
    TVector3 tp = bz->GetTargetPosition();    
    targetxy->Fill(tp.X(),tp.Y());
    targetxz->Fill(tp.X(),tp.Z());
    targetyz->Fill(tp.Y(),tp.Z());


    for(int h=0; h<hi->GetMult(); h++){
      HiCARIHitCalc* hit = hi->GetHit(h);
      if(hit->IsBigRIPS())
	continue;
      if(hit->IsTracking())
	continue;
      h_egam_tgam->Fill(hit->GetTime(),hit->GetEnergy());
      h_egamdc_tgam->Fill(hit->GetTime(),hit->GetDCEnergy());
      if(hit->GetPosition().Theta()>0 && hit->GetEnergy() > 10){
	h_egamdc->Fill(hit->GetDCEnergy());
	h_egam_beta->Fill(bz->GetDopplerBeta(), hit->GetEnergy());
	h_egamdc_beta->Fill(bz->GetDopplerBeta(), hit->GetDCEnergy());
	h_egamdc_theta->Fill(hit->GetPosition().Theta()*180./3.1415, hit->GetDCEnergy());
	h_egamdc_summary->Fill(4*hit->GetCluster()+hit->GetCrystal(), hit->GetDCEnergy());
	if(hit->GetMaxSegment()>-1)
	  h_egamdc_segsummary->Fill((4*hit->GetCluster()+hit->GetCrystal())*6+hit->GetMaxSegment(), hit->GetDCEnergy());
	h_tgam_summary->Fill(4*hit->GetCluster()+hit->GetCrystal(), hit->GetTime());
	if(hit->GetEnergy()>2000)
	  h_tgam_summary_HE->Fill(4*hit->GetCluster()+hit->GetCrystal(), hit->GetTime());
	h_tgam_trigbit->Fill(trigbit, hit->GetTime());
	for(int h2=0; h2<hi->GetMult(); h2++){
	  if(h==h2)
	    continue;
	  HiCARIHitCalc* hit2 = hi->GetHit(h2);
	  h_egamdc_egamdc->Fill(hit->GetDCEnergy(),hit2->GetDCEnergy());
	  h_egamdc_egamdc->Fill(hit2->GetDCEnergy(),hit->GetDCEnergy());
	}
      }//valid hit with position
    }//hits

    if(writeTree>0)
      rtr->Fill();
    if(i%10000 == 0){
      double time_end = get_time();
      cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*i)/nentries<<" % done\t"<<(Float_t)i/(time_end - time_start)<<" events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i<<"s to go \r"<<flush;
    }
    if(writeTree>0 && i%1000000 == 0)
      rtr->AutoSave();

  }
  cout << endl;


  cout << "writing to file" << endl;
  cout << endl;
  ofile->cd();
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
  if(writeTree>0)
    rtr->Write("",TObject::kOverwrite);
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
