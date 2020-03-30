#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLine.h"
#include "TH2F.h"

#include "/home/gamma20/HiCARI/inc/Trace.hh"
#include "/home/gamma20/HiCARI/inc/HiCARI.hh"
char* filename = "/home/gamma20/rootfiles/merge_BR0018_Hi0431.root";
void TimeStampDiff(int firstevt=0, int lastevt=-1){
  TFile *f = new TFile(filename);
  TTree* tr = (TTree*)f->Get("tr");
  unsigned long long int brTS;
  unsigned long long int hiTS;
  unsigned long long int prev_brTS = 0;
  unsigned long long int prev_hiTS = 0;
  unsigned long long int prev_TS = 0;

  tr->SetBranchAddress("brTS",&brTS);
  tr->SetBranchAddress("hiTS",&hiTS);
  double range = 2e7;
  TH1F *tdiff = new TH1F("tdiff","tdiff", 1000,-500,500);
  TH1F *hi_tolast = new TH1F("hi_tolast","hi_tolast", 10000,0,range);
  TH1F *hi_tolasthi = new TH1F("hi_tolasthi","hi_tolasthi", 10000,0,range);
  TH1F *br_tolast = new TH1F("br_tolast","br_tolast", 10000,0,range);
  TH1F *br_tolastbr = new TH1F("br_tolastbr","br_tolastbr", 10000,0,range);
  if(lastevt<0)
    lastevt = tr->GetEntries();
  cout << "prev_TS\tprev_brTS\tprev_hiTS\tbrTS\thiTS" << endl; 
  for(int i=firstevt; i<lastevt-1; i++){ 
    Int_t status = tr->GetEvent(i);
    //cout << prev_TS << "\t" << prev_brTS << "\t" << prev_hiTS << "\t" << brTS << "\t" << hiTS << endl; 
    if(hiTS>0){
      hi_tolast->Fill(hiTS-prev_TS);
      hi_tolasthi->Fill(hiTS-prev_hiTS);
      prev_hiTS = hiTS;
    }
    if(brTS>0){
      br_tolast->Fill(brTS-prev_TS);
      br_tolastbr->Fill(brTS-prev_brTS);
      prev_brTS = brTS;
    }
    if(hiTS>0 && brTS>0){
      tdiff->Fill(hiTS*1.0-brTS*1.0);
      if(hiTS > brTS)
 	prev_TS = hiTS;
      else
	prev_TS = brTS;
    }
    if(i%10000==0)
      cout << i-firstevt <<" / " << lastevt-firstevt << " events,\t " << 100.0*(i-firstevt)/(lastevt-firstevt) << " % done " << endl; 
  }//events
  TCanvas *c = new TCanvas("c","c",1200,1200);
  c->Divide(1,5);
  c->cd(1);
  tdiff->Draw();
  c->cd(2);
  hi_tolast->Draw();
  c->cd(3);
  hi_tolasthi->Draw();
  //hi_tolasthi->SetLineColor(2);
  c->cd(4);
  br_tolast->Draw();
  c->cd(5);
  br_tolastbr->Draw();
  //br_tolastbr->SetLineColor(2);
  //tdiff->Draw();
}
void TimeStampDiffMode3(int firstevt=0, int lastevt=-1){
  TFile *f = new TFile(filename);
  TTree* tr = (TTree*)f->Get("build");
  Mode3Event *mode3 = new Mode3Event;
  tr->SetBranchAddress("mode3Event",&mode3);
  TH1F *tdiff = new TH1F("tdiff","tdiff", 1000,0,1000000);
  if(lastevt<0)
    lastevt = tr->GetEntries();
  for(int i=firstevt; i<lastevt-1; i++){ 
    Int_t status = tr->GetEvent(i);
    if(mode3->GetMult()>0){
      long long int thisTS = mode3->GetHit(0)->GetTS();
      long long int thisTrTS = mode3->GetHit(0)->GetTrace(0)->GetTS();
      status = tr->GetEvent(i+1);
      if(mode3->GetMult()>0){
	long long int nextTS = mode3->GetHit(0)->GetTS();
	long long int nextTrTS = mode3->GetHit(0)->GetTrace(0)->GetTS();
	//cout << thisTS <<"\t"<< nextTS <<"\t"<< nextTS-thisTS << ",\ten = " << mode3->GetHit(0)->GetTrace(0)->GetEnergy() << ",\tts = " << mode3->GetHit(0)->GetTrace(0)->GetTS() << endl;
	cout << thisTS <<"\t"<< nextTS <<"\t"<< nextTS-thisTS << ",\t " << thisTrTS <<"\t"<< nextTrTS << "\t" << nextTrTS-thisTrTS << endl;
	tdiff->Fill(nextTS-thisTS);
      }
    }
  }//events
  tdiff->Draw();
}
