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
char* filename = "~/rootfiles/run0301.root";
int verbose = 0;
int frontBL = 70;
Double_t flinear(Double_t *x, Double_t *par);
void Baseline(int hole, int cry, int slot, int firstevt=0, int lastevt=-1){
  TFile *f = new TFile(filename);
  TTree* tr = (TTree*)f->Get("build");
  Mode3Event *mode3 = new Mode3Event;
  tr->SetBranchAddress("mode3Event",&mode3);
  TCanvas *c = new TCanvas("c","c",900,400);
  c->Divide(2,1);
  TH2F *traces = new TH2F("traces","traces", 200,0,200,1200,-10000,2000);
  TH2F *smokep = new TH2F("smokep","smokep", 2000,0,5e5,1000,0,1000);
  if(lastevt<0)
    lastevt = tr->GetEntries();
  int ctr = 0;
  for(int i=firstevt; i<lastevt; i++){ 
    Int_t status = tr->GetEvent(i);
    for(int e=0;e<mode3->GetMult();e++){
      for(int t=0;t<mode3->GetHit(e)->GetMult();t++){
	Trace *trace = mode3->GetHit(e)->GetTrace(t);
	if (trace==NULL || trace->GetLength()<1){
	  cout << "bad trace, aborting" << endl;
	  continue;
	}
	if(trace->GetHole()!=hole || trace->GetCrystal()!=cry || trace->GetSlot()!=slot || trace->GetChn()!=9)
	  continue;
	if(verbose){
	  cout << "Trace " << ctr++ << " Length " << trace->GetLength() << 
	    "\tEnergy " << trace->GetEnergy() <<
	    "\tBoard " << trace->GetBoard() <<
	    "\tChannel " << trace->GetChn() <<
	    "\tHole " << trace->GetHole() <<
	    "\tCrystal " << trace->GetCrystal();
	  if(trace->GetChn()==9)
	    cout << " <- CC " << endl;
	  else if(trace->GetEnergy()>1000)
	    cout << " <- with net energy " << endl;
	  else
	    cout << endl;
	}//verbose
	double baseline =0;
	for(int j=0;j<trace->GetLength();j++){
	  traces->Fill(j,(int)trace->GetTrace()[j]);
	  if(j<frontBL)
	    baseline+=trace->GetTrace()[j];
	}
	baseline/=frontBL;
	smokep->Fill(trace->GetEnergy(),baseline);
      }//traces
    }//hits
  }//events
  c->cd(1);
  traces->Draw("colz");
  c->cd(2);
  smokep->Draw("colz");
}


void CoreTraces(int firstevt=0, int lastevt=-1, int crystal = -1){
  TFile *f = new TFile(filename);
  TTree* tr = (TTree*)f->Get("build");
  Mode3Event *mode3 = new Mode3Event;
  tr->SetBranchAddress("mode3Event",&mode3);
  TH2F *traces = new TH2F("traces","traces", 200,0,200,1200,-10000,2000);
  if(lastevt<0)
    lastevt = tr->GetEntries();
  for(int i=firstevt; i<lastevt; i++){ 
    Int_t status = tr->GetEvent(i);
    for(int e=0;e<mode3->GetMult();e++){
      Trace *trace = mode3->GetHit(e)->GetCoreTrace();
      if(crystal>-1 && trace->GetCrystal() != crystal)
	continue;
      if (trace==NULL){
	cout << "Null core trace, aborting" << endl;
	continue;
      }
      for(int j=0;j<trace->GetLength();j++){
	traces->Fill(j,(int)trace->GetTrace()[j]);
      }
      
    }
  }//events
  traces->Draw("colz");
}

void ViewCoreTrace(int n){
  TFile *f = new TFile(filename);
  TTree* tr = (TTree*)f->Get("build");
  Mode3Event *mode3 = new Mode3Event;
  tr->SetBranchAddress("mode3Event",&mode3);
  
  Int_t status = tr->GetEvent(n);
  cout << "Event length " << mode3->GetMult() << endl;
  vector<TGraph*> g;
  g.resize(mode3->GetMult());
  for(int e=0;e<mode3->GetMult();e++){
    cout << "Hit length " << mode3->GetHit(e)->GetMult() << endl;
    Trace *trace = mode3->GetHit(e)->GetCoreTrace();
    if (trace==NULL){
      cout << "Null core trace, aborting" << endl;
      return;
    }
    cout << "Trace " << g.size() << " Length " << trace->GetLength() << 
      "\tEnergy " << trace->GetEnergy() <<
      "\tBoard " << trace->GetBoard() <<
      "\tChannel " << trace->GetChn() <<
      "\tHole " << trace->GetHole() <<
      "\tCrystal " << trace->GetCrystal() << endl;
    
    int data[200];
    int x[200];
    
    for(int i=0;i<trace->GetLength();i++){
      x[i] = i;
      data[i] = (int)trace->GetTrace()[i];
    }
    g[e] = new TGraph(trace->GetLength(),x,data);
    if(e==0){
      g[e]->Draw("APL");
    } else {
      g[e]->Draw("PL");
    }
  }
}

vector<TGraph*> ViewTrace(int n, int e){
  TFile *f = new TFile(filename);
  TTree* tr = (TTree*)f->Get("build");
  Mode3Event *mode3 = new Mode3Event;
  tr->SetBranchAddress("mode3Event",&mode3);
  
  Int_t status = tr->GetEvent(n);
  vector<TGraph*> g;
  TMultiGraph *mg = new TMultiGraph();
  cout << "mult " << mode3->GetMult() << endl;
  if(e>= mode3->GetMult()){
    cout << "hit " << e << " not available, only " << mode3->GetMult() << " hits in event" << endl;
    return vector<TGraph*>();;
  }
  cout << "Hit length " << mode3->GetHit(e)->GetMult() << endl;
  for(int t=0;t<mode3->GetHit(e)->GetMult();t++){
    Trace *trace = mode3->GetHit(e)->GetTrace(t);
    if (trace==NULL || trace->GetLength()<1){
      cout << "bad trace, aborting" << endl;
      continue;
    }
    // cout << "Trace " << g.size() << " Length " << trace->GetLength() << 
    //   "\tEnergy " << trace->GetEnergy() <<
    //   "\tBoard " << trace->GetBoard() <<
    //   "\tChannel " << trace->GetChn() <<
    //   "\tHole " << trace->GetHole() <<
    //   "\tCrystal " << trace->GetCrystal();
    // if(trace->GetChn()==9)
    //   cout << " <- CC " << endl;
    // else if(trace->GetEnergy()>1000)
    //   cout << " <- with net energy " << endl;
    // else
    //   cout << endl;
    
    int data[200];
    int cfd[200];
    int x[200];
    int length = 1;
    int delay = 10;
    float fraction = 0.3;
    float baseline = 0;
    int baselinelength = 50;
    float led = -20;
    bool negative = false;
    int crossing = -1;
    for(int i=0;i<baselinelength;i++){
      baseline+=(int)trace->GetTrace()[i];
    }
    baseline/=baselinelength;
    for(int i=0;i<trace->GetLength();i++){
      x[i] = i;
      data[i] = (int)trace->GetTrace()[i]-baseline;
      if(trace->GetChn()==9){
	cfd[i] = 0;
	if(i>baselinelength){
	  for(int j=0;j<length;j++)
	    cfd[i] += fraction*data[i-j] - data[i-j-delay]; 
	}
	if(data[i]<led){
	  if(negative && cfd[i]>0){
	    cout << "crossing " << i << endl;
	    crossing = i;
	    negative = false;
	  }
	  negative = cfd[i]<0;
	}
      }
    }
    g.push_back(new TGraph(trace->GetLength(),x,data));
    mg->Add(g.back(),"LP");
    if(trace->GetChn()==9){
      g.back()->SetLineColor(2);
      g.back()->SetMarkerColor(2);
      g.push_back(new TGraph(trace->GetLength(),x,cfd));
      mg->Add(g.back(),"LP");
      g.back()->SetLineColor(3);
      g.back()->SetMarkerColor(3);
      if(crossing>baselinelength){
	TF1 *ff = new TF1("f",flinear,crossing-3,crossing+3,2);
	ff->SetParameters(10,0);
	g.back()->Fit(ff,"R");
	cout << -ff->GetParameter(1)/ff->GetParameter(0);
      }
    }
  }
  mg->Draw("a");
  return g;
}
void ViewTrace(int n, int e, int p){
  vector<TGraph*> g = ViewTrace(n,e);
  g[p]->Draw("APL");
}
Double_t flinear(Double_t *x, Double_t *par){
  return x[0]*par[0] + par[1];
}
