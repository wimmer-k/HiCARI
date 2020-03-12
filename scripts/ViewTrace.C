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

#include "/home/wimmer/programs/HrROOT/inc/Trace.hh"
char* filename = "root/run0031.root";
void CoreTraces(int firstevt=0, int lastevt=-1){
  TFile *f = new TFile(filename);
  TTree* tr = (TTree*)f->Get("build");
  Mode3Event *mode3 = new Mode3Event;
  tr->SetBranchAddress("mode3Event",&mode3);
  TH2F *traces = new TH2F("traces","traces", 200,0,200,1100,-10000,1000);
  if(lastevt<0)
    lastevt = tr->GetEntries();
  for(int i=firstevt; i<lastevt; i++){ 
    Int_t status = tr->GetEvent(i);
    for(int e=0;e<mode3->GetMult();e++){
      Trace *trace = mode3->GetHit(e)->GetCoreTrace();
      if (trace==NULL){
	cout << "Null core trace, aborting" << endl;
	continue;
      }
      for(int i=0;i<trace->GetLength();i++){
	traces->Fill(i,(int)trace->GetTrace()[i]);
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
    cout << "Trace " << g.size() << " Length " << trace->GetLength() << 
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
    
    int data[200];
    int x[200];
      
    for(int i=0;i<trace->GetLength();i++){
      x[i] = i;
      data[i] = (int)trace->GetTrace()[i];
    }
    g.push_back(new TGraph(trace->GetLength(),x,data));
    mg->Add(g.back(),"LP");
    // if(g.size()==1){
    //   g.back()->Draw("APL");
    // } else {
    //   g.back()->Draw("PL");
    // }
    if(trace->GetChn()==9)
      g.back()->SetLineColor(2);

  }
  mg->Draw("a");
  return g;
}
void ViewTrace(int n, int e, int p){
  vector<TGraph*> g = ViewTrace(n,e);
  
  g[p]->Draw("APL");

}
