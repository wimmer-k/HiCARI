#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLine.h"

#include "/home/users/analysis/kathrin/code/trunk/Trace.hh"
#include "/home/users/analysis/kathrin/code/trunk/S800.hh"

void ViewCoreTrace(int n){
  TFile *f = new TFile("globalraw.root");
  TTree* tr = (TTree*)f->Get("gtr");
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
    cout << "Trace Lenghth " << trace->GetLength() << endl;
    
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
  TFile *f = new TFile("globalraw.root");
  TTree* tr = (TTree*)f->Get("gtr");
  Mode3Event *mode3 = new Mode3Event;
  tr->SetBranchAddress("mode3Event",&mode3);
  
  Int_t status = tr->GetEvent(n);
  vector<TGraph*> g;
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
    if(g.size()==1){
      g.back()->Draw("APL");
    } else {
      g.back()->Draw("PL");
    }
  }

  return g;
}
void ViewTrace(int n, int e, int p){
  vector<TGraph*> g = ViewTrace(n,e);
  
  g[p]->Draw("APL");

}
void ViewCard29(int n, char* filename = "raw_tw3000.root", int tracenum=0){
  TFile *f = new TFile(filename);
  TTree* tr = (TTree*)f->Get("gtr");
  Mode3Event *mode3 = new Mode3Event;
  S800 *s800 = new S800;
  tr->SetBranchAddress("mode3Event",&mode3);
  tr->SetBranchAddress("s800",&s800);
  
  Int_t status = tr->GetEvent(n);
  TGraph* g;
  Mode3Hit* hit = mode3->GetHit(0);
  if(hit==NULL){
    cout << "bad hit, aborting" << endl;
    return;
  }
  cout << "Hit length " << mode3->GetHit(0)->GetMult() << endl;
  Trace *trace = hit->GetTrace(tracenum);
  if(trace==NULL || trace->GetLength()<1){
    cout << "bad trace, aborting" << endl;
    return;
  }
  if(trace->GetHole()!=31 && false){
    cout << "not card 29! Aborting" << endl;
    return;
  }
    
  cout << "Length " << trace->GetLength() << 
    "\tEnergy " << trace->GetEnergy() <<
    "\tBoard " << trace->GetBoard() <<
    "\tChannel " << trace->GetChn() <<
    "\tHole " << trace->GetHole() <<
    "\tCrystal " << trace->GetCrystal() << endl;
  
  int data[200];
  int x[200];
      
  

  double crossed =-1;
  int thresh = 1000;
  for(int i=0;i<trace->GetLength();i++){
    x[i] = i;
    data[i] = (int)trace->GetTrace()[i];
    if(crossed <0 && data[i]>thresh &&i>0){
      crossed = (i-1) + (float)(thresh-data[i-1])/(float)(data[i]-data[i-1]);
    }
  }
  g = new TGraph(trace->GetLength(),x,data);
  //g->Draw("APL");
  g->Draw("A*L");
  
  cout << "crossed at " << crossed  << endl;
  cout << "pileup bit was " << trace->GetPileUp() << endl;
  cout << "timestamp was " << trace->GetTS() << endl;
  cout << "s800 difference was " << trace->GetTS() - s800->GetTS() << endl;

  TLine *lthresh = new TLine(0,thresh,trace->GetLength(),thresh);
  lthresh->SetLineColor(2);
  lthresh->Draw();
  TLine *lcrossed = new TLine(crossed,0,crossed,5000);
  lcrossed->SetLineColor(3);
  lcrossed->Draw();

  return;
}
