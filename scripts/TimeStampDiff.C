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
void TimeStampDiff(int firstevt=0, int lastevt=-1){
  TFile *f = new TFile("test.root");
  TTree* tr = (TTree*)f->Get("build");
  Mode3Event *mode3 = new Mode3Event;
  tr->SetBranchAddress("mode3Event",&mode3);
  TH1F *tdiff = new TH1F("tdiff","tdiff", 1000,2100000,2120000);
  if(lastevt<0)
    lastevt = tr->GetEntries();
  for(int i=firstevt; i<lastevt-1; i++){ 
    Int_t status = tr->GetEvent(i);
    if(mode3->GetMult()>0){
      long long int thisTS = mode3->GetHit(0)->GetTS();
      status = tr->GetEvent(i+1);
      if(mode3->GetMult()>0){
	long long int nextTS = mode3->GetHit(0)->GetTS();
	//cout << thisTS <<"\t"<< nextTS <<"\t"<< nextTS-thisTS << endl;
	tdiff->Fill(nextTS-thisTS);
      }
    }
  }//events
  tdiff->Draw();
}
