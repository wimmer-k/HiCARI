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
char* filename = "/home/gamma20/rootfiles/run0331.root";
void TimeStampDiff(int firstevt=0, int lastevt=-1){
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
