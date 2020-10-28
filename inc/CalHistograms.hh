#ifndef CAL_HISTOGRAMS_HH__
#define CAL_HISTOGRAMS_HH__

#include <iostream>
#include <iomanip>

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TEnv.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdexcept>

#include "Settings.hh"
#include "Gretina.hh"
#include "HiCARI.hh"

using namespace std;

#define CUTFILE_NAMELENGTH 100

class GretinaCalc;

class CalHistograms {
public:
  CalHistograms(Settings* set = NULL, int nentries=100000){
    Init(set,nentries);
  }
  ~CalHistograms(){
    delete fhlist;
  }
  void Init(Settings* set,int nentries){
    if (set!=NULL){
      fSett = set;
    } else {
      fSett = new Settings;
    }
    fhlist = new TList;
    fnentries = nentries;
    fentry =0;
  }
  void AddEntry(){
    fentry++;
  }
  TList* GetHList(){return fhlist;}
  void Write(){
    fhlist->Sort();
    fhlist->Write();
  }
  void FillHistograms(HiCARICalc* hi);
  void FillHistograms(GretinaCalc* gr);

  void FillI(string name,int bins, double low, double high, double value){
    try{
      fhmap1d.at(name)->Fill(value);
    } catch(const out_of_range& e) {
      TH1I* newHist = new TH1I(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value);
      fhlist->Add(newHist);
      fhmap1d[name] = newHist;
    }
  }
  void FillS(string name,int bins, double low, double high, double value){
    try{
      fhmap1d.at(name)->Fill(value);
    } catch(const out_of_range& e) {
      TH1S* newHist = new TH1S(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value);
      fhlist->Add(newHist);
      fhmap1d[name] = newHist;
    }
  }
  void Fill(string name,int bins, double low, double high, double value, double weight =1){
    try{
      fhmap1d.at(name)->Fill(value,weight);
    } catch(const out_of_range& e) {
      TH1F* newHist = new TH1F(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value,weight);
      fhlist->Add(newHist);
      fhmap1d[name] = newHist;
    }
  }
  void Fill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY, 
	    double weight =1){
    try{
      fhmap2d.at(name)->Fill(valueX,valueY,weight);
    } catch(const out_of_range& e) {
      TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			       binsX,lowX,highX,
			       binsY,lowY,highY);
      newHist->Fill(valueX,valueY,weight);
      fhlist->Add(newHist);
      fhmap2d[name] = newHist;
    }
  }

protected:
  TList* fhlist;
  map<string,TH1*> fhmap1d;
  map<string,TH2*> fhmap2d;
  int fnentries;
  int fentry;
  Settings* fSett;


};

#undef CUTFILE_NAMELENGTH

#endif
