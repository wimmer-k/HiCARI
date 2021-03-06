#ifndef RAW_HISTOGRAMS_HH__
#define RAW_HISTOGRAMS_HH__

#include <iostream>
#include <iomanip>

#include "TObject.h"
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
#include "Trace.hh"

using namespace std;


class RawHistograms {
public:
  RawHistograms(Settings* set = NULL,int nentries=100000){
    fentry = -1;
    fhlist = new TList;
    fnentries = nentries;
    fSett = set;
  }
  ~RawHistograms(){
    delete fhlist;
  }
  
  TList* GetHList(){return fhlist;}
  void Write();

  void FillHistograms(Mode3Event* m3e, HiCARI* hi, Gretina* gr);

  void Fill(string name,int bins, double low, double high, double value){
    try{
      fhmap.at(name)->Fill(value);
    } catch(const out_of_range& e) {
      //cout << "New 1-d histogram: " << name << endl;
      TH1F* newHist = new TH1F(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value);
      fhlist->Add(newHist);
      fhmap[name] = newHist;
    }
  }
  void Fill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY){
    try{
      fhmap.at(name)->Fill(valueX,valueY);
    } catch(const out_of_range& e) {
      //cout << "New 2-d histogram: " << name << endl;
      TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			       binsX,lowX,highX,
			       binsY,lowY,highY);
      newHist->Fill(valueX,valueY);
      fhlist->Add(newHist);
      fhmap[name] = newHist;
    }
  }

protected:
  void FillMode3Histograms(Mode3Event* m3e);
  void FillHiCARIHistograms(HiCARI* hi);
  void FillGretinaHistograms(Gretina* gr);

  TList* fhlist;
  map<string,TH1*> fhmap;
  Settings* fSett;
  int fnentries;
  int fentry;
};

#endif
