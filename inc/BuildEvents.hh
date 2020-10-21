#ifndef __BUILDEVENTS_HH
#define __BUILDEVENTS_HH
#include <vector>
#include <algorithm>

#include "TTree.h"

#include "Beam.hh"
#include "FocalPlane.hh"
#include "HiCARI.hh"
#include "Gretina.hh"
#include "Settings.hh"
#include "Globaldefs.h"
#include "MergeHistograms.hh"

/*!
  A container to keep track of the timestamps and corresponding detectors
*/
struct detector{
  //! timestamp of the detector hit
  unsigned long long int TS;
  //! ID for the detector, 0 BigRIPS, 1 HiCARI, 2 Mode2
  int ID;
};

/*!
  A class for building BigRIPS and WASABI combined events
*/
class BuildEvents {
public:
  //! default constructor
  BuildEvents(){};
  //! constructor                                                                                                                                                             
  BuildEvents(Settings* set);
  //! dummy destructor
  ~BuildEvents(){
  };

  //! Initialize trees
  void Init(TTree* brtr, TTree* hitr, TTree* m2tr);
  //! Set the window for event building
  void SetWindow(unsigned long long int window){fwindow = window;};
  //! Set correlation mode
  void SetCorrelationMode(int mode){fmode = mode;};
  //! Set the verbose level
  void SetVerbose(int verbose){fverbose = verbose;};
  //! Set the last event
  void SetLastEvent(int lastevent){flastevent = lastevent;};
  //! Read one entry from the HiCARI tree
  bool ReadHiCARI();
  //! Read one entry from the BigRIPS tree
  bool ReadBigRIPS();
  //! Read one entry from the Mode2 tree
  bool ReadMode2();
  //! Read one entry from each tree
  bool ReadEach();
  
  //! Merge the data streams
  bool Merge();

  //! Set the last event
  int GetNEvents(){return fHIentries + fBRentries + fM2entries;};
  //! Close the event and write to tree
  void CloseEvent();
  //! Get the merged tree
  TTree* GetTree(){return fmtr;};
  //! Get the histograms
  MergeHistograms* GetHistos(){return fmhist;};

private:
  //! Settings to control merger
  Settings* fSett;
  //! BigRIPS input tree
  TTree* fBRtr;
  //! HiCARI input tree
  TTree* fHItr;
  //! Mode2 input tree
  TTree* fM2tr;
  //! merged output tree
  TTree* fmtr;
  //! verbose level
  int fverbose;
  //! hasBR
  bool fhasBR;
  //! hasHI
  bool fhasHI;
  //! hasM2
  bool fhasM2;

  //! list of detector hits
  vector<detector*> fdetectors;

  //! number of events to be read
  int flastevent;

  //! number of bytes read
  int fnbytes;
  //! current timestamp
  unsigned long long int fcurrentts;
  //! last read BigRIPS timestamp
  unsigned long long int flastBRts;
  //! last read HiCARI timestamp
  unsigned long long int flastHIts;
  //! last read Mode2 timestamp
  unsigned long long int flastM2ts;
  //! timestamp jumped in BigRIPS
  bool fBRtsjump;
  //! timestamp jumped in HiCARI
  bool fHItsjump;
  //! timestamp jumped in Mode2
  bool fM2tsjump;

  //! BigRIPS timestamp
  unsigned long long int fBRts;
  //! bigrips checkADC
  int fcheckADC;
  //! bigrips trigbit
  int ftrigbit;
  //! bigrips data
  Beam* fbeam;
  //! bigrips focal plane information
  FocalPlane* ffp[NFPLANES];
  //! HiCARI timestamp
  unsigned long long int fHIts;
  //! HiCARI data
  HiCARICalc* fhicari;
  //! Mode2 timestamp
  unsigned long long int fM2ts;
  //! Mode2 data
  GretinaCalc* fmode2;

  //! local copy of BigRIPS timestamp
  unsigned long long int flocalBRts;
  //! local copy of trigbit
  int flocaltrigbit;
  //! local copy of checkADC
  int flocalcheckADC;
  //! local copy of bigrips data
  Beam* flocalbeam;
  //! local copy of bigrips focal plane information
  FocalPlane* flocalfp[NFPLANES];
  //! local copy of HiCARI timestamp
  unsigned long long int flocalHIts;
  //! local copy of HiCARI data
  HiCARICalc* flocalhicari;
  //! Mode2 timestamp
  unsigned long long int flocalM2ts;
  //! Mode2 data
  GretinaCalc* flocalmode2;

  //! number of bigrips entries
  double fBRentries;
  //! number of HiCARI entries
  double fHIentries;
  //! number of Mode2 entries
  double fM2entries;
  //! current bigrips entry
  unsigned int fBRentry;
  //! current HiCARI entry
  unsigned int fHIentry;
  //! current Mode2 entry
  unsigned int fM2entry;

  //! time window for eventbuilding
  unsigned long long fwindow;
  //! modus for writing the merged data: 0 all, 1 only particle/gamma (BR and Hi)
  int fmode;

  //! Histograms for merged data
  MergeHistograms* fmhist;
  
  
};

/*!
  A class for sorting the detectors in time
*/
class TSComparer {
public:
  //! comparator
  bool operator() ( detector *lhs, detector *rhs) {
    return (*rhs).TS > (*lhs).TS;
  }
};
#endif
