#ifndef __BUILDEVENTS_HH
#define __BUILDEVENTS_HH
#include <vector>
#include <algorithm>

#include "TTree.h"

#include "Beam.hh"
#include "FocalPlane.hh"
#include "HiCARI.hh"
#include "Globaldefs.h"


/*!
  A container to keep track of the timestamps and corresponding detectors
*/
struct detector{
  //! timestamp of the detector hit
  unsigned long long int TS;
  //! ID for the detector, 0 BigRIPS, 1 HiCARI
  int ID;
};

/*!
  A class for building BigRIPS and WASABI combined events
*/
class BuildEvents {
public:
  //! Default constructor
  BuildEvents(){};
  //! Initialize trees
  void Init(TTree* brtr, TTree* hitr);
  //! Set the window for event building
  void SetWindow(unsigned long long int window){fwindow = window;};
  //! Set coincidence mode
  void SetCoincMode(int mode){fmode = mode;};
  //! Set the verbose level
  void SetVerbose(int verbose){fverbose = verbose;};
  //! Set the last event
  void SetLastEvent(int lastevent){flastevent = lastevent;};
  //! Read one entry from the HiCARI tree
  bool ReadHiCARI();
  //! Read one entry from the BigRIPS tree
  bool ReadBigRIPS();
  //! Read one entry from each tree
  bool ReadEach();
  
  //! Merge the data streams
  bool Merge();

  //! Set the last event
  int GetNEvents(){return fHIentries + fBRentries;};
  //! Close the event and write to tree
  void CloseEvent();
  //! Get the merged tree
  TTree* GetTree(){return fmtr;};
  
private:
  //! BigRIPS input tree
  TTree* fBRtr;
  //! HiCARI input tree
  TTree* fHItr;
  //! merged output tree
  TTree* fmtr;
  //! verbose level
  int fverbose;
  //! hasBR
  bool fhasBR;
  //! hasHI
  bool fhasHI;

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
  //! timestamp jumped in BigRIPS
  bool fBRtsjump;
  //! timestamp jumped in HiCARI
  bool fHItsjump;

  //! BigRIPS timestamp
  unsigned long long int fBRts;
  //! bigrips data
  Beam* fbeam;
  //! bigrips focal plane information
  FocalPlane* ffp[NFPLANES];
  //! HiCARI timestamp
  unsigned long long int fHIts;
  //! HiCARI data
  HiCARICalc* fhicari;

  //! local copy of BigRIPS timestamp
  unsigned long long int flocalBRts;
  //! local copy of bigrips data
  Beam* flocalbeam;
  //! local copy of bigrips focal plane information
  FocalPlane* flocalfp[NFPLANES];
  //! local copy of HiCARI timestamp
  unsigned long long int flocalHIts;
  //! local copy of HiCARI data
  HiCARICalc* flocalhicari;

  //! number of bigrips entries
  double fBRentries;
  //! number of HiCARI entries
  double fHIentries;
  //! current bigrips entry
  unsigned int fBRentry;
  //! current HiCARI entry
  unsigned int fHIentry;

  //! time window for eventbuilding
  unsigned long long fwindow;
  //! modus for writing the merged data: 0 all, 1 only particle/gamma (BR and Hi)
  int fmode;
  
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
