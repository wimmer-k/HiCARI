#ifndef __HFC_H 
#define __HFC_H

#include <iostream>
#include <list>
#include "global.h"

using namespace std;

//
// HFC for preproccessing a file made by GEB
// (Gretina Event Builder)
// GEB does sometimes odd things:
// 1) In any mode events out of time sequence appear
// 2) In Mode 3 (traces) several digitizer channels
//    are packed in one payload with one GEB header
//    instead GEB, a channel, GEB, a channel, GEB, ...
//
// HFC takes events (geb header, data) and write them
// in time-ordered manner. It holds a user-defined
// number of events in memory for performing re-ordering.
// If memory is too small, events which do not fit will
// be dicarded.
//  
// User methods:
//
// constructers
// HFC(int num, FILE* out)
// HFC(int num)
// HFC(FILE* out)
// HFC()
// num is user-defined memory depth
// out is file pointer where output is written (default NULL)
//  
// 
// bool add(gebData aGeb, BYTE* data);
// sending event to HFC. TS from aGeb is 
// used for time ordering. Returns false 
// if event can't be processed as its TS 
// is too 'old'
//
// bool add(long long TS, int type, int length, BYTE* data);
// sending event to HFC. TS used for time ordering.
// Same return value as for other add method.
//
// void flush();
// flushs the events held im memory to file. IMPORTANT: This
// method has to be used before closing output file or 
// destroying an HFC object.  
//
// void printstatus();
// prints statistics 
//

struct gebData
{
  int type;
  int length; // payload in bytes
  long long timestamp;
};


struct HFC_item
{
  gebData geb;
  BYTE*  data;
};



class HFC
{
 private:
  int m_memdepth;
  FILE *m_file;
  int m_evt;

  list<HFC_item*> m_HFClist;
  list<HFC_item*>::iterator m_HFClast_it;

  int m_discarded;
    
  void init(int num, FILE* out);

  // we won't expand list, but discard one
  // item by writing it.
  bool addToFullList(HFC_item* hfc);
  
  // Well, this actually does the writing
  bool writeItem(HFC_item* hfc);

  // An item gets inserted at 'right' 
  // position in list, list gets expanded
  bool insert(HFC_item* hfc);

 public:
  HFC();
  HFC(int num);    //num == #events in memory
  HFC(FILE* out);
  HFC(int num, FILE* out);

  // 'user' method for adding event to list
  // data will be copied, so user can free/
  // change data block she/he is pointing to.
  bool add(gebData aGeb, BYTE* data);
  bool add(long long TS, int type, int length, BYTE* data);

  // 'user' method for writing all list items
  // stored in memory
  void flush();

  //
  void printstatus();
};



#endif
