#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <fstream>

#include "global.h"
#include "HFC.h"

using namespace std;

struct Mode3event
{
  int length;
  int board_id;
  int chn_id;
  int module;
  long long LED_ts;
  int en;
  bool en_sign;
  int pileup;
  long long CFD_ts;
  int CFD_1;
  int CFD_2;
  int trace_len;
  INT16 trace[8192];
};

void swapbytes(char* a, char *b)
{
  char tmp=*a;
  *a=*b;
  *b=tmp;
}

// Mode 3 data is high endian format
void 
HEtoLE(char* cBuf, int bytes)
{
  for(int i=0; i<bytes; i+=2)
    swapbytes((cBuf+i), (cBuf+i+1));
}

void 
Mode3Event(char* cBuf, int len, Mode3event* mode3)
{
  // Welcome in the 16 bit world
  UINT16* wBuf= (UINT16*)cBuf;
  
  // digitizer saves in high endian format
  // we're low endian
  HEtoLE(cBuf, len);

  // length now in units of 16bit words
  len/=2;

  // 1st & 2nd word are 0xaaaa
  if((*wBuf != 0xaaaa) &&
     (*(wBuf+1) != 0xaaaa))
    {
      cerr << "0xAAAA header missing" << endl;
      return;
    }
  wBuf+=2;

  // 3rd word
  // Digitizer reports length in 32bit units
  // we convert in 16bit. Furthermore the length
  // doesn't account for the two 0xAAAA words

  mode3->length = (*wBuf & 0x07ff) * 2 + 2;
  if(mode3->length != len)
    {
      cerr << "inconsistent mode3 buffer length "
	   << "Geb Hdr: " << len << " wrds "
	   << "Mode2: " << mode3->length << " wrds"
	   <<endl;
      return;
    }

  // also board id encoded (=GA ?)
  mode3->board_id = *wBuf >> 11;
  wBuf++;
  
  // 4th word
  mode3->chn_id = *wBuf & 0x000f;
  mode3->module = *wBuf >> 4;
  wBuf++;

  // 5th, 6th and 8th word LED timestamp
  // 5th 'middle', 6th 'lowest', 8th 'highest'

  mode3->LED_ts =
    ((long long) *(wBuf+3)) << 32;
  mode3->LED_ts +=
    ((long long) *(wBuf+0)) << 16;
  mode3->LED_ts +=
    ((long long) *(wBuf+1)) << 0 ;
  wBuf+=2; //point 7th
  
  // 7th is low 16bit energy
  // 10th upper 8bit, 9th bit sign
  mode3->en =
    (int) *(wBuf+3) & 0x00ff;
  mode3->en = mode3->en << 16;
  mode3->en += *wBuf;
  
  mode3->en_sign = *(wBuf+3) & 0x0100;
  mode3->pileup  = (*(wBuf+3) & 0x8000) >> 15;
  wBuf+=2; //point 9th

  // 9th, 11th and 12th word CFD ts
  // 9th 'lower, 11th 'highest', 12th 'middle'
  mode3->CFD_ts =
    ((long long) *(wBuf+2)) << 32;
  mode3->CFD_ts +=
    ((long long) *(wBuf+3)) << 16;
  mode3->CFD_ts +=
    ((long long) *(wBuf+0)) << 0 ;
  wBuf+=4; //point 13th

  // 13th, 14th CFD point 1
  mode3->CFD_1 = 
    (int) *(wBuf+1) << 16;
  mode3->CFD_1 +=
    (int) *wBuf;
  wBuf+=2; //point 15th
  
  // 15th, 16th CFD point 2
  mode3->CFD_2 = 
    (int) *(wBuf+1) << 16;
  mode3->CFD_2 +=
    (int) *wBuf;
  wBuf+=2; //point 17th

  // wBuf points at 1st trace element now
  mode3->trace_len = 
    mode3->length - 16; 

  for(int i=0; i<mode3->trace_len/2; i++)
    {
#define OFFSET 512;
      mode3->trace[2*i+1]  =-(INT16)(*(wBuf+1)) + OFFSET;
      mode3->trace[2*i+0]=-(INT16)(*(wBuf+0)) + OFFSET;
      wBuf+=2;
    }

  //cerr << hex
  //   << " LED: 0x" << mode3->LED_ts
  //   << " CFD: 0x" << mode3->CFD_ts
  //   << dec << endl;
  //
  

}

void BrowseData(gebData header)
{
  cerr 
    << "type:" << header.type 
    << " len: " << header.length
    << " ts: 0x" << hex << header.timestamp << dec
    << endl;
}


int HFC_mode3(BYTE* cBuf, HFC* hfc_list)
// return value: processed data in bytes
{
  long long mode3_ts;
  int mode3_len;
 
  
  // 15th and 16th byte is ts high (and endian)
  mode3_ts = 
    ((long long) cBuf[14]) << 40;

  mode3_ts += 
    ((long long) cBuf[15]) << 32;

  // 9th, 10th ts middle 
  mode3_ts += 
    ((long long) cBuf[8]) << 24;

  mode3_ts += 
    ((long long) cBuf[9]) << 16;

  // 11th, 12th ts 'low'
  mode3_ts += 
    ((long long) cBuf[10]) << 8;

  mode3_ts += 
    ((long long) cBuf[11]);


  //
  // 5th, 6th is length, in 32bit units
  mode3_len = 
    ((int)(cBuf[4])) << 8;
  mode3_len +=
    ((int)cBuf[5]);

  mode3_len &= 0x7ff;

  mode3_len *= 4; // convert in bytes

  hfc_list->add(mode3_ts, 2, mode3_len+4, cBuf);

  return (mode3_len + 4); // 0xaaaa 0xaaaa not counted in mode3_len
}



int main(int argc, char** argv)
{
  if(argc==1)
    {
      cerr << argv[0] << " file:"
	   << endl
	   << "brings GRETINA Mode3 event file"
	   << endl
	   << "in proper sequence"
	   << endl;
      exit(0);
    }

  bool pipeflag = false;
  FILE* out;
  if (argc>1 && !strcmp(argv[1],"-p"))
    {
      argc--;
      argv++;
      out = stdout;
      pipeflag = true;
    } else {
    out = fopen64("HFC.dat", "wb");
  }
  
  FILE* in;
  if (argc>1){
    string filename=argv[1];
    in =  fopen64(filename.c_str(), "rb");
    if(!in)
      {
	cerr << "HFC: cannot open file " << filename << endl;
	return 1;
      }
  } else {
    in = stdin;
  }

  
  long long totread;
  int read;
  int EvtCount=0;
  BYTE cBuf[8*16382];
  gebData aGeb;
  HFC hfc_list(50*8192, out);
  // 972: strange mode 2 with mem depth 40*8192 needed
  //

  bool success=true;
	    
  while(fread(&aGeb, sizeof(struct gebData), 1, in)==1)
    {
      read = fread(cBuf, sizeof(char), aGeb.length, in);
      totread+=read+sizeof(struct gebData);
      if(read!=aGeb.length)
	{
	  cerr << aGeb.length << " bytes expected but"
	       << read << " bytes read. Bailing out"
	       << endl;
	  cerr.flush();
	  break;
	}

      EvtCount++;
		
      if((EvtCount % 200)==0 && !pipeflag)
	{
	  cerr << "Event " << EvtCount
	       << " read:" << read
	       << " total read:" << totread/1000000
	       << " Mb \r";
	  cerr.flush();
	}
#if(0)
      cerr << "Event:" << EvtCount 
	   << " #data:" << read
	   << " geb:" << sizeof(struct gebData)
	   << " total bytes read: " << totread
	   << " (0x" << hex << totread << dec << ")"
	   << endl;
#endif	
	
      switch(aGeb.type)
	{
	case 1: // Mode 2, so far always GEB, 1 event, GEB, 1 event,....
	case 6: // S800 scaler etc. before Jan 2013
	case 10: // S800 scaler etc. after Jan 2013
 
	case 8: // card 29

	case 9:  // S800_PHYSICSDATA
	case 11: // GRETINA_G4SIM simulated gamma information
  
	case 5: // S800 event
	  if(!hfc_list.add(aGeb, cBuf) && success)
	    {
	      success = false;
	      cerr << "HFC: adding event in HFC failed with aGeb.type = "
		   << aGeb.type << ", length = " << aGeb.length
		   << endl;
	    }
	  break;

	case 2: // Mode 3 (raw)
	  {
	    int nBytes = aGeb.length;
	    BYTE *data;
	    data=cBuf;
	    while(nBytes)
	      {
		int nread;
		nread = HFC_mode3(data, &hfc_list);
			  
		data+=nread;
		nBytes -=nread;
			  
		if(nBytes < 0)
		  {
		    cerr << "HFC: nBytes negative!!"
			 << endl;
		    exit(0);
		    break;
		  }
	      }
	  }
	  break;
		    
	case 3: // Mode 1
	  cerr << "HFC: Not implemented Mode 1 (GEB type 3)" << endl;
	  break;

	default:
	  cerr << "HFC: Unknown packet type " << aGeb.type
	       << "HFC: skipping that one" << endl;
	}
		
    }

  cerr << "HFC: calling flush" << endl; cerr.flush();
  hfc_list.flush();
  hfc_list.printstatus();

  cerr << "HFC: closing files" << endl; cerr.flush();
  fclose(in);
  fclose(out);
  cerr << "HFC: done" << endl; cerr.flush(); 

  return 0;
}

	 
