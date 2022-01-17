#ifndef __SIMDEFS_H
#define __SIMDEFS_H

#define HICARI_ID                7         //hicari identifier
#define HICARI_TAG              0xf005ba11 //type identifier
#define GAMMA_G4SIM_ID          11     /* ID for simulated gamma data */
#define ZERO_PHYSDATA_ID        13     /* ID for zero deg data */
#define ZEROPHYSDATA_TYPETAG    0x0de90de9

#ifdef BACKGROUND
#define MAX_SIM_GAMMAS 100       /* max. simulated gammas per event */
#elif HIGHMULT
#define MAX_SIM_GAMMAS 40       /* max. simulated gammas per event */
#else
#define MAX_SIM_GAMMAS 10       /* max. simulated gammas per event */
#endif

#define FIRSTMB 201
#define MBCLUST 6
#define MBCRYST 3
#define MBSEGS  6
#define CLOVERS 4
#define CLCRYST 4
#define CLSEGS  8


#if MBCRYST > CLCRYST
#define CRYST MBCRYST
#else
#define CRYST CLCRYST
#endif

#if MBSEGS > CLSEGS
#define SEGS MBSEGS
#else
#define SEGS CLSEGS
#endif

struct sim_seg{
  int seg_id;
  float e;
};

struct sim_clust{
  int type;
  int crystal_id;
  int num;
  float tot_e;
  sim_seg seg[MBSEGS];
};

typedef struct g4sim_emitted_gamma{
  float e;
  float x, y, z;
  float phi, theta;
  float beta;
} EG;

typedef struct g4sim_abcd1234 {
  int type;          /* defined as abcd1234 */
  int num;           /* # of emitted gammas */
  int full;          /* is full energy */
  EG gammas[MAX_SIM_GAMMAS];
} G4SIM_EGS;

typedef struct ZD_physicsdata {
  int type;    /* defined 0deg0deg for indicating this version */
  float ata;    /* dispersive angle        */
  float bta;    /* non-dispersive angle    */
  float xta;    /* dispersive position     */
  float yta;    /* non-dispersive position */
  float betata; /* beta velocity           */
} ZD_PHYSICSDATA;

struct simbroken{
  int Detector;
  int Crystal;
  Short_t Segment;
  bool broken;
};

struct simresolution{
  // resolution = A*sqrt(1 + B*Energy)                                                                                                     
  int Detector;
  int Crystal;
  double A;
  double B;
  double C;
};

struct simthreshold{
  // threshold = 0.5*( 1 + tanh( (Energy - E) / dE ) )                                                                                     
  int Detector;
  int Crystal;
  double E;
  double dE;
};

#endif
