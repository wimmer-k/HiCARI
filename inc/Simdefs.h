#ifndef __SIMDEFS_H
#define __SIMDEFS_H

#define GAMMA_G4SIM_ID          11     /* ID for simulated gamma data */
#define ZERO_PHYSDATA_ID        13     /* ID for zero deg data */
#define MINOS_PHYSDATA_ID       17     /* ID for minos data */
#define ZEROPHYSDATA_TYPETAG    0x0de90de9
#define MINOSXYZDATA_TYPETAG    0xaffec0c0
#ifdef USELISA
#define LISATARDATA_TAG        0x10ca7e17 //type identifier
#define LISA_PHYSDATA_ID       18     /* ID for multi data */
#define MAXTARGETS              20
#endif
#ifdef BACKGROUND
#define MAX_SIM_GAMMAS 100       /* max. simulated gammas per event */
#elif HIGHMULT
#define MAX_SIM_GAMMAS 40       /* max. simulated gammas per event */
#else
#define MAX_SIM_GAMMAS 10       /* max. simulated gammas per event */
#endif

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

typedef struct MINOS_data {
  int type;    /* defined affeec0c0 for indicating this version */
  float betare;    /* beta velocity at reaction  */
  float x, y, z;   /* reaction vertex  */
} MINOS_DATA;

#ifdef USELISA
typedef struct LISA_data {
  int type;    /* defined affec0c0 for indicating this version */
  int ntargets;
  float deltaE[MAXTARGETS];    /* energy losses */
  float x, y, z;   /* reaction vertex  */
  int reaction;   /* which target made the reaction */
} LISA_DATA;
#endif

#endif
