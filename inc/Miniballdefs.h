#define MINIBALL_ID  7 //miniball identifier
#define MINIBALL_TAG 0xf005ba11 //type identifier


#define FIRSTMB 201
#define MBCLUST 8
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
