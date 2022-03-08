/*******************************************************************************
**
**  FILE:    sitsoil.h
**
**  AUTHOR: Kendrick Killian  March 17 2014
**   This strips the soil layer and sitepar structures from the original
**   soilwater.h header file so they can be made available to the fortran
**   routines. To do so, these definitions MUST be converted to fortran syntax.
**   Making this a seperate file facilitates the conversion to the Fortran 2003
**   fortran_iso_C_binding.
**
**   include this in soilwater.h
**
**   Currently we use code to convert these to Fortran common blocks, let the
**   fortran allocate the memory then access the memory in C using as C extern
**   structures. Doing this changes the allocation and pointer statements but
**   nothing else in the C side.
**
**   This is plan B. The original intent was to pass the struct locations to
**   Fortran. For now, the syntax needed to pass struct locations escapes me
**   in the layers of indirection and C vs fortran pointers.
**
********************************************************************************/

#ifndef  MXSWLYR
 #include "swconst.h"
#endif


/* typedef struct {
    int verbose;
    int debug;
  } FLAG_S, *FLAG_SPT;
*/


typedef struct
{
  double swc[MXSWLYR];
  double swcmin[MXSWLYR];
/* integer values put here to allow 64bit integers*/
  int    numlyrs;
  int    nelyrs;
  int    ubnd[CMXLYR];
  int    lbnd[CMXLYR];
/* values input directly drom the soil.in file */
  float  dpthmn[MXSWLYR];
  float  dpthmx[MXSWLYR];
  float  lyblkd[MXSWLYR];
  float  fieldc[MXSWLYR];
  float  wiltpt[MXSWLYR];
  float  ecoeff[MXSWLYR];
  float  tcoeff[MXSWLYR];
  float  sandfrac[MXSWLYR];
  float  clayfrac[MXSWLYR];
  float  orgfrac[MXSWLYR];
  float  satcond[MXSWLYR];
  float  lyrpH[MXSWLYR];
/* calculated values */
  float  width[MXSWLYR];
  float  depth[MXSWLYR];
  float  swclimit[MXSWLYR]; /* swclimit(numlyrs) = wiltpt(numlyrs) - deltamin */
  float  swcfc[MXSWLYR];    /* swcfc[ilyr]  = fieldc[ilyr] * width[ilyr]; */
  float  swcwp[MXSWLYR];    /* swcwp[ilyr]  = wiltpt[ilyr] * width[ilyr]; */
  float  minpot[MXSWLYR];
  float  wfps[MXSWLYR];
  float  thetas[MXSWLYR];
  float  thetas_bd[MXSWLYR];
  float  psis[MXSWLYR];
  float  b[MXSWLYR];
  float  sumecoeff;
} LAYERPAR_S, *LAYERPAR_SPT;


typedef struct
{
  double dDO_fc;            /* WFPS normalized diffusivity in aggregate soil at field capacity */
  double dDO_wp;            /* WFPS normalized diffusivity in aggregate soil at wilting point */
  int    usexdrvrs;
  int    texture;
  int    drainlag;
  int    jdayStart;
  int    jdayEnd;
  int    SnowFlag;          /* snow insulation effect on soil surface temp: 0 = not insulating, 1 = insulating */
  float  avgsand;           /* average sand to 10 cm. used for texture and methane */
  float  sublimscale;
  float  reflec;
  float  albedo;
  float  dmpflux;
  float  hours_rain;
  float  hpotdeep;
  float  ksatdeep;
  float  rlatitude;
/*  float  cldcov[NMONTH+1]; */
  float  sradadj[NMONTH];   /* solar radiation adjustment for cloud cover & transmission coefficent */
  float  dmp;
  float  Ncoeff;
/* changes for WFPS effect on N2O flux Cindy K 7 Aug 2013 */
  float  N2Oadjust_fc;      /* Maximum proportion of nitrified N lost as N2O .. field capacity */
  float  N2Oadjust_wp;      /* minimum proportion of nitrified N lost as N2O .. wilting point */
  float  MaxNitAmt;         /* maximum daily nitrification amount (gN/m^2) */
  float  netmn_to_no3;      /* fraction of new net mineralization that goes to NO3 (0.0-1.0) */
  float  wfpsdnitadj;       /* adjustment on inflection point for water filled pore space */
  float  N2N2Oadj;          /* multiplier on N2/N2O ratio */
/* variables for site specific daylength */
  float  elevation; /* elevation elevation, meters Brookings SDE */
  float  sitslp;    /* sitslp    site slope, degrees    sitslp rather than SLOPE constant*/
  float  aspect;    /* aspect    site aspect, degrees */
  float  ehoriz;    /* ehoriz    site east horizon, degrees */
  float  whoriz;    /* whoriz    site west horizon, degrees */
/* variables for methanogenesis model */
  int    floodN2delay;   /* TEST INPUT: days to transition to flooded N2/N2O ratio */
  float  flood_N2toN2O;    /* N2/N2O ratio for flooded state (100.0) (-1 disable) */
  float  CO2_to_CH4;       /* fraction of heterotrophic soil respiration to CH4 */
/* float  C6H12O6_to_CH4;  methane_production bug fix; constant from reaction dynamics not input KLK 10Mar15*/
  float  methzr;           /* methane fraction emitted via bubbles at 0 root mass*/
  float  frexud;           /* fraction of root production to root exudates (default)*/
  float  mrtblm;           /* root biomass that starts to reduce methane bubble formation (default)*/
  float  frCH4emit;        /* fraction of methane emitted by rice plants  KLK 27Jul15*/
  float  frac_to_exudates; /* fraction of root production to root exudates */
  float  Aeh;              /* methane differential coefficient (Aeh) */
  float  Deh;              /* methane differential coefficient (Deh) */
  float  Beh_flood;        /* methane Eh lower-limit during flooding */
  float  Beh_drain;        /* methane Eh upper-limit during drainage */
  float  zero_root_frac;   /* methane fraction emitted via bubbles at 0 root mass */
  float  ch4rootlim;       /* root biomass that starts to reduce methane bubble formation (default) */
} SITEPAR_S, *SITEPAR_SPT;
