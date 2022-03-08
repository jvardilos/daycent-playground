/*----------------------
   define the maximum freeze/thaw pulse length
   DEFAULT  0 days; no freeze/thaw pulse
   --TPL
     tplen     Thaw N maximum Freeze Thaw Period in days
               NOTE: regardless of trigger, this value is required fo any pulse.

   thaw Pulse snow pack trigger
   DEFAULT  0 FDD; no freeze/thaw pulse
   --PST       Thaw Pulse Length
     tnsad;    snow depth required to set (arm) the pulse
     tnstd;    snow depth below which the pulse is started
               if <0 this is the amount of snow melt which triggers the pulse

   thaw Pulse Freezing Degree Day pulse width limiting and/or trigger.
   pulse length varies between 1 and tplen days as FDD varies between  plmfdd and plxfdd
   If this is specified with a snow trigger it modifies the pulse width
   With no trigger, a thaw in ptrlyr will trigger a pulse.
   DEFAULT  0 FDD; no freeze/thaw pulse shaping or trigger
   --PFT
     plmfdd   minimum FDD needed for a freeze/thaw pulse
              if <0 the a pulse will NOT clear the FDD sum
     plxfdd   FDD needed for a full length pulse
     ptrlyr   soil layer thaw which can trigger a pulse; default 1;
              if <0 the layer does not need to thaw, the temperature must be just be >0

   thaw Pulse denitrification x-inflection point shift
   If FDD limits, dpsfdm and dpsfdx, are given the inflection will linearly change
   from 0 to tnpsif as the FDD sum changes between the limits.
   --PDS
     tnpsif   maximum denitrification x-inflection point shift during pulse
     dpsfdm   Minimum FDD limit for a Denitrification Pore Space shift
     dpsfdx   FDD limit for a full Denitrification Pore Space shift

   saturation fraction during freeze/thaw pulse.
   DEFAULT  no modification of saturation effect
   During thaw pulse this specifies the fraction of the soil subject to the
   flooded N2toN2O ratio
   -—PSF
      sn2r     thaw pulse saturated fraction

   remove respiration restrain above the input layer during pulse
   DEFAULT  ON
   --ntrrd
     tndrrl    deepest layer to remove respiration restraint
               -1 remove restrain for water > field capacity

   change decomposition rate during freeze/thaw pulse
   DEFAULT  2 parameters required
            not specified no change in decomposition
   --TDP
     decrate   multiplier on the nominal decomposition rate
     <decpool> a comma delimited list of the pools to modify
             1  surface metabolic
             2  surface structural
             3  surface som2
             4  surface wood
             5  soil metabolic
             6  soil structural
             7  soil som2 and som3
             8  soil wood

----------------------*/

/*routines needed in denitrify */
  /* set Water Filled Pore Space x_inflection shift based on freezing degree days */
double thwplse(int jday);
  /* retrieve denitrification respiration restraint layers and pulse saturated fraction*/
void   denitrrpuls(int jday, int *thpulsrr, float *thpsat);

/* do the daily pulse book keeping; called from watrflow */
void chkthwpulse(float snowpack, float time, int jday, float melt, float soiltavg[MXSWLYR]);

/* other package routines */
void setdrrl(int *tndrrli);    /* set denitrification respiration restraint layer */
void setsatfr(float *tnsatfi); /* N2N2O ratio saturated fraction */
void sthwplsd(int *tpleni);    /* set thaw pulse maximum length */;
void ssnowtrig(float tnsadi, float tnstdi); /*sets the snow pulse trigger values */
void sfddtrig(float plmfddi, float plxfddi, int ptrlyri); /*freezing degree day pulse trigger values */
void stdpsf(double tnpsifi); /* denitrification Pore Space inflection Factor */
void stdpsfdd(float dpsfdmi, float dpsfdxi); /* denitrification pore space shift freezing degree day limits */
/*
void   thawdecp(int *jday, float anerdcmp[])
void sethawdec(float ithwdcmp[]);
*/
