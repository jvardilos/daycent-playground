/*  swconst.h */

#define SEC_PER_DAY 86400
#define SEC_PER_HOUR 3600
#define HOURS_PER_DAY 24
#define NDAY 366               /* dimension for daily arrays */
#define NMONTH 12              /* # of months in a year */
#define MXSWLYR 21             /* Max # of soil layers (soil water model) */
#define CMXLYR 10              /* Max # of century soil layers */
#define MAXSTLYR 200           /* Max # of soil temperature layers */
#define NTDEPTHS 4             /* Max # of soil regions */
#define VALMISSING 99          /* Used for value missing from histweth file */
#define BAR2CM 1024            /* 1 bar = 1024 cm H2O */
#define PARTDENS 2.65          /* Particle Density (g/cm3) */
/* On Jan 17, 2013, Keough,Cynthia wrote: change definition from, -1.0 to -0.5
   From: Ernie Marx [mailto:erniemarx@gmail.com]
   Keith,
   The DayCent problem that resulted in the loss of denitrification during
   spring thaw has been identified and fixed. This required a code change
   in the model. The problem was related to the impact of litter on soil temperature.

   This issue had a large effect on denitrification when it was triggered.
   The problem appears when there is snow melting on frozen ground and
   litter carbon is >250 gC/m2. Annual N2O emissions could be reduced by
   as much as 50% in these situations.
*/
#define FRZSOIL -0.5           /* temperature at which soil is considered frozen (deg C) */
#define MIN_FRZN_COND 0.00005  /* minimum hydraulic conductivity in frozen */
                               /* soil (cm/sec) */

#define max(a,b)        (((a) > (b)) ? (a) : (b))
#define min(a,b)        (((a) < (b)) ? (a) : (b))

#ifndef PI
  #define PI 3.14159265358979323846264338327950288   /* pi */
#endif

#define FNSOIL "soils.in"
#define FNSITE "sitepar.in"

#define CONVLAI 80    /* biomass needed to produce an LAI of 1 (g/m**2) */
