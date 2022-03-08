/*
!               Copyright 2017 Colorado State University
!                        All Rights Reserved

! ****************************************************************************
!
!   FILE:      thawpulse.c
!
!   MODULE:    float thplse()
!
!   PURPOSE:   Modify the x_inflection for some period after thaw.
!
!   HISTORY:
!     04/10/2017 (KLK)
!
!   INPUTS:
!     atemp - average temperature for the day (degrees C)
!     jday  - day number
!
!   GLOBAL VARIABLES:
!     None
!
!   PUBLIC FUNCTIONs:
!       double thwplse     returns the factor to shift x_inflection
!              ssnowtrig  sets the internal state
!              chkthwpulse changes internal state depending on snowpack
!
!   package variables:
!       static int    daycnt     = 0;  number of days remaining in pulse
!       static int    datstdepth = 0;  date snowpack was last above the trigger level
!       static int    tplen    = 0;    maximum days for the pulse
!       static float  tnsad    = 5;    minimum snow pack needed to arm pulse
!       static float  tnstd    = 0;    snow pack depth needed to trigger pulse (below x or x loss)
!       static double tnpsif = 1.0;    factor to shift x_inflection
!    NOTES: static gives  persistence AND restricts their scope to this file
!           the default factors disable any effect
!
!   INPUTS:
!   define the maximum freeze/thaw pulse length
!   DEFAULT  0 days; no freeze/thaw pulse
!   --TPL
!     tplen     Thaw N maximum Freeze Thaw Period in days
!               NOTE: regardless of trigger, this value is required fo any pulse.
!
!   thaw Pulse snow pack trigger
!   DEFAULT  0 FDD; no freeze/thaw pulse
!   --PST       Thaw Pulse Length
!     tnsad;    snow depth required to set (arm) the pulse
!     tnstd;    snow depth below which the pulse is started
!               if <0 this is the amount of snow melt which triggers the pulse
!
!   thaw Pulse Freezing Degree Day pulse width limiting and/or trigger.
!   pulse length varies between 1 and tplen days as FDD varies between  plmfdd and plxfdd
!   If this is specified with a snow trigger it modifies the pulse width
!   With no trigger, a thaw in ptrlyr will trigger a pulse.
!   DEFAULT  0 FDD; no freeze/thaw pulse shaping or trigger
!   --PFT
!     plmfdd   minimum FDD needed for a freeze/thaw pulse
!              if <0 the a pulse will NOT clear the FDD sum
!     plxfdd   FDD needed for a full length pulse
!     ptrlyr   soil layer thaw which can trigger a pulse; default 1;
!              if <0 the layer does not need to thaw, the temperature must be just be >0
!
!   thaw Pulse denitrification x-inflection point shift
!   If FDD limits, dpsfdm and dpsfdx, are given the inflection will linearly change
!   from 0 to tnpsif as the FDD sum changes between the limits.
!   --PDS
!     tnpsif   maximum denitrification x-inflection point shift during pulse
!     dpsfdm   Minimum FDD limit for a Denitrification Pore Space shift
!     dpsfdx   FDD limit for a full Denitrification Pore Space shift
!
!   saturation fraction during freeze/thaw pulse.
!   DEFAULT  no modification of saturation effect
!   During thaw pulse this specifies the fraction of the soil subject to the
!   flooded N2toN2O ratio
!   -—PSF
!      sn2r     thaw pulse saturated fraction
!
!   remove respiration restrain above the input layer during pulse
!   DEFAULT  ON
!   --ntrrd
!     tndrrl    deepest layer to remove respiration restraint
!               -1 remove restrain for water > field capacity
!
!   change decomposition rate during freeze/thaw pulse
!   DEFAULT  2 parameters required
!            not specified no change in decomposition
!   --TDP
!     decrate   multiplier on the nominal decomposition rate
!     <decpool> a comma delimited list of the pools to modify
!             1  surface metabolic
!             2  surface structural
!             3  surface som2
!             4  surface wood
!             5  soil metabolic
!             6  soil structural
!             7  soil som2 and som3
!             8  soil wood
!
!   INPUT EXAMPLES
!    COMMAND LINE
!      --TPL -6 --PFT -200,285,-1 --PDS 1.0,200,285 --PST 5,-0.6 --ntrrd 3
!    fix.100 input format
!      # freeze pulse
!      # -6             TPLEN       # Thaw N Freeze Thaw Period in days
!      #  200.0         PLMFDD      # Pulse Length Minimum Freezing Degree Day
!      #  285.0         PLXFDD      # Pulse Length maXimum Freezing Degree Day
!      # 1              PTRLYR      # soil layer to sample for trigger
!      #   0.8          TNPSIF      # Pore Space x_inflection shift Factor
!      # 200.0          DPSFDM      # denitrification pore space frozen degree minimum
!      # 285.0          DPSFDX      # denitrification pore space frozen degree maximum
!      #   5.           TNSAD       # minimum Snow Arm Depth
!      #  -0.6          TNSTD       # snow pack needed to trigger pulse (<x or x loss)
!      # 3              TNDRRL      # remove denitrification respiration restraint above layer
!      #! 0              TNSATF      # pulse N2N2O saturated fraction
!      #! 10.           DECRAT(1,2,4,5)  # (removed) decomposition rate during pulse; indexes are effected pools
!
!
! ***************************************************************************/

/* Import needed interfaces: */
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <swconst.h>
#include <thawpulse.h>

/* package status variables: */

  static int    tndrrl     = 0;  /* denitrification respiration restraint off */

  static int    daycnt     = 0;  /* number of days remaining in pulse */
  static int    tplen      = 0;  /* Thaw N Freeze Thaw Period in days */

  static int    datstdepth = 0;  /* date snowpack was last above the trigger level */
  static float  tnsad      = 0;  /* minimum Snow Arm Depth freeze/thaw N release */
  static float  tnstd      = 0;  /* snow pack below which freeze/thaw N release starts */

  static float  rfstf     = 1;   /* Require FDD tigger surface freeze */
  static float  ftsfrz    = 0;   /* FDD trigger, surface frozen */
  static float  plreset   = 1;   /* Pulse Length FDD reset on trigger */
  static float  plmfdd    = 0;   /* Pulse Length Minimum Freezing Degree Day */
  static float  plxfdd    = 0;   /* Pulse Length maXimum Freezing Degree Day */
  static int    ptrlyr    = 1;   /* soil layer to sample for trigger */

  static int    sbsoilwd  = 0;   /* subsoil warm days, days since thaw */
  static float  frzdd     = 0;   /* running sum of frozen degree days */
  static float  tnsatf    = 0;   /* pulse N2N2O saturated fraction */

  static double tnpsif    = 1.0; /* Pore Space x_inflection shift Factor */
  static float  dpsfdm    = 0;   /* denitrification pore space frozen degree minimum */
  static float  dpsfdx    = 0;   /* denitrification pore space frozen degree maximum */

  static float  thawh2o[2] = {0.0,0.0}; /* previous days snowmelt */
  static int    thwptr     =0;          /* index into thaw water history */

  /* debug prints
       0 no print, 1 state transitions, 2 daily status, 3 pulse effects
     increment by setting the pulse length, tplen, <0; can be repeated
  */
  static int  debug = 0;

/* Decomposition effect status variables (removed)
  static int    thawdec    = 0;   Thaw decomposition on
  decomposition multiplier during freeze/thaw
  static float  thwdcmp[8] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,};
*/

/* routines to implement the decomposition rate changes (removed)
  void thawdecp(int *jday, float anerdcmp[]) {
    int ii;
if(debug>2) {printf("TNP; thawdecp day= %d   daycnt= %d %d\n",*jday, thawdec,daycnt);}
  if(thawdec==1  &&  daycnt > 0) {
      for(ii=0; ii<8; ii++) {anerdcmp[ii] *= thwdcmp[ii];}
      return 1;
    }
  }
  void sethawdec(float ithwdcmp[]) {
    int ii;
    thawdec = 1;
    for(ii=0; ii<8; ii++) {thwdcmp[ii] = ithwdcmp[ii];}
    printf("TNP; freeze thaw decomposition rate:");
    for(ii=0; ii<8; ii++) {printf("%.3f  ",thwdcmp[ii]);}
    printf("\n");
  }
*/

  void denitrrpuls(int jday, int *thpulsrr, float *thpsat) {
    *thpulsrr=  0;
    *thpsat =  -1;
    if(daycnt > 0) {
      *thpulsrr=  tndrrl;
      *thpsat =  tnsatf;
     if(debug>2) {printf("TNP; denitrrpuls day= %d   daycnt= %d  %d  %.3f\n",jday,daycnt, *thpulsrr, *thpsat);}
    }
    return;
  }
  void setdrrl(int *tndrrli) {
    tndrrl = *tndrrli;
    printf("TNP; remove denitrification respiration restraint above layer %d\n",tndrrl);
  }
  void setsatfr(float *tnsatfi) {
    tnsatf = *tnsatfi;
    printf("TNP; N2N2O ratio saturated fraction %.3f\n",tnsatf);
  }


  double thwplse(int jday)
  {
    /* set the shift based on freezing degree days */
    double plsxshft = 1.0; /* default factor; no pulse or less than lower limit */
    if(daycnt > 0 && frzdd > dpsfdm) {
      plsxshft = (frzdd >= dpsfdx)? tnpsif :                /* maximum response */
                 tnpsif * (frzdd-dpsfdm)/(dpsfdx-dpsfdm);   /* linear increase with FDD */
    }
    if(debug>2) {printf("TNP; thwplse  jday= %d  daycnt= %d plsxshft %.3f   FDD %.2f (%.2f, %.2f)\n",jday, daycnt,plsxshft, frzdd,dpsfdm,dpsfdx);}
    return plsxshft;
  }


  void sthwplsd(int *tpleni) { /* Set THaW PuLSe Duration */
     /*
        This routine sets the package values:
        tplen;  ! Thaw Pulse LENgth
        Jan 2018 modification removes the requirement fo
     */
     tplen  = abs(*tpleni); /* save without debug flag */
printf("TNP; duration: tplen= %d",tplen);
     if(*tpleni < 0) {debug++;printf("  logging pulses\n");} else {printf("\n");}
  }

  void ssnowtrig(float tnsadi, float tnstdi) {
     /*
        This routine sets the snow level pulse trigger values:
        tnsad;  ! tnsad Thaw N Snow Total Depth
        tnstd;  ! tnstd Thaw N Thaw Snow Depth
     */
     tnsad  = fabsf(tnsadi);
     tnstd  = tnstdi;
printf("TNP; snow triggers: arm tnsad= %.3f,",tnsad);
if(tnstd >=0){printf(" start below tnstd= %.3f snow\n",tnstd);}
else         {printf(" start with  tnstd= %.3f snow melt\n",fabs(tnstd));}
  }

  void sfddtrig(float plmfddi, float plxfddi, int ptrlyri) {
     /*
        This routine sets the freezing degree day pulse trigger values:
        tnmnfd;  ! Thaw N Minimum Freezing Degree days for pulse
        tnmxfd;  ! Thaw N MaXimum Freezing Degree days for pulse
     */
     if(plmfddi <0) {plreset = 0;}
     plmfdd  = fabsf(plmfddi);
     plxfdd  = plxfddi;
     if(ptrlyri <0) {rfstf = 0;}
     ptrlyr = abs(ptrlyri) -1;
     if(ptrlyr<0  || ptrlyr>MXSWLYR) {ptrlyr = 0;}
printf("TNP; freezing degree day pulse length plmfdd= %.3f, plxfdd= %.3f triggering on lyr %d\n",plmfdd,plxfdd,1+ptrlyr);
if(plreset<1)  {printf("TNP; freezing degree day sum NOT RESET on pulse\n");}
if(! rfstf)    {printf("TNP; FDD trigger does not require frozen soil\n");}
  }

  void stdpsf(double tnpsifi) {
     /*
        This routine sets the denitrification Pore Space inflection Factor:
        tnpsif;  ! tniwf Thaw N Pore Space inflection Factor
     */

     tnpsif  = tnpsifi;
printf("TNP; denitrification pore space inflection tnpsif= %.3f\n",tnpsif);
  }


  void stdpsfdd(float dpsfdmi, float dpsfdxi) { /* set denitrification pore space shift freezing degree day limits */
     /*
        This routine sets the package values:
        dpsfdm; ! denitrification pore space frozen degree minimum
        dpsfdx; ! denitrification pore space frozen degree maximum
     */
     dpsfdm = fabsf(dpsfdmi);
     dpsfdx = fabsf(dpsfdxi);
     frzdd  = 0;

printf("TNP; frozen degree day denitrification pore space shift limits dpsfdm= %.3f, dpsfdx= %.3f\n", dpsfdm,dpsfdx);
  }


  void chkthwpulse(float snowpack, float time, int jday, float melt, float soiltavg[MXSWLYR]) {
      /*
         do the daily pulse book keeping; called from watrflow
      */
    float snowmelt;
    int soilthw, spring;
    double trigtemp;

    if(tplen == 0) {return;} /* skip this entire routine if no pulse length */

    /* temperature samples */
    trigtemp = soiltavg[ptrlyr]>0;             /* definition of temperature trigger (2-5 cm) */
    soilthw = soiltavg[1]>0 && soiltavg[2]>0;  /* definition of soil thaw */
    if(soilthw==0 && ftsfrz==0) {ftsfrz = 1; /* FDD trigger; flag that soil has frozen */
      if(debug && tnsad == 0) {printf("TNP; set FDD trig      day= %d\n",jday);}
    }
    if(soilthw){sbsoilwd++;}                   /* number of days since sub soil thaw*/
          else {sbsoilwd=0;}                   /* subsoil frozen */

    /* spring
       snow: trigger armed, soilt > 0 and no more snow
       FDD:  minimum FDD>0, sum >0 and subsoil not frozen for 2 weeks
    */
    spring = ((datstdepth > 0  &&  soiltavg[0] > 0  &&  snowpack == 0) || /* spring clear of an untriggered snow pulse armed */
              (plmfdd != 0     &&  frzdd > 0        &&  sbsoilwd >= 14));

    /* terminate pulse and clear accumulators as needed;
       the pulse completed yesterday so do this before evaluating today's conditions. */
    if(daycnt == 1 || spring) { /* pulse termination */
      /* pulse complete; clear control valriables */
      if(daycnt >= 1) {
        if(debug) {printf("TNP; ending            time,day= %.2f/%d   frzdd= %7.2f arm date %d\n",time,jday,frzdd,datstdepth);}
        if(plreset>0) {frzdd = 0;} /* reset freezing day sum; Only for FDD trigger (tnsad==0  &&  plmfdd>0 && plreset>0)*/
      }

      if(spring) {
        if(debug>0 && (datstdepth > 0 || (frzdd > plmfdd && plmfdd>0))) {
           printf("TNP; clearing triggers time,day= %.2f/%d  frzdd %.3f arm date %d srfctemp %.3f  snow %.3f  sbsoilwd %d\n",
                   time,jday, frzdd, datstdepth, soiltavg[0], snowpack,sbsoilwd);}
        frzdd = 0;       /* clear FDD sum*/
        ftsfrz = 0;      /* spring resets soil frozen flag */
      }
      daycnt = 0;      /* stop  the pulse */
      datstdepth = 0;  /* clear the snow trigger */
    }

    else if(daycnt >0) {
      daycnt--; /* decrement the daycount */

      /* abort pulse but don't reset the triggers */
      if((snowpack > tnsad  && tnstd > 0) || /*if snow trigger and snow over trigger*/
         (tnsad == 0  &&  plmfdd > 0  && soilthw ==0)) { /*FDD trigger (not snow and FDD width set) and soil is frozen*/
        daycnt = 0;
        if(debug>0) {printf("TNP; refreezing        day= %d\n",jday);}
      }
    }


    snowmelt        = melt + thawh2o[0] + thawh2o[1];
    thwptr          = (thwptr+1)%2;  /* increment the pointer for saving the melt  */
    thawh2o[thwptr] = melt;
    if(soiltavg[0]<=0) {frzdd += fabsf(soiltavg[0]);} /* add to the freezing days */

    if(debug>1) {printf("TNP chkthwpulse day= %3d snow %5.2f trgtemp  %7.2f surftemp, %7.2f frzdd %7.2f  melt %.3f %.3f arm: datstdepth %d frzdd %d  thaw sbsoilwd %d pulse day %d\n",
                        jday,snowpack, soiltavg[ptrlyr], soiltavg[0], frzdd,melt,snowmelt, datstdepth, frzdd > plmfdd, sbsoilwd, daycnt);}

    /* can we arm a snow pulse */
    if(snowpack > tnsad  &&  tnsad > 0  && datstdepth == 0) {
      datstdepth = jday; /* arm the pulse */
      if(debug>0) {printf("TNP; snow threshold day= %d, snowpack= %.3f, frzdd= %.3f\n",
      jday,snowpack,frzdd);}
    }

    /* start the pulse count */
    if(daycnt == 0     &&     /*1 if it isn't running (don't restart) */
       frzdd >= plmfdd &&     /*2 FDD pulse width > 0 NOTE: this leaves the pulse armed! */
       ((datstdepth > 0  &&   /* snow pulse armed */
              (snowpack <= tnstd   || /* snow less than trigger level   OR*/
               (tnstd < 0 && snowmelt > -tnstd))) || /* snow melt exceeds (negative) trigger level */
       (tnsad == 0  &&  plmfdd > 0  &&  trigtemp > 0  && (ftsfrz  || ! rfstf)))  /* a FDD trigger but not snow, and surface unfrozen */
      ) {

      daycnt = (plxfdd > plmfdd  &&  frzdd < plxfdd)? /* check for a freezing days adjustment to length */
               lroundf(tplen*((fabsf(frzdd) - plmfdd)/(plxfdd - plmfdd))): /* linearly interpolate between limits */
               tplen;  /* else use a full length pulse */

      /* clear FDD trigger surface frozen flag at START it can reset during the pulse*/
      ftsfrz = 0;

      if(daycnt>0  &&  debug>0) {
        printf("TNP; starting time,day= %.2f/%d for %2d days,  conditions  fdd %6.2f ",time,jday,daycnt,frzdd);
          if(tnsad == 0  &&  plmfdd > 0  &&  trigtemp > 0) {printf(" FDD trigtemp %.3f\n", soiltavg[ptrlyr]);}
          else if(tnstd < 0 && snowmelt > -tnstd)          {printf(" snowmelt     %.3f > %.3f\n", snowmelt, -tnstd);}
          else if(snowpack <= tnstd)                       {printf(" snow depth   %.3f < %.3f\n", snowpack, tnstd);}
      }
    }
  }
