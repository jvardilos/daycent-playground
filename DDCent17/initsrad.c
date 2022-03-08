
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsrad.c
**
**  FUNCTION:  void initsrad()
**
**  PURPOSE:   Initialize the solar radiation submodel.  Build a 366-day
**             arrays of ttmax0, potential radiation and day length values for
**             use in the solar radiation calculations.
**
**  AUTHOR:    Adapted from Peter Thorton's code extracted from the Sipnet model
**             (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
** UPDATES:
**    KLK Apr  2018
**             - Added the c_shwave routine to file; Use the stored solar radiation
**             instead of redoing the calculation
**             -  Use the closed integral for the flat plate insolation
**             (2*(ahou  * sinegeom + sahou * cosegeom) as coded in shwave) instead
**             of the numerical approximation.
**             -  Renamed solar constant and recast inverse relative Earth-Sun distance
**             for readability
**             -  changed sample points in integration loop; eliminates the round-off
**             sensitivity at limits.
**
**  INPUTS:
**    aspect      - site aspect (degrees)
**    ehoriz      - site east horizon (degrees)
**    elevation   - site elevation (meters)
**    latitude    - site latitude (decimal degrees)
**    sitslp      - site slope (degrees)
**    whoriz      - site west horizon (degrees)
**    daylength[] - length of day light (hours)
**
**  GLOBAL VARIABLES:
**    DAYSOFF   - julian day offset of winter solstice (11.25)
**    eps       - minimum declination, radians
**    PI        - pi (3.14159265)
**    RADPERDAY - radians of Earth orbit per julian day (0.017214)
**    RADPERDEG - radians per degree (0.01745329)
**    SECPERRAD - seconds per radian of hour angle (13750.9871)
**    SRADDT    - timestep for radiation routine, seconds (600.0)
**    TBASE     - max inst. trans., 0m, nadir, dry atm, dimensionless (0.870)
**
**  EXTERNAL VARIABLES:
**    ctrl                     - program control variables
**    ctrl->indewpt            - dewpoint temperature input flag, 0=NO, 1=YES
**    params                   - site parameters
**    params->base_elev        - site elevation, meters
**    params->site_asp         - site aspect, degrees, (0=N,90=E,180=S,270=W)
**    params->site_ehoriz      - site east horizon, degrees
**    params->site_elev        - site elevation, meters
**    params->site_lat         - site latitude, decimal degrees, (- for south)
**    params->site_slp         - site slope, degrees
**    params->site_whoriz      - site west horizon, degrees
**    yeardata                 - year long data ararys of ttmax0, potential
**                               radiation, and day length
**    yeardata->daylength[]    - length of day light (seconds)
**    yeardata->flat_potrad[]  - average daylight flux density on a flat surface
**    yeardata->slope_potrad[] - average daylight flux density on sloped surface
**    yeardata->ttmax0[]       - maximum daily total transmittance
**
**  LOCAL VARIABLES:
**    am               - optical air mass
**    ami              - index for accessing value from optical airmass by
**                       degrees array, optam
**    bsg1             - first term in beam-slope geometry calculation
**    bsg2             - second term in beam-slope geometry calculation
**    bsg3             - thrid term in beam-slope geometry calculation
**    cbsa             - cosine of beam-slope angle
**    cosasp           - cosine of site aspect
**    cosdecl          - cosine of declination
**    cosegeom         - cosine of site latitude * cosine of declination
**    cosh             - cosine of hour angle
**    coshss           - daylight hours, -1 = 24 hours of daylight
**                                        1 =  0 hours of daylight
**    coslat           - cosine of site latitude
**    cosslp           - cosine of site slope
**    coszeh           - cosine of zenith angle for east horizon
**    coszwh           - cosine of zenith angle for west horizon
**    cza              - cosine of solar zenith angle
**    decl             - declination
**    dh               - hour-angle increment
**    dir_beam_topa    - extraterrestrial radiation perpendicular to beam
**    dir_flat_topa    - potential flat surface radiation, top of atmosphere
**    dt               - timestep
**    hh               - hour-angle loop index
**    ii               - loop control variable
**    hss              - hour angle at sunset (radians)
**    lat              - Site latitude (decimal degrees)
**    optam[]          - optical airmass by degrees
**    pratio           - pressure ratio
**    sc               - solar constant as a function of yearday (W/m^2)
**    sinasp           - sine of site aspect
**    sindecl          - sine of declination
**    sinegeom         - sine of site latitude * sine of declination
**    sinh             - sine of hour angle
**    sinlat           - sine of site latitude
**    sinslp           - sine of site slope
**    sum_flat_potrad  - total potential radiation on a flat surface for ideal
**                       horizons
**    sum_slope_potrad - sun between east and west horizons, and direct on
**                       slope
**    sum_trans        - daily total transmittance
**    t1               - first term in pressure ratio calculation
**    t2               - second term in pressure ratio calculation
**    trans1           - initial transmittance corrected for elevation
**
**  OUTPUTS:
**    ctrl->indewpt            - dewpoint temperature input flag, 0=NO, 1=YES
**    daylenght[]              - length of day light (hours)
**    params->site_asp         - site aspect, degrees, (0=N,90=E,180=S,270=W)
**    params->site_ehoriz      - site east horizon, degrees
**    params->site_elev        - site elevation, meters
**    params->site_lat         - site latitude, decimal degrees, (- for south)
**    params->site_slp         - site slope, degrees
**    params->site_whoriz      - site west horizon, degrees
**    yeardata->daylength[]    - length of day light (seconds)
**    yeardata->flat_potrad[]  - average daylight flux density on a flat surface
**    yeardata->slope_potrad[] - average daylight flux density on sloped surface
**    yeardata->ttmax0[]       - maximum daily total transmittance
**
**  CALLED BY:
**    initsw()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calcSrad.h"
#include "mssg.h"
#include "swconst.h"

CONTROL_S Control;                   /* program control variables */
CONTROL_SPT ctrl = &Control;
PARAMETER_S Parameters;              /* site parameters */
PARAMETER_SPT params = &Parameters;
DATA_S Site_Data;                    /* site climate data */
DATA_SPT sitedata = &Site_Data;
YEARDATA_S Year_Data;                /* 366-day arrays of ttmax0, potential */
YEARDATA_SPT yeardata = &Year_Data;  /* radiation, and day length */

    void initsrad(double elevation, double latitude, double sitslp,
                  double aspect, double ehoriz, double whoriz,
                  float daylength[NDAY])
    {

      int    ami, ii;
      double t1,t2;
      double pratio;
      double lat,coslat,sinlat,dt,dh,hh;
      double cosslp,sinslp,cosasp,sinasp;
      double cza,cbsa,coszeh,coszwh;
      double decl,cosdecl,sindecl,cosegeom,sinegeom,coshss,sinhss,hss;
      double bsg1,bsg2,bsg3;
      double sc,dir_beam_topa;
      double sum_flat_potrad, sum_slope_potrad, sum_trans;
      double cosh,sinh;
      double dir_flat_topa,am;
      double trans1;
      /* optical airmass by degrees
      optam ~= exp(0.000255563x3 -0.002171847x2 + 0.061897769x 1.058277298)*/
      double optam[21] = {2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07, 4.37,
                          4.72, 5.12, 5.60, 6.18, 6.88, 7.77, 8.90, 10.39,
                          12.44, 15.36, 19.79, 26.96, 30.00};

      /* Initial the variables in the control structure */
      ctrl->indewpt = 0;  /* Dewpoint temperature input? (0=NO, 1=YES) */

      /* Initialize site parameters */
      params->site_lat = latitude;    /* Site latitude, degrees */
      params->site_elev = elevation;  /* Site elevation, meters */
      params->site_slp = sitslp;      /* Site slope, degrees */
      params->site_asp = aspect;      /* Site aspect, degrees */
      params->site_ehoriz = ehoriz;   /* Site east horizon, degrees */
      params->site_whoriz = whoriz;   /* Site west horizon, degrees */

      /* do basic error checking on parameters */
      /* check dewpoint input flag */
      if ((ctrl->indewpt < 0) || (ctrl->indewpt > 1)) {
        printf("WARNING: input dewpoint flag should be 0 or 1: assuming");
        printf(" dewpoint input ...\n");
        ctrl->indewpt = 1;
      }

      /* error checking for site parameters */
      if ((params->site_lat < -90.0) || (params->site_lat > 90.0)) {
        abortmssg("site latitude must be in the range -90.0 - 90.0  degrees");
      }
      if (params->site_elev > 5000) {
        sprintf(msgbfr, "site elev (%.1lf m) > 5000: use meters",params->site_elev);
        abortmssg(msgbfr);
      }
      if (params->site_slp > 60.0) {
        sprintf(msgbfr, "site slope (%.1lf deg) > 60:",params->site_slp);
        abortmssg(msgbfr);
      }
      if (params->site_slp < 0) {
        abortmssg("site slope must be >= 0.0");
      }
      if ((params->site_asp < 0.0) || (params->site_asp > 360.0)) {
        abortmssg("site aspect must be in the range 0.0 - 360.0 degrees");
      }
      if ((params->site_ehoriz < 0.0) || (params->site_ehoriz > 180.0)) {
        abortmssg("site east horizon must be in the range 0.0 - 180.0 degrees");
      }
      if ((params->site_whoriz < 0.0) || (params->site_whoriz > 180.0)) {
        abortmssg("site west horizon must be in the range 0.0 - 180.0 degrees");
      }

      /* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
      t1 = 1.0 - (LR_STD * params->site_elev)/T_STD;
      t2 = G_STD / (LR_STD * (R/MA));
      pratio = pow(t1,t2);

      /* STEP (2) correct initial transmittance for elevation */
      trans1 = pow(TBASE,pratio);

      /* STEP (3) build 366-day array of ttmax0, potential rad, and */
      /* daylength */
      /* precalculate the transcendentals */
      lat = params->site_lat;
      /* check for (+/-) 90 degrees latitude, throws off daylength calc */
      lat *= RADPERDEG;
      if (lat > 1.5707) lat = 1.5707;
      if (lat < -1.5707) lat = -1.5707;
      coslat = cos(lat);
      sinlat = sin(lat);

      cosslp = cos(params->site_slp * RADPERDEG);
      sinslp = sin(params->site_slp * RADPERDEG);
      cosasp = cos(params->site_asp * RADPERDEG);
      sinasp = sin(params->site_asp * RADPERDEG);
      /* cosine of zenith angle for east and west horizons */
      coszeh = cos(1.570796 - (params->site_ehoriz * RADPERDEG));
      coszwh = cos(1.570796 - (params->site_whoriz * RADPERDEG));

      /* sub-daily time and angular increment information */
      dt = SRADDT;                /* set timestep */
      dh = dt / SECPERRAD;        /* calculate hour-angle step */

      /* begin loop through yeardays */
      for (ii=0 ; ii<365 ; ii++) {
        /* decl = asin(sin(eps)* cos(0.01720289*(jday + 10) + EOEcc*sin(0.01720289*(jday - 2))))
            ASCE approximates this as
            decl =  eps*sinf(twoPI*jday/365-1.39)       ~  23.4393*sinf(twoPI*(jday + 284)/365.25]
                 =  eps*sinf(twoPI*(jday-80.75)/365.25) = -eps*cosf(twoPI*(jday-10.50)/365.25)
        */
        decl = eps * sin(RADPERDAY*(ii-DNSS));
        cosdecl = cos(decl);
        sindecl = sin(decl);

        /* do some precalculations for beam-slope geometry (bsg) */
        bsg1 = -sinslp * sinasp * cosdecl;
        bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl;
        bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl;

        /* calculate daylength as a function of lat and decl */
        /* sunset hour angle cos(thetaz) = cos(lat)⋅cos(decl)⋅cosω+sin(lat)⋅sin(decl) where thetaz=0) */
        cosegeom = coslat * cosdecl;
        sinegeom = sinlat * sindecl;
        coshss   = -(sinegeom) / cosegeom;
        sinhss   = (fabs(coshss) >= 1.0)? 0.0f: sqrt(1.0 - coshss*coshss);
        hss      = atan2(sinhss,coshss);

        /* daylength (seconds) */
        yeardata->daylength[ii] = 2.0 * hss * SECPERRAD;
        daylength[ii] = (float)yeardata->daylength[ii] / SEC_PER_HOUR;

        /* solar constant with day's relative Earth-Sun distance (W/m^2)*/
        sc = ScW * (1.0f + EOEcc*cos(RADPERDAY*ii)); /* + 45.6*cos(RADPERDAY*jday)*/

        sum_flat_potrad = sc * SpRad * 2*(hss * sinegeom + sinhss * cosegeom);

        /* extraterrestrial radiation perpendicular to beam, for a timestep */
        dir_beam_topa = sc * dt;

        sum_trans = 0.0;
        sum_slope_potrad = 0.0;

        /* begin sub-daily hour-angle loop, from -hss to hss
           loop modified to sample the center of the interval KLK
        */
        for (hh=-hss+dh/2 ; hh<hss ; hh+=dh) {
          /* precalculate cos and sin of hour angle */
          cosh = cos(hh);
          sinh = sin(hh);

          /* calculate cosine of solar zenith angle */
          cza = cosegeom * cosh + sinegeom;

          /* calculate cosine of beam-slope angle */
          cbsa = sinh * bsg1 + cosh * bsg2 + bsg3;

          /* check if sun is above a flat horizon */
          if (cza > 0.0) {
            /* sun is above the ideal (flat) horizon, do the flat-surface calculations */
            /* to determine daily total transmittance, potential radiation for later */
            /* calculations of diffuse radiation */

            /* potential radiation for this time period, flat surface, top of atmosphere */
            dir_flat_topa = dir_beam_topa * cza;

            /* determine optical air mass */
            am = 1.0/(cza + 0.0000001);
            if (am > 2.9) {
              ami = (int)(acos(cza)/RADPERDEG) - 69;
              if (ami < 0) ami = 0;
              if (ami > 20) ami = 20;
              am = optam[ami];
            }

            /* correct potential radiation with instantaneous transmittance through optical
               air mass to get daily total transmittance */
            sum_trans += pow(trans1,am) * dir_flat_topa;

            /* keep track of whether this time step contributes to component 1 */
            /* (direct on slope) */
            if (((hh<0.0 && cza>coszeh) || (hh>=0.0 && cza>coszwh)) && cbsa>0.0) {
              /* component 1: sun between east and west horizons, and direct on slope. */
              sum_slope_potrad += dir_beam_topa * cbsa;
            }

          } /* end if sun above ideal horizon */
        } /* end of sub-daily hour-angle loop */

        /* calculate maximum daily total transmittance and daylight average */
        /* flux density for a flat surface and the slope */
        if (yeardata->daylength[ii]) {
           yeardata->ttmax0[ii] = sum_trans / sum_flat_potrad;
           yeardata->flat_potrad[ii] =
             sum_flat_potrad / yeardata->daylength[ii];
           yeardata->slope_potrad[ii] =
             sum_slope_potrad / yeardata->daylength[ii];
        } else {
           yeardata->ttmax0[ii] = 0.0;
           yeardata->flat_potrad[ii] = 0.0;
           yeardata->slope_potrad[ii] = 0.0;
        }
      } /* end of ii=365 days loop */

      /* force yearday 366 = yearday 365 */
      yeardata->ttmax0[365] = yeardata->ttmax0[364];
      yeardata->flat_potrad[365] = yeardata->flat_potrad[364];
      yeardata->slope_potrad[365] = yeardata->slope_potrad[364];
      yeardata->daylength[365] = yeardata->daylength[364];
      daylength[365] = (float)yeardata->daylength[365] / SEC_PER_HOUR;
    }


/*****************************************************************************
**
**
**  FILE:      c_shwave.c
**
**  FUNCTION:  float c_shwave()
**
**  PURPOSE:   return the solar radiation [MJ m-2 d-1] outside the atmosphere
**
** UPDATES:
**    KLK Apr  2018
**             Recoded to return the data calculated by initsrad. No use having this
**             done in two places.
**    KLK July 2017
**             Code upgrade:
**             Improved comments and variable traceability to Yao's reference
**             Everybody uses different declination parameters/approximations.
**             Settle on values from Jean Meeus: Astronomical Algorithms; those used
**             by the US Naval Observatory web site.
**    YZ  Oct 2013
**             Change the equations to ASCE (2005) method.
**    KLK May 2011
**             removed transcof. No point in multiplying then dividing by the
**             same value. This removes the list dependence on month but left
**             it in the argument list.
**             Reordered the sahou, cahou calculation to avoid recalculation of
**             the same results since cahou is a component of sahou.
**             Use float trigonometry calls instead of doing repeated
**             conversions between float and double and back.
**
**  INPUTS:
**    jday      - current Julian day (1-366)
**    month     - current month (1-12)
**    rlatitude - latitude of the site (in radians)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ws              - sunset hour angle
**    decl            - declination (radians)
**    tsirad          - total solar irradiance (ly/day) (
**                      was solrad changed to avoid confusion with input variable klk 2014
**
**  OUTPUTS:
**    exsolrad - short wave solar radiation outside the atmosphere
**
**  CALLED BY:
**    snowCent()
**    snowmodel()
**
**  CALLS:
**    None
**
**  Notes:
**   langley (Ly) used to measure solar radiation (or insolation).  Definition
**      1 thermochemical calorie / cm2  =  41840 J/m²  =  11.622 watt-hours/m².
**
**   Solar Constant:    1.367 kW/m2  = 4.9212 MJ/h/m2   = 1.96 langleys (Ly)/minute
**
**   Exposure day: 24/PI
**    Sun exposure for a plate rotating once per day
**    2*Integrate[Cos[t*2*Pi/24],{t, 0, 6}
**    or
**    2*Integrate[
**      Cos[decl] Cos[rlat] Cos[t*2*Pi/24] + Sin[decl] Sin[rlat],
**      {t, 0, ArcCos[-Tan[decl]*Tan[rlat]*12/Pi }]
**    where rlat or decl = 0
**
*****************************************************************************/
    float c_shwave(int month, float rlatitude, int jday)
    {

      double  exsolrad;

      /* Determine the short wave radiation outside the atmosphere */
      exsolrad  = (float)((yeardata->flat_potrad[jday] * yeardata->daylength[jday])/WtL);

/* This is the old calculation with the updated orbital/solar constants
      double  decl, dr, ahou, sahou, cahou, cosegeom, sinegeom;
      /! decl = asin(sin(eps)* cos(0.01720289*(jday + 10) + EOEcc*sin(0.01720289*(jday - 2))))
          ASCE approximates this as
          decl =  23.4393*sinf(twoPI*jday/365-1.39)    ~  23.4393*sinf(twoPI*(jday + 284)/365]
               =  23.4393*sinf(twoPI*(jday-1.39*365/twoPI)/365)          ( 80.75)
               = -23.4393*cosf(twoPI*(jday+(1.39/twoPI - 0.25)*365)/365) (-10.50)
     !/
      decl = eps * sinf((jday-DNSS)* RADPERDAY);

      dr = 1.0f + EOEcc*cos(jday*RADPERDAY); /! inverse relative Earth-Sun distance !/


      /! sunset hour angle cos(thetaz) = cos(lat)⋅cos(decl)⋅cosω+sin(lat)⋅sin(decl) where thetaz=0) !/
      cosegeom = cos(rlatitude) * cos(decl);
      sinegeom = sin(rlatitude) * sin(decl);
      cahou = -sinegeom / cosegeom; /! this is cos(ahou), cosine of the sunset angle hour !/
      sahou = (cahou*cahou >= 1.0)? 0.0f: sqrt(1.0 - cahou*cahou);
      ahou  = atan2(sahou,cahou);

      exsolrad= ScW * dr*SpRad * 2*(ahou  * sinegeom + sahou * cosegeom)/WtL;

     printf("c_shwave day %3d exsolrads %f exsolrad %f diff %f\n",jday, exsolrads, exsolrad, exsolrads/exsolrad-1);
*/

      return((float)exsolrad);
    }
