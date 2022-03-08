/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "soilwater.h"
#include "mssg.h"

#ifndef min
  #define min(a,b) (((a) < (b)) ? (a) : (b))
  #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

   /* Intialize variables */
    const float dx = 2.0;  /* the width of each soil temperature layer */
    const float maxdepth = 300.f;
    const int nd  = 150;   /* maxdepth/dx; */
    static float temrat = 0;
    /* temperature floating roundoff; set to make the dummy difference stable = 2^-18 ~ 3.8e-6*/
    static int nds = 2;    /* start with 2 heat flow time steps per day */

/*****************************************************************************
**
**  FILE:      therm.c
**
**  FUNCTION:  void therm()
**
**  PURPOSE:   Calculate thermal diffusivity.
**
**  HISTORY:
**    Modified by for use with the current daily soil water model
**    Melannie Hartman
**    9/94
**
**  INPUTS:
**    aclay[]   - fraction of clay in soil layer, 0.0 - 1.0
**    aorg[]    - fraction of organic mater in soil layer, 0.0 - 1.0
**    asand[]   - fraction of sand in soil layer, 0.0 - 1.0
**    bulkden[] - bulk density by layer (g/cm3)
**    depth[]   - the distance from the surface to the middle of the soil
**                layer (cm)
**    dx        - thickness of each soil temperature model layer (cm)
**    fieldc[]  - volumetric water content at field capacity for
**                layer (cm H2O/cm of soil)
**    nd        - the number of soil temperature layers
**    numlyrs   - total number of layers in the soil water model soil profile
**    stmtemp[] - the soil temperature in soil temperature model layer (deg C)
**    swc[]     - soil water content in soil layer (cm H2O)
**    tmax      - maximum air temperature for the day (deg C - 2m)
**    tmin      - minimum air temperature for the day (deg C - 2m)
**    width[]   - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    MXSWLYR     - maximum number of soil water model layers
**    MAXSTLYR    - maximum number of 2 centimeter layers for the soil
**                  temperature model (200)
**
**  LOCAL VARIABLES:
**    bden        - bulk density of current soil temperature layer (g/cm3)
**    clay        - fraction of clay in current soil temperature layer,
**                  0.0 - 1.0
**    d1          - depth of the soil water model layer that the top of the
**                  current soil temperature layer occupies
**    d2          - depth of the soil water model layer that the bottom of the
**                  current soil temperature layer occupies
**    dtlyr       - the center of the current soil temperature layer, the
**                  distance from the surface to the middle of the soil
**                  temperature layer  (cm)
**    fci         - volumetric water content at field capacity for current
**                  soil temperature layer (cm H2O/cm of soil)
**    found       - flag, 0 = current soil temperature layer to soil water
**                  model layer correspondence not found, 1 = current soil
**                  temperature layer to soil water model layer correspondence
**                  found
**    hc_dry_soil - heat capacity of dry soil (cal/gm-deg C)
**    hc_h2o      - heat capacity of water (cal/gm-deg C)
**    hcs         - heat capacity of soil (cal/g-deg C)
**    i1          - index of the soil water model layer that the top of the
**                  current soil temperature layer occupies
**    i2          - index of the soil water model layer that the bottom of the
**                  current soil temperature layer occupies
**    ii, jj      - loop control variables
**    org         - fraction of organic mater in current soil temperature
**                  layer, 0.0 - 1.0
**    p           - proportion of current soil temperature model layer that
**                  occupies the bottom soil water model layer (i2)
**    sand        - fraction of sand in current soil temperature layer,
**                  0.0 - 1.0
**    silt        - fraction of silt in current soil temperature layer,
**                  0.0 - 1.0
**    tcair       - thermal conductivity of air
**    tch2o       - thermal conductivity of water
**    tcs         - thermal conductivity of soil
**    tcsat       - thermal conductivity of air saturated with water vapor
**    tcvap       - thermal conductivity of water vapor
**    tds         - thermal diffusivity of the soil
**    topors      - soil pores not occupied by clay, organic matter, sand,
**                  or silt
**    wh2o        - fraction by weight of moisture in soil (frac)
**    vclay       - fraction of soil volume made up of clay
**    vh2oc       - fraction of soil volume made up of water
**    vmuck       - fraction of soil volume made up of organic matter
**    vsand       - fraction of soil volume made up of sand
**    vsilt       - fraction of soil volume made up of silt
**    xair        - intermediate parameter pertaining to the thermal
**                  conductivity of air in the soil temperature model layer
**    xclay       - intermediate parameter pertaining to the thermal
**                  conductivity of clay in the soil temperature model layer
**    xmuck       - intermediate parameter pertaining to the thermal
**                  conductivity of organic matter in the soil temperature
**                  model layer
**    xsand       - intermediate parameter pertaining to the thermal
**                  conductivity of sand in the soil temperature model layer
**    xsilt       - intermediate parameter pertaining to the thermal
**                  conductivity of silt in the soil temperature model layer
**    xgair       - intermediate parameter pertaining to the thermal
**                  conductivity of air in the soil temperature model layer
**
**  OUTPUTS:
**    tdif[] - thermal diffusivity of the soil by soil temperature model layer
**
**  CALLED BY:
**    soiltemp()
**
**  CALLS:
**    None
**
*****************************************************************************/

    void therm(int numlyrs, float width[MXSWLYR], float depth[MXSWLYR],
               float bulkden[MXSWLYR], float fieldc[MXSWLYR],
               double swc[MXSWLYR], int nd, float stmtemp[MAXSTLYR],
               float tdif[MAXSTLYR], float asand[MXSWLYR], float aclay[MXSWLYR],
               float aorg[MXSWLYR], float dx)
    {
      int   ii, i1, i2, jj;
      int   found;
      float org, sand, silt, clay, topors;
      float xsand, xclay, xmuck, xair, xsilt;
      float vsand, vsilt, vclay, vmuck, vh2oc;
      float dtlyr, d1, d2, bden, p, tch2o, tcsat;
      float tcair, tcvap, fci, xgair, tds, tcs, hcs;
      float wh2o;

      const float  hc_h2o = 1.0f;        /* cal/g/K  or cal/cm3/k */
      const float  hc_dry_soil = 0.20f;  /* cal/g/K */

      jj=0;
      for(ii=0; ii<nd; ii++) {
        /* dtlyr is at the center of soil temperature layer ii */
        dtlyr=(float)(ii)*dx + dx/2;
/*        printf("ii = %1d,  dtlyr = %5.2f", ii, dtlyr);  */
        found = 0;

        /* Determine which soil layer contains the current soil  */
        /* temperature layer.                                    */

        if(dtlyr <= depth[0]) {
          found = 1;

          bden = bulkden[0];
          fci  = fieldc[0];
          vh2oc= (float)swc[0]/width[0];
          sand = asand[0];
          clay = aclay[0];
          org  = aorg[0];

        } else if(dtlyr >= depth[numlyrs-1]) {
          i1=numlyrs-1;
          found = 1;

          bden = bulkden[i1];
          fci  = fieldc[i1];
          vh2oc= (float)swc[i1]/width[i1];
          sand = asand[i1];
          clay = aclay[i1];
          org  = aorg[i1];

        } else {
          /* Search for i1 and i2 */
          while ((depth[jj] > dtlyr) || (dtlyr >= depth[jj+1])) {jj++;}

          i1=jj;
          i2=jj+1;
          d1=depth[i1];
          d2=depth[i2];

/*          printf("\nFOUND: "); */
          found = 1;
          /* if p > 0.5 then dtlyr is farther from the center of layer i1 */
          /* than it is to the center of layer i2 */
          p=(dtlyr-d1)/(d2-d1);

          bden = bulkden[i1] + p*(bulkden[i2] - bulkden[i1]);
          fci  = fieldc[i1] + p*(fieldc[i2] - fieldc[i1]);
          vh2oc= p*(float)swc[i2]/width[i2]+(1.0-p)*(float)swc[i1]/width[i1];
          /* vh2oc= (float)((float)swc[i1]/width[i1] + p*((float)swc[i2]/width[i2] - (float)swc[i1]/width[i1])); */
          sand = asand[i1] + p*(asand[i2] - asand[i1]);
          clay = aclay[i1] + p*(aclay[i2] - aclay[i1]);
          org  = aorg[i1] + p*(aorg[i2] - aorg[i1]);

        } /* end search for i1 i2 */

        if (!found) {
          abortmssg("soil layer indices, i1,i2, not valid in therm");
        }

        wh2o = vh2oc/bden; /* cm/cm * cm3/g  */

        /* New silt terms by Bill Parton 9/94  */
        silt = 1.0f - sand - clay - org;
        vsilt = bden*silt/2.65f;
        vsand = bden*sand/2.65f;
        vclay = bden*clay/2.65f;
        vmuck = bden*org/1.30f;
        topors=1.0f-vsand-vclay-vmuck-vsilt;

        /* calculate heat capacity of soil (cal/g-deg C) */
        /* hcs calculation modified by Bill Parton 9/94 */
/*        hcs=(0.46*(vsand+vclay))+(0.6*vmuck)+vh2oc */
/*        hcs=(0.20*(vsand+vclay+vsilt))+(0.30*vmuck)+vh2oc; */

        hcs = hc_h2o*wh2o + hc_dry_soil*(1-wh2o);

        /* calculate the thermal conductivity of air saturated with water */
        /* vapor and water at temperature of layer */

        tch2o=1.33f+(0.0044f*stmtemp[ii]);
        tcair=0.058f+(0.000017f*stmtemp[ii]);
        tcvap=0.052f*(float)exp(0.058*(double)stmtemp[ii]);
        tcsat=tcair+(tcvap*vh2oc/fci);
        if(vh2oc >= fci) {
          tcsat=tcair+tcvap;
        }

        /* calculate an intermediate parameter pertaining to the */
        /* thermal conductivity of each soil constituent */

        xsand=((2.0f/(1.0f+(((20.4f/tch2o)-1.0f)*0.125f)))+
               (1.0f/(1.0f+(((20.4f/tch2o)-1.0f)*0.75f))))/3.0f;
        xclay=((2.0f/(1.0f+(((7.0f/tch2o)-1.0f)*0.125f)))+
               (1.0f/(1.0f+(((7.0f/tch2o)-1.0f)*0.75f))))/3.0f;
        xsilt=((2.0f/(1.0f+(((7.0f/tch2o)-1.0f)*0.125f)))+
               (1.0f/(1.0f+(((7.0f/tch2o)-1.0f)*0.75f))))/3.0f;
        xmuck=((2.0f/(1.0f+(((0.6f/tch2o)-1.0f)*0.5f)))+1.0f)/3.0f;
        xgair=0.0f+(0.333f*vh2oc/fci);
        if(xgair >= 0.3333) {
          xgair=0.3333f;
        }
        xair=((2.0f/(1.0f+(((tcsat/tch2o)-1.0f)*xgair)))+
              (1.0f/(1.0f+(((tcsat/tch2o)-1.0f)*(1.0f-(2.0f*xgair))))))/3.0f;

        /* calculate the thermal conductivity of soil */

        tcs=((vh2oc*tch2o)+(xsand*vsand*20.4f)+(xclay*vclay*7.0f)+
             (xsilt*vsilt*7.0f)+(xmuck*vmuck*0.6f)+
             (xair*(topors-vh2oc)*tcsat))/(vh2oc+(xsand*vsand)+
             (xclay*vclay)+(xsilt*vsilt)+(xmuck*vmuck)+
             (xair*(topors-vh2oc)))/1000.0f;

        /* calculate the thermal diffusivity of the soil */

        tds=tcs/hcs;

        tdif[ii]=tds;
        if (tds < 0.004) {
          tdif[ii]=0.004f;
        } else if (tds > 0.009) {
          tdif[ii]=0.009f;
        }

/*        tdif[ii] = 0.004f; */

        /* If there is a potential melting of the snow, reduce the thermal */
        /* conductivity.  Bill Parton - 12/5/95 */

      } /* for ii */

      return;
    }


    void tempflow(int jday, int nstep, int nd, float temrat, float srfctemp, float stmtemp[MAXSTLYR], float tdif[MAXSTLYR])
   {
       int pn, kk;
       float dummy1, dtemp;
       const double  temproff = pow(2.0,-18);
       char buffr[132];

       for (pn=0; pn < nstep; pn++) {
         /* Calculate change of temperature for each depth (dtemp)
            soil layer depth is 2cm;  time step = 1 day = 24 hrs*3600 sec = 86400 sec. */
         stmtemp[0] = srfctemp;      /* soil surface temperature was calculated in surftemp */
         dummy1=stmtemp[0]-2.0f*stmtemp[1]+stmtemp[2];
         /* define zero temp below measurement roundoff.
            Prevents underflows by differencing roundoff when temperatures are near zero */
         dtemp = (fabsf(dummy1) < temproff)? 0.: temrat * tdif[0] * dummy1;

/*  printf("  stmtemp(%1d) = %7.2f, tdif = %f, temrat = %f, dummy1 = %e, dtemp = %e\n",0,stmtemp[0],tdif[0],temrat * tdif[0],dummy1,dtemp); */
         for (kk=1; kk<nd; kk++) {
           dummy1= (stmtemp[kk] - stmtemp[kk+1] + (stmtemp[kk+2] - stmtemp[kk+1])) + dtemp;
           stmtemp[kk]=stmtemp[kk]+dtemp;
           if (temrat * tdif[kk] > 1) { /* (fabsf(stmtemp[kk]) > 50) */
             sprintf(buffr, "Unstable heat flow in soiltemp - day = %1d layer %d; temperature change fraction %7.2f >1\n",jday, kk, temrat * tdif[kk]);
             abortmssg(buffr);
             /*printf("Problem in soiltemp - day = %1d, stmtemp(%1d) = %7.2f\n",jday, kk, stmtemp[kk]);*/
           }
           dtemp = (fabsf(dummy1) < temproff)? 0.: temrat * tdif[kk] * dummy1;
/*  printf("  stmtemp(%1d) = %7.2f, tdif = %f, temrat = %f, %f, %f, %f, dummy1 = %e, dtemp = %e\n",
   kk,stmtemp[kk],tdif[kk],temrat * tdif[kk],stmtemp[kk], stmtemp[kk+1], stmtemp[kk+2], dummy1,dtemp); */
         }
       }
   }


    void tempinit(float tmn2m[12], float tmx2m[12], double swc[MXSWLYR], float daylength[366])
    {
      extern SITEPAR_SPT sitepar;
      extern SOIL_SPT soil;

      int   ii;
      float Tavean, Tmaxan, Tminan, Tamp;


      /* Intialize variables */
      /* from The Agricultural Production Systems Simulator, soil temperature model
         http://www.apsim.info/Documentation/Model,CropandSoil/SoilModulesDocumentation/SoilTemp.aspx
         "A solution for the heat flow equation assuming a sinusoidal upper boundary
          condition, which might for example be the annual cycle in air temperature, is
          (Carslaw and Jaeger, 1959),
          t(z,t) = Tavean + Tamp*exp(-z/Zd)*sin(omega*t-z/zd)
           zd = sqrt(2*lambda/(omega*C))
      */

      Tavean = tmn2m[0] + tmx2m[0];
      Tmaxan = tmx2m[0];
      Tminan = tmn2m[0];
      for(ii=1; ii<12; ii++) {
        Tavean += ((tmn2m[ii] < 0.0)? 0.0: tmn2m[ii]) +
                  ((tmx2m[ii] < 0.0)? 0.0: tmx2m[ii]);
        if(tmx2m[ii] > Tmaxan) {Tmaxan = tmx2m[ii];}
        if(tmn2m[ii] < Tminan) {Tminan = tmn2m[ii];}
      }
      Tamp = Tmaxan - Tminan;
      Tavean /= 24.0f;
printf("Tminan= %f  Tmaxan= %f  Tamp= %f\n Setting soil profile to Tavean= %f\n",Tminan,Tmaxan,Tamp,Tavean);

      /* set the soil to average temperature */
      for(ii=0; ii<MAXSTLYR ;ii++) {
        soil->stmtemp[ii] = Tavean;
      }
      /* Lambda is the average soil conductivity
      lambda =*/

      return;
    }




    void tempinity(float tmn2m[12], float tmx2m[12], double swc[MXSWLYR],
                   float daylength[366])
   {
     extern LAYERPAR_SPT layers;
     extern SITEPAR_SPT sitepar;
     extern SOIL_SPT soil;

     float Tavean, Tmaxan, Tminan, Tamp;
     int   ii, mnth;

     float tdif[MAXSTLYR];
     float srfctemp, tmns_mlt, tmxs_mlt;


     /* Intialize variables */
     /* from The Agricultural Production Systems Simulator, soil temperature model
        http://www.apsim.info/Documentation/Model,CropandSoil/SoilModulesDocumentation/SoilTemp.aspx
        "A solution for the heat flow equation assuming a sinusoidal upper boundary
         condition, which might for example be the annual cycle in air temperature, is
         (Carslaw and Jaeger, 1959),
         t(z,t) = Tavean + Tamp*exp(-z/Zd)*sin(omega*t-z/zd)
          zd = sqrt(2*lambda/(omega*C))
     */

     Tavean = tmn2m[0] + tmx2m[0];
     Tmaxan = tmx2m[0];
     Tminan = tmn2m[0];
     for(ii=1; ii<12; ii++) {
       Tavean += tmn2m[ii] + tmx2m[ii];
       if(tmx2m[ii] > Tmaxan) {Tmaxan = tmx2m[ii];}
       if(tmn2m[ii] < Tminan) {Tminan = tmn2m[ii];}
     }
     Tamp = Tmaxan - Tminan;
     Tavean /= 24.0f;

     /* set the soil to average temperature */
     for(ii=0; ii<MAXSTLYR ;ii++) {
       soil->stmtemp[ii] = Tavean;
     }
     /* Lambda is the average soil conductivity
     lambda =*/

     for (mnth=0; mnth < 12; mnth++) {

       if (daylength[mnth*30+15] < 12.0f) {
         tmns_mlt = ((12.0f - daylength[mnth*30+15]) * 3.0f + 12.0f) / 24.0f;
       } else {
         tmns_mlt = ((12.0f - daylength[mnth*30+15]) * 1.2f + 12.0f) / 24.0f;
       }

       tmns_mlt = max(0.05f, min(0.95f, tmns_mlt));
       tmxs_mlt = 0.5f - tmns_mlt + 0.5f;
       srfctemp = tmxs_mlt * tmx2m[mnth] + tmns_mlt * tmn2m[mnth]; /* temperature assuming no snow */

       /* Calculate thermal diffusivity */
       therm(layers->numlyrs, layers->width, layers->depth, layers->lyblkd,
              layers->fieldc, layers->swc, nd, soil->stmtemp, tdif,
           layers->sandfrac, layers->clayfrac, layers->orgfrac, dx);
/*      therm(numlyrs, width, depth, lyblkd, fieldc, swc, nd, stmtemp, tdif,
            sand, clay, org, dx);*/

       temrat=sitepar->dmp * (float)(SEC_PER_DAY/nds)/(dx*dx);
       while ((temrat * tdif[0] > 1.0f || temrat * tdif[nd-1] > 1.0f) && nds < 11) {
         nds++;
         if(nds == 1 || nds == 11) {nds++;}
         temrat=sitepar->dmp * SEC_PER_DAY/(float)nds/(dx*dx);
       }

       tempflow(mnth*30+15, 30*nds, nd, temrat, srfctemp, soil->stmtemp, tdif);

     }

     return;
    }


/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      surftemp.c
**
**  FUNCTION:  void surftemp()
**
**  PURPOSE:   Compute the minimum, maximum, and average soil surface
**             temperatures
**
**  HISTORY:
**    Rewritten from Fortran to C to enable this routine to be called from
**    the watrflow subroutine.  Moved the code that is calculating the
**    insulation effects of snow on surface soil temperature from the soiltemp
**    routine for inclusion in this subroutine.
**    Cindy Keough
**    01/2010
**
**  INPUTS:
**    aglivb        - above ground live biomass (g/m2)
**    daylength     - amount of daylight hours (1..24)
**    diurnal_range - the difference between the maximum and minimum surface
**                    soil temperatures (deg C)
**    elitst        - effect of litter on soil temperature relative
**                    to live and standing dead biomass (fix.100)
**    litrcrbn      - amount of carbon in the surface litter (sum of som1c(1),
**                    som2c(1), strucc(1), and metabc(1)
**    pmntmp        - effect of biomass on minimum surface temperature
**                    (fix.100)
**    pmxbio        - maximum biomass for soil temperature calculations
**                    (fix.100)
**    pmxtmp        - effect of biomass on maximum surface temperature
**                    (fix.100)
**    sfclit        - above ground litter biomass (g/m2)
**    SnowFlag      - effect of snow on soil surface temperature calc
**                    -1 = minimum temperature calculated using snow temp buffering at
**                         snow base during run.
**                         for comparison ONLY; run level retention is a bug.
**                     0 = no snow insulation.
**                     1 = snow insulation around basetemp; may exhibit freeze thaw
**                         basetemp = min(littrC/750.0f - 2.0f, 0.0f)
**                         srfctemp = basetemp + tave * 0.3 * max(1.0 - 0.15*snowpack, 0.0)
**                     2 = snow temp buffering from Parton paper; no
**                         tave > 0 basetemp
**                         tave < 0 basetemp + tave * 0.3 * max(1.0 - 0.15*snowpack, 0.0)
**                     3 = low temperature snow insulation (#2) with basetemp defined
**                         as minimum of basetemp or snowfall temperature;
**                         may exhibit freeze thaw
**                     4 = constant at basetemp
**                     5 = minimum temperature calculated using snow temp buffering at
**                         snow base under current snowpack.
**    snowpack      - current snowpack (equiv. cm H2O)
**    srfctemp      - average soil surface temperature (degrees C)
**    stamt         - amount of warming applied to soil surface temperature
**                    (degrees C)
**    stdead        - above ground dead biomass (g/m2)
**    ststart       - start time for soil warming simulation
**    stsys         - flag to indicate if soil warming will be simulated
**                      0 = no, 1 = yes
**    time          - current simulation time
**    tmax          - maximum air temperature (deg C)
**    tmin          - minimum air temperature (deg C)
**    tmns          - minimum soil surface temperature (degrees C)
**    tmxs          - maximum soil surface temperature (degrees C)
**    woodb         - wood biomass, fine branch + large wood (g/m^2)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    basetemp  - impact of surface litter carbon on mineral soil surface
**                temperature under snow
**    bio       - leaf biomass (g/m^2)
**    sfaltmp   - surface temp under snow, for latching models this is the minimum;
**                large value flags no snow
**    basetemp  - temp at base of snow layer
**    snowmult  - effect of snowdepth on surface temperature.  It is
**                smallest with deep snow and multiples average air
**                temperature.  The more snow, the closer soil surface
**                temperature is to freezing (frac).
**    tmns_leaf - minimum temperature with leaf shading (degrees C)
**    tmns_mlt  - fraction used to compute weighted mimimum surface soil
**                temperature (0.5 - 0.95)
**    tmns_wood - minimum temperature with wood shading (degrees C)
**    tmxs_leaf - maximum temperature with leaf shading (degrees C)
**    tmxs_mlt  - fraction used to compute weighted maximum surface soil
**                temperature (1.0 - tmns_mlt)
**    tmxs_wood - maximum temperature with wood shading (degrees C)
**    woodbio   - wood biomass (g/m^2)
**
**  OUTPUTS:
**    diurnal_range - the difference between the maximum and minimum surface
**                    soil temperatures (deg C)
**    srfctemp      - average soil surface temperature (degrees C)
**    tmns          - minimum soil surface temperature (degrees C)
**    tmxs          - maximum soil surface temperature (degrees C)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

    void surftemp(float elitst, float pmxtmp, float pmntmp, float pmxbio,
                  float tmax, float tmin, float *tmxs, float *tmns,
                  float *srfctemp, float daylength, float aglivb,
                  float stdead, float sfclit, float woodb, float stsys,
                  float ststart, float stamt, float time, float snowpack,
                  float *diurnal_range, float litrcrbn, int SnowFlag)
    {

      float bio, woodbio;
      float tmns_mlt, tmxs_mlt;
      float snowmult;
      float basetemp;
      static float sfaltmp = 256.0f;

      if(sfaltmp == 256.0f) {
        sfaltmp = 128.0f;
        if     (SnowFlag ==-1) {printf("SnowFlag = -1; run level minimum insulated surface temp\n");}
        else if(SnowFlag == 0) {printf("SnowFlag = 0;  no snow insulation\n");}
        else if(SnowFlag == 1) {printf("SnowFlag = 1;  snow insulation, around basetemp\n");}
        else if(SnowFlag == 2) {printf("SnowFlag = 2;  snow insulation with basetemp max; Parton paper\n");}
        else if(SnowFlag == 3) {printf("SnowFlag = 3;  low temp snow insulation, around the minimum of basetemp or snowfall temperature\n");}
        else if(SnowFlag == 4) {printf("SnowFlag = 4;  constant; basetemp\n");}
        else if(SnowFlag == 5) {printf("SnowFlag = 5;  snow cover minimum insulated temp\n");}
        else {abortmssg("unknown SnowFlag in surftemp\n");}
      }

      /* Leaf biomass */
      bio = aglivb + stdead + elitst * sfclit;
      bio = min(bio, pmxbio);

      /* Wood biomass */
      woodbio = min(woodb, 5000.0f);

      /* Maximum temperature with leaf shading
         tmxs_leaf = tmax+(25.4f/(1.0f+18.0f*(float)exp(-0.20*tmax))) *
                     (expf(pmxtmp*bio)-0.13f);
         Maximum temperature with wood shading
         tmxs_wood = tmax+(25.4f/(1.0f+18.0f*(float)exp(-0.20*tmax))) *
                     (expf(pmxtmp*0.1f*woodbio)-0.13f);
         *tmxs = min(tmxs_leaf, tmxs_wood);

         tmns_wood = tmin+pmntmp*0.1f*woodbio-1.78f;  Minimum temperature wood shading
         tmns_leaf = tmin+pmntmp*bio-1.78f;           Minimum temperature leaf shading
         *tmns = max(tmns_leaf, tmns_wood);
      */
      *tmxs = tmax+(25.4f/(1.0f+18.0f*expf(-0.20f*tmax))) *
              (min(expf(pmxtmp*bio), expf(pmxtmp*0.1f*woodbio))-0.13f);
      *tmns = tmin + pmntmp*max(bio, 0.1f*woodbio) - 1.78f;

      /* Average surface temperature */
      /* Use day length to compute the values for the multipliers on */
      /* minimum and maximum soil surface temperature, cak - 04/25/03 */
      /* Weigh surface temperature more towards tmns in the winter between */
      /* the fall and spring equinox when nights are longer than days. */
      if (daylength < 12.0f) {
        tmns_mlt = ((12.0f - daylength) * 3.0f + 12.0f) / 24.0f;
      } else {
        tmns_mlt = ((12.0f - daylength) * 1.2f + 12.0f) / 24.0f;
      }
      tmns_mlt = min(0.95f, tmns_mlt);
      tmns_mlt = max(0.05f, tmns_mlt);
      tmxs_mlt = 0.5f - tmns_mlt + 0.5f;
      *srfctemp = tmxs_mlt * *tmxs + tmns_mlt * *tmns; /* temperature assuming no snow */

      /* Surface temperature without snow */
      if (snowpack <= 0.00000001 || SnowFlag == 0) {
        *diurnal_range = *tmxs - *tmns;
        if(SnowFlag > 1) {sfaltmp = 128.0f;} /* under water, no snow fall temp */
      }

      /* Surface temperature with insulation effects of snow 11/30/95 (Bill Parton). */
      else {

        /* Snow base temperature
           impact of surface litter carbon on mineral soil surface temperature
           On Jan 17, 2013, at 2:56 PM, Keough,Cynthia wrote:
           Modify the constants in the equation that is calculating effect of surface
           litter on snow soil surface temperature.
           Code before:
           basetemp = min((2.0f - (-2.0f)) / (1000.0f - 0.0f) *
                         (litrcrbn - 1000.0f) + 2.0f, 0.0f);
                    = min(4.0f * (litrcrbn/ 1000.0f - 1.0f) + 2.0f, 0.0f);

           Revised equation
           basetemp = min(4.0f * (litrcrbn/3000.0f - 1.0f) + 2.0f, 0.0f);
           NOTE basetemp ~ -2.0C for normal litter litrcrbn KLK 6/2016
        */
        basetemp = min(-2.0f + litrcrbn/750.0f, 0.0f); /* Simplified algebra */

        if(sfaltmp == 128.0f) {sfaltmp = *srfctemp;} /* save temperature at snowfall */

        if(SnowFlag == 3) {
          basetemp = min(basetemp, sfaltmp); /* low temperature model base temperature */
        }

        /* these are the pieces of partons snow temperature fit.
           The If splits these for constant and reactive the sub models, 3 and 4 */
        if ((SnowFlag != 4 && tmin + tmax < 0.0) || SnowFlag == 1 || SnowFlag == 3) {
          /* average air temperature is below freezing */
          snowmult = max(1.0f - 0.15f * snowpack , 0.0f);
          *srfctemp = basetemp + (0.3f * (tmin + tmax) / 2.0f) * snowmult;
          /* *srfctemp = min(*srfctemp, 0.0); */
          *diurnal_range = 0.3f * (tmax - tmin) * snowmult;
          if (*diurnal_range / 2.0 + *srfctemp > 0.0) {
            *diurnal_range = basetemp * *srfctemp;
          }
        }

        else if (SnowFlag == 4 || tmin + tmax >= 0.0) {
          /* if there is snow, and average air temperature gets above
             freezing, average soil surface temperature stays at freezing
             NOTE: applies to SnowFlag -1,1,5;  0,1,3 have been handled so omit them from if */
          *srfctemp = basetemp;
          *diurnal_range = 0.3f * (tmax - tmin) / 2.0f;
          if (*diurnal_range / 2.0 + *srfctemp > 0.0) {
            *diurnal_range = basetemp * *srfctemp;
          }
        }

        /* code to hole srfctemp less than a season or run minimum temperature
           ratchets down the minimum temperature seen under the snow layer
           The original coding (SnowFlag == -1) didn't reset the minimum and it continued
           to decrease over THE ENTIRE RUN!
           Obviously A BUG, snow insulation doesn't carry past snowmelt. KLK 6/2016
        */
        if(SnowFlag == -1  || SnowFlag == 5) {
          /* If the soil has frozen prior to snow accumulating on the site
             the insulating effect of the snow will hold the low frozen soil
             temperature until the snow melts, cak 09/29/2011 */
          if (sfaltmp < *srfctemp) {*srfctemp = sfaltmp;}
          sfaltmp = *srfctemp; /* ratchet the previous down */
        }

        /* printf("Snowpack SnowFlag= %d,  tave= %f,  srfctemp= %f, savmintmp= %f, basetemp= %f, snowpack= %f\n",
                SnowFlag, (tmin + tmax)/2., *srfctemp, sfaltmp, basetemp, snowpack); */
      }

      /* Added code to allow warming of the soil surface temperature, */
      /* cak - 07/02/03 */
      if (stsys > 0 && time >= ststart) {
        *srfctemp = *srfctemp + stamt;
      }

      return;
    }


/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      soiltemp.c
**
**  FUNCTION:  void soiltemp()
**
**  PURPOSE:   This subroutine calculates the daily average, maximum, and
**             minimum soil temperature (t[ii][j]) for a specified number of
**             depths (nd).  The inputs include the maximum and minimum air
**             temperature (tmax and tmin), the standing plant biomass
**             (biomass), the soil temperature at the bottom of the soil
**             profile, and the thermal diffusivity of the soil. The model is
**             described in a paper by Parton(1981).
**
**  HISTORY:
**    Rewrite of base temperatue simulation KLK 8/8/2016
**     - run temperature flow to 3 meters
**     - assume bottom boundary condition is yearly average air temperature at 3m
**     - disable the measured sine response for the Pawnee profile.
**
**    Decrease temperature flow time step to 12 hours or less to prevent instabilities.
**    For certain soils, setting the dmp factor large enough to model to 20 cm causes
**    heat flow to diverge, predicting temperature changes larger than the difference,
**    and causing the code to crash. Decreasing the time step allows temperatures to
**    respond to the energy flow and damps this instability. KLK 16 Aug 2016
**
**    Modified for use with the Trace Gas Model    Melannie Hartman    4/97
**
**  INPUTS:
**    lyblkd[]      - bulk density by layer (g/cm3)
**    clay[]        - the fraction of clay in soil layer
**    depth[]       - the distance from the surface to the middle of the soil
**                    layer (cm)
**    diurnal_range - the difference between the maximum and minimum surface
**                    soil temperatures (deg C)
**    fieldc[]      - volumetric water content at field capacity for
**                    layer (cm H2O/cm of soil)
**    jday          - current julian day (1..366)
**    numlyrs       - total number of layers in the soil water model soil
**                    profile
**    org[]         - the fraction of organic matter in soil layer
**    sand[]        - the fraction of sand in soil layer
**    soiltavewk    - average soil temperature in the second soils.in soil
**                    layer over the previous 7 days (deg C)
**    soiltavg[]    - average soil temperature by layer (deg C)
**    soiltmax[]    - maximum soil temperature by layer (deg C)
**    soiltmin[]    - minimum soil temperature by layer (deg C)
**    srfctemp      - soil surface temperature (deg C)
**    stmtemp[]     - the average soil temperature of the soil temperature
**                    model layers
**    swc[]         - soil water content of the soil layer (cm H2O)
**    tmax          - maximum air temperature for the day (deg C - 2m)
**    tmin          - minimum air temperature for the day (deg C - 2m)
**    tmns          - minimum soil surface temperature (deg C)
**    tmxs          - maximum soil surface temperature (deg C)
**    width[]       - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    MXSWLYR     - maximum number of soil water model layers
**    MAXSTLYR    - maximum number of 2 centimeter layers for the soil
**                  temperature model (200)
**    PI          - pi (3.14159265)
**    SEC_PER_DAY - number of seconds in a day (86400)
**
**  EXTERNAL VARIABLES:
**    sitepar->dmp    - time step correction factor, relates to how fast the
**                      heat gets into/out of the soil
**
**  LOCAL VARIABLES:
**    a, b, c, d     - intermediate variables for calculations
**    arrayindx      - position in the circular avgtemparray array
**    asilt[]        - the fraction of silt in soil layer (0.0 - 1.0)
**    avgtemparray[] - array used to store the average soil temperatures in
**                     the 2nd soil layer over the previous 7 days
**    avtd           - the average thermal diffusivity used to calculate
**                     maximum and minimum soil temperature as a function of
**                     depth
**    dbot           - depth at bottom of the soil profile depth(numlyrs) + 2 cm (cm)
**    deltat         - the change in soil temperature between today and
**                     yesterday (deg C)
**    diff1          - intermediate variable for calculations
**    differ         - the difference between the maximum and minimum surface
**                     soil temperatures (deg C)
**    dtemp          - daily change of temperature (deg c/day) in the previous soil layer
**    dummy1         - intermediate variables for calculations
**    dx             - the depth interval for the soil temperature model soil
**                     layers (cm)
**    ii, kk, ll     - loop control variables
**    k1, k2         - intermediate variables for calculations
**    maxdepth       - the depth of the bottom layer in the soil profile (cm)
**    m1, m2         - intermediate variables for calculations
**    nd             - number of soil temperature layers (calculated by model)
**    prevstavg[]    - soiltavg[] from the previous day (deg C)
**    soilenrgy      - total energy absorbed/released from the soil (cal)
**                     Negative soilenergy represents heat going into the
**                     soil.  Positive soilenergy represents heat going to
**                     warm the atmosphere.
**    startofrun     - flag to indicate the start of the run
**    t[][]          - the average ([][0]),maximum ([][1]),and minimum ([][2])
**                     soil temperature for the soil layer (deg C)
**    tdif[]         - thermal diffusivity of the soil by soil temperature
**                     model layer
**    temrat         - intermediate variable for calculations
**    vclay          - fraction of soil volume made up of clay
**    vh2oc          - fraction of soil volume made up of water
**    vmuck          - fraction of soil volume made up of organic matter
**    volheat[]      - volume heat of soil (cal/cm3/deg)*(cm3*deg)=cal
**    vsand          - fraction of soil volume made up of sand
**    vsilt          - fraction of soil volume made up of silt
**
**  OUTPUTS:
**   soiltavewk - average soil temperature in the second soils.in soil layer
**                over the previous 7 days (deg C)
**   soiltavg[] - average soil temperature by layer (deg C)
**   soiltmax[] - maximum soil temperature by layer (deg C)
**   soiltmin[] - minimum soil temperature by layer (deg C)
**   stmtemp[]  - the average soil temperature of the soil temperature model
**                layers
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    therm() - calculate thermal diffusivity
**
*****************************************************************************/

    void soiltemp(int jday, float tmin, float tmax, float depth[MXSWLYR],
                  float width[MXSWLYR], float fieldc[MXSWLYR],
                  float sand[MXSWLYR], float clay[MXSWLYR], float org[MXSWLYR],
                  float lyblkd[MXSWLYR], double swc[MXSWLYR], int numlyrs,
                  float soiltavg[MXSWLYR], float soiltmin[MXSWLYR],
                  float soiltmax[MXSWLYR], float stmtemp[MAXSTLYR], float tmns,
                  float tmxs, float *soiltavewk, float srfctemp,
                  float diurnal_range)
    {
      extern SITEPAR_SPT sitepar;

      float tdif[MAXSTLYR], t[MXSWLYR][3];
      int   ii;
      float differ, diff1, deltat;
      /* float a, b, c, d; */
      float avtd;
      float vsand, vsilt, vclay, vh2oc, vmuck;
      float asilt[MXSWLYR];
      float volheat[MXSWLYR];
      float prevstavg[MXSWLYR];
      int   k1, k2;
      float m1, m2;
      float soilenrgy;

      static int   arrayindx;
      static float avgtemparray[7];
      static int   startofrun = 1;

      /* Intialize variables */

      soilenrgy = 0.0f;

      for (ii=0; ii<MXSWLYR; ii++) {
        prevstavg[ii] = soiltavg[ii];
      }

      /* Calculate the number of soil layers nd
         NOTE: ndd not used
               strange depth calculation the 5.0f probably should be dx
               then the rest of the calculation is the hard way to get the wrong number.
      maxdepth=0.0f;
      for(ii=0; ii<numlyrs; ii++) {
        maxdepth=maxdepth+width[ii];
      }
      dbot=maxdepth+5.0f;
      nd=(int)(dbot/dx -1.0f);
      ndd=nd+1;
 */

      /* Calculate thermal diffusivity */
      therm(numlyrs, width, depth, lyblkd, fieldc, swc, nd, stmtemp, tdif,
            sand, clay, org, dx);
      temrat=sitepar->dmp * (float)(SEC_PER_DAY/nds)/(dx*dx);
      while ((temrat * tdif[0] > 1.0f || temrat * tdif[nd-1] > 1.0f) && nds < 11) {
        nds++;
        if(nds == 1 || nds == 11) {nds++;}
printf("adjusting temperature flow time step %d  temp ratio surface %f 3m %f\n",nds, temrat * tdif[0], temrat * tdif[nd-1]);
        temrat=sitepar->dmp * SEC_PER_DAY/(float)nds/(dx*dx);
      }

       tempflow(jday, nds, nd, temrat, srfctemp, stmtemp, tdif);

      /* Specify the soil temperature at depth(numlyrs) + 2cm
         set to sin function, using max, min and lag time of the minimum at profile depth,
         parameters based on Pawnee 1971-79 ave monthly soil temps at 183 cm.
         removed sine curve temperature boundary condition in favor of a constant T at 3m.
         KLK 8/8/2016
      */

      /* Calculate the average,maximum and minimum soil temperature(t[ii,j])
         for the iith soil depth.  The depth is defined by depth[ii].  The
         average soil temperature is calculated by linearly interpolating
         between soil temperatures calculated at 2cm intervals.  Based on
         the fourier heat transfer equation. */
      avtd=0.00215f;
      for(ii=0; ii<numlyrs; ii++) {
        if(depth[ii] == 0.) {
          t[ii][1]=tmxs;
          t[ii][2]=tmns;
          t[ii][0]=(tmxs+tmns)/2.0f;
        } else {
          k1 = (int)max((depth[ii]-dx/2.0f)/dx, 0.0f);
          k2 = (int)((depth[ii]+dx/2.0f)/dx);
          m1 = (dx*k2-(depth[ii]-dx/2.0f))/dx;
          m2 = (depth[ii]+dx/2.0f - dx*k2)/dx;
          t[ii][0] = (m1*stmtemp[k1] + m2*stmtemp[k2]);

          /* Calculate the maximum and minimum soil temperature at depth
             depth[ii].use eqn presented by parton(1983),which is a function
             of the average thermal diffusivity in the top 15cm, and the
             soil depth(depth[ii]),and the diurnal variation at the soil
             surface. */
          differ = diurnal_range;
          diff1=-depth[ii]* (float)pow((0.00005/(double)avtd),0.5);
          if (diff1 < -60.0) {
            diff1 = -60.0f;
          }
          t[ii][1]=t[ii][0]+differ*(float)exp((double)diff1)/2.0f;
          t[ii][2]=t[ii][0]-differ*(float)exp((double)diff1)/2.0f;
        }
      }

      /* Save today's values of soil temperature for use tomorrow */
      for(ii=0; ii<numlyrs; ii++) {

        /* New volume heat calculation by Bill Parton 9/94  */
        /* The soil volume is 100cm*100cm*width[ii]      */
        asilt[ii] = 1.0f-sand[ii]-clay[ii]-org[ii];
        vsilt = lyblkd[ii]*asilt[ii]/2.65f;
        vsand=lyblkd[ii]*sand[ii]/2.65f;
        vclay=lyblkd[ii]*clay[ii]/2.65f;
        vmuck=lyblkd[ii]*org[ii]/1.30f;
        vh2oc = (float)swc[ii]/width[ii];

        /* If the temperature drops in the soil, energy is given of back
           to the atmosphere, and is therefore positive */
        deltat = prevstavg[ii] - t[ii][0] ;

        volheat[ii] = (0.20f*(vsand+vclay+vsilt)*2.65f + 0.30f*vmuck*1.30f +
                       vh2oc) * deltat*(width[ii]*100*100);

        soiltavg[ii] = t[ii][0];
        soiltmax[ii] = t[ii][1];
        soiltmin[ii] = t[ii][2];

        /* Negative soilenergy represents heat going into the soil
           Positive soilenergy represents heat going to warm the atmosphere */
        soilenrgy = soilenrgy + volheat[ii];
      }

      /* If it is the start of the run initialize all of the values in the */
      /* avgtemparray array using the soil temperature value calculated on */
      /* the first day of the run, cak - 04/27/04 */
      if (startofrun) {
        for (arrayindx = 0; arrayindx < 7; arrayindx++) {
          avgtemparray[arrayindx] = soiltavg[1];
        }
        arrayindx = arrayindx % 7;
        startofrun = 0;
      }

      /* Compute 7 day average soil temperature for the 2nd layer. Store the
         last 7, 2nd layer temperatures in a circular array, cak - 04/27/04 */
      avgtemparray[arrayindx] = soiltavg[1];
      arrayindx = arrayindx + 1;
      if (arrayindx > 6) {
          arrayindx = 0;
      }

      /* Code added to compute average soil temperature for the 2nd layer
         over the previous 7 days, this value is used in the scheduling of
         the start of plant growth when using the growing degree day
         implementation, cak - 04/27/04 */
      *soiltavewk = 0.0f;
      for (ii = 0; ii < 7; ii++) {
        *soiltavewk = *soiltavewk + avgtemparray[ii];
      }
      *soiltavewk = *soiltavewk / 7.0f;

      return;
    }
