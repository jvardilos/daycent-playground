
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      denitrify.c
**
**  FUNCTION:  void denitrify()
**
**  PURPOSE:   Calculate N2O flux and N2:N2O ratio due to denitrification.
**
**  DESCRIPTION:
**    Denitrification:  Calculate N2 and N2O gases flux, reduce soil nitrate
**    by amount of total nitrogen gas flux.
**
**    Total Nitrogen losses during denitrification are controlled by soil
**    water content, available soil Carbon, and soil nitrate (NO3-).
**
**    The fraction of Nitrogen losses as N2 or N2O are controlled by soil
**    water content and soil nitrate (NO3-).
**
**  INPUTS:
**    basef      - the fraction of N leached into baseflow
**    critflow   - critical water flow for leaching of minerals (cm H2O)
**    frlechd[]  - daily leaching fraction for each minerl, frlechd[0] for N
**    inorglch   - N from organic leaching of stream flow (base flow + storm
**                 flow) (g/m2)
**    jday       - current Julian day (1-366)
**    newCO2     - the amount of microbial respiration that has occurred
**                 in the past day (g/m2)
**    newNO3     - the amount of NO3 to be added to the soil column (g/m2)
**    nitrate[]  - total nitrate in soil mineral pool, distributed among soil
**                 water model layers (gN/m2)
**    stdbulkd   - a standard bulk density for the soil based on
**                 sand/silt/clay percentages (g/cm^3)
**    stdfieldc  - a standard field capacity for the soil based on
**                 sand/silt/clay percentages (volumetric fraction)
**    stream[]   - Century's stream output (stream flow & amount leached into
**                 stream flow)
**    stormf     - the fraction of N leached into stormflow
**                 (baseflow_N + stormflow_N = stream(2))
**    wfluxout[] - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**
**  GLOBAL VARIABLES:
**    MXSWLYR  - maximum number of soil water model layers (21)
**    PARTDENS - particle density (g/cm3)
**
**  EXTERNAL VARIABLES:
**    layers               - soil water soil layer structure
**    layers->lyblkd[]     - bulk density by layer (g/cm3)
**    layers->dpthmn[]     - tops of soil layers (depth from surface in cm)
**    layers->dpthmx[]     - bottoms of soil layers (depth from surface in cm)
**    layers->fieldc[]     - volumetric water content at field capacity for layer
**                           (cm H2O/cm of soil)
**    layers->numlyrs      - total layers in the soil water model profile
**    layers->wfps[]       - water-filled pore space by layer (fraction 0.0-1.0)
**                          (fraction of porespace filled with water)
**    layers->width[]      - the thickness of soil water model layers (cm)
**    sitepar->N2N2Oadj    - multiplier on N2/N2O ratio
**    sitepar->wfpsdnitadj - adjustment on inflection point for water filled
**                           pore space effect on denitrification curve
**    surflood             - surface flooded SW > FC
**
**  LOCAL VARIABLES:
**    cDno3           - nitrate parameters to Parton-Innis functions
**                       = {9.23f, 1.556f, 76.91f, 0.00222f};
**    a               - intermediate variable for calculations
**    CM_per_METER    - 1 meter = 100 cm
**    co2_correction  - correction factor when water filled pore space has
**                      reached a threshold value
**    dD0_fc          - the normalized diffusivity in aggregate soil media, at
**                      a standard field capacity (0-1)
**    debug           - flag to set debugging mode, 0 = off, 1 = on
**    Dtotflux        - total soil (N2 + N2O) denitrification flux
**                      (ppm N/day)
**    efscltef        - effective soil cultivation effect.
**                      This is the C flow weighted clteff for correction of
**                      the CO2 concentration when used as a surrogate for O2
**                      availability. This reduces the CO2 produced by tillage
**                      reflecting enhanced gas permeability from the tillage.
**                      Modification recommended S Del Grosso KLK 14Aug2012
**    excess          - amount of N flux calculated in excess of what can
**                      actually be removed from soil (gN/m2/day)
**    fDco2           - denitrification flux due to CO2 concentation
**                      (ppm N/day)
**    fDno3           - denitrificaton flux due to soil nitrate (ppm N/day)
**    fDwfps          - the effect of wfps on Denitrification (0-1)
**    fluxout         - actual total trace gas flux (gN/m2/day)
**    fRno3_co2       - effect of nitrate and CO2 on the ratio N2:N2O
**    fRwfps          - the effect of wfps on the ratio N2:N2O
**    grams_per_ug    - grams per microgram
**    grams_soil[]    - the grams of soil in a layer in a 1m x 1m section
**    ilyr            - current layer from which N2O flux is occurring
**    k[]             - proportion of newNO3 to distribute to each layer
**    k1              - intermediate variable for calulations
**    ksum            - k summed over all layers
**    M               - intermediate variable for calulations
**    min_nitrate     - mimimum nitrate concentration required in a soil layer
**                      for trace gas calculation (ppm N)
**    min_nitrate_end - minimum allowable nitrate per layer at end of day (ppm N)
**    n2frac          - fraction of nitrate flux that goes to N2 (0.0 - 1.0)
**    n2ofrac         - fraction of nitrate flux that goes to N2O (0.0 - 1.0)
**    nitratePPM[]    - soil nitrate (NO3-) concentration (ppm)
**    ntotflux[]      - soil (N2 + N2O) denitrification flux, by layer
**                      (gN/m2)
**    Rn2n2o          - ratio of N2 to N2O
**    ug_per_gram     - micrograms per gram (ppm)
**    wfps_fc[]       - water filled pore space with layer at field capacity (0-1)
**    WFPS_threshold  - treshold value for water filled pore space (0-1)
**    x_inflection    - inflection point on x-axis
**
**  OUTPUTS:
**    co2PPM[]  - soil carbon dioxide (CO2) concentration by layer (ppm)
**    dN2lyr[]  - N2 flux from Denitrification by layer (gN/m2/day)
**    dN2Olyr[] - N2O flux from Denitrification by layer (gN/m2/day)
**    Dn2flux   - N2 flux from Denitrification (gN/m2/day)
**    Dn2oflux  - N2O flux from Denitrification (gN/m2/day)
**    inorglch  - inorganic N leaching (gN/m2)
**    newNO3    - the amount of NO3 to be added to the soil column (g/m2)
**    nitrate[] - total nitrate in soil mineral pool, distributed among soil
**                water model layers (gN/m2)
**    stream[]  - Century's stream output (stream flow & amount leached into
**                stream flow)
**
**  CALLED BY:
**    trace_gas_model()
**
**  CALLS:
**    diffusiv()       - estimate normalized diffusivity in soils
**    f_arctangent()   - arctangent function
**    leachdly()       - daily nitrate leaching thru soil structure
**    wfps()           - compute the daily water-filled pore space by layer
**
*****************************************************************************/

#include <math.h>

#include "n2o_model.h"
#include "n2oparam.h"
#include "soilwater.h"
#include "thawpulse.h"

#ifndef PI
  #define  PI 3.14159265358979323846264338327950288
#endif
    static int debug = 0;

    void denitrify(float *newCO2, double *newNO3, double nitrate[MXSWLYR],
                   float wfluxout[MXSWLYR], float *critflow, float frlechd[],
                   float stream[], float *basef, float *stormf,
                   double *inorglch, double *Dn2oflux, double *Dn2flux,
                   float stdfieldc, float stdbulkd, float *efscltef,
                   float co2PPM[], double dN2lyr[], double dN2Olyr[], int *jday,
                   int *surflood)
    {
      int    denitrRR; /* respiration restraint on denit */
      int    ilyr;
      float  a;
      float  fDno3;
      float  fDco2;
      float  fDwfps;
      float  Dtotflux;
      float  fRno3_co2;
      float  fRwfps;
      double ntotflux[MXSWLYR];
      float  grams_soil[MXSWLYR];
      double nitratePPM[MXSWLYR];
      float  wfps_fc[MXSWLYR];
      float  co2_correction[MXSWLYR];
      double Rn2n2o;
      double n2ofrac, n2frac;
      double excess;
      double fluxout;
      double k[MXSWLYR], ksum;
      float  k1;
      float  dD0_fc;
      double x_inflection, x_inflbase;
      float  WFPS_threshold;

      /* variables for flooded N2N2O ratio */
      static int lyrfldcnt[MXSWLYR];    /* count of number of days flooded by layer */
      static int initlfc =0;     /* layer count initialized */
      static int floodcnt =0;    /* count of number of days surface flooded (depricated) */
      int thpulsrr = 0;          /* respiration restraint during thaw pulse */
      float N2Osatf; /* saturated fraction */
      float pulssatf = -1; /* saturated fraction during thaw pulse */

      /* delay in transitioning from aerobic to anaerobic n2/n2o ratio
         Per Parton's suggestion set at 1 week KLK 5 June 2015
      */

      /* parameters */
      double min_nitrate = 0.1;
      double min_nitrate_end = 0.05;
      static double ug_per_gram = 1.0E6;
      static double grams_per_ug = 1.0E-6;
      static float  CM_per_METER = 100.0;  /*    float srfc_fc, srfc_swc; */

      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;

      /* Initializations */
      pulssatf = -1;
      for (ilyr=0; ilyr < MXSWLYR; ilyr++) {
        co2PPM[ilyr] = 0.0f;
        dN2lyr[ilyr] = 0.0;
        dN2Olyr[ilyr] = 0.0;
      }
      if(initlfc==0) {
        for (ilyr=0; ilyr < MXSWLYR; ilyr++) {lyrfldcnt[ilyr] = 0.0;}
        initlfc = 1;
      }

      *Dn2oflux = 0.0;
      *Dn2flux = 0.0;

      ksum = 0.0;
      for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
/*        k[ilyr] = exp(K*layers->dpthmn[ilyr]) -  exp(K*layers->dpthmx[ilyr]); */
        /* Use rooting density from soils.in file to weight CO2 by depth, */
        /* cak - 04/07/04 */
        k[ilyr] = max(layers->tcoeff[ilyr], 0.001);
        ksum += k[ilyr];
        ntotflux[ilyr] = 0.0;
      }

      for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
        k[ilyr] /= ksum;
        nitrate[ilyr] += k[ilyr] * (*newNO3);
      }

      /* newNO3 used up after being distributed thru the soil profile */
      *newNO3 = 0;
      wfps(layers);

      /* Mineral leaching */
      leachdly(wfluxout, layers->numlyrs, nitrate, *critflow, frlechd, stream,
               *basef, *stormf, inorglch);

      /* layer independent x_inflection adjustments include an after snow thaw shift KLK 11/04/2017 */
      x_inflbase = (double)sitepar->wfpsdnitadj * thwplse(*jday);

      denitrrpuls(*jday, &thpulsrr, &pulssatf); /* denitrification pulse description */

      /* is today in the window to suppress CO2 effect on denitrification, cak - 02/19/04 */
      denitrRR = (*jday >= sitepar->jdayStart && *jday <= sitepar->jdayEnd);

      /* update flooded day count based on surface flooding */
      /* set surface flooding since dailymoist doesn't calculate it in this version */
      *surflood = (layers->swc[0] + layers->swc[1] + layers->swc[2]) >
          (layers->fieldc[0] * layers->width[0] +
           layers->fieldc[1] * layers->width[1] +
           layers->fieldc[2] * layers->width[2])? 1: 0;
      floodcnt= (*surflood != 0)? floodcnt+1: 0;
 if(debug) {printf("\ndenitrrpuls = %d %d %f    floodcnt %d  flood_N2toN2O %f  floodN2delay %d\n",
 *jday, thpulsrr, pulssatf, floodcnt, sitepar->flood_N2toN2O, sitepar->floodN2delay);}

      /* Denitrification by layer */
      for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
        /* Convert nitrate (gN/m2) to nitratePPM (ppm N) */
        grams_soil[ilyr] = layers->lyblkd[ilyr] * layers->width[ilyr] *
                           CM_per_METER * CM_per_METER;
        nitratePPM[ilyr] = nitrate[ilyr] / grams_soil[ilyr] * ug_per_gram;
        if (nitratePPM[ilyr] < min_nitrate) {
          if (debug>1) {
            fprintf(stdout, "\nCANNOT DENITRIFY layer %1d, nitrate[%1d] = %f\n",
                    ilyr, ilyr, nitrate[ilyr]);
          }
          continue;
        }
        /*
          On Wed, 15 Aug 2012 22:52:13, Keith Paustian wrote:
          we simply ‘ratio’ out the CO2 effect from tillage in the ‘corrected co2’
          variable in the denitrification subroutine, is the best solution. I don’t
          think the simulation of decomposition from tillage events per se causes
          greater N2O from either CO2 effect on O2 or on electron donors. We can
          still get impacts of tillage on N2O from a host of indirect effects.
           Keith
        */
        co2PPM[ilyr] = (float)k[ilyr] * (*newCO2 / *efscltef) /
                       grams_soil[ilyr] * (float)ug_per_gram;
        wfps_fc[ilyr] = layers->fieldc[ilyr] /
                        (1.0f - layers->lyblkd[ilyr]/(float)PARTDENS);
/*        dD0_fc = diffusiv(&stdfieldc, &stdbulkd, &wfps_fc[ilyr]); */
        /* dD0 calc changed 6/20/00 -mdh */
        dD0_fc = diffusiv(&layers->fieldc[ilyr], &layers->lyblkd[ilyr], &wfps_fc[ilyr]);
        WFPS_threshold = (dD0_fc >= 0.15) ? 0.80f : (dD0_fc*250 + 43)/100;
        if (layers->wfps[ilyr] <= WFPS_threshold) {
          co2_correction[ilyr] =  co2PPM[ilyr];
        } else {
          a = (dD0_fc >= 0.15) ? 0.004f : (-0.1f * dD0_fc + 0.019f);
          co2_correction[ilyr] = co2PPM[ilyr] * (1.0f + a *
                                 (layers->wfps[ilyr] - WFPS_threshold)*100);
        }

        /* Compute the Nitrate effect on Denitrification */
        /* Changed NO3 effect on denitrification based on the paper
           "General model for N2O and N2 gas emissions from soils due to
           denitrification", Del Grosso et. al, GBC, 12/00,  -mdh 5/16/00 */

        /* fDno3 = 1.15 * pow(nitratePPM[ilyr], 0.57); */
        /* cDno3[3] = {9.23f, 1.556f, 76.91f, 0.00222f};  was just A */
        fDno3 = max(0.0f, f_arctangent((float)nitratePPM[ilyr], cDno3));

        /* Compute the Carbon Dioxide effect on Denitrification (fDco2, ppm N) */
        /* CO2 effect on denitrification based on Del Grosso GBC paper (above) -mdh 5/16/00 */

        /* The CO2 effect should use the corrected CO2 concentration, cak - 07/31/02 */
        /*      fDco2 = 0.1 * pow(co2PPM[ilyr], 1.3); */
        fDco2 = max(0.0f,
                    (float)((0.1 * pow((double)co2_correction[ilyr], 1.3)) - min_nitrate));

        /* Compute wfps effect on denitrification, (fDwfps, 0-1) */
        /* Changed wfps effect on denitrification based on Del Grosso GBC paper (above)
          -mdh 5/16/00 */

        /*CO2 as a surragate for O2 demand, anaerobic */
        /* The x_inflection calculation should take into account the
           corrected CO2 concentration, cak - 07/31/02 */
        /*  x_inflection = (float)(9.0 - (0.145 - 1.25*min(0.113, dD0_fc))); */
        x_inflection = (9.0 - (0.145 - 1.25*min(0.113, dD0_fc)) * (double)co2_correction[ilyr]) * x_inflbase;

        /* Changed fDwfps calculation - cak - 9/18/00
          fDwfps = 0.5 + (atan(0.6*PI*(0.1*layers->wfps[ilyr]*100-x_inflection))) / PI; */
        fDwfps = fmaxf(0.0f, (float)(0.45 +
                 atan(0.6*PI*(10.0*(double)layers->wfps[ilyr]- x_inflection)) / PI));

        if (debug>1) {
          fprintf(stdout, "\ndenitrify fDwfps= %6.4f   fDno3= %6.4f   fDco2= %6.4f %6.4f %6.4f %6.4f\n",
                 fDwfps,fDno3,fDco2, (0.1 * pow((double)co2_correction[ilyr], 1.3)) - min_nitrate, co2_correction[ilyr], min_nitrate);
        }

        /* Nitrate effect on N2 to N2O ratio. */
          /* if   a flooded ratio is specified
           AND  the volumetric soil water content in the top 3 layers has been greater
                than field capacity for floodN2delay days (nominally a week)
           AND  If we are in a pulse, the saturation is < 100%
        */
        lyrfldcnt[ilyr] = (layers->swc[ilyr]                           - layers->swcfc[ilyr] >
                          (layers->thetas_bd[ilyr]*layers->width[ilyr] - layers->swcfc[ilyr])/128)? lyrfldcnt[ilyr]+1: 0.0;

        /* Compute the N fluxes (N2 + N2O) for the current layer, ppm N */
        Dtotflux = (fDno3 < fDco2 ||
                    (denitrRR == 1 || ilyr<thpulsrr || (thpulsrr ==-1 && lyrfldcnt[ilyr]>0)))? fDno3: fDco2;

        /* Minimum value for potential denitrification in top 2 soil layers */
        /* ppm N, 9/18/00 -cindyk */
        if (ilyr < 2) {Dtotflux = max(0.066f, Dtotflux);}

        /* Account for water limitation */
        if (debug>1) {printf("Dtotflux = %8.4f (%8.4f)  lyrfldcnt %3d  fDno3 %f fDco2 %f fDwfps %f  dRRoff %1d\n",
        Dtotflux*fDwfps, Dtotflux,lyrfldcnt[ilyr],fDno3,fDco2,fDwfps,
        (denitrRR == 1 || ilyr<thpulsrr || (thpulsrr ==-1 && lyrfldcnt[ilyr]>0)));}
        Dtotflux *= fDwfps;

        /* if there is a non-zero flooded n2/n20 ratio, determine if the layer (>0) or
           profile, surface, (<0) has been satureated for the correct time, >floodN2delay
           OTHERWISE the ratio is 0 (unsaturated)
        */
        N2Osatf =((sitepar->flood_N2toN2O < 0  &&  floodcnt        > sitepar->floodN2delay) ||
                  (sitepar->flood_N2toN2O > 0  &&  lyrfldcnt[ilyr] > sitepar->floodN2delay))? 1.0: 0.;
        /* if saturated (N2Osatf==1) and in a thaw pulse (pulssatf != -1) then use the pulse factor */
        if(N2Osatf==1.0  && pulssatf >= 0) {N2Osatf = pulssatf;}

        /* calculate RN2N20 by mixing the unsaturated and saturated ratios*/
        Rn2n2o = (N2Osatf > 0)? fabsf(sitepar->flood_N2toN2O) * N2Osatf: 0.; /* a fraction produces flooded N2 to N20 ratio */
        if(N2Osatf < 1.) {  /* add in the normal calculation if not fully saturated */

          /* Maximum N2/N2O ratio soil respiration function
             the NO3 and CO2 effect on the N2/N2O ratio based on "General model
             for N2O and N2 gas emissions from soils due to denitrification",
             Del Grosso et. al, GBC, 12/00,
              -mdh 5/16/00
             fRno3_co2 estimates the ratio as a function of electron donor to substrate
             -mdh 5/17/00
          */
          k1 = (float)max(1.5, 38.4 - 350 * dD0_fc);
          /* fRno3_co2 = k1 * (float)max(0.16,
                                exp(-0.8 * nitratePPM[ilyr]/(double)co2PPM[ilyr]));
             limit exponent, 1.0e-7 < nitratePPM/co2PPM < 2.29, to prevent underflow  KLK 21/7/2015
          */
#define  MinfRno3_co2 log(0.16)/(-0.8)
          fRno3_co2 = k1 * (((float)nitratePPM[ilyr] > co2PPM[ilyr]*MinfRno3_co2)? 0.16: /* larger exponents give < 0.16 */
                                 expf(-0.8f * (float)nitratePPM[ilyr]/co2PPM[ilyr]));    /* do the exp */

          /* WFPS effect on the N2/N2O Ratio
             wfps effect on the N2/N2O ratio based on  Del Grosso GBC paper (above)
            -mdh 5/16/00 */
          fRwfps = max(0.1f, 0.015f * layers->wfps[ilyr]*100 - 0.32f);

          /* Compute unsaturated N2:N2O Ratio and that to the total */
          Rn2n2o += max(0.1, fRno3_co2 * fRwfps * sitepar->N2N2Oadj)*(1.0-N2Osatf);

          if (debug>1) {
            printf("fRwfps nitratePPM co2PPM fRno3_co2 k1 dD0   %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
                    fRwfps,nitratePPM[ilyr],co2PPM[ilyr],fRno3_co2,k1,dD0_fc);
          }
        }
        if (debug) {
          fprintf(stdout, "%3d  ilyr %2d  Ratio N2:N2O = %10.4f  N2Osatf %5.2f pulssatf %5.2f   floodcount %d\n", *jday, ilyr, Rn2n2o, N2Osatf, pulssatf, (sitepar->flood_N2toN2O < 0)?floodcnt:lyrfldcnt[ilyr]);
        }

        /* Compute N2 and N2O flux by layer (Denitrification) */
        /* convert ppm N to gN/m^2 */
        ntotflux[ilyr] = Dtotflux * grams_soil[ilyr] * grams_per_ug;
        dN2Olyr[ilyr] = ntotflux[ilyr] / (Rn2n2o + 1.0);
        *Dn2oflux += dN2Olyr[ilyr];
        *Dn2flux += ntotflux[ilyr] - dN2Olyr[ilyr];
        dN2lyr[ilyr] = ntotflux[ilyr] - dN2Olyr[ilyr];
        if (debug) {
          printf("lyr  n2o  n2 %3d  %6.4f  %6.4f\n", ilyr, dN2Olyr[ilyr], dN2lyr[ilyr]);
        }
      } /* for ilyr */

      if (debug) {
        printf("In denitrify, Dn2oflux = %12.10f   Dn2flux = %12.10f\n", *Dn2oflux, *Dn2flux);
      }

      if (*Dn2oflux < 1.0E-25) {*Dn2oflux = 0.0;}
      if (*Dn2flux  < 1.0E-25) {*Dn2flux = 0.0;}

      if (debug) {
        printf("In denitrify(1), Dn2oflux = %12.10f   Dn2flux = %12.10f\n", *Dn2oflux, *Dn2flux);
      }

      /* Reduce nitrate in soil by the amount of N2-N N2O-N that is lost */
      /* Do not let nitrate in any layer go below min_nitrate */
      if (*Dn2oflux + *Dn2flux > 1.0E-30) {
        n2ofrac = *Dn2oflux/(*Dn2oflux + *Dn2flux);
        n2frac  = *Dn2flux/(*Dn2oflux + *Dn2flux);
        excess  = 0.0;

        for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
          if (nitratePPM[ilyr] < min_nitrate) {
            /* No trace gas flux from this layer */
            excess += ntotflux[ilyr];
            if (debug) {
              printf("First IF check in loop, excess = %12.10f\n", excess);
            }
          } else if ((nitrate[ilyr] - ntotflux[ilyr]) >
                      (min_nitrate_end * grams_soil[ilyr] * grams_per_ug)) {
            /* Remove N in calculated trace gas flux from the layer */
            nitrate[ilyr] -= ntotflux[ilyr];
            if (debug) {
              printf("Second IF check in loop, nitrate[%1d] = %12.10f", ilyr, nitrate[ilyr]);
              printf("    nitratePPM[%1d] = %12.10f\n", ilyr, nitratePPM[ilyr]);
            }
          } else {
            /* Reduce trace gas flux in layer so soil N won't fall below */
            /* the minimum value */
            fluxout = (nitratePPM[ilyr] - min_nitrate_end) *
                      grams_soil[ilyr] * grams_per_ug;
            excess += (ntotflux[ilyr] - fluxout);
            nitrate[ilyr] = min_nitrate_end * grams_soil[ilyr] * grams_per_ug;
            if (debug) {
              printf("Third IF check in loop, excess = %12.10f\n",excess);
              printf("nitrate[%1d] = %12.10f", ilyr, nitrate[ilyr]);
              printf("     nitratePPM[%1d] = %12.10f\n", ilyr, nitratePPM[ilyr]);
            }
          }
        } /* for ilyr */

        *Dn2oflux -= n2ofrac * excess;
        *Dn2flux  -= n2frac * excess;

        if (debug) {
          printf("  n2ofrac = %12f   n2frac = %12f  excess = %f\n", n2ofrac, n2frac, excess);
          printf("  n2o_excess = %12.10f   n2_excess = %12.10f\n", n2ofrac * excess, n2frac * excess);
        }
      } else {
        *Dn2oflux = 0.0;
        *Dn2flux = 0.0;
      }

      if (debug) {
        fprintf(stdout, "leaving denitrify(2), Dn2oflux = %12.10f   Dn2flux = %12.10f\n", *Dn2oflux, *Dn2flux);
      }

      return;
    }
