
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */
/*****************************************************************************
**
**  FILE:      nitrify.c
**
**  FUNCTION:  void nitrify()
**
**  PURPOSE:   Nitrification - Produces N2O gas.
**             The N2O flux during nitrification is controlled by soil
**             temperature, soil water content, and soil NH4 level.
**
**  INPUTS:
**    ammonium - total ammonium in soil mineral pool (gN/m2)
**    maxt     - long term average maximum monthly air temperature of the
**               hottest month (Celsius)
**    nreduce  - reduction factor on nitrification rates due to nitrification
**               inhibitors added with fertilizer (0-1)
**    pHscale  - optional scalar on soil pH
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    layers           - soil water soil layer structure
**    layers->lyblkd[] - bulk density by layer (g/cm3)
**    layers->lyrpH[]  - pH of soil layer
**    layers->wfps[]   - water-filled pore space by layer (fraction 0.0-1.0)
**                      (fraction of a porespace that is filled with water)
**    layers->width[]  - the thickness of soil water model layers (cm)
**
**  LOCAL VARIABLES:
**    cNnh4[]          - ammonium concentration parameters for Parton-Innis function
**                       = { 1.0f, -0.0105f, 0.0f, 0.0f}
**    cNph[]           - ph parameters for Parton-Innis function
**                       = { 5.0f,  0.56f,   1.0f, 0.45f}
**    cNsoilt[]        - temperature effect parameters for Parton-Innis function
**                       = {35.0f, -5.0f,    4.5f, 7.0f}
**    cNwfps[]         - Water Filled Pore Space parameters exponential
**                       = {30.0f, -9.0f}
**    a, b, c, d       - intermediate variables for calculations
**    absoluteMaxRate  - maximum amount of nitrogen that can be nitrified in
**                       one day (gN/m^2)
**    avgstemp         - weighted average of the average soil temperature in
**                       the second and third soil layer used when calculating
**                       the temperature effect on nitrification (degrees C)
**    avgwfps          - average wfps in top 15 cm (0-1)
**    base_flux        - equivalent to 0.1 gN/ha/day
**    base1            - intermediate base variable for calculations
**    base2            - intermediate base variable for calculations
**    debug            - flag to set debugging mode, 0 = off, 1 = on
**    e1, e2           - intermediate exponent variables for calculations
**    fNph             - pH effect on nitrification
**    fNsoilt          - effect of soil temperature on nitrification (0-1)
**    fNwfps           - effect of wfps on nitrification (0-1)
**    MaxRate          - maximum fraction of ammonium that goes to NO3 during
**                       nitrification (0-1)
**    min_ammonium     - ammonium will not be allowed to go below min_ammonium (gN/m^2)
**
**  OUTPUTS:
**    ammonium  - total ammonium in soil mineral pool (gN/m2)
**    nh4_2_no3 - the amount of NH4 that is converted to NO3 due to
**                nitrification
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    f_arctangent()          - arctangent function
**    f_exponential()         - exponential function
**    f_gen_poisson_density() - generalized poisson density function
**    wfps()                  - compute daily water-filled pore space by layer
**
*****************************************************************************/

#include <math.h>
#include <stdlib.h>

#include "soilwater.h"
#include "n2o_model.h"
#include "n2oparam.h"


    void nitrify(double *ammonium, double *nh4_2_no3, float *maxt,
                 float *nreduce, float *pHscale)
    {
      static int debug = 0;
      int ilyr;
      static const double MaxRate = 0.15;
      static const double base_flux = 0.1/10000.0;/* The base_flux is equivalent to 0.1 gN/ha/day */
      float fNsoilt;
      float a;
      float fNwfps;
      float fNph;
      double min_ammonium = 0.03;
      float abiotic;
      float  rel_wc, avg_rel_wc, avgstemp;
      double absoluteMaxRate;

      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;
      extern SOIL_SPT soil;

      *nh4_2_no3 = 0.0;

      if (*ammonium < min_ammonium) {
        if (debug) {
          fprintf(stdout, "CANNOT NITRIFY, ammonium too small %8.6f\n",*ammonium);
        }
        goto RET;
      }

      /*Compute the Ammonium effect on Nitrification
       **Collected and commented these unused calculations KLK Jan 2015
             Oct 29, 2014, at 10:13 AM, Keough,Cynthia <Cindy.Keough@colostate.edu> wrote:
             I found a nitrify.c source code file from 2001 that has an equation that
             has been commented out that uses the fNnh4 variable but I don't have
             any notes on when or why the modification was made.

         float fNnh4      - soil ammonium effect on nitrification
         float grams_soil - soil density, grams/m2, in top 15 centimeters
         float nh4_conc   - soil ammonium (NH4+) concentration (ppm)
         float cNnh4[4] = {1.0f, -0.0105f, 0.0f, 0.0f}; - f_exponential function parameters

        Convert ammonium (g/m2) to nh4_conc (ppm) assuming ammonium is in top 15 cm
         NOTE: Nitrification should occur as a function of depth rather than
               assuming that the top 3 soil layers will be 15 cm
      grams_soil = (layers->lyblkd[0]*layers->width[0] +
                    layers->lyblkd[1]*layers->width[1] +
                    layers->lyblkd[2]*layers->width[2])*100*100;
      nh4_conc = (float)*ammonium/grams_soil*1.0E6f;

          Ammonium effect on Nitrification
      fNnh4 = 1.0f - f_exponential(nh4_conc, cNnh4);

      if (debug > 1) {
           fprintf(stdout, "ammonium = %10.4f,  nh4_conc = %10.4f,  fNnh4 =%6.4f\n",
                   *ammonium, nh4_conc, fNnh4);
      }
    */

      /* Compute the effect of soil water on Nitrification (0-1).
         Use relative water content when the soil is drier than field capacity.
         When soil is wetter use water filled pore space.  cak - 06/16/04 */

      /* Compute relative water content in 2nd and 3rd soil layers, cak - 08/19/04 */
      avg_rel_wc = 0.f;
      for (ilyr = 1; ilyr < 3; ilyr ++) {
        rel_wc = ((float)layers->swc[ilyr]/(layers->width[ilyr]) -
                        layers->swclimit[ilyr]) /
                        (layers->fieldc[ilyr] - layers->swclimit[ilyr]);
        if      (rel_wc < 0.0) {rel_wc = 0.0f;}
        else if (rel_wc > 1.0) {rel_wc = 1.0f;}
        avg_rel_wc += rel_wc * layers->width[ilyr];
      }
      avg_rel_wc /= (layers->width[1] + layers->width[2]);

      if (avg_rel_wc < 1.0) {
        /*equation was fNwfps = (float)(1.0/(1.0 + 30.0 * exp(-9.0 * avg_rel_wc))); */
        /* float fNwfps = 1.0f/(1.0f + cNwfps[0] * expf(cNwfps[1] * avg_rel_wc)); */
        fNwfps = (float)(1.0/(1.0 + (double)cNwfps[0] * exp((double)cNwfps[1] * avg_rel_wc)));
      } else {
        /* average water filled pore space in 2nd and 3rd soil layers, cak - 08/19/04
        Corrected a discontinuity at field capacity stemming from the different units
        used to define the linear decrease. The equation SHOULD BE entirely in WFPS,
        fNwfps = (1.0f - avgwfps) / (1.0f - avgwfpsFC); not a mixture of SW and WFPS
        KLK Nov 2014

        Code was:
        wfps(layers);
        avgwfps = (layers->wfps[1]*layers->width[1] +
                   layers->wfps[2]*layers->width[2]) /
                  (layers->width[1] + layers->width[2]);
        avgfc = (layers->fieldc[1]*layers->width[1] +
                 layers->fieldc[2]*layers->width[2]) /
                (layers->width[1] + layers->width[2]);
        fNwfps = (0.0f - 1.0f) / (1.0f - avgfc) * (avgwfps - 1.0f) + 0.0f;
               = (1.0f - avgwfps) / (1.0f - avgfc);
        Should be:        fNwfps = (1.0f - avgwfps) / (1.0f - avgwfpsFC);

        porsp[i]    = (1.0f - lyblkd[i] / (float)PARTDENS)*width[i] = thetas_bd[i]/0.95
        wfps[i]     = swc[i] / layers->width[i]/(1.0f - lyblkd[i] / PARTDENS) = swc[i]/ porsp[i]
        wfpsFC[i]   = swcfc[i]/ porsp[i] = fieldc[ilyr] * width[i]/ porsp[i]
        wfpsSAT[i]  = swcsat[i]/0.95/ porsp[i] = thetas_bd[ilyr]/0.95 * width[ilyr]/ porsp[i] = 1.0
                      NOTE: the 0.95 fraction is in the thetas_bd definition in initlyrs

        The WFPS fraction by layer:
        fNwfps[i]   = (1.0f - wfps[i]) / (1.0f - wfpsFC[i])
                    = (wfpsSAT[i] - wfps[i]) / (wfpsSAT[i] - wfpsFC[i])
                    = (swcsat[i]/0.95 - swc[i]) / (swcsat[i]/0.95 - swcfc[i])
*/
        fNwfps = (layers->width[1] * ((layers->thetas_bd[1] * layers->width[1] - layers->swc[1]) /
                                      (layers->thetas_bd[1] * layers->width[1] - layers->swcfc[1])) +
                 (layers->width[2] * ((layers->thetas_bd[2] * layers->width[2] - layers->swc[2]) /
                                      (layers->thetas_bd[2] * layers->width[2] - layers->swcfc[2])))) /
                      (layers->width[1] + layers->width[2]);
      }

      /* Compute the soil temperature effect on Nitrification
      A[0] = 35.0f;  A[1] = -5.0f;  A[2] = 4.5f;  A[3] = 7.0f; */
      /* Nitrification rates were too low at low soil temperatures, shift the curve
         so that the nitrification rates are effectively higher for cooler sites,
         this change does not affect sites with hot temperatures, cak - 11/25/03 */
      /* cNsoilt  was  A[4] = {35.0f, -5.0f, 4.5f, 7.0f}; */
      avgstemp = (soil->soiltavg[1] * layers->width[1] + soil->soiltavg[2] * layers->width[2]) /
                 (layers->width[1] + layers->width[2]);
      if (*maxt >= cNsoilt[0]) {
        a = cNsoilt[0];     /* save the input maximum temperature (35)*/
        cNsoilt[0] = *maxt; /* replace first element with maxt*/
        fNsoilt = f_gen_poisson_density(avgstemp,cNsoilt);
        cNsoilt[0] = a;     /* restore the first element (35)*/
      } else {
        fNsoilt = f_gen_poisson_density(avgstemp+(cNsoilt[0]-*maxt),cNsoilt);
      }

      /* Compute pH effect on nitrification
       cNph  was float A[4] = {5.0f, 0.56f, 1.0f, 0.45f}; */
      fNph = f_arctangent(layers->lyrpH[1]*(*pHscale), cNph);

      /* Compute amount of ammonium that goes to nitrate during nitrification */
      if (debug > 1) {
        fprintf(stdout, "cNph    =%6.4f %6.4f %6.4f %6.4f\n",
           cNph[0],cNph[1],cNph[2],cNph[3]);
        fprintf(stdout, "cNsoilt =%6.4f %6.4f %6.4f %6.4f\n",
           cNsoilt[0],cNsoilt[1],cNsoilt[2],cNsoilt[3]);
        fprintf(stdout, "cNwfps  =%6.4f %6.4f\n",
           cNwfps[0],cNwfps[1]);
        fprintf(stdout, "cDno3   =%6.4f %6.4f %6.4f %8.6f\n",
           cDno3[0],cDno3[1],cDno3[2],cDno3[3]);
        fprintf(stdout, "fNwfps  =%6.4f  fNsoilt=%6.4f  Ncoeff=%6.4f  fNph=%6.4f\n",
           fNwfps, fNsoilt, sitepar->Ncoeff, fNph);
        debug = 0;
      }

      abiotic = max(fNwfps * fNsoilt, fabsf(sitepar->Ncoeff));
/*
On Apr 19, 2012, at 9:01 AM, delgro wrote:
Our best guess for now is that netmn_to_no3 should be 0.2 and for the
absoluteMaxRate lets use 1.0 instead of 0.4 in the min equation
      absoluteMaxRate = min(0.4, *ammonium * MaxRate);
      use  MaxNitAmt for the limit
*/
      absoluteMaxRate = min(sitepar->MaxNitAmt, *ammonium * MaxRate);
      *nh4_2_no3 = absoluteMaxRate * fNph * abiotic * *nreduce + base_flux;

      if ((*ammonium - *nh4_2_no3) > min_ammonium) {
        *ammonium -= *nh4_2_no3;
      } else {
        *nh4_2_no3 = min(*nh4_2_no3, *ammonium - min_ammonium);
        *ammonium = min_ammonium;
      }

RET:  return;
    }
