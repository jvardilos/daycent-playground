
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      methane_production.c
**
**  FUNCTION:  void methane_production()
**
**  PURPOSE:   Calculate methane production.
**
**  DESCRIPTION:
**    Methane production is a function of soil texture, soil temperature,
**    litter and SOM decomposition rates, aboveground live carbon, and soil
**    redox potential.
**
**  REFERENCE:
**    Predicting methanogenesis from rice paddies using the DAYCENT ecosystem model.
**    Kun Cheng, Stephen M. Ogle, William J. Parton, Genxing Pan
**    Ecological Modelling 261–262 (2013) 19–31
**
**    A semi-empirical model of methane emission from flooded rice paddy
**    soils.  Y. Huang, R.L. Sass, and F.M. Fisher, Jr.  1998.  Global Change
**    Biology 4: 247-268.
**
**    Modeling methane emission from rice paddies with various agricultrual
**    pratices.  Y. Huang, W. Zhang, X. Zheng, J. Li, and Y. Yu.  Journal of
**    Geophysical Research.  2004.  19: DO8113
**
**    Modeling methane emissions from rice paddies.  M. Cao, J.B. Dent, and
**    O.W. Heal.  1995.  Global Biogeochemical Cycles.  Vol 9, No. 2, pages
**    183-195.
**
**  INPUTS:
**    Com        - sum of CO2 losses from heterotrophic decomposition of
**                 metabc(1), metabc(2), strucc(1), strucc(2), som1c(1),
**                 som1c(2), som2c(1), som2c(2), and som3c (g/m^2)
**    Eh         - effect of water management on soil redox potential (Eh)
**                 (-250.0 mv - 300.0 mv)
**    prev_bgprd - previous day's fine root carbon production (g/m^2)
**    watertable - 1 = simulate water table (flooding stage),
**                 0 = no water table (drainage stage)
**    watrflag   - 1 = water added to the system automatically to bring the
**                     soil water content to saturation
**                 0 = rain and irrigation events used to add water to the
**                     soil potentially bringing it to saturation
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    layers                  - soil water soil layer structure
**    layers->sandfrac[]      - sand fraction in soil layer, 0.0 - 1.0
**    layers->width[]         - the thickness of soil water model layers (cm)
**    sitepar                 - site specific parameters structure for soil water model
**    sitepar->Aeh            - differential coefficient for calculation of soil Eh
**    sitepar->Beh_drain      - upper-limit value of Eh during drainage course (mv)
**    sitepar->Beh_flood      - low-limit value for Eh during flooding course (mv)
**    sitepar->Deh            - differential coefficient for calculation of soil Eh
**    soil                    - soil temperature structure
**    soil->soiltavg[]        - average soil temperature of layer (degrees C)
**
**  LOCAL VARIABLES:
**    aglivc_biomass - aboveground live carbon converted from grams of carbon
**                     to grams of biomass
**    avgsand        - average fraction sand content in top 3 soil layers
**    Q10            - temperature coefficient for methane emission
**
**  OUTPUTS:
**    avgst_10cm - average soil temperature in top 10 cm of soil profile
**                 (degrees C)
**    CH4_prod   - methane production (gC/m2/day)
**    Cr         - carbohydrates derived from rice plants (g/m^2)
**    Eh         - effect of water management on soil redox potential (Eh)
**                 (-250.0 mv - 300.0 mv)
**    Feh        - reduction factor of effect of soil redox potential (Eh) on
**                 methane production (0.0 - 1.0)
**    SI         - soil texture index for methane production (0.0 - 1.0)
**    TI         - soil temperature index for methane production (0.0 - 1.0)
**
**  CALLED BY:
**    trace_gas_model()
**
**  Changes:
**    bug fix KLK 10 March 2015
**    - corrected the value C6H12O6_to_CH4. The Cheng paper uses the same value for
**      equation 1 and 2 but were coded as multiplicative inverses.
**    - moved the Com weighting from tcmodel to this routine (ref 1 eq 1)
**    - C6H12O6_to_CH4 was is a local parameter. Cheng (ref 1) calculates the value
**      of 0.5 from reaction dynamics and not a free parameter.
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"
#include "n2o_model.h"
#include <math.h>
#include <stdio.h>
static float SI = -1.0f;

    void methane_production(float *prev_bgprd, float *Com, float *avgst_10cm,
                            float *TI, float *Cr, float *Eh,
                            float *Feh, float *CH4_prod, int *watertable,
                            int *watrflag)
    {
      float Q10 = 3.0f;

      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;
      extern SOIL_SPT soil;

      #define C6H12O6_to_CH4 0.5   /* mole weight CH4 produced from C6H12O6 consumed in
                                      anaerobic carbohydrate fermentation  CO2 */

      /* Soil texture index for methane production
      From: Keough,Cynthia
      Sent: Tuesday, June 17, 2014 11:36 AM
      Subject: RE: soil texture index
       At one point the SI value was being used in the calculation of Cr
       in the methane_production but a model modification that was implemented
       on 3/22/2012 modified the equation so that the SI value is no longer
       being used.  The VI variable was also present in the old equation.
       -- Removed SI from call arguments and output KLK &/2014

       -- Average sand content in top 3 soils.in soil layers (10 cm)
      */

      if(SI < 0) {SI = 0.325f + 2.25f * sitepar->avgsand;} /* use fraction instead of percent sand 0.0255 -> 2.25 */

      *Com = C6H12O6_to_CH4 * sitepar->CO2_to_CH4 * SI * (*Com);

      /* Average soil temperature in top 3 soils.in soil layers (10 cm) */
      *avgst_10cm = (soil->soiltavg[0] * layers->width[0] +
                     soil->soiltavg[1] * layers->width[1] +
                     soil->soiltavg[2] * layers->width[2]) /
                    (layers->width[0] + layers->width[1] + layers->width[2]);
      /* Soil temperature index for methane production */
      *TI = min(*avgst_10cm, 30);
      *TI = (float)pow(Q10, ((*TI - 30.0) / 10.0));

      /* Use the previous day's fine root production to calculate */
      /* carbohydrates derived from rice plants */
      *Cr = sitepar->frac_to_exudates * SI * (*prev_bgprd);

      /* Soil redox potential */
      if (*watertable == 1) {
        /* Flooding stage */
        if (*watrflag == 1) {
          /* Water added automatically to keep soil at saturation */
          *Eh = *Eh - sitepar->Deh * (sitepar->Aeh + min(1.0f, *Com)) *
                (*Eh - sitepar->Beh_flood);
        } else {
          /* Water added via rain and irrigation events */
          *Eh = -20.0f;
        }
      } else {
        /* Drainage stage */
        *Eh = *Eh - sitepar->Deh * (sitepar->Aeh + 0.7f) *
              (*Eh - sitepar->Beh_drain);
      }

      /* Reduction factor of effect of soil redox potential on methane */
      /* production */
      if (*Eh < -150.0) {
        /* *Feh = (float)exp(-1.7 * ((150.0 + (-150.0)) / 150.0)); */
        *Feh = 1.0f;
      } else {
        *Feh = (float)exp(-1.7 * ((150.0 + *Eh) / 150.0));
      }

      /* Methane production */
      *CH4_prod = C6H12O6_to_CH4 * *Feh * (*Com + (*TI * (*Cr)));

      return;
    }


/*****************************************************************************
**
**  FILE:      methane_oxidation.c
**
**  FUNCTION:  void methane_oxidation()
**
**  PURPOSE:   Calculate methane oxidation.
**
**  DESCRIPTION:
**    Methane oxidation is a function of soil water content, temperature,
**    porosity, and field capacity.
**
**  REFERENCE:
**    General CH4 Oxidation Model and Comparison of CH4 Oxidation in Natural
**    and Managed Systems.  S.J. Del Grosso, W.J. Parton, A.R. Mosier, D.S.
**    Ojima, C.S. Potter, W. Borken, R. Brumme, K. Butterbach-Bahi, P.M.
**    Crill, K. Dobbie, and K.A. Smith.  2000.  Global Biogeochemical Cycles
**    14:999-1019.
**
**  INPUTS:
**    isdecid - flag, set to 1 for deciduous forest, 0 otherwise
**    isagri  - flag, set to 1 for agricultural system or grassland that has
**              been fertilized or cultivated, 0 otherwise
**
**  GLOBAL VARIABLES:
**    PARTDENS - particle density (g/cm3) (2.65)
**
**  EXTERNAL VARIABLES:
**    layers           - soil water soil layer structure
**    layers->lyblkd[] - bulk density by layer (g/cm3)
**    layers->dpthmn[] - tops of soil layers (depth from surface in cm)
**    layers->dpthmx[] - bottoms of soil layers (depth from surface in cm)
**    layers->fieldc[] - volumetric water content at field capacity for layer
**                       (cm H2O/cm of soil)
**    layers->numlyrs  - total number of layers in the soil water model
**                       soil profile
**    layers->wfps[]   - water-filled pore space by layer (fraction 0.0-1.0)
**                       (fraction of a porespace that is filled with water)
**    layers->width[]  - the thickness of soil water model layers (cm)
**    soil             - soil temperature structure
**    soil->soiltavg[] - average soil temperature of layer (degrees C)
**
**  LOCAL VARIABLES:
**    agri_adjust   - adjustment factor for agricultural soils
**    bulkdensity   - weighted average of bulk density in top 15 cm of soil
**                    (g/cm3)
**    CH4DEPTH      - maximum depth at which CH4 oxidation occurs
**    CH4max        - maximum CH4 oxidation rate (gC/ha/day)
**    Dopt          - ratio of gas diffusivity through soil at Wopt to gas
**                    diffusivity through air
**    fieldcapacity - weighted average field capacity of top 15 cm of soil
**                    (cm H20 / cm of soil)
**    ilyr          - current layer from which methane oxidation is occurring
**    percentlayer  - fraction of the 15 cm that a layer adds to the CH4 15cm
**    soiltemp      - weighted average of soil temperature in top 15 cm of
**                    soil (deg C)
**    soilwater     - soil water content of top 15 cm of soil
**                    (volumetric, cm3 H2O / cm3 soil * 100)
**    temp_adjust   - adjustment factor for soil temperature
**    watr_adjust   - adjustment factor for water limitation
**    wfps          - weighted average of water filled pore space in top 15 cm
**                    of soil (fraction 0.0-1.0)
**    wfps_adjust   - adjustment factor for water filled pore space
**    Wmax          - maximum water content
**    Wmin          - minimum water content
**    Wopt          - optimum water content
**
**  OUTPUTS:
**    CH4_oxid - methane oxidation (gC/ha/day)
**
**  CALLED BY:
**    trace_gas_model()
**
**  CALLS:
**    diffusiv() - estimate normalized diffusivity in soils
**
*****************************************************************************/

#define CH4DEPTH 15.0f

    void methane_oxidation(double *CH4_oxid, int *isdecid, int *isagri)
    {
      int    ilyr;
      float  bulkdensity;
      float  fieldcapacity;
      double soiltemp;
      double soilwater;
      double wfps;
      double CH4max;
      float  Dopt;
      double Wmin;
      double Wmax;
      double  Wopt;
      double agri_adjust;
      double temp_adjust;
      double watr_adjust;
      double wfps_adjust;
      float  percentlayer;
      float  temp;

      extern LAYERPAR_SPT layers;
      extern SOIL_SPT soil;

      /* Compute a weighted average for soil temperature, field capacity, */
      /* bulk density, water filled pore space, and volumetric soil water */
      /* content in top 15 cm of soil profile */
      bulkdensity = 0.0f;
      fieldcapacity = 0.0f;
      soiltemp = 0.0;
      soilwater = 0.0;
      wfps = 0.0;
      for (ilyr = 0; (ilyr < layers->numlyrs && layers->dpthmn[ilyr] < CH4DEPTH); ilyr++) {
        /* the width is either the width or:
          layers->width[ilyr] * percentlayer =
          layers->width[ilyr] * (CH4DEPTH - layers->dpthmn[ilyr]) / layers->width[ilyr] */
        percentlayer = (layers->dpthmx[ilyr] <= CH4DEPTH)? layers->width[ilyr] / CH4DEPTH:
                                             (CH4DEPTH - layers->dpthmn[ilyr]) / CH4DEPTH;
        bulkdensity += layers->lyblkd[ilyr] * percentlayer;
        fieldcapacity += layers->fieldc[ilyr] * percentlayer;
        soiltemp += soil->soiltavg[ilyr] * percentlayer;
        soilwater += layers->wfps[ilyr] * percentlayer;
        wfps += layers->wfps[ilyr] * percentlayer;
      }
      /* Convert from water filled pore space to volumetric water
         multiply soilwater by 10 instead of 100 then again by 0.1
         each time it is used. Saves 5 multiplications, none of which
         are exact in binary   KLK 7Mar12 */
      soilwater = 10.0 * soilwater * (1.0 - (bulkdensity / PARTDENS));
/*      soilwater *= 100.0;*/

      /* CH4 oxidation for a deciduous system */
      if (*isdecid) {
        CH4max = 40.0 - 18.3 * bulkdensity;
        temp_adjust = 0.0209 * soiltemp + 0.845;
        /* Use bounded value for wfps_adjust if wfps falls below a critical */
        /* value, cak - 11/12/02 */
        if (wfps <= 0.05) {
          wfps_adjust = 0.1;
        } else {
          wfps_adjust = pow((10.0 * wfps - 0.5) / (1.84 - 0.5), 0.13);
          wfps_adjust *= pow((10.0 * wfps - 55) / (1.84 - 55),
                             (0.13 * (55 - 1.84)) / (1.84 - 0.5));
          wfps_adjust = max(0.1, wfps_adjust);
        }
        *CH4_oxid = CH4max * wfps_adjust * temp_adjust;

      } else {
        /* CH4 oxidation for a grassland/coniferous/tropical system */
        Wmin = 3.0 * fieldcapacity - 0.28;
        Wopt = 6.3f * fieldcapacity - 0.58f;
        Wmax = 10.6 * fieldcapacity + 1.9;
        temp = Wopt * 0.1f / (1.0f - (bulkdensity / (float)PARTDENS));
        Dopt = diffusiv(&fieldcapacity, &bulkdensity, &temp);
        CH4max = 53.8 * Dopt + 0.58;
        /* Use lower bound when Wopt <= Wmin, fieldcapacity <= 0.0909091.
           The equation blows up for Wopt <= Wmin. The base of the first power function
           goes negative generating NaN, (negative numbers to fractional powers is
           undefined over reals [Note: the equivalent ((-a)^2)^0.2 is defined but still
           strongly divergent]) The second term is defined but strongly divergent
           (the reciprocal of a value<1 to a large exponent.) KLK 7Mar12*/
        if ((soilwater <= Wmin) || soilwater >= Wmax || Wopt <= Wmin) {
          watr_adjust = 0.1;
        } else {
          watr_adjust = pow((soilwater - Wmin) / (Wopt - Wmin), 0.4) *
                        pow((soilwater - Wmax) / (Wopt - Wmax),
                             (0.4 * (Wmax - Wopt)) / (Wopt - Wmin));
          watr_adjust = max(0.1, watr_adjust);
        }
        if (*isagri) {
          if (Dopt < 0.15) {
            agri_adjust = 0.9;
          } else if (Dopt > 0.28) {
            agri_adjust = 0.28;
          } else {
            agri_adjust = -4.6 * Dopt + 1.6;
          }
        } else {
          agri_adjust = 1.0;
        }
        temp_adjust = (soiltemp * max(0.11, Dopt) * 0.095) + 0.9;
        /*printf("CH4 %5.2f %5.2f %5.2f %5.2f\n", CH4max, watr_adjust, temp_adjust, agri_adjust);*/
        *CH4_oxid = CH4max * watr_adjust * temp_adjust * agri_adjust;
      }

      return;
    }


/*****************************************************************************
**
**  FILE:      methane_emission.c
**
**  FUNCTION:  void methane_emission()
**
**  PURPOSE:   Calculate methane emission via plants and bubbles.
**
**  DESCRIPTION:
**    Methane emission is a function of aboveground plant size and root
**    biomass.
**
**  REFERENCE:
**    A semi-empirical model of methane emission from flooded rice paddy
**    soils.  Y. Huang, R.L. Sass, and F.M. Fisher, Jr.  1998.  Global Change
**    Biology 4: 247-268.
**
**    Modeling methane emission from rice paddies with various agricultrual
**    pratices.  Y. Huang, W. Zhang, X. Zheng, J. Li, and Y. Yu.  Journal of
**    Geophysical Research.  2004.  19: DO8113
**
**  HISTORY:
**    On Jul 1, 2014, at 9:46 AM, Keough,Cynthia <Cindy.Keough@colostate.edu> wrote:
**    Steve William's testing of methane emission submodel showed that when
**    root biomass was very small that the methane emission via bubbles was larger
**    than the methane produced.  I have attached a copy of the modified subroutine
**    and notes on the modifications that I made to the source code.
**
**  INPUTS:
**    aglivc         - aboveground live carbon (g/m^2)
**    avgst_10cm     - average soil temperature in top 10 cm of soil profile
**                     (degrees C)
**    bglivc         - belowground live carbon (g/m^2)
**    CH4_Ep         - methane emitted via plants (g/m^2)
**    CH4_Ebl        - methane emitted via bubbles (g/m^2)
**    CH4_prod       - methane production (gC/m2/day)
**    tmxbio         - theoretical maximum biomass where no more methane is
**                     emitted by the plant (g biomass/m^2)
**    zero_root_frac - fraction of CH4 emitted via bubbles when there is no
**                     root biomass (0.0 - 1.0)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    agbiomass       - aboveground live biomass (g biomass/m^2)
**    bgbiomass       - belowground live biomass (g biomass/m^2)
**    CH4_prod_remain - methane production minus methane emitted via plants (g/m^2)
**    Fp              - fraction of produced methane emitted via plants (0.0 - 0.55)
**    P0              - criterion when bubbles occur (g/m^2)
**
**  OUTPUTS:
**    CH4_Ep  - methane emitted via plants (g/m^2)
**    CH4_Ebl - methane emitted via bubbles (g/m^2)
**
**  CALLED BY:
**    trace_gas_model()
**
**  Changes:
**    - replaced 0.55, with frCH4emit the maximum fraction of methane
**      emitted by rice KLK 27 Jul 2015
**
**  CALLS:
**    None
**
*****************************************************************************/

    void methane_emission(float *aglivc, float *tmxbio, float *bglivc,
                          float *avgst_10cm, float *CH4_prod, float *CH4_Ep,
                          float *CH4_Ebl, float zero_root_frac, float frCH4emit,
                          float ch4rootlim)
    {
      float P0 = 0.002f;
      float Fp;
      float CH4_prod_remain;

      /* printf("methane_emission aglivc %18.15f tmxbio %18.15f\n",*aglivc,*tmxbio);  fflush(stdout); */
      /* Methane emission via plants
        was: Fp = 0.55f * (float)pow(max(1.0 - (*aglivc * 2.5f / *tmxbio),0.0), 0.25);
      */
      if (*aglivc > 0.0  &&  *tmxbio > 0.0) {
        Fp = frCH4emit * (float)pow(max(1.0 - (*aglivc * 2.5f / *tmxbio),0.0), 0.25);
        *CH4_Ep = Fp * *CH4_prod;
      } else {
        *CH4_Ep = 0.0f;
      }
        /* printf("methane_emission   CH4_Ep= %6.4f\n",*CH4_Ep);  fflush(stdout); */

      /* Methane emission via bubbles */
      CH4_prod_remain = max(*CH4_prod - *CH4_Ep - P0, 0.0f);
        /* printf("methane_emission   CH4_prod_remain= %6.4f\n",CH4_prod_remain);  fflush(stdout); */

      *CH4_Ebl = 0.0f;
      if (*avgst_10cm > 1.0  &&  CH4_prod_remain > 0.0) {
        /* CH4 emitted via bubbles when there is no root biomass at optimum temperature*/
        *CH4_Ebl = zero_root_frac * CH4_prod_remain;

        if(*avgst_10cm < 2.718282f) {*CH4_Ebl *= (float)log(*avgst_10cm);} /* temperature effect */

        /* decrease CH4 bubbles above a user defined root biomass threshold (was 1.0) */
        if ((*bglivc * 2.5f) > ch4rootlim) {*CH4_Ebl *= (ch4rootlim / (*bglivc * 2.5f));}

        /* make sure we don't exceed the total potential */
        if (*CH4_Ebl > CH4_prod_remain) {*CH4_Ebl = CH4_prod_remain;}

        /* printf("methane_emission   CH4_Ebl=%6.4f zero_root_frac= %6.4f CH4_prod_remain= %6.4f avgst= %6.4f  bglivB=  %6.4f ch4rootlim= %6.4f\n",
           *CH4_Ebl, zero_root_frac, CH4_prod_remain, *avgst_10cm, (*bglivc * 2.5f), ch4rootlim);  fflush(stdout); */
      }

      return;
    }
