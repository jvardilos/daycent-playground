
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      tgmodel.c
**
**  FUNCTION:  void trace_gas_model()
**
**  PURPOSE:   Trace Gas Model - calculates daily N2, N2O, NO fluxes from the
**             soil.
**
**  HISTORY:
**     New module routine, pulled from dailymoist. -mdh 6/20/00
**
**     fixed a N leak into the system.             -KLK April 2016
**     - immobilization now can set the ammonium and surface nitrate
**        pools negative. This mirrors the behavior of minerl(1) used
**        in the decomposition model and is REQUIRED to signal that
**        model to limit immobilization.
**     - recoded the immobilization use of nitrate in the soil profile.
**       This calculation probably should be reconsidered since it is not
**       tied to any mineral flow process or the depth of the soil SOM
**        a  prevent negative nitrate in layers below the surface.
**        b  no longer access nitrate in deep storage layer
**     - disable NOflux and Dn2oflux flux calculations when potential_NOflux
**       is less thane or equal to zero. Also disable NOflux block
**       when the ammonium is negative. The later is required to
**       keep the code from resetting an ammonium deficit.
**     - added a block for nitrate from nitrification to "re-pay" negative
**       nitrate layers values BEFORE denitrify distributes the excess.
**       NOTE: distributing nitrate using root distribution should be
**          reconsidered since most decomposition occurs at the surface.
**
**
**  INPUTS:
**    afiel      - field capacity in the top 15cm (vol.frac)
**    ammonium   - ammonium concentration (gN/m^2)
**    avgst_10cm - average soil temperature in top 10 cm of soil profile
**                 (degrees C)
**    avgwfps    - avg wfps in top 3 soil layers (~top 10 cm) 0-1  15cm?
**    basef      - fraction of base flow
**    bglivc     - belowground live carbon (g/m^2)
**    tmxbio     - theoretical maximum biomass where no more methane is
**                 emitted by the plant (g biomass/m^2)
**    btolai     - biomass to LAI conversion factor for trees
**    bulkd      - Century soil bulk density (g/cm^3)
**    CH4_Ep     - methane emitted via plants (g/m^2)
**    CH4_Ebl    - methane emitted via bubbles (g/m^2)
**    CH4_oxid   - methane oxidation (g/m^2)
**    CH4_prod   - methane production (g/m^2)
**    clay       - fraction of clay (0.0 - 1.0)
**    co2PPM[]   - CO2 concentration by layer (ppm)
**    Com        - sum of CO2 losses from heterotrophic decomposition of
**                 metabc(1), metabc(2), strucc(1), strucc(2), som1c(1),
**                 som1c(2), som2c(1), som2c(2), and som3c (g/m^2)
**    Cr         - carbohydrates derived from rice plants (g/m^2)
**    critflow   - amount of flow between lyrs to leach (cm H2O)
**    crpstore   - crop storage, NO absorped by canopy is added to this pool
**    Dn2flux    - denitrification N2 flux (gN/m^2/day)
**    Dn2oflux   - denitrification N2O flux (gN/m^2/day)
**    dN2lyr[]   - N2 flux from Denitrification by layer (gN/m2/day)
**    dN2Olyr[]  - N2O flux from Denitrification by layer (gN/m2/day)
**    Eh         - effect of water management on soil redox potential (Eh)
**                 (-250.0 mv - 300.0 mv)
**    Feh        - reduction factor of effect of soil redox potential (Eh) on
**                 methane production (0.0 - 1.0)
**    forstore   - forest storage, NO absorped by canopy is added to this pool
**    frlechd[]  - leaching fractions (?)
**    inorglch   - inorganic N leaching (gN/m^2) (?)
**    isdecid    - flag, set to 1 for deciduous forest, 0 otherwise
**    isagri     - flag, set to 1 for agricultural system or grassland that
**                 has been fertilized or cultivated, 0 otherwise
**    jday       - current day of year (1-366)
**    maxt       - max of montly max temperatures (deg C)
**    newCO2     - amount of soil respiration (gC/m^2/day)
**    newminrl   - net mineralization (gN/m^2/day)
**    nit_amt    - gross nitrification (gN/m^2/day)
**    nitrate[]  - nitrate concentration by lyr (gN/m^2)
**    Nn2oflux   - nitrification N2O flux (gN/m^2/day)
**    NOflux     - NO flux (gN/m^2/day)
**    nreduce    - reduction factor on nitrification rates due to
**                 nitrification inhibitors added with fertilizer (0-1)
**    pHscale    - optional scalar on soil pH
**    ppt        - daily precip (cm)
**    prev_bgprd - previous day's fine root carbon production (g/m^2)
**    rleavc     - leaf carbon (g/m^2)
**    sand       - fraction of sand (0.0 - 1.0)
**    silt       - fraction of silt (0.0 - 1.0)
**    snow       - snow cover (cm SWE)
**    stormf     - fraction of storm flow
**    stream[]   - stream flow (?)
**    texture    - soil texture index (see swconst.h)
**    TI         - soil temperature index for methane production (0.0 - 1.0)
**    time       - current time in years
**    watertable - flag, 1 = simulate water table, 0 = no water table
**    watrflag   - 1 = water added to the system automatically to bring the
**                     soil water content to saturation
**                 0 = rain and irrigation events used to add water to the
**                     soil potentially bringing it to saturation
**    wfluxout[] - amount of flow between lyrs (cm H2O)
**
**  GLOBAL VARIABLES:
**    MXSWLYR - maximum number of soil water model layers (21)
**    PI      - pi (3.14159265)
**
**  EXTERNAL VARIABLES:
**    sitepar               - site parameter structure
**    sitepar->dDO_fc       - the normalized diffusivity in aggregate soil
**                            media at field capacity (calculated by initlyrs)
**    sitepar->dDO_wp       - the normalized diffusivity in aggregate soil
**                            media at wilting point (calculated by initlyrs)
**    sitepar->N2Oadjust_fc - maximum proportion of nitrified N lost as N2O at
**                            field capacity
**    sitepar->N2Oadjust_wp - minimum proportion of nitrified N lost as N2O at
**                            wilting point
**    sitepar->netmn_to_no3 - fraction of net mineralization going to NO3 (0.0-1.0)
**
**  LOCAL VARIABLES:
**    avg_vswc         - average soil water content in top 3 soil layers
**                       (volumetric)
**    canopy_reduction - reduction factor for NO absorbed by canopy (0-1)
**    dDO              - normalized diffusivity in aggregate soil media (0-1)
**    debug            - flag to set debugging mode, 0 = off, 1 = on
**    grass_lai        - amount of LAI in grass/crop leaf canopy (m^2/m^2)
**    ilyr             - current layer in the soil profile
**    surflood         - soil saturation flag, (not used yet) 1 saturate / 0 not
**    krainNO          - increase of NO due to moisture and rain >=1.0
**    newNH4           - new NH4 produced (gN/m^2/day)
**    newNO3           - new NO3 produced (gN/m^2/day)
**    NH4_to_NO        - ammonium converted to NO to reach potential NO flux
**                       (gN/m^2/day)
**    nh4_to_no3       - amount of NH4 converted to NO3 by nitrification
**                       (gN/m^2/day)
**    NO_N2O_ratio     - NO/N2O ratio <= 1.0
**    NOabsorp         - NO absorbed by canopy (gN/m^2/day)
**    Nloss            - the major change in minerl(N) needed to evaluate N balance
**    npool_sum        - sum of N in nitrate and ammonium pools (gN/m^2)
**    potential_NOflux - maximum possible NO flux based on NO/N2O and krain
**                       (gN/m^2/day)
**    stdbulkd         - standard bulk density based on soil texture (g/cm^3)
**    stdfieldc        - standard field capacity based on soil texture
**                       (vol.frac)
**    total_lai        - total LAI at site (m^2/m^2)
**    tree_lai         - amount of LAI in tree leaf canopy (m^2/m^2)
**    turnovfrac       - fraction of new NO3 that goes to N2O
**
**  OUTPUTS:
**    ammonium     - ammonium concentration (gN/m^2)
**    avgst_10cm   - average soil temperature in top 10 cm of soil profile
**                   (degrees C)
**    CH4_Ep       - methane emitted via plants (g/m^2)
**    CH4_Ebl      - methane emitted via bubbles (g/m^2)
**    CH4_oxid     - methane oxidation (g/m^2)
**    CH4_prod     - methane production (g/m^2)
**    co2PPM[]     - CO2 concentration by layer (ppm)
**    Cr           - carbohydrates derived from rice plants (g/m^2)
**    crpstore     - crop storage, NO absorped by canopy is added to this pool
**    dN2lyr[]     - N2 flux from Denitrification by layer (gN/m2/day)
**    dN2Olyr[]    - N2O flux from Denitrification by layer (gN/m2/day)
**    Dn2flux      - denitrification N2 flux (gN/m^2/day)
**    Dn2oflux     - denitrification N2O flux (gN/m^2/day)
**    Eh           - effect of water management on soil redox potential (Eh)
**                   (-250.0 mv - 300.0 mv)
**    Feh          - reduction factor of effect of soil redox potential (Eh) on
**                   methane production (0.0 - 1.0)
**    forstore     - forest storage, NO absorped by canopy is added to this pool
**    inorglch     - inorganic N leaching (gN/m^2) (?)
**    nitamt       - gross nitrification (gN/m^2/day)
**    nitrate[]    - nitrate concentration by lyr (gN/m^2)
**    Nn2oflux     - nitrification N2O flux (gN/m^2/day)
**    NOflux       - NO flux (gN/m^2/day)
**    SI           - soil texture index for methane production (0.0 - 1.0)
**    stream[]     - stream flow (?)
**    TI           - soil temperature index for methane production (0.0 - 1.0)
**
**  CALLED BY:
**    dailymoist() - FORTRAN routine
**
**  CALLS:
**    denitrify()   - calculate N2O flux and N2:N2O ratio due to
**                    denitrification
**    diffusiv()    - estimate normalized diffusivity in soils
**    getsoilprop() - determine soil texture class based on percent of sand,
**                    silt, clay
**    nitrify()     - nitrification, produces N2O gas
**    nox_pulse()   - increase NO due to moisture and rain
**
*****************************************************************************/

#include "soilwater.h"
#include "n2o_model.h"
#include "swconst.h"
#include <math.h>

    void trace_gas_model(double *newminrl, double *ammonium, double nitrate[],
                         int *texture, float *sand, float *silt, float *clay,
                         float *afiel, float *bulkd, float *maxt, float *ppt,
                         float *snow, float *avgwfps, float *stormf,
                         float *basef, float frlechd[], float stream[],
                         double *inorglch, float *critflow, float wfluxout[],
                         float *newCO2, float co2PPM[], float *efscltef,
                         float *time, double *NOflux, double *Nn2oflux,
                         double *Dn2oflux, double *Dn2flux, double *CH4_oxid,
                         int *isdecid, int *isagri, float *aglivc, float *rleavc,
                         float *btolai, float *crpstore, float *forstore,
                         double *nit_amt, float *nreduce, int *cyear, int *jday,
                         float *pHscale, double dN2lyr[], double dN2Olyr[],
                         float *esrsnkN, float minerl[], int *nlayer,
                         float *prev_bgprd, float *Com,                                /* methane */
                         float *avgst_10cm, float *TI, float *Cr, float *Eh,           /* methane */
                         float *Feh, float *CH4_prod, float *CH4_Ep, float *CH4_Ebl,   /* methane */
                         int *watertable, int *watrflag, float *bglivc, float *tmxbio) /* methane */

    {

      /* Local Variables */

      static int    debug = 0;
      int    ilyr;
      int    surflood =0;
      double turnovfrac =  0.012; /* old value stops compiler warning */
      double newNH4;
      double newNO3;
      double nh4_to_no3;
      double krainNO;
      double potential_NOflux;
      double dDO;
      float  stdbulkd;
      float  stdfieldc;
      double NO_N2O_ratio;
      double NH4_to_NO;
      double npool_sum, Nloss;
      double imobfrac;
      float canopy_reduction;
      double NOabsorp;
      float total_lai;
      float grass_lai;
      float tree_lai;

      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;

      *Nn2oflux = 0.0;
      *NOflux = 0.0;
      *Dn2oflux = 0.0;
      *Dn2flux = 0.0;

      /* Compute fraction of new mineralization that is converted to NH4 */
      /* and NO3 */

      if (debug) {
        printf("trace_gas_model newminrl = %6.4f\n", *newminrl); fflush(stdout);
      }

      Nloss = 0.0;
      if (*newminrl <= 0.0) {

        /* Immobilization */
        /* Distribute N loss proportionally between ammonium and nitrate   */
        /* layers.  There is no check that these N pools won't go negative */
        /* once immobilization is accounted for.  It is assumed that the   */
        /* immobilization calculation by the decomp model is moderated by  */
        /* the supply of minerl N.                                         */
        if(nitrate[0] > 0) {Nloss = *newminrl;}
        imobfrac = 0.0;

        /* ignore nitrate in the deep storage layer */
        npool_sum = 0.0;
        for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {
          npool_sum += (nitrate[ilyr] > 0.0) ? nitrate[ilyr] : 0.0;
        }
        npool_sum += (*ammonium > 0.0) ? *ammonium : 0.0;

         /* N available to immobilize; withdraw nitrate from layers
            Don't mine subsurface layers if the surface nitrate is already in deficit
         */
        if(npool_sum > 0 && nitrate[0] > 0) { /* No available N or we already in deficit */
          /* N available to immobilize; calculate fraction available for withdrawals */
          imobfrac   = *newminrl/ npool_sum; /* fraction of N available to withdraw */
          imobfrac   = (imobfrac <= -1)? 1./128.: 1+imobfrac; /* limit to 99.5% of available N */
          *newminrl -= npool_sum * (imobfrac -1.0); /* unfilled immobilization */

          if (*ammonium > 0.0) {
            *ammonium *= imobfrac;  /* += *newminrl * (*ammonium/npool_sum); */
          }
          npool_sum = *ammonium;
          for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {
            if (nitrate[ilyr] > 0.0) {
              nitrate[ilyr] *= imobfrac;  /* += *newminrl * (nitrate[ilyr]/npool_sum); */
              npool_sum += nitrate[ilyr];
            }
          }
        }

        if(*newminrl < 0.0) {
          /* there wasn't enough N, allow the the surface to go negative */
          *ammonium  += *newminrl * (1.0 - (double)sitepar->netmn_to_no3); /* * (0.5 - sitepar->netmn_to_no3 - 0.5);*/
        /*  *newminrl  *= (double)sitepar->netmn_to_no3; */
          /* allocate to the surface the same way as updateN */
          nitrate[0] += *newminrl * (double)sitepar->netmn_to_no3 * 0.07;
          nitrate[1] += *newminrl * (double)sitepar->netmn_to_no3 * 0.20;
          nitrate[2] += *newminrl * (double)sitepar->netmn_to_no3 * 0.73;
          *newminrl   = 0.0;
        }

        newNH4 = 0.0;
        newNO3 = 0.0;

      } else {
        /* Mineralization */
        if(nitrate[0] < 0) {Nloss = *newminrl;}
        newNO3 = *newminrl * (double)sitepar->netmn_to_no3;
        newNH4 = *newminrl - newNO3; /* * (0.5 - (double)sitepar->netmn_to_no3 - 0.5);*/
      }

      if (debug) {
        printf("newNH4 = %6.4f    newNO3 = %6.4f\n", newNH4, newNO3); fflush(stdout);
      }

      *ammonium += newNH4;

      /* Compute the amount of NH4 that is converted to NO3 due to nitrification */

      nitrify(ammonium, &nh4_to_no3, maxt, nreduce, pHscale);
      *nit_amt = nh4_to_no3;
      newNO3  += nh4_to_no3;

      if (debug) {
        printf("texture = %1d\n", sitepar->texture);
        printf("nh4_to_no3 = %6.4f\n", nh4_to_no3);
        printf("maxt = %6.4f\n", *maxt);
        fflush(stdout);
      }

      /* Compute fraction of new NO3 that is converted to N2O and NO */

      krainNO = nox_pulse(ppt, snow);

      getsoilprop(sand, silt, clay, &stdbulkd, &stdfieldc, &sitepar->texture);
      if (debug) {printf("getsoilprop(%6.4f %6.4f %6.4f %6.4f %6.4f)\n",
                  *sand, *silt, *clay, stdbulkd, stdfieldc);  fflush(stdout);}

      /* Use standard field capacity and bulk density according */
      /* to the soil class in the texture triangle -mdh 10/26/99 */
      /*      dDO = diffusiv(afiel(1), bulkd, *avgwfps) */
      /* No, change back to soils.in field capacity and bulk density.  -mdh 6/20/00 */
      dDO = diffusiv(afiel, bulkd, avgwfps); /* not  diffusiv(&stdfieldc, &stdbulkd, avgwfps); */
      if (debug) {printf("dDO = %6.4f = diffusiv(%6.4f %6.4f); newNO3 %6.4f %6.4f\n",
                  dDO, *afiel, *bulkd, newNO3, nh4_to_no3);  fflush(stdout);}

      if (newNO3 > 1.0E-30) {
/*      was *Nn2oflux = newNO3 * 0.02 * sitepar->N2Oadjust;

        turnovfrac = (sitepar->N2Oadjust_wp - sitepar->N2Oadjust_fc) /
                     (sitepar->dDO_wp - sitepar->dDO_fc) *
                     (dDO - sitepar->dDO_wp) + sitepar->N2Oadjust_wp;
        turnovfrac = max(turnovfrac, sitepar->N2Oadjust_wp); //omit  turnovfrac>=0
        turnovfrac = min(turnovfrac, sitepar->N2Oadjust_fc);

        reorder the conditional and the equation to make it obvious that it adds
        a fraction of the input adjustment range based on the current diffusivity
        turnovfrac = N2Oadjust_wp +
                     (N2Oadjust_fc - N2Oadjust_wp) * (dDO - dDO_wp)/(dDO_fc - dDO_wp);
*/
        if     (dDO <= sitepar->dDO_wp) {turnovfrac = sitepar->N2Oadjust_wp;}
        else if(dDO >= sitepar->dDO_fc) {turnovfrac = sitepar->N2Oadjust_fc;}
        else {turnovfrac = sitepar->N2Oadjust_wp +
                           (sitepar->N2Oadjust_fc - sitepar->N2Oadjust_wp) *
                  (dDO - sitepar->dDO_wp)/(sitepar->dDO_fc - sitepar->dDO_wp);
        }
        *Nn2oflux = newNO3 * turnovfrac;
        newNO3 -= *Nn2oflux;

        /* Another update to NO flux calculation -mdh 10/26/99
           NO_N2O_ratio = 15.23 + (35.45*atan(0.676*PI*(10*dDO-1.86)))/PI; */
        NO_N2O_ratio = 8.0 + (18.0*atan(0.75*PI*(10*dDO-1.86)))/PI;
        /* If this is an agricultural system adjust the NO to N2O ratio */
        /* cak - 01/28/03 */
        if (*isagri) {NO_N2O_ratio *= 0.5;}  /* *= 0.2 commented */
        potential_NOflux = NO_N2O_ratio * *Nn2oflux * krainNO;

        if (potential_NOflux <= newNO3) {
          *NOflux = potential_NOflux;
          newNO3 -= *NOflux;
        } else {
          /* take N out of ammonimum to get max NOflux possible */
          NH4_to_NO = min(*ammonium, (potential_NOflux-newNO3));
          *NOflux = newNO3 + NH4_to_NO;
          *ammonium -= NH4_to_NO;
          newNO3 = 0;
        }

        if (*NOflux < 1.0E-30) {
          *NOflux = 0.0;
        }

      } else {
        NO_N2O_ratio = 0.0;
      }

      /* =============================================================
      ... preferentially add nitrate to negative surface layers in the profile
      */
      if (newNO3 > 1.0E-30 && (nitrate[0] < 0 || nitrate[1] < 0 || nitrate[2] < 0)) {
        npool_sum = 0;
        for (ilyr=3; ilyr >= 0; ilyr--) {
          if(nitrate[ilyr] <0  &&  newNO3 > 1.0E-30) {
            if(nitrate[ilyr] + newNO3 > 0.0) {
              npool_sum += nitrate[ilyr];
              newNO3 += nitrate[ilyr];
              nitrate[ilyr] = 0.0;
            } else {
              npool_sum += newNO3;
              nitrate[ilyr] += newNO3;
              newNO3 = 0.0;
              break;
            }
          }
        }
      }

      /* Compute the denitrification N2O (Dn2oflux) and N2 (Dn2flux) fluxes */
      denitrify(newCO2, &newNO3, nitrate, wfluxout, critflow, frlechd, stream,
                basef, stormf, inorglch, Dn2oflux, Dn2flux, stdfieldc,
                stdbulkd, efscltef, co2PPM, dN2lyr, dN2Olyr, jday, &surflood);

      /* Now compute NOflux from denitrification (new calculation */
      /* -mdh 6/1/00 */
/*      potential_NOflux = NO_N2O_ratio * *Dn2oflux * krainNO; */
      /* For denitrification, krainNO is >= 1.0 -mdh 6/22/00 */

      potential_NOflux = NO_N2O_ratio * *Dn2oflux * min(1.0, krainNO);

      if (potential_NOflux <= *ammonium) {
        /* Take all N out of ammonimum pool */
        *NOflux += potential_NOflux;
        *ammonium -= potential_NOflux;
      }
      else if (potential_NOflux > 0) {
        /* skip for 0 potential flux; N values less than lower available limits */
        /* Take N out of available ammonium, then convert some Dn2oflux to NOflux */
        if (*ammonium > 0) { /* don't mess with negative ammonium */
          *NOflux += *ammonium;
          potential_NOflux -= *ammonium;
          *ammonium = 0.0;
        }
        if (potential_NOflux <= *Dn2oflux) {
          *NOflux += potential_NOflux;
          *Dn2oflux -= potential_NOflux;
        }
      }

      /* Compute the amount of the soil NO flux that is absorped by the */
      /* canopy, cak - 09/23/03 */
      grass_lai = (*aglivc * 2.5f) / CONVLAI;
      tree_lai = (*rleavc * 2.5f) * *btolai;
      total_lai = grass_lai + tree_lai;
      if (debug) {printf("total_lai = %6.4f %6.4f %6.4f\n", total_lai, grass_lai, tree_lai);  fflush(stdout);}

      NOabsorp = 0.0;
      if (total_lai > 0.0) {
        /* canopy_reduction appears to be the reabsorbed fraction.
           This equation is a parabola with the minimum about 8.28. It reduces
           absorption for LAI above that so limit LAI to 8.0. The problem seems
           to be the simple biomass to LAI ratio. 200 bu/acre corn is
           aglivC > 1100 and an LAI > 34. Obviously its not all leaves. */
          canopy_reduction = (total_lai > 8.0)? 0.4428f:
                (float)(0.0077 * total_lai*total_lai + -0.13 * total_lai + 0.99);

        /* We need to retain the soil flux value */
        /* NOsoil = *NOflux;  *NOflux *= canopy_reduction;  NOabsorp = NOsoil - *NOflux;
           NOabsorp = NOsoil - *NOflux = *NOflux * (1- canopy_reduction)
           KLK 2Nov12  algebraically reduce the above steps to:
        */
        NOabsorp = *NOflux * (1.0f - canopy_reduction);

        /* NO absorped by canopy goes to crop storage and forest storage */
        *crpstore += (float)(NOabsorp * (grass_lai / total_lai));
        *forstore += (float)(NOabsorp * (tree_lai /  total_lai));

        *NOflux -= NOabsorp; /* reduce the NOflux by absorption
                                NOTE: not sure why but the 2011 version and the
                                      methane version does NOT reduce the NOflux */
      }

      if (*NOflux   < 1.0E-30) {*NOflux = 0.0;}
      if (*Nn2oflux < 1.0E-30) {*Nn2oflux = 0.0;}
      if (*Dn2oflux < 1.0E-30) {*Dn2oflux = 0.0;}
      if (*Dn2flux  < 1.0E-30) {*Dn2flux = 0.0;}

      *CH4_oxid = 0.0f;
      *CH4_prod = 0.0f;
      *CH4_Ep = 0.0f;
      *CH4_Ebl = 0.0f;

      /* Calculate methane oxidation */
      methane_oxidation(CH4_oxid, isdecid, isagri);
      if (debug) {printf("methane_oxidation = %6.4f\n", *CH4_oxid);  fflush(stdout);}

      /* Calculate methane production */
      methane_production(prev_bgprd, Com, avgst_10cm, TI, Cr, Eh, Feh,
                         CH4_prod, watertable, watrflag);
      if (debug) {printf("methane_production CH4_prod = %18.15f\n",*CH4_prod);  fflush(stdout);}

      /* Calculate methane emission via plants and bubbles */
      if (*CH4_prod > 0.0) {
        methane_emission(aglivc, tmxbio, bglivc, avgst_10cm, CH4_prod, CH4_Ep,
                         CH4_Ebl, sitepar->zero_root_frac, sitepar->frCH4emit,
                         sitepar->ch4rootlim);
       if (debug) {printf("methane_emission CH4_prod= %12.7f CH4_Ebl= %12.7f  %12.7f %12.7f %12.7f\n",
          *CH4_prod, *CH4_Ebl, *aglivc, *tmxbio, *bglivc);  fflush(stdout);}
      }

      /* remove fluxes from minerl  and return escaped N to esrsnk(N) */
      *esrsnkN += (*Nn2oflux + *Dn2oflux + *Dn2flux + *NOflux);
      minerl[0] -= (*Nn2oflux + *Dn2oflux + *Dn2flux + (*NOflux + NOabsorp));
      if(Nloss <= 0.0) {Nloss -= (*Nn2oflux + *Dn2oflux + *Dn2flux + (*NOflux + NOabsorp));}

      bal_npool(jday, cyear, nlayer, minerl, ammonium, nitrate, inorglch, Nloss);

      if (debug) {printf("tgmodel return\n");  fflush(stdout);}
      return;
    }
