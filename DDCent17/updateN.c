
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      updateN.c
**
**  FUNCTION:  void update_npool()
**
**  PURPOSE:   This subroutine updates the ammonium and nitrate pools
**             whenever a soil minerl layer is updated.
**
**  INPUTS:
**    amt       - amount of minerl N added or subtracted from minerl(clyr,N)
**                (g/m2) amt < 0 is a loss from minerl(clyr,N),
**                amt > 0 is a gain to minerl(clyr,N)
**    clyr      - Century soil mineral layer (1..nlayer) to which amt
**                has been added.  Since this function will be called from
**                FORTRAN, subtract 1 from this index
**    frac_nh4  - the fraction of amt to be added to the ammonium pool
**    frac_no3  - the fraction of amt to be added to the nitrate pool
**    nameid[]  - name of calling subroutine
**    nitrate[] - total nitrate in soil mineral pool, distributed among soil
**                water model layers (gN/m2)
**    subname   - name of calling routine, space delimited max 11 characters
**                NO longer modified: copy to a local string to prevent termination from
**                causing a memory crash KLK2018/08/18
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    layers          - soil water soil layer structure
**    layers->lbnd[]  - the index of the lower soil water model layer which
**                      corresponds to clyr in Century
**    layers->ubnd[]  - the index of the upper soil water model layer which
**                      corresponds to layer clyr in Century
**    layers->width[] - the thickness of soil water model layers (cm)
**
**  LOCAL VARIABLES:
**    cum_no3       - accumulator for amount added to nitrate pool (g/m2)
**    cum_tot       - accumulator for amount added to ammonium pool and
**                    nitrate pool (g/m2)
**    debug         - flag to set debugging mode, 0 = off, 1 = on
**    ilyr          - current layer in the soil profile
**    nh4amt        - the amount to be added to the ammonium pool (g/m2)
**    no3amt        - the amount to be added to the nitrate pool (g/m2)
**    no3_clyr_psum - total amount of nitrate in Century layer clyr (g/m2)
**    callr         - name of calling routine for debug prints
**    tdepth        - the depth of Century layer clyr (cm)
**
**  OUTPUTS:
**    ammonium  - total ammonium in soil mineral pool (gN/m2)
**    nitrate[] - total nitrate in soil mineral pool, distributed among soil
**                water model layers (gN/m2)
**
**  CALLED BY:
**    detiv()
**    firrtn()
**    grem()
**    growth()
**    partit()
**    pschem()
**    simsom()
**    trees()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "mssg.h"
#include "soilwater.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
   int warncnt = 0;

    void update_npool(int *clyr, float *amt, double *frac_nh4,
                      double *frac_no3, double *ammonium, double nitrate[],
                      char *subname, long int strlen)
    {

      int    debug = 0;
      int    ilyr;
      double nh4amt, no3amt;
      double tdepth;
      double cum_no3;
      double cum_tot;
      double no3_clyr_psum;

      char callr[12] = {'\0'}; /* set array to initial all NUL bytes */

      extern LAYERPAR_SPT layers;

      if (fabs(*amt) < 1.0E-30) {
        return;
      }

      /* Initializations */
      ilyr = 0;
      while(subname[ilyr] != ' '  &&  ilyr <= 10) {callr[ilyr] = subname[ilyr]; ilyr++;}
      callr[ilyr] = '\0';

      tdepth = 0.0;
      cum_no3 = 0.0;
      cum_tot = 0.0;
      no3_clyr_psum = 0.0;

      nh4amt = *frac_nh4 * (double)(*amt);
      no3amt = *frac_no3 * (double)(*amt);
      if (debug) {
        fprintf(stdout, "\nUPDATE %s: clyr=%1d: amt = %12.10f\n", callr,
                *clyr, *amt);
        fprintf(stdout, "UPDATE %s: no3amt = %12.10f\t frac_no3 = %12.10f\n",
                callr, no3amt, *frac_no3);
        fprintf(stdout, "UPDATE %s: nh4amt = %12.10f\t frac_nh4 = %12.10f\n",
                callr, nh4amt, *frac_nh4);
      }

      *ammonium += nh4amt;

      if (*frac_no3 > 0.0) {
        if ((*clyr == 1) && (no3amt >= 0.0)) {
          /* ADD nitrate to top 15 cm */
          nitrate[0] += no3amt*0.07;
          nitrate[1] += no3amt*0.20;
          nitrate[2] += no3amt*0.73;
          cum_no3 = no3amt;
        } else {
          for(ilyr = layers->ubnd[*clyr-1]; ilyr <= layers->lbnd[*clyr-1];
              ilyr++) {
            tdepth += layers->width[ilyr];
            if (debug > 1) {
              fprintf(stdout, "nitrate[%1d] = %12.10f\n", ilyr,
                      nitrate[ilyr]);
            }
            if (nitrate[ilyr] > 0.0) {
              no3_clyr_psum += nitrate[ilyr];
            }
          }
          if (debug) {
            fprintf(stdout, "no3_clyr_psum = %12.10f\n", no3_clyr_psum);
          }
          for(ilyr = layers->ubnd[*clyr-1]; ilyr <= layers->lbnd[*clyr-1];
              ilyr++) {
            /* allocate a demand against negative nitrate as if both
               were positive (no3amt >= 0.0)
              KLK 5Nov2011*/
            if (no3amt >= 0.0  ||  (no3_clyr_psum == 0. && no3amt < 0.)) {
              /* ADD nitrate */
              if (debug) {fprintf(stdout, "Add to layer %1d\n", ilyr);}
              nitrate[ilyr] += ((double)(layers->width[ilyr])/tdepth)*no3amt;
              cum_no3 += ((double)(layers->width[ilyr])/tdepth)*no3amt;
            } else {
              /* ADD or REMOVE nitrate */
              if (nitrate[ilyr] > 0.0) {
                if (debug) {
                  fprintf(stdout, "   %1d: %12.10f + %12.10f = ", ilyr,
                          nitrate[ilyr],(nitrate[ilyr]/no3_clyr_psum)*no3amt);
                }
                cum_no3 += (nitrate[ilyr]/no3_clyr_psum) * no3amt;
                nitrate[ilyr] += (nitrate[ilyr]/no3_clyr_psum) * no3amt;
                if (debug) {fprintf(stdout, "%12.10f\n", nitrate[ilyr]);}
              }
            }
          }
        }
      }
      cum_tot = cum_no3 + nh4amt;
      if (debug) {
        fprintf(stdout, "In update_npool: nh4amt = %12.10f\n", nh4amt);
        fprintf(stdout, "In update_npool: cum_no3 = %12.10f\n", cum_no3);
        fprintf(stdout, "In update_npool: cum_tot = %12.10f\n", cum_tot);
      }

      /*changed this condition to only check cum_no3 when it can be set.
        cum_no3 will NOT be set when nitrate < 0 and amt < 0
         if (fabs(cum_no3-no3amt) > 1.0E-9) {*/
      if (fabs(cum_no3-no3amt) > 1.0E-9) {
        /* if(no3_clyr_psum >0. || no3amt > 0.0) {*/
        sprintf(msgbfr, "NO3 loss in update_npool: %s.  no3amt = %12.10f,  cum_no3 = %12.10f",
                callr, no3amt, cum_no3);
         abortmssg(msgbfr);
      }
      if(no3_clyr_psum==0. && no3amt < -1.0E-7) {
        if(warncnt++ < 10) {
          fprintf(stdout,
          "Warning: No nitrate available update_npool: %s.  no3amt = %12.10f",
              callr, no3amt);
          if(warncnt == 10) {fprintf(stdout, "         further messages suppressed");}
          fprintf(stdout, "\n");
        }
      }

      return;
    }
