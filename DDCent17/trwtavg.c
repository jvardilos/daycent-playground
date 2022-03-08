
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      trwtavg.c
**
**  FUNCTION:  float trwtavg()
**
**  PURPOSE:   Compute weighted average of soil water potential to be
**             used for transpiration calculations.
**
**  AUTHOR:    Susan Chaffee  4/30/92
**
**  REWRITE:   Melannie Hartman  9/20/93 - 9/20/93
**
**  HISTORY:
**    11/09/95 (MDH) Added PILPS modifications.
**    05/19/08 (CAK) Modify subroutine to remove the code that is using
**                   shallow, intermediate, deep, and very deep soil depths to
**                   calculate a weighted average value for calculating
**                   transpiration.  Replace this code with a routine that
**                   will return the soil water potential of the wettest soil
**                   layer within the plant rooting zone to be used for
**                   calculating transpiration.
**    07/25/11 (CAK) Return a 0.0 to 1.0 multiplier instead of the soil water
**                   potential of the wettest soil layer within the plant
**                   rooting zone.
**     6/24/14 (KLK) removed the unused callname string to avoid compiler warning.
**                   Moved swp calculation, exponential, outside the max loop.
**                   max swp = swp max rwcf
**
**  INPUTS:
**    flags          - structure containing debugging flags
**    flags->debug   - flag to set debugging mode, 0 = off, 1 = on
**    layers         - soil water soil layer structure
**    layers->lbnd[] - the index of the lower soil water model layer which
**                     corresponds clyr in Century
**    layers->swc[]  - soil water content by layer (cm H2O)
**    nlaypg         - number of layers in plant rooting zone
**
**  GLOBAL VARIABLES:
**
**  LOCAL VARIABLES:
**    callname   - call name for subroutine
**    ilyr       - current layer in the soil profile
**    swp        - soil water potential accumulator
**    rwcf       - relative water content
**    swptemp    - intermediate variable in calculation of soil water potential
**
**  OUTPUTS:
**    swp - soil water potential
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
**
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "soilwater.h"

    float trwtavg(int *nlaypg, LAYERPAR_SPT layers, FLAG_SPT flags)
    {
      int ilyr;
      float swp;  /* soil water potential accumulator */
      double rwcf;
      #define  callname = "trwtavg"

      if (flags->debug > 2) {
        printf("Entering function trwtavg\n");
      }

      /* Find the wettest layer in the plant rooting zone */
      rwcf = 0.0;
      for(ilyr=0; ilyr <= layers->lbnd[*nlaypg-1]; ilyr++) {
        rwcf = max((layers->swc[ilyr]/layers->width[ilyr] - layers->swclimit[ilyr]) /
                   (layers->fieldc[ilyr] - layers->swclimit[ilyr]),
                   rwcf);
     }
     swp = 1.0f - expf(-3.0f * (float)rwcf);
     /*swp = max(1.0f - (float)exp(-3.0 * rwcf), swp);  ANSI compatible double EXP*/

      if (flags->debug > 2) {
        printf("Exiting function trwtavg\n");
      }

      return(swp);
    }
