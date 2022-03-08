
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      swpotentl.c
**
**  FUNCTION:  float swpotentl()
**
**  PURPOSE:   Calculate the soil water potential of a soil layer, given its
**             soil water content.
**
**  AUTHOR:    Susan Chaffee    April 2, 1992
**
**  REWRITE:   Melannie Hartman  9/20/93 - 10/6/93
**
**  HISTORY:
**    09/01/92 (SLC) If swc comes in as zero, set swpotentl to upperbnd.
**                   (Previously, we flagged this as an error, and set
**                   swpotentl to zero).
**
**  DESCRIPTION:
**    The equation and its coefficients are based on a paper by Cosby,
**    Hornberger, Clapp, Ginn in WATER RESOURCES RESEARCH June 1984.
**    Moisture retention data was fit to the power function
**      soil water potential = psis*(theta/thetas)**(-b)
**
**  COMMENT:
**    See the routine "watreqn" for a description of how the variables
**    psis, b, and thetas are initialized.
**
**  INPUTS:
**    callname         - call name for subroutine which called swpotentl
**                       function
**    ilyr             - current layer in the soil profile
**    layers           - soil water soil layer structure
**    layers->b[]      - slope of retention curve
**    layers->psis[]   - "saturation" matric potential of "ilyr" (cm H2O ?)
**    layers->thetas[] - volumetric soil water content at saturation for layer
**                       (% volume)
**    layers->width[]  - the thickness of soil water model layers (cm)
**    swc              - soil water content of the current layer (cm H2O)
**
**  LOCAL VARIABLES:
**    base     - base value for power function
**    expon    - exponent value for power function
**    theta    - volumetric soil water content * 100
**    upperbnd - upper bound for soil water potential (bars)
**
**  GLOBAL VARIABLES:
**    BAR2CM  - conversion factor for bars to centimeters H2O (1024)
**              (1 bar = 1024 cm H2O)
**
**  OUTPUTS:
**    swptnl - soil water potential of the current layer (bars)
**
**  CALLED BY:
**    h2oflux()
**    initdaily()
**    potbse()
**    soiltransp()
**    trwtavg()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mssg.h"
#include "soilwater.h"

    float swpotentl(double swc, int ilyr, LAYERPAR_SPT layers,
                    char callname[])
    {

      float upperbnd = 80.0f;
      float theta;
      float swptnl;
      double base;
      double expon;

      /* get the soil water content of the current layer, needed to */
      /* compute soil water potential. */

      if (swc > 0) {
        theta = ((float)swc / layers->width[ilyr])*100;
        base =  (double)(theta / layers->thetas[ilyr]);
        expon = (double)(layers->b[ilyr]);
        swptnl = (layers->psis[ilyr] / (float)pow(base,expon)) / BAR2CM;
      } else {
        swptnl = upperbnd;
        printf("Potential problem (%s) with swpotentl", callname);
        printf("\tswc[%1d] = %10.8f\n", ilyr, swc);
        fprintf(stderr, "Potential problem (%s) with swpotentl, ", callname);
        fprintf(stderr, "\tswc[%1d] = %10.8f\n", ilyr, swc);
      }

/*      printf("theta = %7.5f  swc = %6.4f  swptnl = %7.4f base = %5.3f  ",
             theta, swc, swptnl, base); */

      if (swptnl > upperbnd) {
        swptnl = upperbnd;
      }

      return(swptnl);
    }


/*****************************************************************************
**
**  FILE:      watreqn.c
**
**  FUNCTION:  void watreqn()
**
**  PURPOSE:   Compute parameters to be used in equations to compute water
**             content or matric potential based on the regression equations
**             by Cosby, Hornberger, Clapp, Ginn (see WATER RESOURCES
**             RESEARCH, Vol.20, No.6, pp.682-690, June 1984).
**             These are a set of regression equations based on soil texture
**             composition as a function of water content and matric
**             potential.
**
**  DESCRIPTION:
**    The regression equation can be used to compute either soil water
**    content:
**         soil water content=(thetas*(psis/(bars*barconv))**b)*0.01
**  OR
**    used to compute matric potential:
**         matric potential = (psis/(theta/thetas)**b)/BARCONV
**
**    where:    BARCONV = 1024 cm water per bar
**              bars    = matric potential at which to compute soil water
**              theta   = volumetric soil water content * 100.
**
**  AUTHOR:  Susan Chaffee    March 23, 1992
**
**  REWRITE: Melannie Hartman  9/10/93 - 9/29/93
**           KLK 28Sept11  diverted err print to output routine
**           KLK 14Aug15   removed unused binverse
**
**  INPUTS:
**    clay  - amount of clay in the soil (fraction, 0.0 - 1.0)
**    sand  - amount of sand in the soil (fraction, 0.0 - 1.0)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    err - error check, used to check for divide by zero
**
**  OUTPUTS:
**    b        - slope of the retention curve
**    psis     - "saturation" matric potential (cm H2O ?)
**    thetas   - volumetric soil water content at saturation (% volume)
**
**  CALLED BY:
**    initlyrs()
**
**  CALLS:
**    None
**
*****************************************************************************/


    void watreqn(float sand, float clay, float *thetas, float *psis, float *b)
    {

      float err = 0.00001f;

      *thetas = -14.2f*sand - 3.7f*clay + 50.5f;

      *psis = (float)pow(10.0, (double)(-1.58f*sand - 0.63f*clay + 2.17f));

      *b = -0.3f*sand + 15.7f*clay + 3.10f;

      if (*b < err) {
        abortmssg("possible division by zero in function 'watreqn()'");
      }

      return;
    }