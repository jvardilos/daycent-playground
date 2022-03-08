
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      fracbslos.c
**
**  FUNCTION:  void fracbslos()
**
**  PURPOSE:   Calculate fraction of water loss from bare soil
**             evaporation and transpiration
**
**  REWRITE:   Melannie Hartman	9/22/93 - 9/29/93
**
**  HISTORY:
**    11/11/2017 KLK  prevent underflow by setting fbse = 0   when fbst == 1
**    11/11/2017 KLK  changed calculation to match the form used by the FAO Irrigation
**                    and Drainage Paper, No 56, Crop Evapotranspiration  ch 9 eq(97) pg 186
**                    This changes the fraction so transpiration is zero at 0 biomass
**                    and soil evaporation has a minimum of 5% even with infinit canopy.
**                    This is about half of the FAO 10% for their reference canopy.
**    4/30/92  (SLC)
**
**  INPUTS:
**    blivelai - live biomass leaf area index
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    bsepar1 - parameter used to compute transpiration fraction of water loss
**    bsepar2 - parameter used to compute transpiration fraction of water loss
**    bsemax  - maximum bare soil evaporation fraction  (depricated)
**    minbse  - minimum bare soil evaporation fraction
**
**  OUTPUTS:
**    fbse - fraction of water loss from bare soil evaporation
**    fbst - fraction of water loss from bare soil transpiration
**
**    NOTE: fbse = expf(-blivelai*bsepar1) + bsepar2;
**          fbse + fbst = 1.0;
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
****************************************************************************/

#include <math.h>

    void fracbslos(float *fbse, float *fbst, float blivelai)
    {

/* Initializations */
#define  explim -19.f   /* smaller exponents round below 1 ULP causing no roundoff difference*/
#define  bsepar1 1.5f

      *fbst = -blivelai*bsepar1; /* temporary storage for the exponent */

/*    original equation
#define  bsepar2  0.0f
#define  bsemax   0.995f
      *fbse = fmin(bsemax, (*fbst < explim)? bsepar2: ((*fbst>=0)? 1: expf(*fbst)) + bsepar2);
*/

/* use the FAO ch 9 eq(97) pg 186 equation  cast into the *fbse form
      *fbst= (*fbst>=0)? 0: (*fbst < explim)? (1.0-minbse): (1.0-minbse)*(1-expf(*fbst));
*/
#define  minbse 0.05f
      *fbse = (*fbst>=0)? 1.0: (*fbst < explim)? minbse: (1.0-minbse)*expf(*fbst) + minbse;

      *fbst = 1 - *fbse;

      return;
    }
