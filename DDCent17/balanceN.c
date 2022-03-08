
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      balanceN.c
**
**  FUNCTION:  void bal_npool()
**
**  PURPOSE:   This routine updates the minerl[:,N] distribution after the
**             trace_gas_model calculates the nitrogen losses and the leaching.
**             It also checks the total N in the pools for N balance errors.
**
**  INPUTS:
**    ammonium  - total ammonium in soil mineral pool (gN/m2)
**    inorglch  - N from organic leaching of stream flow (base flow + storm
**                flow) (g/m2)
**    minerl[]  - soil minerl pool, 1..nlayer (g/m2)
**                Minerl is a 2-D array in FORTRAN.  FORTRAN stores multi-
**                dimensioned arrays in column-major order.  In this routine
**                it is a 1-D array, where minerl(clyr+1,iel+1) in FORTRAN is
**                equivalent to minerl[(iel)*CMXLYR+(clyr)] in this routine.
**    nitrate[] - total nitrate in soil mineral pool, distributed among soil
**                water model layers (gN/m2)
**    nlayer    - number of layers in Century soil profile
**
**  GLOBAL VARIABLES:
**    CMXLYR    - maximum number of Century soil layers (10)
**    MXSWLYR   - maximum number of soil water model layers (21)
**
**  EXTERNAL VARIABLES:
**    layers          - soil water soil layer structure
**    layers->lbnd[]  - the index of the lower soil water model layer which
**                      corresponds to clyr in Century
**    layers->numlyrs - total number of layers in the soil water model profile
**    layers->ubnd[]  - the index of the upper soil water model layer which
**                      corresponds to layer clyr in Century
**    Nloss           - the major change in minerl(N). Used to evaluate N balance
**  LOCAL VARIABLES:
**    clyr       - Century soil mineral layer (1..nlayer)
**    iel        - current element, set to 0 for nitrogen
**    ilyr       - current layer in the soil profile
**    minerl_sum - minerl[], summed over soil layers
**    npool_sum  - nitrate[], summed over soil layers
**
**  OUTPUTS:
**    minerl[] - soil minerl pool, 1..nlayer (g/m2)
**               Minerl is a 2-D array in FORTRAN.  FORTRAN stores multi-
**               dimensioned arrays in column-major order.  In this routine
**               it is a 1-D array, where minerl(clyr+1,iel+1) in FORTRAN is
**               equivalent to minerl[(iel)*CMXLYR+(clyr)] in this
**               routine.
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    None
**
**  CHANGES:
**   April 2016 KLK
**    recoded N balance print statements so there is 1 set of prints
**    with different error messages.
**    -  Set warning level to about 1.5e-5 of the total mineral N
**    -  Left Error message at 0.1 g/m2
**    minerl print now shows both the old and new values.
**    The error can be routed to STDERR but removed the profile output to STDERR.

*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "soilwater.h"

    static int    debug   = 0;
    static int    errprnt = 0;

    void printnpool(int *doy, int *cyear, int *nlayer, float minerl[], double *newminrl, double *ammonium, double nitrate[])
    {
      int    clyr, ilyr;
      const int    iel = 0;           /* Nitrogen */
      double npool_sum, minerl_sum;

      extern LAYERPAR_SPT layers;

      /* Sum ammonium and nitrate pools */
      npool_sum = nitrate[0];
      for (ilyr=1; ilyr < MXSWLYR; ilyr ++) {
        npool_sum += nitrate[ilyr];
      }
      npool_sum += *ammonium;

      minerl_sum = 0;
      for (clyr=0; clyr < CMXLYR; clyr++) {
        minerl_sum += minerl[iel*CMXLYR + clyr];
      }

 fprintf(stdout, "%d day %d npool_sum = %.10f  minerl_sum = %.10f  npool-minrl = %.4g    newminrl = %.4g\n",
          *cyear, *doy, npool_sum, minerl_sum, npool_sum - minerl_sum, *newminrl);
      return;
    }



    void bal_npool(int *doy, int *cyear, int *nlayer, float minerl[], double *ammonium,
                   double nitrate[], double *inorglch, double Nloss)
    {
      int    clyr, ilyr;
      const int    iel = 0;           /* Nitrogen */
      float  oldminrl[CMXLYR];
      double npool_sum, minerl_sum,  nmax;
      const float balwarn = pow(2.0, -17); /* 2-17 or < 7.6e-6 */

      extern LAYERPAR_SPT layers;

      /* Sum ammonium and nitrate pools */
      npool_sum = nitrate[0];
      for (ilyr=1; ilyr < MXSWLYR; ilyr ++) {
        npool_sum += nitrate[ilyr];
      }
      npool_sum += *ammonium;

      minerl_sum =0.0;
      nmax  = fabs(Nloss); /* keep a maximum entry in the sum */
      if(fabs(*ammonium) > nmax) {nmax = fabs(*ammonium);}
      for (clyr=0; clyr < CMXLYR; clyr++) {
        if(fabsf(minerl[clyr]) > nmax) {nmax  = fabsf(minerl[clyr]);}
        oldminrl[clyr] = minerl[clyr];
        minerl_sum += minerl[iel*CMXLYR + clyr];
        minerl[iel*CMXLYR + clyr] = 0.0f;
      }
      if(fabs(minerl_sum) > nmax) {nmax = fabs(minerl_sum);}
      minerl_sum -=*inorglch; /* start sum with the leaching loss */

      /* Sum ammonium and nitrate pools and compare to sum of minerl pool */

if(debug >1) {
 fprintf(stdout, "%d day %d npool_sum = %.10f  minerl_sum = %.10f  npool-minrl = %.4g,   Nloss= %.4g  nmax=%.4g   err = %.4g\n",
          *cyear, *doy, npool_sum, minerl_sum, npool_sum - minerl_sum, Nloss, nmax,
          fabs(npool_sum - minerl_sum)/nmax);
}

      minerl[iel*CMXLYR + 0] = (float)*ammonium;
      /* minerl(1,N) = ammonium */
      for (clyr=0; clyr < *nlayer; clyr++) {
        for(ilyr = layers->ubnd[clyr]; ilyr <= layers->lbnd[clyr]; ilyr++) {
          minerl[iel*CMXLYR + clyr] += (float)nitrate[ilyr];
        }
      }
      minerl[iel*CMXLYR + *nlayer] = (float)nitrate[layers->numlyrs];
      /* compare ammonium and nitrate pools to the original minerl pool sum*/

      if (fabs(npool_sum- minerl_sum) > nmax*balwarn) {
      /* balance beyond threshold. do we need a warning */

        if (fabs(npool_sum - minerl_sum) > 0.1) { /* severe error! at one time this was fatal */
          fprintf(stdout, "ERROR! N balance failure %d/%d  npool= %13.10f  minerl= %.8f  npool-minrl= %.5g, err= %.3E\n",
                *doy, *cyear, npool_sum, minerl_sum, npool_sum - minerl_sum,
                (npool_sum - minerl_sum)/nmax);
          if (errprnt) {fprintf(stdout, "ERROR! N balance failure %d/%d\n",*doy,*cyear);}
        }

        else if (fabs(npool_sum) < 0.001  ||
               (fabs(npool_sum- minerl_sum) > nmax/32768.0 )) /* 3e-5 */ {
          fprintf(stdout, "N balance warning %d/%d  npool= %13.10f  minerl= %.8f  npool-minrl= %.5g, err= %.3E\n",
                *doy, *cyear, npool_sum, minerl_sum, npool_sum - minerl_sum,
                (npool_sum - minerl_sum)/nmax);
          if (errprnt) {fprintf(stdout, "N balance warning %d/%d.\n",*doy,*cyear);}
        }

        if(debug) {
          fprintf(stdout, "ammonium   = %13.10f\n", *ammonium);
          for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {
            fprintf(stdout, "nitrate[%1d] = %13.10f\n", ilyr, nitrate[ilyr]);
          }
          for (clyr=0; clyr <= *nlayer; clyr++) {
            fprintf(stdout, "minerl(%1d, N)    old %13.10f      new  %13.10f\n", clyr+1,
                       oldminrl[clyr], minerl[iel*CMXLYR + clyr]);
          }
        }
      }

      return;
    }
