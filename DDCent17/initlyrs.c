
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initlyrs.c
**
**  FUNCTION:  void initlyrs()
**
**  PURPOSE:   Read in layers of the soil structure for the site.
**
**  AUTHOR:    Susan Chaffee    March 12, 1992
**
**  REWRITE:   Melannie Hatman    9/7/93 - 9/29/93
**
**  HISTORY:
**   Added WFPS initialization calculations; KLK 18Sept13
**   The texture from sand fraction using existing profile is now set; KLK 18Sept13
**   Added a routine to set the soil layers from the site file; KLK
**
**  INPUTS:
**    flags          - structure containing debugging flags
**    flags->debug   - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose - flag to set verbose debugging mode, 0 = off, 1 = on
**    layers         - soil water soil layer structure
**
**  GLOBAL VARIABLES:
**    COARSE     - designates a coarse, sandy soil, texture (1)
**    FINE       - designates a fine soil texture (3)
**    INPTSTRLEN - maximum length of input file line (120)
**    MXSWLYR    - maximum number of soil water model layers (21)
**    MEDIUM     - designates a medium, loamy soil, texture (2)
**
**  LOCAL VARIABLES:
**    sitepar->avgsand  - weighted average for sand fraction in the top 3 soil layers
**    deltamin - minimum volumetric soil water content below wilting point for
**               soil layer
**    fp_in    - pointer to input file
**    ilyr     - current layer in the soil profile
**
**  OUTPUTS:
**  The following values are read from the file soilname per layer.
**    layers->lyblkd[]   - bulk density by layer (g/cm3)
**    layers->clayfrac[] - clay fraction in soil layer, 0.0 - 1.0
**    layers->dpthmn[]   - tops of soil layers (depth from surface in cm)
**    layers->dpthmx[]   - bottoms of soil layers (depth from surface in cm)
**    layers->ecoeff[]   - bare-soil evaporation water absorption coefficients
**                         by layer (ND)
**    layers->fieldc[]   - volumetric water content at field capacity for
**                         layer (cm H2O/cm of soil)
**    layers->orgfrac[]  - organic matter in soil layer, fraction 0.0 - 1.0
**    layers->lyrpH[]    - pH of soil layer
**    layers->sandfrac[] - sand fraction in soil layer, 0.0 - 1.0
**    layers->satcond[]  - saturated hydraulic conductivity by layer (cm/sec)
**    layers->tcoeff[]   - transpiration water absoption coefficients by layer
**                         (ND)
**    layers->wiltpt[]   - volumetric water content at wilting point for layer
**                         (cm H2O/cm of soil)
**
**  Calculated:
**    layers->depth[]     - the distance from the surface to the middle of the
**                          soil layer (cm)
**    layers->nelyrs      - number of layers to consider in evaporation
**    layers->numlyrs     - total number of layers in the soil water model
**                          soil profile
**    layers->sumecoeff   - sum of evaporation coefficients
**    layers->swcfc[]     - volumetric soil water content at field capacity
**                          for layer (cm H2O/cm of soil)
**    layers->swclimit[]  - minimum volumetric soil water content of a layer,
**                          fraction 0.0 - 1.0
**    layers->swcwp[]     - volumetric soil water content at wilting point for
**                          layer (cm H2O)
**    layers->thetas_bd[] - volumetric soil water content at saturation by
**                          layer computed using bulk density (% volume)
**    layers->width[]     - the thickness of soil water model layers (cm)
**    sitepar->texture    - texture classification for trace gas model
**                          (1 = coarse, 2 = medium, 3 = fine)
**
**  The next four parameters are used to compute soil water content or matric
**  potential as a function of soil texture per layer.  See the routine
**  "watreqn" for a description.
**    layers->b[]        - slope of retention curve
**    layers->psis[]     - the saturation matric potential by layer (cm H2O ?)
**    layers->thetas[]   - volumetric soil water content at saturation for
**                         layer (% volume)
**
**  CALLED BY:
**    detiv()
**
**  CALLS:
**    watreqn() - compute parameters used later to calcuate soil water
**                content or matric potential.
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mssg.h"
#include "soilwater.h"
#include "n2o_model.h"

static int read=0; /* set if input data was read from the site file */

/*    void initlyrs(char *soilname, SITEPAR_SPT sitepar, LAYERPAR_SPT layers,
                  FLAG_SPT flags) */
    void initlyrs(SITEPAR_SPT sitepar, LAYERPAR_SPT layers, FLAG_SPT flags)
    {
      int   ilyr;
      float porosity;
      float siltf;
      float osiltf = 0.0f;
      float avg_bd, avg_fc, avg_wp;
      float wfps_fc, wfps_wp;

      if (flags->debug) {printf("Entering function initlyrs\n");}

      /* initialize counters and accumulators */
      layers->nelyrs = 0;
      layers->sumecoeff = 0.0f;

      for(ilyr=0; ilyr<layers->numlyrs; ilyr++) {

        layers->width[ilyr] = layers->dpthmx[ilyr] - layers->dpthmn[ilyr];
        layers->depth[ilyr] = (layers->dpthmn[ilyr] + layers->dpthmx[ilyr])/2.0f;

        /* Set up parameters to be used later for computing matric */
        /* potential or soil water content based on the regression */
        /* equations in the subroutine watreqn. */

        if ((layers->sandfrac[ilyr] <= 0) && (layers->clayfrac[ilyr] <= 0)) {
          abortmssg("soil desription  sandfrac or clayfrac <0");
        }

        /* Use a higher precision subtraction for roundoff then check that the silt
           fraction still does not go negative. Check only if the silt changes,
           so we don't get multiple warnings on the same structure    klk Nov 11 */
        siltf = 0.5f - (layers->sandfrac[ilyr] + layers->clayfrac[ilyr]) + 0.5f;
        if(siltf < -0.001f) {
          abortmssg("soil description  sandfrac + clayfrac >1 (siltf <0)");
        }
        else if(siltf < 0.0f) {
          /* roundoff error subtract double the overage from the largest component.
             This leaves a small positive silt fraction */
          if(layers->sandfrac[ilyr] > layers->clayfrac[ilyr]) {
            layers->sandfrac[ilyr] += (siltf+siltf);
          } else {layers->clayfrac[ilyr] += (siltf+siltf);}
          /* print a warning if the first time this value is changed */
          if(siltf != osiltf) {fprintf(stderr, "adjusting  sandf %2.6f  clayf %2.6f\n",
                 layers->sandfrac[ilyr],layers->clayfrac[ilyr]);}
        }
        osiltf = siltf;

        /* Stop if field capacity <= wilting point */
        if (layers->fieldc[ilyr] <= layers->wiltpt[ilyr]) { /* not swcwp which is wiltpt * width  KLK 3/2018*/
           sprintf(msgbfr,"soil layer %d  field capacity %f <= wilting point %f", ilyr, layers->fieldc[ilyr], layers->swcwp[ilyr]);
           abortmssg(msgbfr);
        }

        porosity = 0.9999f - layers->lyblkd[ilyr] / PARTDENS;
        /* Calculation of thetas_bd added 9/18/00. -cindyk
           save as fraction since it is always used that way  KLK 28 Apr 2015 */
        layers->thetas_bd[ilyr] = 0.95*porosity;
        /* check that field capacity is less than saturation. KLK 14 Nov 2016
           the inversion can occur in the Saxton equations with high (>70%) clay soils.
        */
        if (layers->fieldc[ilyr] > layers->thetas_bd[ilyr]) {
          if(ilyr==0  ||  layers->lyblkd[ilyr] != layers->lyblkd[ilyr-1]) {
            printf("Warning: lowering %s field capacity %2.5f below %2.5f in layer %2d and below",
                   FNSOIL,layers->fieldc[ilyr],porosity,ilyr);
            printf(" to match Bulk density %2.5f\n",layers->lyblkd[ilyr]);
          }
          layers->fieldc[ilyr] = porosity*(243.0F/256.0f); /* about 5% below porosity */

          /* now place saturation between field capacity and porosity */
          layers->thetas_bd[ilyr] = (porosity + layers->fieldc[ilyr])/2.0;
          if(fabsf(porosity-layers->fieldc[ilyr]) > 0.05f * layers->fieldc[ilyr]) {
            sprintf(msgbfr,"porosity-fieldc > 0.05f * fieldc for layer %1d", ilyr);
            abortmssg(msgbfr);
          }
        }

        /* Set up parameters to be used later for computing matric */
        /* potential or soil water content based on the regression */
        /* equations in the subroutine watreqn. */

        watreqn(layers->sandfrac[ilyr], layers->clayfrac[ilyr],
                &layers->thetas[ilyr], &layers->psis[ilyr],
                &layers->b[ilyr]);

        /* Convert volumetric soil water content to  */
        /* soil water content per layer (in cm H2O) */

        layers->swcfc[ilyr]  = layers->fieldc[ilyr] * layers->width[ilyr];
        layers->swcwp[ilyr]  = layers->wiltpt[ilyr] * layers->width[ilyr];

        /* Set up evaporation parameters to be used later.  This is done */
        /* here to save having to do this every time step in watrflow. */

        if (layers->ecoeff[ilyr] > 0.0) {
          layers->nelyrs = ilyr+1;
          layers->sumecoeff += layers->ecoeff[ilyr]*layers->width[ilyr];
        }

        if (flags->debug) {
          printf("%3d\t%5.1f\t%5.1f\t%4.2f\t%6.4f\t%6.4f\t",
                 ilyr, layers->dpthmn[ilyr], layers->dpthmx[ilyr],
                 layers->lyblkd[ilyr], layers->fieldc[ilyr],
                 layers->wiltpt[ilyr]);
          printf("%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%7.5f\t%4.2f\n",
                 layers->ecoeff[ilyr], layers->tcoeff[ilyr],
                 layers->sandfrac[ilyr], layers->clayfrac[ilyr],
                 layers->orgfrac[ilyr], layers->swclimit[ilyr],
                 layers->satcond[ilyr], layers->lyrpH[ilyr]);
          fflush(stdout);
        }
      } /* while */

      /* Set the TEXTURE parameter based on the sand content in the top 3 */
      /* layers of the soil profile.  The TEXTURE parameter is used for */
      /* computing the soil moisture effect in the nitrify and calcdefac */
      /* subroutines.  If these routines are changed to use something other */
      /* than the top 3 soil layers for computing the average wfps value */
      /* make the appropirate change in this calculation for setting the */
      /* TEXTURE parameter.  CAK - 05/31/01) */
      sitepar->avgsand = (layers->sandfrac[0]*layers->width[0] +
                          layers->sandfrac[1]*layers->width[1] +
                          layers->sandfrac[2]*layers->width[2]) /
                (layers->width[0] + layers->width[1] + layers->width[2]);
      if      (sitepar->avgsand > 0.7) {sitepar->texture = COARSE;}
      else if (sitepar->avgsand < 0.3) {sitepar->texture = FINE;}
      else {sitepar->texture = MEDIUM;}

      /* Set the aggregate soil normalized diffusivity at field capacity and */
      /* wilting point for use in the trace_gas_model subroutine cak - 07/2013 */
      avg_fc  = (layers->swcfc[0] + layers->swcfc[1] + layers->swcfc[2] ) /
                (layers->width[0] + layers->width[1] + layers->width[2]);
      avg_wp  = (layers->swcwp[0] + layers->swcwp[1] + layers->swcwp[2] ) /
                (layers->width[0] + layers->width[1] + layers->width[2]);
      avg_bd  = (layers->lyblkd[0] * layers->width[0] +
                 layers->lyblkd[1] * layers->width[1] +
                 layers->lyblkd[2] * layers->width[2]) /
                (layers->width[0] + layers->width[1] + layers->width[2]);

      wfps_fc = avg_fc / (1.0f - avg_bd/(float)PARTDENS);
      wfps_wp = avg_wp / (1.0f - avg_bd/(float)PARTDENS);
      sitepar->dDO_fc  = diffusiv(&avg_fc, &avg_bd, &wfps_fc);
      sitepar->dDO_wp  = diffusiv(&avg_fc, &avg_bd, &wfps_wp);


      if (flags->verbose) {
        printf("\nSoil Structure: \n");
        printf("lyr\tdpthmn\tdpthmx\tbulkd\tfieldc\twiltpt\tecoeff\t");
        printf("tcoeff\tsand\tclay\torg\tswclim\tsatcond\tpH\n");
        for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
          printf("%3d\t%5.1f\t%5.1f\t%4.2f\t%6.4f\t%6.4f\t",
                 ilyr, layers->dpthmn[ilyr], layers->dpthmx[ilyr],
                 layers->lyblkd[ilyr], layers->fieldc[ilyr],
                 layers->wiltpt[ilyr]);
          printf("%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%7.5f\t%4.2f\n",
                 layers->ecoeff[ilyr], layers->tcoeff[ilyr],
                 layers->sandfrac[ilyr], layers->clayfrac[ilyr],
                 layers->orgfrac[ilyr], layers->swclimit[ilyr],
                 layers->satcond[ilyr], layers->lyrpH[ilyr]);
        }

        printf("\nCalculated soil parameters: \n");
        printf("lyr\twidth\tdepth\tthetas\tpsis\tb\n");
        for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
          printf("%3d\t%5.1f\t%5.1f\t%7.3f\t%7.3f\t%6.3f\n", ilyr,
                 layers->width[ilyr], layers->depth[ilyr],
                 layers->thetas[ilyr], layers->psis[ilyr], layers->b[ilyr]);
        }
        printf("nelyrs = %4d\n", layers->nelyrs);
      }

      if (flags->debug) {
        printf("Exiting function initlyrs\n");
      }

      read = 1;
      return;
    }


    void inputlyrs(SITEPAR_SPT sitepar, LAYERPAR_SPT layers, FLAG_SPT flags)
/*****************************************************************************
**
**  FUNCTION:  void inputlyrs()
**
**  PURPOSE:   reads the soil.in file.
**
**  Modifications:    Kendrick Killian    Aug, 2012
**    This was the first section of inputlyrs
**
**  SAMPLE DATA FILE:
**
**   0.0   1.0  1.5  0.211  0.142  0.8  0.01  0.7  0.15  0.01  0.062  0.00116  6.3
**   1.0   4.0  1.5  0.211  0.142  0.2  0.12  0.7  0.15  0.01  0.062  0.00116  6.3
**   4.0  15.0  1.5  0.211  0.142  0.0  0.32  0.7  0.15  0.01  0.042  0.00116  6.3
**  15.0  30.0  1.5  0.211  0.142  0.0  0.28  0.7  0.15  0.01  0.012  0.00116  6.3
**  30.0  45.0  1.5  0.211  0.142  0.0  0.17  0.7  0.15  0.01  0.000  0.00116  6.3
**  45.0  60.0  1.5  0.211  0.142  0.0  0.06  0.7  0.15  0.01  0.000  0.00116  6.3
**  60.0  75.0  1.5  0.211  0.142  0.0  0.02  0.7  0.15  0.01  0.000  0.00116  6.3
**  75.0  90.0  1.5  0.211  0.142  0.0  0.02  0.7  0.15  0.01  0.000  0.00116  6.3
**
**  Column  1 - Minimum depth of soil layer (cm)
**  Column  2 - Maximum depth of soil layer (cm)
**  Column  3 - Bulk density of soil layer (g/cm^3)
**  Column  4 - Field capacity of soil layer, volumetric
**  Column  5 - Wilting point of soil layer, volumetric
**  Column  6 - Evaporation coefficient for soil layer
**  Column  7 - Percentage of roots in soil layer, these values must sum to
**              1.0
**  Column  8 - Fraction of sand in soil layer, 0.0 - 1.0
**  Column  9 - Fraction of clay in soil layer, 0.0 - 1.0
**  Column 10 - Organic matter in soil layer, fraction 0.0 - 1.0
**  Column 11 - Minimum volumetric soil water content below wilting point for
**              soil layer, soil water content will not be allowed to drop
**              below this value
**  Column 12 - Saturated hydraulic conductivity of soil layer in centimeters
**              per second
**  Column 13 - pH of soil layer
**
**  AUTHOR:    Kendrick Killian    March 9, 2012
*****************************************************************************/

    {
      int   ilyr;
      FILE *fp_in;
      float deltamin;

      if (flags->debug) {printf("Entering function inputlyrs\n");}

      /* skip the file read if this we this from site file*/
      if(read ==0) {

        /* initialize counters and accumulators */
        layers->nelyrs = 0;

        if((fp_in = fopen(FNSOIL, "r")) == NULL) {
          sprintf(msgbfr, "Cannot open file %s", FNSOIL);
          abortmssg(msgbfr);
        }

        ilyr=0;
        while(fscanf(fp_in, "%f %f %f %f %f %f %f %f %f %f %f %f %f",
              &layers->dpthmn[ilyr], &layers->dpthmx[ilyr],
              &layers->lyblkd[ilyr], &layers->fieldc[ilyr],
              &layers->wiltpt[ilyr], &layers->ecoeff[ilyr],
              &layers->tcoeff[ilyr], &layers->sandfrac[ilyr],
              &layers->clayfrac[ilyr], &layers->orgfrac[ilyr],
              &deltamin, &layers->satcond[ilyr],
              &layers->lyrpH[ilyr]) != EOF) {

          layers->swclimit[ilyr] = layers->wiltpt[ilyr] - deltamin;
          if (layers->swclimit[ilyr] < 0.0) {
            sprintf(msgbfr, "\nswclimit[%1d] is negative.  Check soils.in",ilyr);
            abortmssg(msgbfr);
          }

          if (ilyr >= MXSWLYR-1) {
            sprintf(msgbfr, "Number of soil layers exceeds MXSWLYR in swconst.h.");
            abortmssg(msgbfr);
          }


          ilyr++;
        } /* while */

        layers->numlyrs=ilyr;

        fclose(fp_in);
      }

      if (flags->debug) {printf("Exiting function inputlyrs\n");}

      read = 1;
      return;
    }
