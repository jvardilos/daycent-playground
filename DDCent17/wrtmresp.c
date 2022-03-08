
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtmresp.c
**
**  FUNCTION:  void wrtmresp()
**
**  PURPOSE:   Write out the maintenance respiration values. 
**
**  AUTHOR:    Cindy Keough 02/02
** 
**  INPUTS:
**    carbostg11       - unlabeled C in carbohydrate storage for grass/crop
**                       system (gC/m^2)
**    carbostg12       - labeled C in carbohydrate storage for grass/crop
**                       system (gC/m^2)
**    carbostg21       - unlabeled C in carbohydrate storage for forest system
**                       (gC/m^2)
**    carbostg22       - labeled C in carbohydrate storage for forest system
**                       (gC/m^2)
**    cmrspflux1       - amount of daily maintenance respiration flux from
**                       aboveground grass/crop material that flows from the
**                       grass/crop carbohydrate storage pool
**                       (carbostg(1,*)) to the C source/sink pool (csrsnk)
**                       (gC/m^2)
**    cmrspflux2       - amount of daily maintenance respiration flux from
**                       juvenile belowground grass/crop material that flows
**                       from the grass/crop carbohydrate storage pool
**                       (carbostg(1,*)) to the C source/sink pool (csrsnk)
**                       (gC/m^2)
**    cmrspflux3       - amount of daily maintenance respiration flux from
**                       mature belowground grass/crop material that flows
**                       from the grass/crop carbohydrate storage pool
**                       (carbostg(1,*)) to the C source/sink pool (csrsnk)
**                       (gC/m^2)
**    curday           - the day of the year (1..366)
**    fmrspflux1       - amount of daily maintenance respiration flux from
**                       live leaf material that flows from the tree
**                       carbohydrate storage pool (carbostg(2,*)) to the C
**                       source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux2       - amount of daily maintenance respiration flux from
**                       live juvenile fine root material that flows from the
**                       tree carbohydrate storage pool (carbostg(2,*)) to
**                       the C source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux6       - amount of daily maintenance respiration flux from
**                       live mature fine root material that flows from the
**                       tree carbohydrate storage pool (carbostg(2,*)) to
**                       the C source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux3       - amount of daily maintenance respiration flux from
**                       live fine branch material that flows from the tree
**                       carbohydrate storage pool (carbostg(2,*)) to
**                       the C source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux4       - amount of daily maintenance respiration flux from
**                       live large wood material that flows from the tree
**                       carbohydrate storage pool (carbostg(2,*)) to
**                       the C source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux5       - amount of daily maintenance respiration flux from
**                       live coarse root material that flows from the tree
**                       carbohydrate storage pool (carbostg(2,*)) to
**                       the C source/sink pool (csrsnk) (gC/m^2)
**    mcprd1           - daily NPP for shoots for grass/crop system (gC/m^2)
**    mcprd2           - daily NPP for juvenile roots for grass/crop system
**                       (gC/m^2)
**    mcprd3           - daily NPP for mature roots for grass/crop system
**                       (gC/m^2)
**    mfprd1           - daily NPP for live leaves for tree system (gC/m^2)
**    mfprd2           - daily NPP for live juvenile fine roots for tree
**                       system (gC/m^2)
**    mfprd6           - daily NPP for live mature fine roots for tree system
**                       (gC/m^2)
**    mfprd3           - daily NPP for live fine branches for tree system
**                       (gC/m^2)
**    mfprd4           - daily NPP for live large wood for tree system
**                       (gC/m^2)
**    mfprd5           - daily NPP for live coarse roots for tree system
**                       (gC/m^2)
**    mrspann1         - accumulator for annual maintenance respiration for
**                       grass/crop (gC/m^2)
**    mrspann2         - accumulator for annual maintenance respiration for
**                       tree (gC/m^2)
**    mrspflux1        - daily maintenance respiration flow from storage pool
**                       (carbostg(1,*) to C source/sink for grass/crop
**                       system (gC/m^2)
**    mrspflux2        - daily maintenance respiration flow from storage pool
**                       (carbostg(2,*) to C source/sink for tree system
**                       (gC/m^2)
**    mrspTempEffect11 - temperature effect on maintenance respiration for
**                       aboveground crop/grass components
**    mrspTempEffect12 - temperature effect on maintenance respiration for
**                       belowground crop/grass components
**    mrspTempEffect21 - temperature effect on maintenance respiration for
**                       leaves, fine branch, and large wood forest components
**    mrspTempEffect22 - temperature effect on maintenance respiration for
**                       fine root and coarse root forest components
**    mrspWaterEffect1 - water effect on maintenance respiration for
**                       crop/grass system
**    mrspWaterEffect2 - water effect on maintenance respiration for forest
**                       system
**    tavedly          - mean air temperature over production period (deg C)
**    time             - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files              - structure containing information about output files
**    files->fp_mresp    - file pointer to mresp.out output file
**    files->write_mresp - flag to indicate if mresp.out output file should be
**                         created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**     None
**
**  CALLED BY:
**     simsom()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void wrtmresp(float *time, int *curday, float *mrspflux1, float *mrspflux2,
                  float *cmrspflux1, float *cmrspflux2, float *cmrspflux3,
                  float *fmrspflux1, float *fmrspflux2, float *fmrspflux6,
                  float *fmrspflux3, float *fmrspflux4, float *fmrspflux5,
                  float *mcprd1, float *mcprd2, float *mcprd3, float *mfprd1,
                  float *mfprd2, float *mfprd6, float *mfprd3, float *mfprd4,
                  float *mfprd5, float *carbostg11, float *carbostg12,
                  float *carbostg21, float *carbostg22, float *mrspann1,
                  float *mrspann2, float *tavedly, float *mrspTempEffect11,
                  float *mrspTempEffect12, float *mrspWaterEffect1,
                  float *mrspTempEffect21, float *mrspTempEffect22,
                  float *mrspWaterEffect2)
    {
      extern FILES_SPT files;

      if (!files->write_mresp) {
        return;
      }

      fprintf(files->fp_mresp, "%6.2f  %2d  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *time, *curday, *mrspflux1, *mrspflux2, *cmrspflux1,
              *cmrspflux2);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *cmrspflux3, *fmrspflux1, *fmrspflux2, *fmrspflux6,
              *fmrspflux3);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ",
             *fmrspflux4, *fmrspflux5, *mcprd1, *mcprd2, *mcprd3);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ", 
              *mfprd1, *mfprd2, *mfprd6, *mfprd3, *mfprd4);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *mfprd5, *carbostg11, *carbostg12, *carbostg21, *carbostg22);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  ", *mrspann1,
              *mrspann2, *tavedly);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  ",
              *mrspTempEffect11, *mrspTempEffect12, *mrspWaterEffect1);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f\n",
              *mrspTempEffect21, *mrspTempEffect22, *mrspWaterEffect2);
      return;
    }
