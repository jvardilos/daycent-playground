
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtgrsp.c
**
**  FUNCTION:  void wrtgrsp()
**
**  PURPOSE:   Write out the growth respiration values. 
**
**  AUTHOR:    Cindy Keough 01/2007
** 
**  INPUTS:
**    carbostg11  - unlabeled C in carbohydrate storage for grass/crop system
**                  (gC/m^2)
**    carbostg12  - labeled C in carbohydrate storage for grass/crop system
**                  (gC/m^2)
**    carbostg21  - unlabeled C in carbohydrate storage for forest system
**                  (gC/m^2)
**    carbostg22  - labeled C in carbohydrate storage for forest system
**                  (gC/m^2)
**    cgrspflux1  - amount of daily growth respiration flux from aboveground
**                  grass/crop material that is blown off into the atmosphere
**                  during plant carbon production (gC/m^2)
**    cgrspflux2  - amount of daily growth respiration flux from juvenile
**                  belowground grass/crop material that is blown off into the
**                  atmosphere during plant carbon production (gC/m^2)
**    cgrspflux3  - amount of daily growth respiration flux from mature
**                  belowground grass/crop material that is blown off into the
**                  atmosphere during plant carbon production (gC/m^2)
**    curday      - the day of the year (1..366)
**    fgrspflux1  - amount of daily growth respiration loss from live leaf
**                  material that is blown off into the atmosphere during
**                  plant carbon production (gC/m^2)
**    fgrspflux2  - amount of daily growth respiration loss from live
**                  juvenile fine root material that is blown off into the
**                  atmosphere during plant carbon production (gC/m^2)
**    fgrspflux6  - amount of daily growth respiration loss from live
**                  mature fine root material that is blown off into the
**                  atmosphere during plant carbon production (gC/m^2)
**    fgrspflux3  - amount of daily growth respiration loss from live fine
**                  branch material that is blown off into the atmosphere
**                  during plant carbon production (gC/m^2)
**    fgrspflux4  - amount of daily growth respiration loss from live large
**                  wood material that is blown off into the atmosphere
**                  during plant carbon production (gC/m^2)
**    fgrspflux5  - amount of daily growth respiration loss from live coarse
**                  root material that is blown off into the atmosphere
**                  during plant carbon production (gC/m^2)
**    grspann1    - accumulator for annual growth respiration for grass/crop
**                  (gC/m^2)
**    grspann2    - accumulator for annual growth respiration for tree
**                  (gC/m^2)
**    grspflux1   - daily growth respiration flow from storage pool 
**                  (carbostg(1,*) to C source/sink for grass/crop system
**                  (gC/m^2)
**    grspflux2   - daily growth respiration flow from storage pool
**                  (carbostg(2,*) to C source/sink for tree system (gC/m^2)
**    mcprd1      - daily NPP for shoots for grass/crop system (gC/m^2)
**    mcprd2      - daily NPP for juvenile roots for grass/crop system (gC/m^2)
**    mcprd3      - daily NPP for mature roots for grass/crop system (gC/m^2)
**    mfprd1      - daily NPP for live leaves for tree system (gC/m^2)
**    mfprd2      - daily NPP for live juvenile fine roots for tree system
**                  (gC/m^2)
**    mfprd6      - daily NPP for live mature fine roots for tree system
**                  (gC/m^2)
**    mfprd3      - daily NPP for live fine branches for tree system (gC/m^2)
**    mfprd4      - daily NPP for live large wood for tree system (gC/m^2)
**    mfprd5      - daily NPP for live coarse roots for tree system (gC/m^2)
**    time        - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files              - structure containing information about output files
**    files->fp_gresp    - file pointer to gresp.out output file
**    files->write_gresp - flag to indicate if gresp.out output file should be
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

    void wrtgresp(float *time, int *curday, float *grspflux1, float *grspflux2,
                  float *cgrspflux1, float *cgrspflux2, float *cgrspflux3,
                  float *fgrspflux1, float *fgrspflux2, float *fgrspflux6,
                  float *fgrspflux3, float *fgrspflux4, float *fgrspflux5,
                  float *mcprd1, float *mcprd2, float *mcprd3, float *mfprd1,
                  float *mfprd2, float *mfprd6, float *mfprd3, float *mfprd4,
                  float *mfprd5, float *carbostg11, float *carbostg12,
                  float *carbostg21, float *carbostg22, float *grspann1,
                  float *grspann2)
    {
      extern FILES_SPT files;

      if (!files->write_gresp) {
        return;
      }

      fprintf(files->fp_gresp, "%6.2f  %2d  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *time, *curday, *grspflux1, *grspflux2, *cgrspflux1,
              *cgrspflux2);
      fprintf(files->fp_gresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *cgrspflux3, *fgrspflux1, *fgrspflux2, *fgrspflux6,
              *fgrspflux3);
      fprintf(files->fp_gresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *fgrspflux4, *fgrspflux5, *mcprd1, *mcprd2, *mcprd3);
      fprintf(files->fp_gresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *mfprd1, *mfprd2, *mfprd6, *mfprd3, *mfprd4);
      fprintf(files->fp_gresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *mfprd5, *carbostg11, *carbostg12, *carbostg21, *carbostg22);
      fprintf(files->fp_gresp, "%12.4f  %12.4f\n", *grspann1, *grspann2);

      return;
    }
