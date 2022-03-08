/*****************************************************************************
**
**  FILE:      wrtdels.c
**
**  FUNCTION:  void wrtdels()
**
**  PURPOSE:   This function writes out daily delta 13C/14C values
**
**  INPUTS:
**    dcrtjresp   - daily growth and maintenance respiration from crop/grass
**                  juvenile fine root pool (gC/m^2)
**    dcrtmresp   - daily growth and maintenance respiration from crop/grass
**                  mature fine root pool (gC/m^2)
**    ddccarbostg - delta 13C/14C value for crop/grass carbohydrate storage
**                  pool
**    ddcgresp    - daily delta 13C/14C value for crop/grass growth
**                  respiration
**    ddcmresp    - daily delta 13C/14C value for crop/grass maintenance
**                  respiration
**    ddeloe      - daily delta 13C/14C value for heterotrophic respiration
**                  for the OE layer (surface som2c)
**    ddeloi      - daily delta 13C/14C value for heterotrophic respiration
**                  for the OI layer (surface metabolic, structural, and
**                  som1c)
**    ddfcarbostg - delta 13C/14C value for forest carbohydrate storage pool
**    ddfgresp    - daily delta 13C/14C value for forest growth respiration
**    ddfmresp    - daily delta 13C/14C value for forest maintenance
**                  respiration
**    ddhetresp   - daily delta 13/14C value for heterotrophic respiration
**    ddsmnrl     - daily delta 13C/14C value for heterotrophic respiration
**                  for the mineral soil (soil metabolic, structural, som1c,
**                  som2c, and som3c)
**    ddsoilresp  - daily delta 13/14C value for soil respiration
**                  (hetrotrophic + root autotrophic)
**    ddsrfclit   - daily delta 13C/14C value for heterotrophic respiration
**                  for the surface litter (surface metabolic, structural,
**                  som1c, and som2c)
**    dfrtcresp   - daily growth and maintenance respiration from forest
**                  coarse root pool (gC/m^2)
**    dfrtjresp   - daily growth and maintenance respiration from forest
**                  juvenile fine root pool (gC/m^2)
**    dfrtmresp   - daily growth and maintenance respiration from forest
**                  mature fine root pool (gC/m^2)
**    dgresp      - daily growth respiration (gC/m^2)
**    dhresp      - daily heterotrophic respiration (gC/m^2)
**    dsmnrlrsp   - daily heterotrophic respiration from mineral soil (gC/m^2)
**    dmresp      - daily maintenance respiration (gC/m^2)
**    doeresp     - daily heterotrophic respiration from OE layer (gC/m^2)
**    doiresp     - daily heterotrophic respiration from OI layer (gC/m^2)
**    dslitrsp    - daily heterotrophic respiration from surface litter
**                  (gC/m^2)
**    dsresp      - daily soil respiration (heterotropic + root autotrophic)
**                  (gC/m^2)
**    jday        - current day of year (1..366)
**    time        - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_dels    - file pointer to dels.out output file
**    files->write_dels - flag to indicate if dels.out output file should
**                        be created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    dailymoist()
**
*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "soilwater.h"

    void wrtdels(float *time, int *jday, float *ddeloi, float *ddeloe,
                 float *ddsrfclit, float *ddsmnrl, float *ddhetresp,
                 float *ddsoilresp, float *ddcmresp, float *ddfmresp,
                 float *ddcgresp, float *ddfgresp, float *ddccarbostg,
                 float *ddfcarbostg, float *doiresp, float *doeresp,
                 float *dslitrsp, float *dsmnrlrsp, float *dhresp,
                 float *dcrtjresp, float *dcrtmresp, float *dfrtjresp,
                 float *dfrtmresp, float *dfrtcresp, float *dsresp,
                 float *dmresp, float *dgresp)
    {
      extern FILES_SPT files;

      if (!files->write_dels) {
        goto ex;
      }

      fprintf(files->fp_dels, "%8.2f  %4d  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *time, *jday, *ddeloi, *ddeloe, *ddsrfclit, *ddsmnrl);
      fprintf(files->fp_dels, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *ddhetresp, *ddsoilresp, *ddcmresp, *ddfmresp, *ddcgresp);
      fprintf(files->fp_dels, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *ddfgresp, *ddccarbostg, *ddfcarbostg, *doiresp, *doeresp);
      fprintf(files->fp_dels, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *dslitrsp, *dsmnrlrsp, *dhresp, *dcrtjresp, *dcrtmresp);
      fprintf(files->fp_dels, "%10.4f  %10.4f  %10.4f  %10.4f  ",
              *dfrtjresp, *dfrtmresp, *dfrtcresp, *dsresp);
      fprintf(files->fp_dels, "%10.4f  %10.4f\n", *dmresp, *dgresp);

ex:   return;
    }
