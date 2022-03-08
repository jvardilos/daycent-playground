
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      pevapdly.c
**
**  FUNCTION:  float pevapdly()
**
**  PURPOSE:   Century's pevap function converted to daily timestep to
**             calculate the potential evapotranspiration rate.
**
**  HISTORY:
**    This is equation 52, chapter 3 of the report and Ra is in mm h2O/day.
**    Changed constant so it converts shwave langleys to cm water.   KLK 04/29/14
**    Calculate PET using the FAO Penman-Monteith equation, cak - 04/07/03
**    Reference:  http://www.fao.org/docrep/X0490E/x0490e08.htm
**
**  INPUTS:
**    jday    - Julian day (1-366)
**    month   - current month (1-12)
**    sitlat  - latitude (degrees)
**    tmax    - maximum air temperature for the day (deg C - 2m)
**    tmin    - minimum air temperature for the day (deg C - 2m)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    avgmth[] - average monthly temperature (degrees C)
**    e        - potential evaoptranspiration rate (mm/day)
**    elev     - elevation
**    highest  - highest average monthly temperature (degrees C)
**    ii       - loop control variable
**    lowest   - lowest average monthly temperature (degrees C)
**    ra       - temperature range between highest average month temperature
**               and lowest average month temperature
**    t        - intermediate variable for calculations
**    td       - intermediate variable for calculations
**    tm       - intermediate variable for calculations
**    tr       - temperature range for day, tmax - tmin
**
**  OUTPUTS:
**    pet - potential evapotranspiration rate (cm/day)
**
**  CALLED BY:
**    calcpet()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"
#include <math.h>

    float pevapdly(float tmin, float tmax, float sitlat, int jday, int month)
    {
      float pet, rlatitude, trange, tmean;

      float const1 = 0.0023f;
      float const2 = 17.8f;
      float langleys2cm = 585.4f; /* 58.54 langley/mm then Convert mm to cm */
      /* should be
         langley2mm = 41868*0.408/1,000,000 = 0.01708 = 1/58.540661
         equivalent evaporation [mm day-1] = 0.408 x Radiation [MJ m-2 day-1]
        1 Langley = 0.041868 MJ/m2 = 41868 J/m2/1,000,000
      */

      rlatitude = sitlat * (float)(PI/180.0);
      trange = tmax-tmin;
      tmean = (tmax+tmin)/2.0f;
      pet = const1 / langleys2cm * (tmean + const2) * (float)sqrt(trange) *
            c_shwave(month, rlatitude, jday);

      if (pet < 0.01) {
        pet = 0.01f;
      }

      return(pet);
    }
