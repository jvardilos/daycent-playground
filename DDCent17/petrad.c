
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      petrad.c
**
**  FUNCTION:  float petrad()
**
**  PURPOSE:   Calculate the potential evapotranspiration rate using
**             Penmans equation (1948).
**
**  REWRITE:  Melannie Hartman  9/21/93 - 10/6/93
**
**  NOTES:    Mathematical function substitution:
**              atan2(y,x) is unchanged - arctan(y/x) in the range [-PI, PI]
**
**  HISTORY:
**    2/03/2104 KLK simplified and commented the UNUSED the hour angle calculation
**                  removed an if that always uses plant albedo for Jan 1
**                  removed unused formal parameters, jday, month, latitude
**                  changed cldcov    - average cloud cover for the month (%, 1..100)
**                     to clrsky an alias for the sradadj input
**    4/30/92  (SLC)
**
**  INPUTS:
**    albedo    - fraction of light reflected by snow
**    avgtemp   - average air temperature for the day (deg C)
**    clrsky      - relative amount of clear sky (decimal fraction)
**    jday      - current julian day (1..366)
**    month     - current month of the year (1..12)
**    reflec    - fraction of light reflected by vegetation
**    rhumid    - average relative humidity for the day (% 1..100)
**    rlatitude - latitude of the site (in radians)
**    snowpack  - current snowpack (equiv. cm H2O)
**    solrad    - total incoming shortwave radiation (langleys/day)
**    windsp    - average daily windspeed at 2 meters (mph)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    arads       - intermediate variable for calculations
**    fhumid      - vapor pressure (mm Hg)
**    ftemp       - intermediate variable for calculations
**    kelvin      - avgtemp in kelvin degrees
**    KELVIN_CONV - conversion factor for converting degrees C to degrees K
**    netreflec   - net reflectivity (use either reflec or albedo)
**    par1, par2  - parameters in computation of pet
**    shwave      - short wave solar radiation (mm H2O)
**
**  OUTPUTS:
**    pet - potential evapotranspiration rate (cm/day)
**
**  CALLED BY:
**    calcpet()
**
**  CALLS:
**    svapor() - calculate the saturation vapor pressure of water for air
**               temperature
**
*****************************************************************************/

#include <math.h>
#include "soilwater.h"

#define KELVIN_CONV 273

    float petrad(int jday, float rlatitude,
                 float windsp, float rhumid, float avgtemp, float clrsky, 
                 float reflec, float solrad, float albedo, float snowpack)
    {
      float pet;      /* return value */
      /* float declin, ahou; */
      float shwave,kelvin,arads;
      float fhumid,ftemp;
      float par1, par2;
      float netreflec;

      /* calculate the short wave solar radiation */
      /* on a clear day using a equation presented by sellers(1965) */

      /* if there is snow on the ground, use snow albedo as */
      /* reflectivity, else use vegetation albedo.  Suggested */
      /* by Bill Parton 9/15/94 */

      netreflec = reflec;
      /*if (jday > 1) { no reason day 1 reflectivity should ignore snow KLK*/
        if (snowpack > 0) {
          netreflec = albedo;
        }


      /* hour angle calculations  Commented since ahou is never used    KLK
         Original FORTRAN statement for par1:
         par1=sqrt(1.-(-tan(rlatitude)*tan(declin))**2)

        declin = 0.401426f*(float)sin(2*PI*(jday-77)/365);
        ahou = (-tan((double)rlatitude)*tan(declin));
        ahou = atan2(sqrt(1-pow(ahou,2)),ahou);
        if(ahou < 0) {ahou = 0.0;}
      */

      /* shwave = mm of water equivalent for total incoming radiation */
      /* 0.0168 converts langleys/day to mm H2O (-mdh 9/94) */
      /* 10*(1/597) = 0.0168: mult by 10 to get mm, divide by latent heat */
      /* of vaporization (cal/cm3) at 0 deg C */

      shwave = solrad*0.0168f;

      kelvin = avgtemp + KELVIN_CONV;
      arads = svapor(avgtemp) * (3010.21f)/(kelvin*kelvin);

   /* clrsky = 1 - cldcov/100; */
      fhumid = rhumid * svapor(avgtemp)/100;

      ftemp = (avgtemp + KELVIN_CONV)*0.01f;
      ftemp = (float)pow((double)ftemp, 4) * (0.201f);

      par1 = (0.35f)*(svapor(avgtemp)-fhumid) * (1+0.0098f*(windsp*24));

/*      par2 = shwave*(1-reflec)*(0.18 + 0.55*clrsky) -
             ftemp*(0.56 - 0.092*sqrt(fhumid))*(0.10 + 0.90*clrsky); */

      par2=shwave*(1-netreflec)-ftemp*(0.56f-0.092f*
           (float)sqrt((double)fhumid))*(0.10f+0.90f*clrsky);

      pet = ((arads*par2 + 0.27f*par1)/(arads + 0.27f))/10;

      if(pet < 0.03) {
        pet = 0.03f;
      }

      return(pet);
    }
