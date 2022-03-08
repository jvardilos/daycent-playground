
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:     calcSrad.h
**
**  PURPOSE:  Header file for calcSrad.c.
**
**  AUTHOR:   Adapted from Peter Thorton's code extracted from the Sipnet
**            model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**            Cindy Keough - 04/2009
**
*****************************************************************************/

/***********************
**                    **
**  GLOBAL CONSTANTS  **
**                    **
***********************/
#define PI        3.14159265358979323846264338327950288 /* pi */
#define RADPERDAY 0.01720242383895848453641420059290625 /* radians of Earth orbit per julian day = 2*PI/365.25*/
#define RADPERDEG 0.01745329251994329576923690768488613 /* radians per degree */
#define SECPERRAD 13750.98708313975701043155715538524   /* seconds per radian of hour angle */

#define DNSS       80.06        /* day number of spring solstice */
#define EOEcc      0.0334       /* earth orbit eccentricity */

#define ABASE      -6.1e-5      /* vapor pressure effect on transmittance  (1/Pa) */
#define B0          0.013       /* radiation parameter (dim) */
#define B1          0.201       /* radiation parameter (dim) */
#define B2          0.185       /* radiation parameter (dim) */
#define C           1.5         /* radiation parameter (dim) */
#define CP       1010.0         /* specific heat of air (J kg-1 K-1) */
#define DIF_ALB     0.6         /* diffuse albedo for horizon correction  (dim) */
#define EPS         0.62196351  /* unitless ratio of molec weights (MW/MA) */
#define G_STD       9.80665     /* standard gravitational accel. (m s-2) */
#define LR_STD      0.0065      /* standard temperature lapse rate (-K m-1) */
#define MA         28.9644e-3   /* molecular weight of air (kg mol-1) */
#define P_STD  101325.0         /* standard pressure at 0.0 m elevation (Pa) */
#define R           8.3144598   /* (48) gas law constant (m3 Pa mol-1 K-1) */
#define RAIN_SCALAR 0.75        /* correction to trans. for rain day (dim) */
#define SNOW_TCRIT -6.0         /* critical temperature for snowmelt (deg C) */
#define SNOW_TRATE  0.042       /* snowmelt rate (cm/degC/day) */
#define SRADDT    600.0         /* timestep for radiation routine (seconds) */
#define T_STD     288.15        /* standard temp at 0.0 m elevation (K) */
#define TBASE       0.870       /* max inst. trans., 0m, nadir, dry atm (dim) */
#define TDAYCOEF    0.45        /* daylight air temperature coefficient (dim) */


      /*solar constant 4.914 MJ m-2/hr for a day; (from NASA mynasadata.larc.nasa.gov/glossary/solar-constant-2/)
          to convert to Langleys/day /(4.184 J/cal/ * 10000 cm2/m2 *10^-6 MJ/J);
          MJ /m/hr (12h/Pi rad) * (0.01 m2/cm2 J/MJ * 4.184 cal/j) * 2 (integration constant)
      */
      #define ScW 1365.0 /* solar constant Watts */
      static double  WtL = 0.04184e6; /* Watt/Langley */
      static double  SpRad = 3600 * (12/PI); /* S/Rad; */

      /* according to Jean Meeus: Astronomical Algorithms
          eps = the obliquity of the ecliptic is:
              23 deg 26' 21.448" -46.815" T -0.00059"T*T -0.001813"T*T*T
              23.4393 - 0.0130042 T - 1.66667*10^-7 T^2 - 5.02778*10^-7 T^3
               where T=(JD - 2451545.0)/36525 the number of Centuries since Jan 1 2000
          Forget the procession and use   23.43929  radians = 0.4090927848297817
      */
      static double eps = PI/180 *(23. + (26. + 21.448/60)/60);


/****************************
**                         **
**  STRUCTURE DEFINITIONS  **
**                         **
****************************/
typedef struct {
  double site_lat;       /* site latitude, dec. degrees (- for south) */
  double site_elev;      /* site elevation, meters */
  double site_slp;       /* site slope, degrees */
  double site_asp;       /* site aspect, degrees */
  double site_ehoriz;    /* site east horizon, degrees */
  double site_whoriz;    /* site west horizon, degrees */
} PARAMETER_S, *PARAMETER_SPT;

typedef struct {
  double s_tdew;         /* site dewpoint temperature value (degrees C) */
  double s_tmax;         /* site tmax value */
  double s_tmin;         /* site tmin value */
  double s_tday;         /* site daylight temperature value */
  double s_prcp;         /* site prcp value */
  double s_srad;         /* site shortwave radiation value */
  double s_dayl;         /* site daylength value */
  double s_swe;          /* site snow water equivalent value (cm) */
} DATA_S, *DATA_SPT;

typedef struct {
  double ttmax0[366];
  double flat_potrad[366];
  double slope_potrad[366];
  double daylength[366];
} YEARDATA_S, *YEARDATA_SPT;

typedef struct {
  int indewpt;           /* input dewpoint temperature flag (0=NO, 1=YES) */
} CONTROL_S, *CONTROL_SPT;


/***************************
**                        **
**  FUNCTION DEFINITIONS  **
**                        **
***************************/
double atm_pres(double elev);

double calc_pet(double rad, double ta, double pa, double dayl);

void calc_srad_humidity(CONTROL_SPT ctrl, PARAMETER_SPT params,
                        DATA_SPT sitedata, YEARDATA_SPT yeardata, int yday);

void calc_srad_humidity_iterative(CONTROL_SPT ctrl, PARAMETER_SPT params,
                                  DATA_SPT sitedata, YEARDATA_SPT yeardata,
                                  int yday);

void snowpack(DATA_SPT sitedata);
