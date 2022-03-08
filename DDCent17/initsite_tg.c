
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsite_tg.c
**
**  FUNCTION:  void initsite()
**
**  PURPOSE:   Read in site specific parameters (from sitepar.in)
**
**  AUTHOR:    Susan Chaffee    March 10, 1992
**
**  HISTORY:
**    24/03/14 KLK Rewrite
**             require all lines to include the variable name.
**              - input values are assigned by variable name NOT input order
**              - file input is order independent.
**              - all unlabeled ot unrecognized lines are echoed with a warning.
**             Common practice was to label the values anyway. This allows legacy
**             files to be parsed with error checking and without assigning values
**             to the wrong location.
**             Put the entire file in a read loop.
**              1 Read inputs the character record.
**              2 scan for a real value. Error report if nothing parses.
**              3 If then else block searches for EACH variable name and
**                assigns values to permanent location, rescanning for integer
**                ore multiple reads as required.
**             NEW ROUTINE defaultsitpar to set sitepar value defaults.
**                required to set initial values. Read routine does not check
**                that all values input.
**             REMOVED methane parameters. Moved to FIX.100
**
**    24/03/14 KLK
**             eliminated cldcov    - average cloud cover for the month (%, 1..100)
**              in preference to sradadj input
**    REWRITE:   Melannie Hartman  9/9/93 - 9/23/93
**    8/13/92 (SLC) - Change the way the parameter for minimum soil water
**                    content is used.  No longer a function of wilting point,
**                    now it simply gives the percent of water at each layer.
**    5/08/08 (CAK) - Add checks for end of file, if the end of the file is
**                    reached before all of the variables have been read and
**                    initialized this will result in a fatal error.
**   11/21/13 (KLK) - For consistency, use abortmssg for the missing variable aborts.
**   12/13 (KLK)    - Added methanogenesis variables,.
**
**  INPUTS:
**    flags          - structure containing debugging flags
**    flags->debug   - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose - flag to set verbose debugging mode, 0 = off, 1 = on
**    sitename       - data file name containing site parameters
**    sitepar        - site specific parameters structure for soil water model
**    layers         - soil water soil layer structure
**
**  GLOBAL VARIABLES:
**    INPTSTRLEN - maximum length of input file line (120)
**    NTDEPTHS   - maximum number of soil regions (4)
**
**  LOCAL VARIABLES:
**    errmsg[]   - string containing error message
**    fp_in      - pointer to input file
**    imo        - current month (1..12)
**    inptline[] - line read from input file
**    m          - month
**    C6H12O6_to_CH4 - converted to frCH4emit KLK 10Mar2015
**
**  OUTPUTS:
**    Read from file sitename (changed to read order):
**    sitepar->usexdrvrs     - -1 = determine using weather file format
**                              0 = use air temperature to drive PET rates
**                              1 = use extra weather drivers (solrad, rhumid,
**                                  windsp) for PET calculation
**                              2 = read srad from weather data file
**                              3 = use extra drivers for both PET and
**                                  read srad from weather data file
**                              4 = use air temperature and read EVI
**    sitepar->sublimscale   - multiplier to scale sublimation
**    sitepar->reflec        - fraction of light reflected by vegetation
**    sitepar->albedo        - fraction of light reflected by snow
**    sitepar->dmpflux       - damping factor for soil water flux (in h2oflux)
**    sitepar->hours_rain    - the duration of the rainfall/snowmelt event (hours)
**    sitepar->drainlag      - number of days that soil drainage should lag
**                             behind rain/irrigation/melt events (1-5)
**    sitepar->hpotdeep      - hydraulic water potential of deep storage layer, the more
**                             negative the number the dryer the soil layer (units?)
**    sitepar->ksatdeep      - saturated hydraulic conductivity of deep storage (cm/sec)
**    sitepar->dmp           - damping factor for calculating soil temperature by layer
**    sitepar->Ncoeff        - minimum water/temperature limitation coefficient
**                             for nitrification
**    sitepar->jdayStart     - the Julian day to start the turning off of the
**                             restriction of the CO2 effect on denitrification
**    sitepar->jdayEnd       - the Julian day to end the turning off of the
**                             restriction of the CO2 effect on denitrification
**    sitepar->N2Oadjust_fc  - maximum proportion of nitrified N lost as N2O, occurs
**                             when soil water content is at field capacity  8/2013
**    sitepar->N2Oadjust_wp  - minimum proportion of nitrified N lost as N2O, occurs
**                             when soil water content is at wilting point  8/2013
**    sitepar->MaxNitAmt     - maximum daily nitrification amount (gN/m^2)        8/2013
**    sitepar->SnowFlag      - effect of snow on soil surface temperature calc
**                             1 = insulating effect, 0 = no insulating effect
**    sitepar->netmn_to_no3  - fraction of new net mineralization going to NO3 (0.0-1.0)
**    sitepar->wfpsdnitadj   - adjustment on inflection point for water filled pore space
**                             effect on denitrification curve    8/2013
**    sitepar->N2N2Oadj      - multiplier on N2/N2O ratio                         8/2013
**    sitepar->CO2_to_CH4    - fraction of CO2 from heterotrophic soil respiration that
**                             is used to produce CH4 (0.0 - 1.0)
**    sitepar->frac_to_exudates - fraction of fine root production that is root exudates
**    sitepar->Aeh           - differential coefficient for calculation of soil Eh
**    sitepar->Deh           - differential coefficient for calculation of soil Eh
**    sitepar->Beh_flood     - low-limit value for Eh during flooding course (mv)
**    sitepar->Beh_drain     - upper-limit value of Eh during drainage course (mv)
**    sitepar->elevation     - site elevation (meters)
**    sitepar->sitslp        - site slope (degrees)
**    sitepar->aspect        - site aspect (degrees)
**    sitepar->ehoriz        - site east horizon, degrees
**    sitepar->whoriz        - site west horizon, degrees
**    sitepar->sradadj[]     - solar radiation adjustment for cloud cover and
**                             transmission coefficient
**
**  EXAMPLE INPUT FILE:
**  0        / usexdrvrs      - 0 = no extra drivers,
**            1 = use extra weather drivers for PET (solrad, rhumid, windsp)
**            2 = read solar radiation
**            3 = use extra drivers for PET and read solar radiation
**            4 - EVI
**  1.0      / sublimscale
**  0.18     / reflec         - vegetation reflectivity (frac) hardwoods = 0.20 spruce = 0.10
**  0.65     / albedo         - snow albedo (frac)
**  0.000008 / dmpflux        - in h2oflux routine (0.000001 = original value)
**  10       / hours_rain     - duration of each rain event
**  0        / drainlag       - # of days between rainfall event and drainage of soil (-1=computed)
**  -200     / hpotdeep       - hydraulic water potential of deep storage layer (units?)
**  0.0003   / ksatdeep       - saturated hydraulic conductivity of deep storage layer (cm/sec)
**  0.003    / dmp               - damping factor for calculating soil temperature by layer
**  0.03     / Ncoeff            - min water/temperature limitation coefficient for nitrify
**  0 0      / jdayStart jdayEnd - turn off respiration restraint on denit between these Julian dates
**  0.01     / N2Oadjust_fc      - maximum proportion of nitrified N lost as N2O @ field capacity
**  0.002    / N2Oadjust_wp      - minimum proportion of nitrified N lost as N2O @ wilting point
**  0.4      / MaxNitAmt         - maximum daily nitrification amount (gN/m^2)
**  1        / SnowFlag          - snow insulation effect on soil surface temp: 0 = not insulating, 1 = insulating
**  0.2      / netmn_to_no3      - fraction of new net mineralization that goes to NO3 (0.0-1.0)
**  1.0      / wfpsdnitadj       - adjustment on inflection point for WFPS effect on denit
**  1.0      / N2N2Oadj          - N2/N2O ratio adjustment coefficient
**  0.5      / CO2_to_CH4        - fraction of CO2 from soil respiration used to produce CH4
**  0.0      / frCH4emit         - fraction of methane emmitted
**  0.45     / frac_to_exudates  - fraction of root production that is root exudates
**  0.23     / Aeh               - differential coefficient (Aeh)
**  0.16     / Deh               - differential coefficient (Deh)
**  -250.0   / Beh_flood         - low-limit value for Eh during flooding course (mv)
**  300.0    / Beh_drain         - upper-limit value of Eh during drainage course (mv)
**  0.7      / zero_root_frac    - fraction CH4 emitted via bubbles when zero root biomass (0.0-1.0)
**  500.0    / elevation         - elevation, meters Brookings SDE
**  0.0      / slope             - site slope, degrees
**  0.0      / aspect            - site aspect, degrees
**  0.0      / ehoriz            - site east horizon, degrees
**  0.0      / whoriz            - site west horizon, degrees
**  1  0.42  / solar radiation adjustment for cloud cover & transmission coeffient
**  2  0.50  / solar radiation adjustment for cloud cover & transmission coeffient
**  3  0.53  / solar radiation adjustment for cloud cover & transmission coeffient
**  4  0.57  / solar radiation adjustment for cloud cover & transmission coeffient
**  5  0.62  / solar radiation adjustment for cloud cover & transmission coeffient
**  6  0.69  / solar radiation adjustment for cloud cover & transmission coeffient
**  7  0.71  / solar radiation adjustment for cloud cover & transmission coeffient
**  8  0.66  / solar radiation adjustment for cloud cover & transmission coeffient
**  9  0.58  / solar radiation adjustment for cloud cover & transmission coeffient
**  10  0.52 / solar radiation adjustment for cloud cover & transmission coeffient
**  11  0.46 / solar radiation adjustment for cloud cover & transmission coeffient
**  12  0.45 / solar radiation adjustment for cloud cover & transmission coeffient
**
**  CALLED BY:
**    detiv()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"
#include "mssg.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int first=1;

    void initsite(SITEPAR_SPT sitepar, LAYERPAR_SPT layers, FLAG_SPT flags)
    {

      int  imo, m, fcnt;
      char inptline[INPTSTRLEN];
      char errmsg[INPTSTRLEN];
      float fval;
      FILE *fp_in;

      if (flags->debug) {
        printf("Entering function initsite\n");fflush(stdout);
      }
      if(first !=1) {return;}

      if ((fp_in = fopen(FNSITE, "r")) == NULL) {
        sprintf(errmsg, "Cannot open file %s\n", FNSITE);
        perror(errmsg);
        abortmssg(errmsg);
      }

      while (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        fcnt = sscanf(inptline, "%f", &fval);
        if(fcnt != 1) {
          printf("WARNING: unable to parse sitepar.in: %s",inptline);
          continue;
        }

        if(strstr(inptline,"usexdrvrs") != NULL) {
          sscanf(inptline, "%d", &sitepar->usexdrvrs);
          if (flags->verbose) {printf("sitepar; usexdrvrs: %1d\n", sitepar->usexdrvrs);}
          if (sitepar->usexdrvrs < -1 || sitepar->usexdrvrs > 4) {
            abortmssg("ERROR:  Invalid usexdrvrs in sitepar.in");
          }
        }

        else if(strstr(inptline,"sublimscale") != NULL) {
          sitepar->sublimscale = fval;
          if (flags->verbose) {
            printf("sitepar; sublimation scalar sublimscale: %f\n", sitepar->sublimscale);
          }
        }

        else if(strstr(inptline,"reflec") != NULL) {
          sitepar->reflec = fval;
          if (flags->verbose) {
            printf("sitepar; vegetation reflectivity, reflec: %f\n", sitepar->reflec);
          }
        }

        else if(strstr(inptline,"albedo") != NULL) {
          sitepar->albedo = fval;
          if (flags->verbose) {
            printf("sitepar; snow albedo, albedo: %f\n", sitepar->albedo);
          }
        }

        else if(strstr(inptline,"dmpflux") != NULL) {
          sitepar->dmpflux = fval;
          if (flags->verbose) {
            printf("sitepar; soil water flux damping factor, dmpflux: %f\n", sitepar->dmpflux);
          }
        }

        else if(strstr(inptline,"hours_rain") != NULL) {
          sitepar->hours_rain = fval;
          if (flags->verbose) {
            printf("sitepar; rainfall duration, hours_rain: %f\n", sitepar->hours_rain);
          }
        }

        /* Allow user to set number of days between rainfall event and */
        /* drainage of soil profile.  If a value of -1 is entered set the */
        /* number of days to drainage based in the soil texture.  Constrain */
        /* the number of days to drainage to be <=5 to prevent numerical */
        /* instabilities in the h2oflux subroutine.  cak - 02/13/04 */
        else if(strstr(inptline,"drainlag") != NULL) {
          sscanf(inptline, "%d", &sitepar->drainlag);
          if (sitepar->drainlag < 0) {
            sitepar->drainlag = sitepar->texture - 1;
          }
          if (sitepar->drainlag > 5) {
            printf("lag period for drainage too long, setting to max value\n");
            sitepar->drainlag = 5;
          }
          if (flags->verbose) {
            printf("drainlag: %d\n", sitepar->drainlag);
          }
        }

        else if(strstr(inptline,"hpotdeep") != NULL) {
          sitepar->hpotdeep = fval;
          if (flags->verbose) {
            printf("sitepar; deep storage hydraulic potential, hpotdeep: %f\n", sitepar->hpotdeep);
          }
        }

        else if(strstr(inptline,"ksatdeep") != NULL) {
          sitepar->ksatdeep = fval;
          if (flags->verbose) {
            printf("sitepar; deep storage saturated hydraulic conductivity, ksatdeep: %f\n", sitepar->ksatdeep);
          }
        }

        /* removed sine temperature boundary condition klk 2016
           The texture parameter is being set in the initlyrs subroutine  CAK - 05/31/01
           texture parameter replaced by the minimum and maximum temperatures for bottom soil layer
        */
        else if(strstr(inptline,"tbotmn") != NULL) {
          printf("sitepar NOTE: ignoring unused soil boundary temperature values: tbotmn, tbotmx\n");
        }

        /* Added dmp parameter to the sitepar.in file, cak - 12/16/02 */
        else if(strstr(inptline,"dmp") != NULL) {
          sitepar->dmp = fval;
          if (flags->verbose) {
            printf("sitepar; soil temperature damping factor, dmp: %f\n", sitepar->dmp);
          }
        }

        /* Added timlag parameter to the sitepar.in file, cak - 04/24/03 */
        else if(strstr(inptline,"timlag") != NULL) {
          printf("sitepar NOTE: ignoring unused soil boundary temperature value: timlag\n");
        }

        /* Added Ncoeff parameter to the sitepar.in file, cak - 04/08/03 */
        else if(strstr(inptline,"Ncoeff") != NULL) {
          if(fval > 1.0f) {fval = 1.0f;} else if (fval < 0.0) {fval = 0.0f;}
          sitepar->Ncoeff = fval;
          if (flags->verbose) {
            printf("sitepar; coefficient for nitrification, Ncoeff: %f\n", sitepar->Ncoeff);
          }
        }

        else if(strstr(inptline,"jdayStart") != NULL) {
          sscanf(inptline, "%d %d", &sitepar->jdayStart, &sitepar->jdayEnd);
          if (flags->verbose) {
            printf("sitepar; denitrification respiration restraint dates, jdayStart/jdayEnd: %d %d\n",
                   sitepar->jdayStart, sitepar->jdayEnd);
          }
        }

        /* nitrification field capacity N2O adjustment factor, cak - 05/09/2013 */
        else if(strstr(inptline,"N2Oadjust_fc") != NULL) {
          sitepar->N2Oadjust_fc = fval;
          if (flags->verbose) {
            printf("sitepar; N2O adjustment factor at fc, N2Oadjust_fc: %f\n", sitepar->N2Oadjust_fc);
          }
        }

        /* nitrification wilting point N2O adjustment factor, cak - 05/09/2013 */
        else if(strstr(inptline,"N2Oadjust_wp") != NULL) {
          sitepar->N2Oadjust_wp = fval;
          if (flags->verbose) {
            printf("sitepar; N2O adjustment factor at fc, N2Oadjust_wp: %f\n", sitepar->N2Oadjust_wp);
          }
        }

        /* coefficient to control the maximum daily nitrification amount - cak 09/21/2011 */
        else if(strstr(inptline,"MaxNitAmt") != NULL) {
          sitepar->MaxNitAmt = fval;
          if (flags->verbose) {
            printf("sitepar; maximum daily nitrification, MaxNitAmt: %f\n", sitepar->MaxNitAmt);
          }
        }

        /* flag to turn off the insulating effect of snow on soil surface temperature - cak 10/07/2011 */
        else if(strstr(inptline,"SnowFlag") != NULL) {
          sscanf(inptline, "%d", &sitepar->SnowFlag);
          if (flags->verbose) {
            printf("sitepar; Snow insulating flag, SnowFlag: %d\n", sitepar->SnowFlag);
          }
        }

        /* coefficient to control the fraction of new net mineralization going to NO3 - cak 01/09/2012 */
        else if(strstr(inptline,"netmn_to_no3") != NULL) {
          sitepar->netmn_to_no3 = fval;
          if (flags->verbose) {
            printf("sitepar; Fraction of net mineralization to NO3, netmn_to_no3: %f\n", sitepar->netmn_to_no3);
          }
        }

       /* inflection point for the water filled pore space effect on denitrification - cak 01/06/2013 */
        else if(strstr(inptline,"wfpsdnitadj") != NULL) {
          sitepar->wfpsdnitadj = fval;
          if (flags->verbose) {
            printf("sitepar; WFPS inflection point adjustment, wfpsdnitadj: %f\n", sitepar->wfpsdnitadj);
          }
        }

        /* Added N2/N2O ratio adjustment coefficient - cak 01/06/2013 */
        else if(strstr(inptline,"N2N2Oadj") != NULL) {
          sitepar->N2N2Oadj = fval;
          if (flags->verbose) {
            printf("sitepar; N2/N2O ratio adjustment, N2N2Oadj: %f\n", sitepar->N2N2Oadj);
          }
        }

        /* Added elevation (meters) to sitepar.in file, cak - 04/15/2009 */
        else if(strstr(inptline,"elevation") != NULL) {
          sitepar->elevation = fval;
          if (flags->verbose) {
            printf("sitepar; Elevation (meters), elevation: %f\n", sitepar->elevation);
          }
        }

        /* Added site slope, aspect, east horizon, and west horizon, cak - 05/08/2009 */
        else if(strstr(inptline,"sitslp") != NULL  || strstr(inptline,"slope") != NULL) {
          sitepar->sitslp = fval;
          if (flags->verbose) {
            printf("sitepar; Site slope (degrees), elevation: %f\n", sitepar->sitslp);
          }
        }

        else if(strstr(inptline,"aspect") != NULL) {
          sitepar->aspect = fval;
          if (flags->verbose) {
            printf("sitepar; Site aspect (degrees), aspect: %f\n", sitepar->aspect);
          }
        }

        else if(strstr(inptline,"ehoriz") != NULL) {
          sitepar->ehoriz = fval;
          if (flags->verbose) {
            printf("sitepar; Site east horizon (degrees), ehoriz: %f\n", sitepar->ehoriz);
          }
        }

        else if(strstr(inptline,"whoriz") != NULL) {
          sitepar->whoriz = fval;
          if (flags->verbose) {
            printf("sitepar; Site west horizon (degrees), whoriz: %f\n", sitepar->whoriz);
          }
        }

        else if(strstr(inptline,"sradadj") != NULL  ||
                strstr(inptline,"solar radiation") != NULL) {
          if(fabsf(fval - 1.0f) < 0.00001) {
            sscanf(inptline, "%d  %f", &m, &sitepar->sradadj[0]);
          } else {sitepar->sradadj[0] = fval;}
          for(imo=1; imo<12; imo++) {
            if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
              sscanf(inptline, "%f", &sitepar->sradadj[imo]);
              if(fabsf(sitepar->sradadj[imo] - (float)(imo+1)) < 0.00001) {
                sscanf(inptline, "%d  %f", &m, &sitepar->sradadj[imo]);
              }
            } else {
              printf("WARNING: Missing Solar radiation transmission, sradadj[%d]\n", imo);
            }
          }
          if (flags->verbose) {
            printf("sitepar; Solar radiation transmission, sradadj:");
            for(imo=0; imo<12; imo++) {printf(" %f", sitepar->sradadj[imo]);}
            printf("\n");
          }
        }

        /* unable to find a variable name */
        else {printf("WARNING: Unknown sitepar record: %s", inptline);fflush(stdout);}
      }

      fclose(fp_in);
      first=0;

      if (flags->debug) {printf("Exiting function initsite\n");fflush(stdout);}
      return;
    }

    void defaultsitpar(SITEPAR_SPT sitepar, LAYERPAR_SPT layers, FLAG_SPT flags)
    {
      /* default values from sitpar.in */
      int  imo;

      /* average solar radiation adjustment for cloud cover & transmission coeffient */
      for(imo=0; imo<=11; imo++) {sitepar->sradadj[imo] = 0.56;}
      sitepar->hours_rain   = 10;       /* duration of each rain event */
      sitepar->elevation    = 0;        /* elevation, meters Brookings SDE */
      sitepar->sitslp       = 0.0;      /* site slope, degrees */
      sitepar->aspect       = 0.0;      /* site aspect, degrees */
      sitepar->ehoriz       = 0.0;      /* site east horizon, degrees */
      sitepar->whoriz       = 0.0;      /* site west horizon, degrees */
      sitepar->sublimscale  = 1;        /* PET to sublimation conversion scaler*/
      sitepar->reflec       = 0.18;     /* vegetation reflectivity (frac) hardwoods = 0.20 spruce = 0.10 */
      sitepar->albedo       = 0.65;     /* snow albedo (frac) */
      sitepar->dmpflux      = 0.000001; /* h2oflux routine (8.E-6?) Parton recommends 10^-6 26 Aug 15 */
      sitepar->hpotdeep     = -200;     /* hydraulic water potential of deep storage layer (units?) */
      sitepar->ksatdeep     = 0.0003;   /* saturated hydraulic conductivity of deep storage layer (cm/sec) */
      sitepar->dmp          = 0.003;    /* damping factor for calculating soil temperature by layer */
      sitepar->Ncoeff       = 0.03;     /* min water/temperature limitation coefficient for nitrify */
      sitepar->jdayStart    = 0;        /* turn off respiration restraint on denit between these Julian dates */
      sitepar->jdayEnd      = 0;
      sitepar->N2Oadjust_fc = 0.012;    /* maximum proportion of nitrified N lost as N2O @ field capacity   1) 0.010  2) 0.015  3) 0.012 */
      sitepar->N2Oadjust_wp = 0.012;    /* minimum proportion of nitrified N lost as N2O @ wilting point    1) 0.002  2) 0.010  3) 0.012 */
      sitepar->MaxNitAmt    = 1.0;      /* maximum daily nitrification amount (gN/m^2) (0.4)  1.0 is value suggested by delgrosso Apr 19, 2012 */
      sitepar->netmn_to_no3 = 0.2;      /* fraction of new net mineralization that goes to NO3 (0.0-1.0) */
      sitepar->wfpsdnitadj  = 1.0;      /* adjustment on inflection point for WFPS effect on denit */
      sitepar->N2N2Oadj     = 1.0;      /* N2/N2O ratio adjustment coefficient */

      return;
    }
