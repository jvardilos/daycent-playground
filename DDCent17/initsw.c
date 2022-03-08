
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsw.c
**
**  FUNCTION:  void initsw()
**
**  PURPOSE:   Initialize the soil water model
**
**  INPUTS:
**    sitlat - latitude (degrees)
**
**  GLOBAL VARIABLES:
**    BAR2CM   - conversion factor for bars to centimeters H2O (1024)
**               (1 bar = 1024 cm H2O)
**    FNSITE   - file name for site specific input parameters (sitepar.in)
**    FNSOIL   - file name for soil layer structure input file (soils.in)
**    MXSWLYR  - maximum number of soil water model layers (21)
**    MAXSTLYR - maximum number of 5 centimeter layers for the soil
**               temperature model (200)
**    PI       - pi (3.14159265)
**
**  EXTERNAL VARIABLES:
**    files                  - structure containing information about output
**                             files
**    flags                  - structure containing debugging flags
**    layers                 - soil water soil layer structure
**    layers->numlyrs        - total number of layers in the soil water model
**                             soil profile
**    layers->swcfc[]        - volumetric soil water content at field capacity
**                             for layer (cm H2O/cm of soil)
**    rwcf                    initial soil water content from site.100
**                             <0 field capacity fraction, >0 rwcf,
**    layers->swclimit[]     - minimum volumetric soil water content of a
**                             layer, fraction 0.0 - 1.0
**    layers->swcwp[]        - volumetric soil water content at wilting point
**                             for layer (cm H2O)
**    layers->width[]        - the thickness of soil water model layers (cm)
**    sitepar                - site specific parameters structure for soil
**                             water model
**    soil                   - soil temperature structure
**
**  LOCAL VARIABLES:
**    callname - call name for subroutine
**    ilyr     - current layer in the soil profile
**    latitude - site latitude (decimal degrees)
**    line[]   - buffer containing line read from input file
**    MAXL     - maximum length of line read from input file
**    fswcinit - flag for initializing swc from rwcf
**
**  OUTPUTS:
**    daylength[]           - length of day light (hours)
**    files->fp_bio         - file pointer to bio.out output file
**    files->fp_cflows      - file pointer to cflows.out output file
**    files->fp_co2         - file pointer to co2.out output file
**    files->fp_dcsip       - file pointer to dc_sip.csv output file
**    files->fp_deadc       - file pointer to deadc.out output file
**    files->fp_dels        - file pointer to dels.out output file
**    files->fp_dN2lyr      - file pointer to dN2lyr.out output file
**    files->fp_dN2Olyr     - file pointer to the dN2Olyr.out output file
**    files->fp_gresp       - file pointer to gresp.out output file
**    files->fp_harv        - file pointer to the harvest.csv output file
**    files->fp_livec       - file pointer to livec.out output file
**    files->fp_methane     - file pointer to methane.out output file
**    files->fp_mresp       - file pointer to mresp.out output file
**    files->fp_outf        - file pointer to outfiles.in input file
**    files->fp_soilc       - file pointer to soilc.out output file
**    files->fp_soiln       - file pointer to soiln.out output file
**    files->fp_soiltavg    - file pointer to soiltavg.out output file
**    files->fp_soiltmax    - file pointer to soiltmax.out output file
**    files->fp_soiltmin    - file pointer to soiltmin.out output file
**    files->fp_stempdx     - file pointer to stemp_dx.out output file
**    files->fp_swc         - file pointer to vswc.out output file
**    files->fp_sysc        - file pointer to sysc.out output file
**    files->fp_tgmonth     - file pointer to the tgmonth.out output file
**    files->fp_wb          - file pointer to watrbal.out output file
**    files->fp_wflux       - file pointer to wflux.out output file
**    files->fp_wfps        - file pointer to wfps.out output file
**    files->fp_yearsum     - file pointer to year_summary.out output file
**    files->fp_yrcflows    - file pointer to year_cflows.out output file
**    files->write_bio      - flag to indicate if bio.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_cflows   - flag to indicate if cflows.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_co2      - flag to indicate if co2.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_dcsip    - flag to indicate if dc_sip.csv output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_deadc    - flag to indicate if deadc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dels     - flag to indicate if dels.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_dN2lyr   - flag to indicate if dN2lyr.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dN2Olyr  - flag to indicate if dN2Olyr.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_gresp    - flag to indicate if gresp.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_harvest  - flag to indicate if harvest.csv output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_livec    - flag to indicate if livec.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_methane  - flag to indicate if methane.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_mresp    - flag to indicate if mresp.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soilc    - flag to indicate if soilc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiln    - flag to indicate if  output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soiltavg - flag to indicate if soiltavg.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmax - flag to indicate if soiltmax.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmin - flag to indicate if soiltmin.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_stempdx  - flag to indicate if stemp_dx.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_swc      - flag to indicate if vswc.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_sysc     - flag to indicate if sysc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_tgmonth  - flag to indicate if tgmonth.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wb       - flag to indicate if watrbal.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wflux    - flag to indicate if wflux.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_wfps     - flag to indicate if wfps.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_yearsum  - flag to indicate if year_summary.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_yrcflow  - flag to indicate if year_cflows.out output file
**                            should be created, 0 = do not create, 1 = create
**    flags->debug          - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose        - flag to set verbose debugging mode, 0 = off,
**                            1 = on
**    layers->minpot[]      - minimum matric potential by layer based on
**                            swcmin (-cm)
**    layers->swc[]         - soil water content by layer (cm H2O)
**    layers->swcmin[]      - lower bound on soil water content by layer
**                            (cm H2O) swc will not be allowed to drop below
**                            this minimum
**    sitepar->rlatitude    - latitude of the site (in radians)
**    soil->soiltavg[]      - average soil temperature of layer (degrees C)
**    soil->soiltmax[]      - maximum soil temperature by layer (degrees C)
**    soil->soiltmin[]      - minimum soil temperature by layer (degrees C)
**    soil->stmtemp[]       - the average soil temperature of the soil
**                            temperature model layers (degrees C)
**    sradadj[]             - solar radiation adjustment for cloud cover and
**                            transmission coeffient
**
**  Changes:
**    KLK 5/4/2015 KLK
**      Since SWC is available to detive in sitsoil, eliminate the float copy swcinit
**      simplified watrbal header to a single printf
**
**  CALLED BY:
**    detiv
**
**  CALLS:
**    initlyrs,soils.in, and initsite, sitepar.in, are being called from detive
**     with other setup routines
**      equvalent sequence would be: inputlyrs();initlyrs(); initsite();
**      with sitepar, layers, flags as arguments to all 3 routines*
**    initsrad()  - initialize the solar radiation submodel
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "soilwater.h"
#include "mssg.h"

#define MAXL 150

static int first=1;

FLAG_S flagstruct;                 FLAG_SPT flags = &flagstruct;
/* LAYERPAR_S lyrstruct;  LAYERPAR_SPT layers = &lyrstruct;
   SITEPAR_S sitestruct;  SITEPAR_SPT sitepar = &sitestruct; */
extern LAYERPAR_S layerparstruct;  LAYERPAR_SPT layers  = &layerparstruct;
extern SITEPAR_S  siteparstruct;   SITEPAR_SPT  sitepar = &siteparstruct;
SOIL_S soilstruct;                 SOIL_SPT soil = &soilstruct;
FILES_S filestruct;                FILES_SPT files = &filestruct;

    void initsw(float *sitlat, float daylength[NDAY],
                float rwcf[CMXLYR], float adep[CMXLYR])
    {

      int  ilyr, clyr;
      static char *callname = "initsw";

      float totd;

      flags->debug = 0;
      flags->verbose = 0;

      /* make the cast of these sitepar float's to double explicit */
      initsrad((double)sitepar->elevation, (double)*sitlat, (double)sitepar->sitslp,
            (double)sitepar->aspect,(double)sitepar->ehoriz, (double)sitepar->whoriz,
            daylength);

      sitepar->rlatitude = *sitlat * (float)(PI/180.0);

      if (flags->verbose) {
        printf("sitlat = %6.2f\n", *sitlat);
        printf("rlatitude = %6.2f\n", sitepar->rlatitude);
        printf("numlyrs = %6d\n", layers->numlyrs);
      }

      if (flags->verbose) printf("%2s  %8s  %8s  %8s\n", "ly", "swc", "swcmin", "minpot");

       /* sum the rwcf. negative values are fraction of field capacity.
         Note for this to work all rwcf elements must be properly initialized
         by the site file read */
       float fswcinit=0;
       for (ilyr=0; ilyr<CMXLYR; ilyr++) {fswcinit += rwcf[ilyr];}
       clyr=0;
       totd=adep[clyr];

       for (ilyr=0; ilyr<layers->numlyrs; ilyr++) {

        if(first ==1) {
          if(layers->swc[ilyr] == 0) {
            if(layers->dpthmx[ilyr] > totd) {clyr++; totd += adep[clyr];}
            if     (rwcf[clyr] >= 1.0f || rwcf[clyr] <= -1.0f) {
                  layers->swc[ilyr] = layers->swcfc[ilyr];}
            else if(fswcinit < 0.0f) {
                  layers->swc[ilyr] = -rwcf[clyr] * layers->swcfc[ilyr];
            }
/*On Aug 7, 2014, at 9:44 AM, Hartman,Melannie wrote:
      I talked to Bill about the rwcf calculation.  We want to keep it as it is,
      using the minimum water content (swclimit) of each layer rather than the
      wilting point: rwcf=(vswc-swclimit)/(fieldc-swclimit).
      changed the input calculation to reflect the swclimit limit for rwcf  KLK 8 2014*/
            else {layers->swc[ilyr] = (rwcf[clyr]*(layers->fieldc[ilyr]-layers->swclimit[ilyr]) +
                                       layers->swclimit[ilyr])*layers->width[ilyr];
            }
          }
        }

        /* Set the lower limit on soil water content and potential. */

        layers->swcmin[ilyr] = min(layers->swclimit[ilyr]*layers->width[ilyr],
                                   layers->swcwp[ilyr]);
        layers->minpot[ilyr] = -swpotentl(layers->swcmin[ilyr],ilyr,layers,
                                          callname)*BAR2CM;

        if (flags->verbose) {
          printf("%2d  %8.4f  %8.4f  %8.4f\n", ilyr, layers->swc[ilyr],
                 layers->swcmin[ilyr], layers->minpot[ilyr]);
        }
      }
      first = 0; /* we have now set the swc */

      for (ilyr=0; ilyr<MXSWLYR; ilyr++) {
        soil->soiltavg[ilyr] = 0.0f;
        soil->soiltmin[ilyr] = 0.0f;
        soil->soiltmax[ilyr] = 0.0f;
      }

      for (ilyr=0; ilyr<MAXSTLYR; ilyr++) {
        soil->stmtemp[ilyr] = 0.0f;
      }

      layers->swc[layers->numlyrs] = 0.0;

      if (flags->debug) {
        printf("Exiting initsw...\n");
      }

      return;
    }


/*****************************************************************************
**
**  FILE:      initsw.c
**
**  FUNCTION:  void outfils()
**
**  PURPOSE:   read the auxiliary output files list
**
**  INPUTS:
**
**  GLOBAL VARIABLES:
**
**  EXTERNAL VARIABLES:
**    files                  - structure containing information about output
**                             files
**    flags                  - structure containing debugging flags
**
**  LOCAL VARIABLES:
**    callname - call name for subroutine
**    line[]   - buffer containing line read from input file
**    lcnt     - count of number of input file lines read
**    wrt      - flag, 0 = do not write to output file,
**               1 = do not write to output file
**
**  OUTPUTS:
**    daylength[]           - length of day light (hours)
**    files->fp_bio         - file pointer to bio.out output file
**    files->fp_cflows      - file pointer to cflows.out output file
**    files->fp_co2         - file pointer to co2.out output file
**    files->fp_dcsip       - file pointer to dc_sip.csv output file
**    files->fp_deadc       - file pointer to deadc.out output file
**    files->fp_dels        - file pointer to dels.out output file
**    files->fp_dN2lyr      - file pointer to dN2lyr.out output file
**    files->fp_dN2Olyr     - file pointer to the dN2Olyr.out output file
**    files->fp_gresp       - file pointer to gresp.out output file
**    files->fp_harv        - file pointer to the harvest.csv output file
**    files->fp_livec       - file pointer to livec.out output file
**    files->fp_methane     - file pointer to methane.out output file
**    files->fp_mresp       - file pointer to mresp.out output file
**    files->fp_outf        - file pointer to outfiles.in input file
**    files->fp_soilc       - file pointer to soilc.out output file
**    files->fp_soiln       - file pointer to soiln.out output file
**    files->fp_soiltavg    - file pointer to soiltavg.out output file
**    files->fp_soiltmax    - file pointer to soiltmax.out output file
**    files->fp_soiltmin    - file pointer to soiltmin.out output file
**    files->fp_stempdx     - file pointer to stemp_dx.out output file
**    files->fp_swc         - file pointer to vswc.out output file
**    files->fp_sysc        - file pointer to sysc.out output file
**    files->fp_tgmonth     - file pointer to the tgmonth.out output file
**    files->fp_wb          - file pointer to watrbal.out output file
**    files->fp_wflux       - file pointer to wflux.out output file
**    files->fp_wfps        - file pointer to wfps.out output file
**    files->fp_yearsum     - file pointer to year_summary.out output file
**    files->fp_yrcflows    - file pointer to year_cflows.out output file
**    files->write_bio      - flag to indicate if bio.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_cflows   - flag to indicate if cflows.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_co2      - flag to indicate if co2.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_dcsip    - flag to indicate if dc_sip.csv output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_deadc    - flag to indicate if deadc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dels     - flag to indicate if dels.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_dN2lyr   - flag to indicate if dN2lyr.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dN2Olyr  - flag to indicate if dN2Olyr.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_gresp    - flag to indicate if gresp.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_harvest  - flag to indicate if harvest.csv output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_livec    - flag to indicate if livec.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_mresp    - flag to indicate if mresp.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soilc    - flag to indicate if soilc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiln    - flag to indicate if  output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soiltavg - flag to indicate if soiltavg.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmax - flag to indicate if soiltmax.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmin - flag to indicate if soiltmin.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_stempdx  - flag to indicate if stemp_dx.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_swc      - flag to indicate if vswc.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_sysc     - flag to indicate if sysc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_tgmonth  - flag to indicate if tgmonth.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wb       - flag to indicate if watrbal.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wflux    - flag to indicate if wflux.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_wfps     - flag to indicate if wfps.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_yearsum  - flag to indicate if year_summary.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_yrcflow  - flag to indicate if year_cflows.out output file
**                            should be created, 0 = do not create, 1 = create
**    flags->debug          - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose        - flag to set verbose debugging mode, 0 = off,
**                            1 = on
**
**           write_methane  - flag to indicate if methane.out output file
**                            should be created, 0 = do not create, 1 = create
**  Changes:
**    KLK 5/4/2015 KLK
**      Pulled the output files into it's own subroutine
**
**  CALLED BY:
**    detiv
**
*****************************************************************************/

    void outfils()
    {

      int wrt, lcnt;
      static char *callname = "outfils";
      char line[MAXL];

      /* variables to switch on the fortran written daily.out, nflux.out
         and summary.out. If the file is to be used, write the header and close.
         dailymoist will append data using the same writes as before*/
      FILE *fp_fortran;
      int  write_daily;
      int  write_nflux;
      int  write_summary;
      int  write_methane;

      if ((files->fp_outf = fopen("outfiles.in", "r")) == NULL) {
        abortmssg("Cannot open file outfiles.in");
      }

      files->write_bio = 0;      {remove("bio.out");}
      files->write_soiln = 0;    {remove("soiln.out");}
      files->write_soiltavg = 0; {remove("soiltavg.out");}
      files->write_soiltmax = 0; {remove("soiltmax.out");}
      files->write_soiltmin = 0; {remove("soiltmin.out");}
      files->write_stempdx = 0;  {remove("stemp_dx.out");}
      files->write_swc = 0;      {remove("vswc.out");}
      files->write_wb = 0;       {remove("watrbal.out");}
      files->write_wfps = 0;     {remove("wfps.out");}
      files->write_co2 = 0;      {remove("co2.out");}
      files->write_wflux = 0;    {remove("wflux.out");}
      files->write_mresp = 0;    {remove("mresp.out");}
      files->write_yearsum = 0;  {remove("year_summary.out");}
      files->write_livec = 0;    {remove("livec.out");}
      files->write_deadc = 0;    {remove("deadc.out");}
      files->write_soilc = 0;    {remove("soilc.out");}
      files->write_sysc = 0;     {remove("sysc.out");}
      files->write_tgmonth = 0;  {remove("tgmonth.out");}
      files->write_dN2lyr = 0;   {remove("dN2lyr.out");}
      files->write_dN2Olyr = 0;  {remove("dN2Olyr.out");}
      files->write_gresp = 0;    {remove("gresp.out");}
      files->write_dels = 0;     {remove("dels.out");}
      files->write_dcsip = 0;    {remove("dc_sip.csv");}
      files->write_harvest = 0;  {remove("harvest.csv");}
      files->write_cflows = 0;   {remove("cflows.out");}
      files->write_yrcflows = 0; {remove("year_cflows.out");}
      write_daily = 0;           {remove("daily.out");}
      write_nflux = 0;           {remove("nflux.out");}
      write_summary = 0;         {remove("summary.out");}
      write_methane = 0;         {remove("methane.out");}

      lcnt = 0;
      while( fgets(line, MAXL, files->fp_outf) != NULL) {
/*        printf("%s", line);*/
        lcnt++;
        wrt=0;
          sscanf(line, "%d", &wrt);
/*          printf("wrt = %d\n", wrt);*/
        if     (lcnt == 1)  {continue;}
             if(lcnt == 2)  {files->write_bio = wrt;}
        else if(lcnt == 3)  {files->write_soiln = wrt;}
        else if(lcnt == 4)  {files->write_soiltavg = wrt;}
        else if(lcnt == 5)  {files->write_soiltmax = wrt;}
        else if(lcnt == 6)  {files->write_soiltmin = wrt;}
        else if(lcnt == 7)  {files->write_stempdx = wrt;}
        else if(lcnt == 8)  {files->write_swc = wrt;}
        else if(lcnt == 9)  {files->write_wb = wrt;}
        else if(lcnt == 10) {files->write_wfps = wrt;}
        else if(lcnt == 11) {files->write_co2 = wrt;}
        else if(lcnt == 12) {files->write_wflux = wrt;}
        else if(lcnt == 13) {files->write_mresp = wrt;}
        else if(lcnt == 14) {files->write_yearsum = wrt;}
        else if(lcnt == 15) {files->write_livec = wrt;}
        else if(lcnt == 16) {files->write_deadc = wrt;}
        else if(lcnt == 17) {files->write_soilc = wrt;}
        else if(lcnt == 18) {files->write_sysc = wrt;}
        else if(lcnt == 19) {files->write_tgmonth = wrt;}
        else if(lcnt == 20) {files->write_dN2lyr = wrt;}
        else if(lcnt == 21) {files->write_dN2Olyr = wrt;}
        else if(lcnt == 22) {files->write_gresp = wrt;}
        else if(lcnt == 23) {files->write_dels = wrt;}
        else if(lcnt == 24) {files->write_dcsip = wrt;}
        else if(lcnt == 25) {files->write_harvest = wrt;}
        else if(lcnt == 26) {files->write_cflows = wrt;}
        else if(lcnt == 27) {files->write_yrcflows = wrt;}
        else if(lcnt == 28) {       write_daily = wrt;}
        else if(lcnt == 29) {       write_nflux = wrt;}
        else if(lcnt == 30) {       write_summary = wrt;}
        else if(lcnt == 31) {       write_methane = wrt;}
        if(wrt != 0) {printf("wrt = %d  line %d:     %s", wrt, lcnt, line);}
      }

      if (files->write_soiltavg) {
        files->fp_soiltavg = fopen("soiltavg.out", "w");
      }
      if (files->write_soiltmax) {
        files->fp_soiltmax = fopen("soiltmax.out", "w");
      }
      if (files->write_soiltmin) {
        files->fp_soiltmin = fopen("soiltmin.out", "w");
      }
      if (files->write_stempdx) {
        files->fp_stempdx = fopen("stemp_dx.out", "w");
      }
      if (files->write_swc) {
        files->fp_swc = fopen("vswc.out", "w");
      }
      if (files->write_wfps) {
        files->fp_wfps = fopen("wfps.out", "w");
      }

      if (files->write_wb) {
        files->fp_wb = fopen("watrbal.out", "w");
        fprintf(files->fp_wb, "#0=ppt+dswc-intrcpt-evap-transp-outflow");
        fprintf(files->fp_wb, " (when sublim = 0)\n");
        fprintf(files->fp_wb, "#0=melt+dswc-intrcpt-evap-transp-outflow");
        fprintf(files->fp_wb, " (when sublim > 0)\n");
        fprintf(files->fp_wb, "#%4s %8s %5s %7s %7s %7s %7s %7s %7s %7s "
                              "%7s %8s %9s %7s %7s %7s\n", "time", "dayofyr", "ppt",
                "accum", "dsnlq", "melt", "intrcpt", "evap", "transp",
                "sublim", "dswc", "outflow", "balance", "snow", "snlq", "runoff");
      }

      if (files->write_soiln) {
        printf("Open soiln.out\n");
        files->fp_soiln = fopen("soiln.out", "w");
        fprintf(files->fp_soiln, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "ammonium", "NO3_ppm[0]", "NO3_ppm[1]",
                "NO3_ppm[2]", "NO3_ppm[3]");
        fprintf(files->fp_soiln, "%12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                "NO3_ppm[4]", "NO3_ppm[5]", "NO3_ppm[6]", "NO3_ppm[7]",
                "NO3_ppm[8]", "NO3_ppm[9]", "etc...");
      }

      if (files->write_co2) {
        printf("Open co2.out\n");
        files->fp_co2 = fopen("co2.out", "w");
        fprintf(files->fp_co2, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "CO2_ppm[0]", "CO2_ppm[1]", "CO2_ppm[2]",
                "CO2_ppm[3]", "CO2_ppm[4]");
        fprintf(files->fp_co2, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "CO2_ppm[5]", "CO2_ppm[6]", "CO2_ppm[7]", "CO2_ppm[8]",
                "CO2_ppm[9]", "etc...");
      }

      if (files->write_wflux) {
        printf("Open wflux.out\n");
        files->fp_wflux = fopen("wflux.out", "w");
        fprintf(files->fp_wflux, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "wflux[0]", "wflux[1]", "wflux[2]",
                "wflux[3]", "wflux[4]");
        fprintf(files->fp_wflux, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "wflux[5]", "wflux[6]", "wflux[7]", "wflux[8]", "wflux[9]",
                "etc...");
      }

      if (files->write_bio) {
        printf("Open bio.out\n");
        files->fp_bio = fopen("bio.out", "w");
        fprintf(files->fp_bio, "%8s  %8s  %10s  %10s  %10s  ",
                "time", "dayofyr", "aglivc", "bglivcj", "bglivcm");
        fprintf(files->fp_bio, "%10s  %10s  %10s  ",
                "aglivn", "bglivnj", "bglivnm");
        fprintf(files->fp_bio, "%10s  %10s  %10s  %10s  %10s  %10s",
                "rleavc", "frootcj", "frootcm", "fbrchc", "rlwodc", "crootc");
        fprintf(files->fp_bio, "%10s  %10s\n",
                "h2ogef(1)", "h2ogef(2)");
      }

      if (files->write_mresp) {
        printf("Open mresp.out\n");
        files->fp_mresp = fopen("mresp.out", "w");
        fprintf(files->fp_mresp, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "mrspflux(1)", "mrspflux(2)", "cmrspflux(1)",
                "cmrspflux(2)", "cmrspflux(3)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "fmrspflux(1)", "fmrspflux(2)", "fmrspflux(6)",
                "fmrspflux(3)", "fmrspflux(4)", "fmrspflux(5)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "mcprd(1)", "mcprd(2)", "mcprd(3)", "mfprd(1)", "mfprd(2)",
                "mfprd(6)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "mfprd(3)", "mfprd(4)", "mfprd(5)", "carbostg(1,1)",
                "carbostg(1,2)", "carbostg(2,1)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s  %12s  ",
                "carbostg(2,2)", "mrspann(1)", "mrspann(2)", "tavedly");
        fprintf(files->fp_mresp, "%12s  %12s  %12s  ",
                "mrspTempEffect(1,1)", "mrspTempEffect(1,2)",
                "mrspWaterEffect(1)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s\n",
                "mrspTempEffect(2,1)", "mrspTempEffect(2,2)",
                "mrspWaterEffect(2)");
      }

      if (files->write_yearsum) {
        printf("Open year_summary.out\n");
        files->fp_yearsum = fopen("year_summary.out", "w");
        fprintf(files->fp_yearsum, "%8s  %12s  %12s  %12s  %12s  %12s",
                "time", "N2Oflux", "NOflux", "N2flux", "CH4_oxid", "NIT");
        fprintf(files->fp_yearsum, "%12s\n", "ANNPPT");
      }

      if (files->write_livec) {
        printf("Open livec.out\n");
        files->fp_livec = fopen("livec.out", "w");
        fprintf(files->fp_livec, "%8s  %8s  %10s  %10s  %10s  ", "time",
                "dayofyr", "aglivc", "bglivcj", "bglivcm");
        fprintf(files->fp_livec, "%10s  %10s  %10s  ",
                "rleavc", "frootcj", "frootcm");
        fprintf(files->fp_livec, "%10s  %10s  %10s  %10s\n",
                "fbrchc", "rlwodc", "crootc", "frnutc");
      }

      if (files->write_deadc) {
        printf("Open deadc.out\n");
        files->fp_deadc = fopen("deadc.out", "w");
        fprintf(files->fp_deadc, "%8s  %8s  %10s  %10s  %10s  %10s  ", "time",
                "dayofyr", "stdedc", "metabc(1)", "strucc(1)", "wood1c");
        fprintf(files->fp_deadc, "%10s  %10s\n",
                "wood2c", "wood3c");
      }

      if (files->write_soilc) {
        printf("Open soilc.out\n");
        files->fp_soilc = fopen("soilc.out", "w");
        fprintf(files->fp_soilc, "%8s  %8s  %10s  %10s  %10s  %10s  ", "time",
                "dayofyr", "metabc(2)", "strucc(2)", "som1c(1)", "som1c(2)");
        fprintf(files->fp_soilc, "%10s  %10s %10s\n",
                "som2c(1)", "som2c(2)", "som3c");
      }

      if (files->write_sysc) {
        printf("Open sysc.out\n");
        files->fp_sysc = fopen("sysc.out", "w");
        fprintf(files->fp_sysc, "%8s  %8s  %10s  %10s  %10s  %10s %10s\n",
                "time", "dayofyr", "livec", "deadc", "soilc", "sysc", "CO2resp");
      }

      if (files->write_tgmonth) {
        printf("Open tgmonth.out\n");
        files->fp_tgmonth = fopen("tgmonth.out", "w");
        fprintf(files->fp_tgmonth, "%8s  %12s  %12s  %12s  %12s  %12s",
                "time", "N2Oflux", "NOflux", "N2flux", "CH4_oxid", "NIT");
        fprintf(files->fp_tgmonth, "%12s\n", "PPT");
      }

      if (files->write_dN2lyr) {
        printf("Open dN2lyr.out\n");
        files->fp_dN2lyr = fopen("dN2lyr.out", "w");
        fprintf(files->fp_dN2lyr, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "dN2_g/m2[0]", "dN2_g/m2[1]", "dN2_g/m2[2]",
                "dN2_g/m2[3]", "dN2_g/m2[4]");
        fprintf(files->fp_dN2lyr, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "dN2_g/m2[5]", "dN2_g/m2[6]", "dN2_g/m2[7]", "dN2_g/m2[8]",
                "dN2_g/m2[9]", "etc...");
      }

      if (files->write_dN2Olyr) {
        printf("Open dN2Olyr.out\n");
        files->fp_dN2Olyr = fopen("dN2Olyr.out", "w");
        fprintf(files->fp_dN2Olyr, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "dN2O_g/m2[0]", "dN2O_g/m2[1]",
                "dN2O_g/m2[2]", "dN2O_g/m2[3]", "dN2O_g/m2[4]");
        fprintf(files->fp_dN2Olyr, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "dN2O_g/m2[5]", "dN2O_g/m2[6]", "dN2O_g/m2[7]",
                "dN2O_g/m2[8]", "dN2O_g/m2[9]", "etc...");
      }

      if (files->write_gresp) {
        printf("Open gresp.out\n");
        files->fp_gresp = fopen("gresp.out", "w");
        fprintf(files->fp_gresp, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "grspflux(1)", "grspflux(2)", "cgrspflux(1)",
                "cgrspflux(2)", "cgrspflux(3)");
        fprintf(files->fp_gresp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "fgrspflux(1)", "fgrspflux(2)", "fgrspflux(6)",
                "fgrspflux(3)", "fgrspflux(4)", "fgrspflux(5)");
        fprintf(files->fp_gresp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "mcprd(1)", "mcprd(2)", "mcprd(3)", "mfprd(1)", "mfprd(2)",
                "mfprd(6)");
        fprintf(files->fp_gresp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "mfprd(3)", "mfprd(4)", "mfprd(5)", "carbostg(1,1)",
                "carbostg(1,2)", "carbostg(2,1)");
        fprintf(files->fp_gresp, "%12s  %12s  %12s\n",
                "carbostg(2,2)", "grspann(1)", "grspann(2)");
      }

      if (files->write_dels) {
        printf("Open dels.out\n");
        files->fp_dels = fopen("dels.out", "w");
        fprintf(files->fp_dels, "%8s  %8s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "deloi", "deloe", "dsrfclit", "dsmnrl");
        fprintf(files->fp_dels, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "dhetresp", "dsoilresp", "dcmresp", "dfmresp", "dcgresp",
                "dfgresp");
        fprintf(files->fp_dels, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "dccarbostg", "dfcarbostg", "oiresp", "oeresp", "slitrsp",
                "sminrlrsp");
        fprintf(files->fp_dels, "%12s  %12s  %12s  %12s  %12s  ",
                "hresp", "crtjresp", "crtmresp", "frtjresp", "frtmresp");
        fprintf(files->fp_dels, "%12s  %12s  %12s  %12s\n",
                "frtcresp", "sresp", "mresp", "gresp");
      }

      if (files->write_dcsip) {
        printf("Open dc_sip.csv\n");
        files->fp_dcsip = fopen("dc_sip.csv", "w");
        fprintf(files->fp_dcsip, "%s,%s,", "time", "dayofyr");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,",
                "trandly", "evapdly", "intrcpt", "sublim", "drain", "runoff");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,%s,",
                "ppt", "accum", "melt", "snow", "snlq", "petdly", "stemp");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,",
                "wc_2cm", "wc_3cm", "wc_5cm", "wc_10cm", "wc_15cm");
        fprintf(files->fp_dcsip, "%s,%s,",
                "wc_30cm", "CO2resp");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,",
                "mcprd(1)", "mcprd(2)", "mcprd(3)", "mfprd(1)", "mfprd(2)",
                "mfprd(6)");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,",
                "mfprd(3)", "mfprd(4)", "mfprd(5)", "NPP", "NEE");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,%s,",
                "aglivc", "bglivcj", "bglivcm", "rleavc", "frootcj",
                "frootcm", "fbrchc");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,",
                "rlwodc", "crootc", "tlai", "stdedc", "wood1c", "wood2c");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,",
                "wood3c", "strucc(1)", "metabc(1)", "strucc(2)", "metabc(2)");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,\n",
                "som1c(1)", "som1c(2)", "som2c(1)", "som2c(2)", "som3c",
                "totsysc");
      }

      if (files->write_harvest) {
        printf("Open harvest.csv\n");
        files->fp_harv = fopen("harvest.csv", "w");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,%s,", "time", "dayofyr",
                "crpval", "agcacc", "bgcjacc", "bgcmacc");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "cgrain", "egrain(N)", "egrain(P)", "egrain(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "crmvst", "ermvst(N)", "ermvst(P)", "ermvst(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "cstraw", "estraw(N)", "estraw(P)", "estraw(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "stdstraw", "estdstraw(N)", "estdstraw(P)", "estdstraw(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "addsdc", "addsde(N)", "addsde(P)", "addsde(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "resid", "reside(N)", "reside(P)", "reside(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,%s,%s,",
                "irrapp", "fertapp(N)", "fertapp(P)", "fertapp(S)",
                "afertapp(N)", "afertapp(P)", "afertapp(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "omadapp", "omaeapp(N)", "omaeapp(P)", "omaeapp(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,%s,%s,%s,",
                "strmac(1)", "strmac(2)", "strmac(3)", "strmac(4)",
                "strmac(5)", "strmac(6)", "strmac(7)", "strmac(8)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "cgracc", "egracc(N)", "egracc(P)", "egracc(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "accrst", "accrste(N)", "accrste(P)", "accrste(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "ctubesj", "etubesj(N)", "etubesj(P)", "etubesj(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "ctubesm", "etubesm(N)", "etubesm(P)", "etubesm(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "srfclittrj", "esrfclittrj(N)", "esrfclittrj(P)", "esrfclittrj(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "soillittrj", "esoillittrj(N)", "esoillittrj(P)", "esoillittrj(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "srfclittrm", "esrfclittrm(N)", "esrfclittrm(P)", "esrfclittrm(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s\n",
                "soillittrm", "esoillittrm(N)", "esoillittrm(P)", "esoillittrm(S)");
      }

      if (files->write_cflows) {
        files->fp_cflows = fopen("cflows.out", "w");
        fprintf(files->fp_cflows, "%8s  %8s  %10s  %10s  %10s  %10s  %10s  ",
                "time", "dayofyr", "som11tosom21", "som12tosom22",
                "som12tosom3", "som21tosom11", "som21tosom22");
        fprintf(files->fp_cflows, "%10s  %10s  %10s  %10s  %10s  ",
                "som22tosom12", "som22tosom3", "som3tosom12", "metc1tosom11",
                "metc2tosom12");
        fprintf(files->fp_cflows, "%10s  %10s  %10s  %10s  ",
                "struc1tosom11", "struc1tosom21", "struc2tosom12",
                "struc2tosom22");
        fprintf(files->fp_cflows, "%10s  %10s  %10s  %10s  ",
                "wood1tosom11", "wood1tosom21", "wood2tosom11",
                "wood2tosom21");
        fprintf(files->fp_cflows, "%10s  %10s\n",
                "wood3tosom12", "wood3tosom22");
      }

      if (files->write_yrcflows) {
        files->fp_yrcflows = fopen("year_cflows.out", "w");
        fprintf(files->fp_yrcflows, "%8s  %10s  %10s  %10s  %10s  %10s  ",
                "time", "asom11tosom21", "asom12tosom22", "asom12tosom3",
                "asom21tosom11", "asom21tosom22");
        fprintf(files->fp_yrcflows, "%10s  %10s  %10s  %10s  %10s  ",
                "asom22tosom12", "asom22tosom3", "asom3tosom12",
                "ametc1tosom11", "ametc2tosom12");
        fprintf(files->fp_yrcflows, "%10s  %10s  %10s  %10s  ",
                "astruc1tosom11", "astruc1tosom21", "astruc2tosom12",
                "astruc2tosom22");
        fprintf(files->fp_yrcflows, "%10s  %10s  %10s  %10s  ",
                "awood1tosom11", "awood1tosom21", "awood2tosom11",
                "awood2tosom21");
        fprintf(files->fp_yrcflows, "%10s  %10s\n",
                "awood3tosom12", "awood3tosom22");
      }

      if (write_daily) {
        fp_fortran = fopen("daily.out", "w");
        fprintf(fp_fortran,
                "%8s  %8s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n",
                "time", "dayofyr", "PET(cm)", "agdefac", "bgdefac",
                "stemp(C)", "snow", "snlq", "thermunits", "srad");
        fclose(fp_fortran);
      }
      if (write_nflux) {
        fp_fortran = fopen("nflux.out", "w");
        fprintf(fp_fortran,
                "%8s  %8s  %10s  %10s  %10s  %10s  %10s  %10s\n",
                "time", "dayofyr", "nit_N2O-N", "dnit_N2O-N", "dnit_N2-N",
                "NO-N", "CUM-N2O(gN/ha)", "CUM-NO(gN/ha)");
        fclose(fp_fortran);
      }
      if (write_summary) {
        fp_fortran = fopen("summary.out", "w");
        fprintf(fp_fortran,
                "%8s  %8s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n",
                "time", "dayofyr", "tmax", "tmin", "ppt",
                "N2Oflux", "NOflux", "CH4oxid", "NIT", "CO2resp");
        fclose(fp_fortran);
      }

      if (write_methane) {
        fp_fortran = fopen("methane.out", "w");
        fprintf(fp_fortran, "%4s   %3s  %10s  %10s  %10s  ",
                "year", "DOY", "aglivc", "bglivcj", "bglivcm");
        fprintf(fp_fortran, "%11s  %11s  %11s  ",
                "prev_mcprd1", "prev_mcprd2", "prev_mcprd3");
        fprintf(fp_fortran, "%10s  %10s  %10s  %10s  %10s  ",
                "COM", "ppt", "irri", "watr2sat", "avgst_10cm");
        fprintf(fp_fortran, "%10s  %10s  %10s  %10s  ",
                "TI", "Cr", "Eh", "Feh");
        fprintf(fp_fortran, "%10s  %10s  %10s  %10s\n",
                "CH4_prod", "CH4_Ep", "CH4_Ebl", "CH4_oxid");
        fclose(fp_fortran);
      }

/*!!
files->fp_snow = fopen("snow.out", "w");
fprintf(files->fp_snow, "%8s %8s %7s %7s %7s %7s %7s", "time", "dayofyr", "tave",
        "rain", "pet", "snow1", "snlq1");
fprintf(files->fp_snow,"%7s %7s %7s %7s %7s %7s", "snow2",
       "snlq2", "sublim", "snow3", "snlq3", "melt");
fprintf(files->fp_snow, "%7s %7s %7s %7s %7s %7s\n", "snow4",
       "snlq4", "snlq5", "pptrem", "petrem", "runoff");
!!*/

      fclose(files->fp_outf);

      if (flags->debug) {
        printf("Exiting outfils...\n");
      }

      return;
    }
