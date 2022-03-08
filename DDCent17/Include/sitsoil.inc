!******************************************************************************
!**
!**  FILE:    sitsoil.h
!**
!**  AUTHOR: Kendrick Killian  March 17 2014
!**   This strips the soil layer and sitepar structures from the original
!**   soilwater.h header file so they can be made available to the fortran
!**   routines. To do so, these definitions MUST be converted to fortran syntax.
!**   Making this a seperate file facilitates the conversion to the Fortran 2003
!**   fortran_iso_C_binding.
!**
!**   include this in soilwater.h
!**
!**   Currently we use code to convert these to Fortran common blocks, let the
!**   fortran allocate the memory then access the memory in C using as C extern
!**   structures. Doing this changes the allocation and pointer statements but
!**   nothing else in the C side.
!**
!**   This is plan B. The original intent was to pass the struct locations to
!**   Fortran. For now, the syntax needed to pass struct locations escapes me
!**   in the layers of indirection and C vs fortran pointers.
!**
!*******************************************************************************



! typedef struct {
!    int verbose;
!    int debug;
!  } FLAG_S, *FLAG_SPT;



!      TYPE, BIND(C) :: layerparstruct
         real(c_double) :: swc(MXSWLYR);
         real(c_double) :: swcmin(MXSWLYR);
! integer values put here to allow 64bit integers
         integer(c_int) :: numlyrs;
         integer(c_int) :: nelyrs;
         integer(c_int) :: ubnd(CMXLYR);
         integer(c_int) :: lbnd(CMXLYR);
! values input directly drom the soil.in file
         real(c_float) :: dpthmn(MXSWLYR);
         real(c_float) :: dpthmx(MXSWLYR);
         real(c_float) :: lyblkd(MXSWLYR);
         real(c_float) :: fieldc(MXSWLYR);
         real(c_float) :: wiltpt(MXSWLYR);
         real(c_float) :: ecoeff(MXSWLYR);
         real(c_float) :: tcoeff(MXSWLYR);
         real(c_float) :: sandfrac(MXSWLYR);
         real(c_float) :: clayfrac(MXSWLYR);
         real(c_float) :: orgfrac(MXSWLYR);
         real(c_float) :: satcond(MXSWLYR);
         real(c_float) :: lyrpH(MXSWLYR);
! calculated values
         real(c_float) :: width(MXSWLYR);
         real(c_float) :: depth(MXSWLYR);
         real(c_float) :: swclimit(MXSWLYR); ! swclimit(numlyrs) = wiltpt(numlyrs) - deltamin
         real(c_float) :: swcfc(MXSWLYR);    ! swcfc[ilyr]  = fieldc[ilyr] * width[ilyr];
         real(c_float) :: swcwp(MXSWLYR);    ! swcwp[ilyr]  = wiltpt[ilyr] * width[ilyr];
         real(c_float) :: minpot(MXSWLYR);
         real(c_float) :: wfps(MXSWLYR);
         real(c_float) :: thetas(MXSWLYR);
         real(c_float) :: thetas_bd(MXSWLYR);
         real(c_float) :: psis(MXSWLYR);
         real(c_float) :: b(MXSWLYR);
         real(c_float) :: sumecoeff;
!      END TYPE layerparstruct 

       common /layerparstruct/ swc, swcmin, numlyrs, nelyrs, ubnd,       &
     &   lbnd, dpthmn, dpthmx, lyblkd, fieldc, wiltpt, ecoeff, tcoeff,   &
     &   sandfrac, clayfrac, orgfrac, satcond, lyrpH, width, depth,      &
     &   swclimit, swcfc, swcwp, minpot, wfps, thetas, thetas_bd, psis,  &
     &   b, sumecoeff
       SAVE /layerparstruct/
       BIND(C) :: /layerparstruct/ 





!      TYPE, BIND(C) :: siteparstruct
         real(c_double) :: dDO_fc;            ! WFPS normalized diffusivity in aggregate soil at field capacity
         real(c_double) :: dDO_wp;            ! WFPS normalized diffusivity in aggregate soil at wilting point
         integer(c_int) :: usexdrvrs;
         integer(c_int) :: texture;
         integer(c_int) :: drainlag;
         integer(c_int) :: jdayStart;
         integer(c_int) :: jdayEnd;
         integer(c_int) :: SnowFlag;          ! snow insulation effect on soil surface temp: 0 = not insulating, 1 = insulating
         real(c_float) :: avgsand;           ! average sand to 10 cm. used for texture and methane
         real(c_float) :: sublimscale;
         real(c_float) :: reflec;
         real(c_float) :: albedo;
         real(c_float) :: dmpflux;
         real(c_float) :: hours_rain;
         real(c_float) :: hpotdeep;
         real(c_float) :: ksatdeep;
         real(c_float) :: rlatitude;
!  float  cldcov(NMONTH+1);
         real(c_float) :: sradadj(NMONTH);   ! solar radiation adjustment for cloud cover & transmission coefficent
         real(c_float) :: dmp;
         real(c_float) :: Ncoeff;
! changes for WFPS effect on N2O flux Cindy K 7 Aug 2013
         real(c_float) :: N2Oadjust_fc;      ! Maximum proportion of nitrified N lost as N2O .. field capacity
         real(c_float) :: N2Oadjust_wp;      ! minimum proportion of nitrified N lost as N2O .. wilting point
         real(c_float) :: MaxNitAmt;         ! maximum daily nitrification amount (gN/m^2)
         real(c_float) :: netmn_to_no3;      ! fraction of new net mineralization that goes to NO3 (0.0-1.0)
         real(c_float) :: wfpsdnitadj;       ! adjustment on inflection point for water filled pore space
         real(c_float) :: N2N2Oadj;          ! multiplier on N2/N2O ratio
! variables for site specific daylength
         real(c_float) :: elevation; ! elevation elevation, meters Brookings SDE
         real(c_float) :: sitslp;    ! sitslp    site slope, degrees    sitslp rather than SLOPE constant
         real(c_float) :: aspect;    ! aspect    site aspect, degrees
         real(c_float) :: ehoriz;    ! ehoriz    site east horizon, degrees
         real(c_float) :: whoriz;    ! whoriz    site west horizon, degrees
! variables for methanogenesis model
         integer(c_int) :: floodN2delay;   ! TEST INPUT: days to transition to flooded N2/N2O ratio
         real(c_float) :: flood_N2toN2O;    ! N2/N2O ratio for flooded state (100.0) (-1 disable)
         real(c_float) :: CO2_to_CH4;       ! fraction of heterotrophic soil respiration to CH4
! float  C6H12O6_to_CH4;  methane_production bug fix; constant from reaction dynamics not input KLK 10Mar15
         real(c_float) :: methzr;           ! methane fraction emitted via bubbles at 0 root mass
         real(c_float) :: frexud;           ! fraction of root production to root exudates (default)
         real(c_float) :: mrtblm;           ! root biomass that starts to reduce methane bubble formation (default)
         real(c_float) :: frCH4emit;        ! fraction of methane emitted by rice plants  KLK 27Jul15
         real(c_float) :: frac_to_exudates; ! fraction of root production to root exudates
         real(c_float) :: Aeh;              ! methane differential coefficient (Aeh)
         real(c_float) :: Deh;              ! methane differential coefficient (Deh)
         real(c_float) :: Beh_flood;        ! methane Eh lower-limit during flooding
         real(c_float) :: Beh_drain;        ! methane Eh upper-limit during drainage
         real(c_float) :: zero_root_frac;   ! methane fraction emitted via bubbles at 0 root mass
         real(c_float) :: ch4rootlim;       ! root biomass that starts to reduce methane bubble formation (default)
!      END TYPE siteparstruct 

       common /siteparstruct/ dDO_fc, dDO_wp, usexdrvrs, texture,        &
     &   drainlag, jdayStart, jdayEnd, SnowFlag, avgsand, sublimscale,   &
     &   reflec, albedo, dmpflux, hours_rain, hpotdeep, ksatdeep,        &
     &   rlatitude, sradadj, dmp, Ncoeff, N2Oadjust_fc, N2Oadjust_wp,    &
     &   MaxNitAmt, netmn_to_no3, wfpsdnitadj, N2N2Oadj, elevation,      &
     &   sitslp, aspect, ehoriz, whoriz, floodN2delay, flood_N2toN2O,    &
     &   CO2_to_CH4, methzr, frexud, mrtblm, frCH4emit,                  &
     &   frac_to_exudates, Aeh, Deh, Beh_flood, Beh_drain,               &
     &   zero_root_frac, ch4rootlim
       SAVE /siteparstruct/
       BIND(C) :: /siteparstruct/ 


