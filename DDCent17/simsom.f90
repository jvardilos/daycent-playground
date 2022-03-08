
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine simsom()
      use iso_c_binding;
      use calflow;
      implicit none

      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'jday.inc'
      include 'ligvar.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

! ... Simulate flow of carbon, nitrogen, phosphorous, and sulfur.
! ... This routine is executed each time step.  It calls the decomposition
! ... submodel and a producer submodel.  It also includes a bunch of
! ... N fixation stuff that needs to be rewritten and put in its own routine.

!     Changes:
!     Sept15: K. Killian
!              moved the tree basal area calculation to potprd where it is used.
!              It now responds daily and does not require the tree be defined at the
!              start of the month.

! ... Added new local variable FSOL.  Added calls to new function FSFUNC
! ... to calculate the amount of mineral P that is in solution.  Added
! ... call to new subroutine PSCHEM, which calculates and schedules the
! ... Phosophorus and Sulfur flows during decomposition.  Previously
! ... this was calculated in the DECOMP routine.  -rm 6/91

! ... Removed decomposition submodel and put in call to dailymoist for
! ... the daily water budget version of Century.  -mdh 8/94

! ... Fortran to C prototype
      INTERFACE

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, ammonium, nitrate, &
          subname) bind(c, name= "update_npool")
          INTEGER          clyr
          REAL             amt
          DOUBLE PRECISION frac_nh4
          DOUBLE PRECISION frac_no3
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          character :: subname(10)  ! character(kind=c_char) :: string(*)
          ! character(kind=c_char) :: subname(10)  ! character(kind=c_char) :: string(*)
        END SUBROUTINE update_npool

!       use iso_c_binding, only: C_CHAR, C_NULL_CHAR
!       interface
!         subroutine print_c(string) bind(C, name="print_C")
!           use iso_c_binding, only: c_char
!           character(kind=c_char) :: string(*)
!         end subroutine print_c
!       end interface
!       call print_c(C_CHAR_"Hello World"//C_NULL_CHAR)


        SUBROUTINE wrtbio(time, curday, aglivc, bglivcj, bglivcm, &
                        aglive, bglivej, bglivem, rleavc, frootcj, &
                        frootcm, fbrchc, rlwodc, crootc, h2ogef1, h2ogef2) &
          bind(c, name= "wrtbio")
          REAL    time
          INTEGER curday
          REAL    aglivc
          REAL    bglivcj
          REAL    bglivcm
          REAL    aglive
          REAL    bglivej
          REAL    bglivem
          REAL    rleavc
          REAL    frootcj
          REAL    frootcm
          REAL    fbrchc
          REAL    rlwodc
          REAL    crootc
          REAL    h2ogef1
          REAL    h2ogef2
        END SUBROUTINE wrtbio

        SUBROUTINE wrtdcsip(time, curday, trandly, evapdly, intrcpt, &
                          sublim, outflow, runoffdly, ppt, saccum,  &
                          melt, snow, snlq, petdly, stemp, wc_2cm,  &
                          wc_3cm, wc_5cm, wc_10cm, wc_15cm, wc_30cm, &
                          CO2resp, mcprd1, mcprd2, mcprd3, mfprd1,  &
                          mfprd2, mfprd6, mfprd3, mfprd4, mfprd5,  &
                          npptot, nee, aglivc, bglivcj, bglivcm,   &
                          rleavc, frootcj, frootcm, fbrchc, rlwodc, &
                          crootc, tlai, stdedc, wood1c, wood2c,   &
                          wood3c, strucc1, metabc1, strucc2, metabc2, &
                          som1c1, som1c2, som2c1, som2c2, som3c, systemc) &
          bind(c, name= "wrtdcsip")
          REAL    time
          INTEGER curday
          REAL    trandly, evapdly, intrcpt
          REAL    sublim, outflow, runoffdly, ppt, saccum
          REAL    melt, snow, snlq, petdly, stemp, wc_2cm
          REAL    wc_3cm, wc_5cm, wc_10cm, wc_15cm, wc_30cm
          REAL    CO2resp, mcprd1, mcprd2, mcprd3
          REAL    mfprd1, mfprd2, mfprd6, mfprd3, mfprd4, mfprd5
          REAL    npptot, nee, aglivc, bglivcj, bglivcm
          REAL    rleavc, frootcj, frootcm, fbrchc, rlwodc
          REAL    crootc, tlai, stdedc, wood1c, wood2c, wood3c
          REAL    strucc1, metabc1, strucc2, metabc2
          REAL    som1c1, som1c2, som2c1, som2c2, som3c
          REAL    systemc
        END SUBROUTINE wrtdcsip

        SUBROUTINE wrtdeadc(time, curday, stdedc, metabc1, strucc1, &
                           wood1c, wood2c, wood3c) &
          bind(c, name= "wrtdeadc")
          REAL    time
          INTEGER curday
          REAL    stdedc
          REAL    metabc1
          REAL    strucc1
          REAL    wood1c
          REAL    wood2c
          REAL    wood3c
        END SUBROUTINE wrtdeadc

      SUBROUTINE wrtdels(time, curday, ddeloi, ddeloe, ddsrfclit, &
                         ddsmnrl, ddhetresp, ddsoilresp, ddcmresp, &
                         ddfmresp, ddcgresp, ddfgresp, ddccarbostg, &
                         ddfcarbostg, doiresp, doeresp, dslitrsp, &
                         dmnrlrsp, dhresp, dcrtjresp, dcrtmresp, &
                         dfrtjresp, dfrtmresp, dfrtcresp, dsresp, &
                         dmresp, dgresp) bind(c, name= "wrtdels")
          REAL    time
          INTEGER curday
          REAL    ddeloi, ddeloe, ddsrfclit
          REAL    ddsmnrl, ddhetresp, ddsoilresp, ddcmresp
          REAL    ddfmresp, ddcgresp, ddfgresp, ddccarbostg
          REAL    ddfcarbostg, doiresp, doeresp, dslitrsp
          REAL    dsmnrlrsp, dhresp, dcrtjresp, dcrtmresp
          REAL    dfrtjresp, dfrtmresp, dfrtcresp, dsresp, dmresp, dgresp
        END SUBROUTINE wrtdels

        SUBROUTINE wrtlivec(time, curday, aglivc, bglivcj, bglivcm, &
                    rleavc, frootcj, frootcm, fbrchc, rlwodc, crootc, frnutc) &
          bind(c, name= "wrtlivec")
          REAL    time
          INTEGER curday
          REAL    aglivc
          REAL    bglivcj
          REAL    bglivcm
          REAL    rleavc
          REAL    frootcj
          REAL    frootcm
          REAL    fbrchc
          REAL    rlwodc
          REAL    crootc
          REAL    frnutc
        END SUBROUTINE wrtlivec

        SUBROUTINE wrtmresp(time, curday,mrspflux1, mrspflux2, &
                  cmrspflux1, cmrspflux2, cmrspflux3, &
                  fmrspflux1, fmrspflux2, fmrspflux6, &
                  fmrspflux3, fmrspflux4, fmrspflux5, &
                  mcprd1, mcprd2, mcprd3, mfprd1, &
                  mfprd2, mfprd6, mfprd3, mfprd4, &
                  mfprd5, carbostg11, carbostg12, &
                  carbostg21, carbostg22, mrspann1, &
                  mrspann2,tavedly, mrspTempEffect11, &
                  mrspTempEffect12, mrspWaterEffect1, &
                  mrspTempEffect21, mrspTempEffect22, &
                  mrspWaterEffect2) &
                   bind(c, name= "wrtmresp")
          REAL    time
          INTEGER curday
          REAL    doiresp
          REAL    mrspflux1
          REAL    mrspflux2
          REAL    cmrspflux1
          REAL    cmrspflux2
          REAL    cmrspflux3
          REAL    fmrspflux1
          REAL    fmrspflux2
          REAL    fmrspflux6
          REAL    fmrspflux3
          REAL    fmrspflux4
          REAL    fmrspflux5
          REAL    mcprd1
          REAL    mcprd2
          REAL    mcprd3
          REAL    mfprd1
          REAL    mfprd2
          REAL    mfprd6
          REAL    mfprd3
          REAL    mfprd4
          REAL    mfprd5
          REAL    carbostg11
          REAL    carbostg12
          REAL    carbostg21
          REAL    carbostg22
          REAL    mrspann1
          REAL    mrspann2
          REAL    tavedly
          REAL    mrspTempEffect11
          REAL    mrspTempEffect12
          REAL    mrspWaterEffect1
          REAL    mrspTempEffect21
          REAL    mrspTempEffect22
          REAL    mrspWaterEffect2
        END SUBROUTINE wrtmresp


        SUBROUTINE wrtsoilc(time, curday, metabc2, strucc2, som1c1, &
                            som1c2, som2c1, som2c2, som3c) &
                   bind(c, name= "wrtsoilc")
          REAL    time
          INTEGER curday
          REAL    metabc2
          REAL    strucc2
          REAL    som1c1
          REAL    som1c2
          REAL    som2c1
          REAL    som2c2
          REAL    som3c
        END SUBROUTINE wrtsoilc

        SUBROUTINE wrtsysc(time, curday, livec, deadc, soilc, sysc, CO2resp) &
          bind(c, name= "wrtsysc")
          REAL    time
          INTEGER curday
          REAL    livec
          REAL    deadc
          REAL    soilc
          REAL    sysc
          REAL    CO2resp
        END SUBROUTINE wrtsysc

      END INTERFACE

! ... Function declarations
      real      fsfunc, delLCout, irrigt
      external  fsfunc, delLCout, irrigt

! ... Local variables
      integer   :: iel, clyr
      integer   :: curday
      integer   :: ostat
      integer,save :: plntd  = 0
      integer,save :: wrtmeth = 0
      logical   :: methopn
      real      :: biof, cmn
      real      :: frlech(MAXIEL), fwdfx, fxbiom
      real      :: satm, sirr, texeff
      real      :: wdfnp, wdfxm
      real      :: irradd, agdefacsum, bgdefacsum, tfrac
!      real     ::  prevgromin, grominwk
      real      :: bgwfunc
      real      :: prevstream(8)
      real      :: tavedly, tavemth, tavewk
      double precision frac_nh4, frac_no3
      real      :: co2val
      real,save :: CO2resp
      real,save :: respsum(2), arspsum(2,2)
      real      :: saccum, evapdly, intrcpt, melt, outflow, petdly
      real      :: runoffdly, sublim, trandly
      real      :: ddeloi, ddeloe, ddsrfclit, ddsmnrl, ddhetresp
      real      :: ddsoilresp
      real      :: oiresp1(2), oiresp2(2), newoiresp(2)
      real      :: oeresp1(2), oeresp2(2), newoeresp(2)
      real      :: slitrsp1(2), slitrsp2(2), newslitrsp(2)
      real      :: mnrlrsp1(2), mnrlrsp2(2), newmnrlrsp(2)
      real      :: hetresp1(2), hetresp2(2), newhetresp(2)
      real      :: ddccarbostg, ddfcarbostg
      real      :: doiresp, doeresp, dslitrsp, dmnrlrsp, dhresp
      real      :: dcrtresp, dcrtjresp, dcrtmresp
      real      :: dfrtresp, dfrtjresp, dfrtmresp, dfrtcresp
      real      :: dsresp, dmresp, dgresp, crplblstg, forlblstg
      real      :: dcmresp, dfmresp, ddcmresp, ddfmresp
      real      :: dcgresp, dfgresp, ddcgresp, ddfgresp
      real      :: npptot, systemc, wc_2cm, wc_3cm, wc_5cm
      real      :: wc_10cm, wc_15cm, wc_30cm, tlai, nee
      real      :: tmns_mlt
      real      :: avgstemp
      real      :: accum(ISOS)
      real      :: asthc, astlc
      character (len=10) :: subname

      ! methane output
      real      :: watr2sat, avgst_10cm, Com
      real      :: prev_mcprd(3)
      real      :: TI, Cr, Feh, CH4_prod, CH4_Ep, CH4_Ebl
      character (len=256) :: buffr
      double precision CH4_oxid
      real,save :: Eh = 300.0


      ! Initialize local variables
      accum  = 0.0
      subname = "simsom"

      if(wrtmeth .eq. 0) then
        wrtmeth = -1
        inquire(file='methane.out', exist=methopn)
        if(methopn) then
          wrtmeth = 73
          open(unit=wrtmeth,file='methane.out', ACCESS='APPEND',iostat=ostat)
          if(ostat .gt. 0) call abortrun('opening output file methane.out')
!          write(70,'(a10,1x,a8,4(1x,a12),1x,a16,5(1x,a12))') 'time',
!     &    'dayofyr','nit_N2O-N','dnit_N2O-N','dnit_N2-N', 'NO-N',
!     &    'CUM-N2O(gN/ha)','CUM-NO(gN/ha)'
        endif
      endif

! ... Initialize accumulators
      call cycle()

! ... N Fixation
! ... Does not take into account the effect of irrigation
      if (nsnfix .eq. 1 .and. nelem .ge. P) then

! ..... Compute mineral N:P ratio for N-Fixation (use suface layer only)
! ..... rnpml1 is used to control soil N-fixation using a regression
! ..... equation based on Kansas data. This ratio is used only if nelem = 2.
! ..... rnpml1 is flagged if either minerl(1,1) or minerl(1,2) is zero.

        rnpml1 = minerl(1,N)/minerl(1,P) * fsfunc(minerl(1,P),pslsrb,sorpmx)

! ..... Wet-dry fixation of nitrogen -- monthly flow
! ..... Atmospheric fixation is split between monthly dry fall and
! ..... wdfnp is the N:P ratio control function for non-symbiotic
! ..... soil fixation.
! ..... Both wdfnp and fxbiom are relative functions
! ..... which vary from 0 to 1.
! ..... wdfnp computed as a negative natural log of rnpml1
! ..... symnfx is the symbiotic N fixation by legumes derived from Cole and
! ..... Heil (1981) using data from Walker et al. 1959.
        if (rnpml1 .eq. 0) then
          wdfnp = 1.
        else
          wdfnp = min(1., ((-alog(rnpml1))/fxnpb)-.2)
        endif

! ..... The definitions of fxmca and fxmcb originally refered to water,
! ..... not biomass. (Alister 11/91)
        biof  = fxmca + (aglivc+stdedc+strucc(SRFC)) * fxmcb
        fxbiom = 1 - biof
        fxbiom = min(1.,fxbiom)
        if (wdfnp .lt. 0 .or. fxbiom .lt. 0 .or. stemp .lt. 7.5) then
          fwdfx = 0.0
        else
          fwdfx = wdfnp * fxbiom
        endif

! ..... Compute N-fixation for the month

! ..... Wet fall depending on the monthly precipitation (wdfxma)
!        wdfxma = wdfxa *  prcurr(month)/precipyr
! ..... Add an optional multiplier on N deposition, cak - 04/05/04
        wdfxma = baseNdep * ((precip(month) * precscalar(month)) / precipyr)
        if (Ninput .eq. 2 .or. Ninput .eq. 3) wdfxma = wdfxma * Nscalar(month)
        wdfxms = fxmxs * fwdfx
        wdfxm  = wdfxma + wdfxms

! ..... Compute annual N-fixation accumulators for atmosphere and soils
        frac_nh4 = 0.5
        frac_no3 = 0.5
        wdfxas = wdfxas + wdfxms
        wdfxaa = wdfxaa + wdfxma
        clyr = 1
        subname = 'simsom1a  '
        call update_npool(clyr, wdfxm, frac_nh4, frac_no3, ammonium, &
                          nitrate, subname)
        call flow(esrsnk(N),minerl(1,N),time,wdfxm)
! ..... nfixac should accumulate SYMBIOTIC N-fixation -mdh 11/16/01
!        nfixac = nfixac+wdfxm

! ... Monthly N-fixation based on annual parameters
      else
! ..... USE PRCURR/precipyr INSTEAD OF DT
!        wdfxms = wdfxs * prcurr(month)/precipyr
!        wdfxma = wdfxa * prcurr(month)/precipyr

        frac_nh4 = 0.5
        frac_no3 = 0.5
        wdfxms = wdfxs * ((precip(month) * precscalar(month))/precipyr)
! ..... Add an optional multiplier on N deposition, cak - 04/05/04
        wdfxma = baseNdep * ((precip(month) * precscalar(month)) / precipyr)
        if (Ninput .eq. 2 .or. Ninput .eq. 3) wdfxma = wdfxma * Nscalar(month)
        wdfxas = wdfxas + wdfxms
        wdfxaa = wdfxaa + wdfxma
        wdfxm = wdfxma + wdfxms
        clyr = 1
        subname = 'simsom2a  '
        call update_npool(clyr, wdfxm, frac_nh4, frac_no3, ammonium, &
                          nitrate, subname)
        call flow(esrsnk(N),minerl(1,N),time,wdfxm)
! ..... nfixac should accumulate SYMBIOTIC N-fixation -mdh 11/16/01
!        nfixac = nfixac + wdfxm
      endif

! ... Accumulate values for annual N deposition output variables,
! ... cak - 04/05/04
      wdfx = wdfx + wdfxma
      wdfxa = wdfxa + wdfxma

! ... Monthly atmospheric S deposition
! ... Note: irract is last month's irrigation at this point -mdh 12/9/96
      if (nelem .eq. S) then
!        satm = satmt * prcurr(month) / precipyr
        satm = satmt * ((precip(month) * precscalar(month)) / precipyr)
        satmac = satmac + satm
        if (irrcnt .ne. 0) then
          sirr = sirri * irract * 0.01
          sirrac = sirrac + sirr
        else
          sirr = 0
        endif
        call flow(esrsnk(S), minerl(SRFC,S), time, satm + sirr)
      endif

! ... BEGIN MONTHLY INITIALIZATION

      prevstream = 0 ! array

!      prevgromin = 0.0

! ... Reset the monthly accumulator
      call mthacc(month, agdefacsum, bgdefacsum)

      CO2resp = 0

! ... BEGIN DAILY LOOP...

      do 200 curday = frstdy(month), lstdy(month)
        doy = curday ! save curday in the timvar common for debugging

! ..... Initialize NPP and growth and maintenance respiration variables
        prev_mcprd = mcprd
        mcprd = 0.0
        mfprd = 0.0
        mrspdyflux  = 0.0
        cmrspdyflux = 0.0
        grspdyflux  = 0.0
        cgrspdyflux = 0.0
        fgrspdyflux = 0.0
        fmrspdyflux = 0.0

! ..... Call schedl to determine scheduling options for this day of the year
        call schedl(curday)

        tfrac = 1.0/dysimo(month)

! ..... Get a day's worth of weather from the weather file
        call getwth(curday, month, cyear, tempmax, tempmin, avgtemp, ppt,   &
                    solrad, rhumid, windsp, tavewk, petdly, fwloss,  &
                    sitlat, snow, tmn2m, tmx2m, tavemth, wthinput,   &
                    precscalar, tmaxscalar, tminscalar, srad, & !, usexdrvrs
                    leapyr, yrmatch, wstrtyr)

! ..... Set the value of tave for output to the *.bin file
        tave = tavemth

        ! Calculate an average daily temperature for use in growth equations.
        ! Making the daily minimum and maximum air temperature multipliers a
        ! function of day length captures the diurnal temperature variation
        if (daylength(curday) .lt. 12.0) then
          tmns_mlt = ((12.0 - daylength(curday)) * 3.0 + 12.0) / 24.0
        else
          tmns_mlt = ((12.0 - daylength(curday)) * 1.2 + 12.0) / 24.0
        endif
        tmns_mlt = max(0.05, min(0.95, tmns_mlt))

! From: Hartman,Melannie Thursday, October 17, 2013 12:03 PM
!  tavedly = tmxs_mlt * tempmax(curday) + tmns_mlt * tempmax(curday)
! This essentially sets tavedly to tempmax.
!  tavedly = tmxs_mlt * tempmax(curday) + tmns_mlt * tempmin(curday)
!          = tempmax - tempmax*tmns + tempmin*tmns = tempmax - (tempmax-tempmin)*tmns
        tavedly = tempmax(curday) - (tempmax(curday)-tempmin(curday)) * tmns_mlt

        irradd = 0
        ! Stop growing season auto-irrigation if nothing is growing
        ! assumes the other growth will will be zero in forsys and crpsys
        if(irrcnt .lt. 0  .and.  (crpgrw.eq.0 .and. forgrw.eq.0)) irrcnt = 0  ! clear the growing season irrigation

        if (irrcnt .ne. 0.) then
          ! get the irrigation if it is not monthly or if it is the week counter trips
          ! ignore  domirri on the monthly since we don't care if we get a double true
          if(.not. domirri  .or.  mod(irrcnt, 7).eq.0) irradd = irrigt(ppt(curday), petdly)
          ! decrement the irrigation count. For any finite period it decreases the time
          !  remaining. For growing season irrigation is the negative of the days in event
          irrcnt = irrcnt - 1;
        endif

        irract = irract + irradd      ! monthly irrigation    (plot1)
        rain   = rain + ppt(curday)     ! monthly precipitation (plot1)
        prcann = prcann + ppt(curday)   ! annual precipitation  (plot1)
        ppt(curday) = ppt(curday) + irradd  ! daily irrigation and precipitation
        annppt = annppt + ppt(curday)    ! yearly irrigation and precipitation

! ..... Set the flags for growth start before other events based on the
! ..... day of the year the event is scheduled to occur
        if (dofone .and. (foneday .eq. curday)) then
          forgrw = 1
          fgrwdys = 0
          call initprod(FORSYS, month)
        endif
        if ((dofrst .and. (frtcindx .lt. 3) .and. (frstday .eq. curday)) .or. &
           (frstschd .and. (frtcindx .eq. 3) .and. (frstday .le. curday) .and. &
           (soiltavewk .ge. tmpgerm))) then
          crpgrw = 1
          cgrwdys = 0
          thermunits = 0.0
          accumdd = .true.
          frstschd = .false.
          frstday = savefrstday
          plntkill = .false.
          call initprod(CRPSYS, month)
        endif
        if ((doplnt .and. (frtcindx .lt. 3) .and. (plntday .eq. curday)) .or.   &
            (plntschd .and. (frtcindx .ge. 4) .and. (plntday .le. curday) .and. &
            (soiltavewk .ge. tmpgerm))) then
          seedl = 1
          plntd = 1
          crpgrw = 1
          cgrwdys = 0
          falprc = 0
          prcfal = 0.0
          thermunits = 0.0
! ....... Do not set the degree day accumulator flag for crop types that
! ....... require a vernalization period, i.e. winter wheat, cak - 06/02/05
          if (frtcindx .ne. 6) accumdd = .true.
          plntschd = .false.
          plntday = saveplntday
          plntcnt = 0
          plntkill = .false.
          grnfill = .false.
          gwstress = 0.0
          grnfldys = 0
          grnhrvt = .false.
          call initprod(CRPSYS, month)
        endif

! ..... Check if days since planting (plntcnt) needs to be updated
        if (plntd .eq. 1 .and. stemp .ge. rtdtmp) then
          plntcnt = plntcnt + 1
        endif

! ..... Reset the reduction factor on nitrification rates due to
! ..... nitrification inhibitors as necessary, cak - 12/01/03
        if ((fertcnt .gt. 0) .and. (fertcnt .lt. (ninhtm * 7.0))) then
          fertcnt = fertcnt + 1
        else
          fertcnt = 0
          nreduce = 1.0
        endif

! ..... Save previous stream values from inorganic and organic leaching.  -mdh 1/97
        do iel = 1, nelem
          prevstream(iel+1) = stream(iel+1)
          prevstream(iel+5) = stream(iel+5)
        enddo
        strm5l = 0.0
        strm5u = 0.0

! ..... Track respiration over the daily timestep  !  both labeled and unlabeled
        autoresp1 = cautoresp + fautoresp
        hetresp1 = mt1c2 + mt2c2 + st1c2 + st2c2 + s11c2 + s12c2 + &
                   s21c2 + s22c2 + s3c2  + wd1c2 + wd2c2 + wd3c2
        mnrlrsp1 = mt2c2 + st2c2 + s12c2 + s22c2 + s3c2 + wd3c2
        oiresp1  = mt1c2 + st1c2 + s11c2
        oeresp1  = s21c2
        slitrsp1 = mt1c2 + st1c2 + s11c2 + s21c2

! ..... Set frlech to leaching fraction vek june90
! ..... Compute normal value for frlech.  Recompute in flood routine
! ..... if flooding occurs.

        texeff = fleach(1) + fleach(2) * sand
        frlech(1:nelem) = texeff * fleach(3:nelem+2)
        if(nelem >= P) frlech(P) = frlech(P) * fsfunc(minerl(SRFC,P), pslsrb, sorpmx)

! ..... Determine co2 effect on transpiration, pass to dailymoist
        if (cursys .eq. SAVSYS) then
          if (aglivc + rleavc .eq. 0.0) then
            co2val = 1.0
          else
            co2val = (co2ctr(CRPSYS)*aglivc + co2ctr(FORSYS)*rleavc) / &
                     (aglivc + rleavc)
          endif
        else
          co2val = co2ctr(cursys)
        endif

        call dailymoist(curday, agdefacsum, bgdefacsum, bgwfunc, &
                        frlech, co2val, saccum, evapdly, intrcpt, melt, &
                        outflow, petdly, runoffdly, sublim, trandly, &
                        avgstemp, watr2sat, prev_mcprd(2) + prev_mcprd(3), &
                        Com, avgst_10cm, TI, Cr, Eh, Feh, &
                        CH4_prod, CH4_Ep, CH4_Ebl, CH4_oxid)

! ..... Compute potential production
        call potprod(trbasl, tavedly, petdly, tfrac, tavemth, curday, srad(curday))

! ..... Volatilization loss of nitrogen as a function of
! ..... gross mineralization
! ..... Compute the grosmin which occured since the last time step
!        grominwk = gromin(1) - prevgromin
!        prevgromin = gromin(1)
        volgm = 0.0
!        volgm = vlossg*gromin(1)
!        volgm = vlossg*grominwk
!        minerl(SRFC,N) = minerl(SRFC,N) - volgm
!        write(*,*) 'volgma (volatilization) = ', volgma
!        esrsnk(N) = esrsnk(N) + volgm

! ..... Soil erosion
        if (doerod .and. (erodday .eq. curday) .or. &
           (erodcnt .gt. 0. .and. erodcnt .lt. 31)) then
          call erosn(psloss*tfrac,bulkd,edepth,enrich,lhzci,lhze,nelem,aceqcnt)
          if (erodday .eq. curday) then
            erodcnt = 0
          endif
          erodcnt = erodcnt + 1
        else
          scloss = 0.0
          erodcnt = 0
        endif

! ..... Fertilization option
! ..... This code has been moved to dailymoist, cak - 04/17/03

        ! moved the omad addition here instead of applying it in both the
        ! crop and tree routines independantly
        ! Allow for two types of organic matter addition.  Litter, for example
        ! wheat straw, is added to the structural and metabolic pools by the
        ! partit subroutine.  Partially decomposed organic matter, for example
        ! compost, is added directly into the slow pool (som2c).  cak - 11/16/2005
        if (doomad) then
          asthc = astgc*OMADscalar(month) * omadtyp
          if(omadtyp .gt. 0) then
            ! add the humas fraction directly som2.  cak - 11/16/2005
            call csched(asthc, astlbl, 1.0, csrsnk(UNLABL), som2ci(SRFC,UNLABL), &
                        csrsnk(LABELD), som2ci(SRFC,LABELD), 1.0, accum)

            cinput = cinput + asthc ! add decomposed omad into Carbon inputs

            ! now do the associated N
            do iel = 1, nelem
              call flow(esrsnk(iel), som2e(SRFC,iel), time, asthc*astrec(iel))
            end do
          endif


          astlc = astgc*OMADscalar(month) - asthc
          if(astlc .gt. 0) then
            ! partition the litter
            call partit(astlc, astrec,1, csrsnk, esrsnk, astlig, astlbl)
          endif

          ! Update OMAD accumulator output variables, cak - 07/13/2006
          asthc  = astgc*OMADscalar(month)
          omadac = omadac + asthc
          omadmth(month) = omadmth(month)  + asthc
          omadtot = omadtot + asthc
          hrvomadc = hrvomadc + asthc

          omadae(1:nelem) = omadae(1:nelem) + asthc * astrec(1:nelem)
          omadmte(month,1:nelem) = omadmte(month,1:nelem) + asthc*astrec(1:nelem)
          omaetot(1:nelem) = omaetot(1:nelem) + asthc*astrec(1:nelem)
          hrvomade(1:nelem) = hrvomade(1:nelem) + asthc*astrec(1:nelem)

          doomad = .FALSE. ! clear the flag now that the deed is done
        endif


        ! Available nutrients
        ! tminrl is the total amount of each element available in mineral form.

        ! Plants can only uptake from a layer with nutrients; don't sum negative values.
        tminrl(1:nelem) = sum(minerl(1:nlayer,1:nelem), 1, (minerl(1:nlayer,1:nelem) .gt. 0.))
        if (nelem .ge. P) tminrl(P) = tminrl(P) * fsfunc(minerl(SRFC,P), pslsrb, sorpmx)

!        write(*,*) 'SIMSOM: available nutrients = ', tminrl(1)

! ..... Compute the fraction of labile (non-sorbed) P in the surface
! ..... layer available to plants
        ! converted to do both crop and tree availability  KLK 30Jan13
        favail(2,:) = max(favail(4,:), min(favail(4,:) + &
                      minerl(SRFC,N)*(favail(5,:) - favail(4,:)) / favail(6,:), &
                      favail(5,:)))

! ..... Add to fallow rain
        if (falprc .eq. 1) then
!          prcfal = prcfal + prcurr(month)
          prcfal = prcfal + ppt(curday)
        endif

! ..... Call the producer submodels

! ..... Determine if daylength is increasing or decreasing
        if (daylength(curday) .lt. dayhrs) then
          hrsinc = .FALSE.
        else if (daylength(curday) .gt. dayhrs) then
          hrsinc = .TRUE.
        endif

! ..... Save the today's day length for use in tomorrow's calculation
        dayhrs = daylength(curday)

! ..... Crop and Forest removal options - moved here from CROP
! ..... and TREES so the mineral pools are not radically changed
! ..... between the growth routines. - rm 7/94

        ! senescence occurs once a month, at end of month.  -mdh 2/4/98
        ! moved dosene calculation to simsom;     KLK 15 Nov 2012
        ! again it controls all senescence actions effecting both death and growth
        if(senmschd .and. ( &
            (frtcindx .lt. 3  .and.  seneday .eq. curday) &             ! Non-GDD model
              .or. &
            (frtcindx .ge. 3 .and. &            ! Growing Degree Day submodel
              ! ((seneday .le. curday .and. thermunits .ge. ddbase) .or. & ! thermal limit
              ! changed (seneday <= curday) to (thermunits > 0) since the former
              ! breaks on a winter crop  KLK 12/04/2015
              ((thermunits .gt. 0. .and. thermunits .ge. ddbase) .or. & ! thermal limit
               accumdd .and. plntkill)) )) then                         ! Killing frost
          dosene = .true.
          senecnt = 1
          senmschd = .false.
        endif

        afrtup(1:nelem) = 0.  ! clear the daily automatic fertilizer
        if (cursys .eq. CRPSYS) then
          call crop(time, bgwfunc, tfrac, tavedly, curday, avgstemp)
          if ((dofire(CRPSYS) .and. (fireday .eq. curday)) .or. &
              (dograz .and. (grazday .eq. curday))) then
            ! if (dofire(CRPSYS)) firecnt = 1
            if (dograz) grazcnt = 1
            ! Burning of aboveground live, standing dead, and litter layer or grazing
            call grem()
          endif

        else if (cursys .eq. FORSYS) then
          call trees(bgwfunc, tavewk, tfrac, tavedly, avgstemp)
          if (dotrem .and. (tremday .eq. curday)) then
            ! Burning of live wood or cutting events
            call frem()
          endif
          if (dofire(FORSYS) .and. (fireday .eq. curday)) then
!            firecnt = 1
            ! Burning of dead wood and litter layer, cak - 08/23/02
            call grem()
          endif

        else if (cursys .eq. SAVSYS) then
          call crop(time, bgwfunc, tfrac, tavedly, curday, avgstemp)
          call trees(bgwfunc, tavewk, tfrac, tavedly, avgstemp)
          if (dotrem .and. (tremday .eq. curday))  call frem() ! tree burn or cutting
          if ((dofire(SAVSYS) .and. (fireday .eq. curday)) .or. &
              (dograz .and. (grazday .eq. curday))) then
!            if (dofire(SAVSYS)) firecnt = 1
            if (dograz)  grazcnt = 1
! ......... Burning of aboveground live, standing dead, litter
! ......... layer, and dead wood or grazing
            call grem()
          endif
        endif

        ! Cultivation
        ! moved cultiv call here so tree's, orchards, can be cultivated to remove
        ! leaf/fruit litter
        if (docult .and. (cultday .eq. curday)) call cultiv(curday, pltlig)
        ! zombie control; do we need to do the same for tree removal?

        if(sum(afrtup) .gt.0 .and. afue .gt. 0) then
          ! schedule the flows that add the Automatic fertilizer to the mineral pool
          afrtup(1:nelem)= afrtup(1:nelem)/afue
          clyr = SRFC
          call update_npool(clyr, afrtup(1), frac_nh4_fert, &
                           frac_no3_fert, ammonium, nitrate, subname)
          do iel = 1, nelem
            call flow(esrsnk(iel),minerl(SRFC,iel),time,afrtup(iel))
            ! Automatic fertilizer added to mineral pool
            fertot(2,iel)      = fertot(2,iel) + afrtup(iel)
            fertac(2,iel)      = fertac(2,iel) + afrtup(iel)
            fertmth(month,iel) = fertmth(month,iel) + afrtup(iel)
          end do
        endif

! ..... Update state variables and accumulators and sum carbon isotopes.
        call flowup(time)
        call sumcar

! ..... Track respiration over the daily timestep
        autoresp2   = cautoresp + fautoresp
        newautoresp = autoresp2 - autoresp1
        hetresp2    = mt1c2 + mt2c2 + st1c2 + st2c2 + s11c2 + s12c2 + &
                      s21c2 + s22c2 + s3c2  + wd1c2 + wd2c2 + wd3c2
        newhetresp = hetresp2 - hetresp1
        dhresp = sum(newhetresp)

        mnrlrsp2 = mt2c2 + st2c2 + s12c2 + s22c2 + s3c2  + wd3c2
        newmnrlrsp = mnrlrsp2 - mnrlrsp1
        dmnrlrsp = sum(newmnrlrsp)

        oiresp2 = mt1c2 + st1c2 + s11c2
        newoiresp = oiresp2 - oiresp1
        doiresp = sum(newoiresp)

        oeresp2 = s21c2
        newoeresp = oeresp2 - oeresp1
        doeresp = sum(newoeresp)

        slitrsp2 = mt1c2 + st1c2 + s11c2 + s21c2
        newslitrsp = slitrsp2 - slitrsp1
        dslitrsp = sum(newslitrsp)

! ..... Harvest may be performed after updating flows.  Put here for
! ..... consistency with the Savanna model - moved calls to flowup,
! ..... sumcar and harvst from CROP routine to here. -rm 7/94
!...... Non growing degree day submodel
        if ((dohrvt .and. (frtcindx .le. 2) .and. (hrvtday .eq. curday)) .or. &
           ! Perennial and non-grain filling annuals growing degree day submodel
           (harvschd .and. &
           ((frtcindx .eq. 3) .or. (frtcindx .eq. 4)) .and. &
           (hrvtday .le. curday) .and. (thermunits .ge. ddbase)) .or. &
           ! Grain filling annual growing degree day submodel
           (harvschd .and. (frtcindx .gt. 4) .and. &
           (hrvtday .le. curday) .and. (grnhrvt)) .or. &
           ! Killing frost, growing degree day submodel
           (harvschd .and. accumdd .and. (frtcindx .ge. 3) .and. plntkill) .or. &
           ! Force a harvest for FRTCINDX 3-5 if the day length is less
           ! than 12 hours and decreasing, cak - 05/19/2008
           (harvschd .and. &
            (((frtcindx .eq. 3) .or. (frtcindx .eq. 4) .or.  &
             (frtcindx .eq. 5)) .and. (dayhrs .lt. 12) .and. &
             accumdd .and. (.not. hrsinc)))) then
          call harvst(month, pltlig, curday)
          plntd = 0
          falprc = 1
          prcfal = 0.0
          harmth = 1
          harvschd = .false.

          ! stop the zombie growth; NO live C after harvest implies a LAST;
          ! copies 6/12/2001 century change KLK 6/07/2018
          if(aglivc .le. 1.0e-5 .and. (bglivcj + bglivcm) .le. 1.0e-5) then
            crpgrw = 0
            dolast = .true.
            lastday = curday
          endif

          ! This is also a LAST event for the grain filling annuals for
          ! the growing degree day implementation
          if (frtcindx .ge. 4) then
            crpgrw = 0
            dolast = .true.
            lastday = curday
            cgrwdys = 0
            thermunits = 0.0
            accumdd = .false.
            plntcnt = 0
            grnfill = .false.
            gwstress = 0.0
            grnfldys = 0
            grnhrvt = .false.
            call inprac(CRPSYS)
          endif
        endif

        ! Forest harvest
        ! This harvest location retains convention that grazing, fire, and forest removal,
        ! occur BEFORE harvest
        ! SHOULD MAKE FIRE CONSISTENT and look at that convention

        ! tree harvest;  do if flag is set and one of these conditions
        ! 1. not growing  or   2. forest last
        if (dothrv .and. (doflst .or. .not. dofngw)) then
          call trehrvst(tfhrvfr, fnhrvc, fnhrve)
          dothrv = .false. ! reset tree harvest flag to indicate it is complete
        endif

! ..... Set the flags for growth end after other events based on the
! ..... day of the year the event is scheduled to occur
        if (doflst .and. (flstday .eq. curday)) then
          forgrw = 0
          fgrwdys = 0
          call inprac(FORSYS)
        endif
        if (dolast .and. (lastday .eq. curday)) then
          crpgrw = 0
          thermunits = 0.0
          accumdd = .false.
          cgrwdys = 0
          plntcnt = 0
          call inprac(CRPSYS)
          ! For an annual plant reset the rooting depth back to 1 on a LAST event
          if ((frtcindx .eq. 2) .or. (frtcindx .ge. 4)) then
            claypg = 1
            ! If this is an annual plant flow the carbohydrate storage pool to the csrsnk
            ! call flow directly to clear stored C;
            ! eliminates a balance bug when input stored C doesn't match cisofr and simpler
            ! than determining correct ratio so csched can invert calculation KLK 2018/10/20
            call flow(carbostg(CRPSYS,UNLABL),csrsnk(UNLABL),time,carbostg(CRPSYS,UNLABL))
            call flow(carbostg(CRPSYS,LABELD),csrsnk(LABELD),time,carbostg(CRPSYS,LABELD))

           ! crpcarbostg = carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) !total stored C
           ! call csched(crpcarbostg,carbostg(CRPSYS,LABELD),crpcarbostg, & !flow the pools, (total, labeled, total) not cisofr
           !             carbostg(CRPSYS,UNLABL),csrsnk(UNLABL), &
           !             carbostg(CRPSYS,LABELD),csrsnk(LABELD), &
           !             1.0,accum)

          endif
        endif

        ! Update state variables and accumulators and sum carbon isotopes.
        call flowup(time)
        call sumcar

        ! Accumulate leached C,N,P,S
        csrsnk(UNLABL) = csrsnk(UNLABL) + strm5u
        csrsnk(LABELD) = csrsnk(LABELD) + strm5l
        stream(5) = stream(5) + strm5u + strm5l
        do iel = 1, nelem
          esrsnk(iel) = esrsnk(iel)+(stream(iel+1)-prevstream(iel+1))+ &
                       (stream(iel+5)-prevstream(iel+5))
        end do

! ..... Volatilization loss as a function of the mineral n which
! ..... remains after uptake by plants

        volex = 0.0
!        if (minerl(SRFC,N) .gt. 0.0) then
!          volex = vlosse*minerl(SRFC,N)*dt
!          volex = vlosse*minerl(SRFC,N)*dt*tfrac
!          minerl(SRFC,N) = minerl(SRFC,N) - volex
!          write(*,*) 'volex (volatilization) = ', volex
!          esrsnk(N) = esrsnk(N) + volex
!        endif

! ..... Volatilization
!        volgma = volgma+volgm
!        volexa = volexa+volex
!        volgac = volgac+volgm
!        voleac = voleac+volex

        if (time .ge. strplt) then
          call wrtbio(time, curday, aglivc, bglivcj, bglivcm, aglive(N), &
                      bglivej(N), bglivem(N), rleavc, frootcj, frootcm, &
                      fbrchc, rlwodc, crootc, h2ogef(1), h2ogef(2))
          call wrtmresp(time, curday, mrspdyflux(CRPSYS), &
                        mrspdyflux(FORSYS), cmrspdyflux(ABOVE), &
                        cmrspdyflux(BELOWJ), cmrspdyflux(BELOWM), &
                        fmrspdyflux(LEAF), fmrspdyflux(FROOTJ), &
                        fmrspdyflux(FROOTM), fmrspdyflux(FBRCH), &
                        fmrspdyflux(LWOOD), fmrspdyflux(CROOT), &
                        mcprd(ABOVE), mcprd(BELOWJ), mcprd(BELOWM), &
                        mfprd(LEAF), mfprd(FROOTJ), mfprd(FROOTM), &
                        mfprd(FBRCH), mfprd(LWOOD), mfprd(CROOT), &
                        carbostg(CRPSYS,UNLABL), carbostg(CRPSYS,LABELD), &
                        carbostg(FORSYS,UNLABL), carbostg(FORSYS,LABELD), &
                        mrspann(CRPSYS), mrspann(FORSYS), tavedly, &
                        mrspTempEffect(CRPSYS,SRFC), mrspTempEffect(CRPSYS,SOIL), &
                        mrspWaterEffect(CRPSYS), &
                        mrspTempEffect(FORSYS,SRFC), mrspTempEffect(FORSYS,SOIL), &
                        mrspWaterEffect(FORSYS))
          call wrtlivec(time, curday, aglivc, bglivcj, bglivcm, rleavc, &
                        frootcj, frootcm, fbrchc, rlwodc, crootc, frnutc)
          call wrtdeadc(time, curday, stdedc, metabc(1), strucc(1), &
                        wood1c, wood2c, wood3c)
          call wrtsoilc(time, curday, metabc(2), strucc(2), som1c(1), &
                        som1c(2), som2c(1), som2c(2), som3c)
          call wrtsysc(time, curday, &
                       aglivc + bglivcj + bglivcm + rleavc + frootcj + &
                       frootcm + fbrchc + rlwodc + crootc, &
                       stdedc + metabc(1) + strucc(1) + wood1c + &
                       wood2c + wood3c, &
                       metabc(2) + strucc(2) + som1c(1) + som1c(2) + &
                       som2c(1) + som2c(2) + som3c, &
                       aglivc + bglivcj + bglivcm + rleavc + frootcj + &
                       frootcm +  fbrchc + rlwodc + crootc + stdedc + &
                       metabc(1) + strucc(1) + wood1c + wood2c + &
                       wood3c + metabc(2) + strucc(2) + som1c(1) + &
                       som1c(2) + som2c(1) + som2c(2) + som3c, &
                       (st1c2(UNLABL) + st1c2(LABELD) + st2c2(UNLABL) + &
                        st2c2(LABELD) + mt1c2(UNLABL) + mt1c2(LABELD) + &
                        mt2c2(UNLABL) + mt2c2(LABELD) + s11c2(UNLABL) + &
                        s11c2(LABELD) + s12c2(UNLABL) + s12c2(LABELD) + &
                        s21c2(UNLABL) + s21c2(LABELD) + s22c2(UNLABL) + &
                        s22c2(LABELD) + s3c2(UNLABL)  + s3c2(LABELD)  + &
                        wd1c2(UNLABL) + wd1c2(LABELD) + wd2c2(UNLABL) + &
                        wd2c2(LABELD) + wd3c2(UNLABL) + &
                        wd3c2(LABELD)) - CO2resp)
          CO2resp = st1c2(UNLABL) + st1c2(LABELD) + st2c2(UNLABL) + &
                    st2c2(LABELD) + mt1c2(UNLABL) + mt1c2(LABELD) + &
                    mt2c2(UNLABL) + mt2c2(LABELD) + s11c2(UNLABL) + &
                    s11c2(LABELD) + s12c2(UNLABL) + s12c2(LABELD) + &
                    s21c2(UNLABL) + s21c2(LABELD) + s22c2(UNLABL) + &
                    s22c2(LABELD) + s3c2(UNLABL)  + s3c2(LABELD)  + &
                    wd1c2(UNLABL) + wd1c2(LABELD) + wd2c2(UNLABL) + &
                    wd2c2(LABELD) + wd3c2(UNLABL) + wd3c2(LABELD)
          call wrtgresp(time, curday, grspdyflux(CRPSYS), &
                        grspdyflux(FORSYS), cgrspdyflux(ABOVE), &
                        cgrspdyflux(BELOWJ), cgrspdyflux(BELOWM), &
                        fgrspdyflux(LEAF), fgrspdyflux(FROOTJ), &
                        fgrspdyflux(FROOTM), fgrspdyflux(FBRCH), &
                        fgrspdyflux(LWOOD), fgrspdyflux(CROOT), &
                        mcprd(ABOVE), mcprd(BELOWJ), mcprd(BELOWM), &
                        mfprd(LEAF), mfprd(FROOTJ), mfprd(FROOTM), &
                        mfprd(FBRCH), mfprd(LWOOD), mfprd(CROOT), &
                        carbostg(CRPSYS,UNLABL), carbostg(CRPSYS,LABELD), &
                        carbostg(FORSYS,UNLABL), carbostg(FORSYS,LABELD), &
                        grspann(CRPSYS), grspann(FORSYS))

! ....... Calculate daily output values for the dels.out file
          dcmresp = sum(cmrspdyflux)
          dfmresp = sum(fmrspdyflux)
          dmresp  = dcmresp + dfmresp
          dcgresp = sum(cgrspdyflux)
          dfgresp = sum(fgrspdyflux)
          dgresp = dcgresp + dfgresp

          dcrtjresp = cmrspdyflux(BELOWJ) + cgrspdyflux(BELOWJ)
          dcrtmresp = cmrspdyflux(BELOWM) + cgrspdyflux(BELOWM)
          dfrtjresp = fmrspdyflux(FROOTJ) + fgrspdyflux(FROOTJ)
          dfrtmresp = fmrspdyflux(FROOTM) + fgrspdyflux(FROOTM)
          dfrtcresp = fmrspdyflux(CROOT)  + fgrspdyflux(CROOT)
          dcrtresp = dcrtjresp + dcrtmresp
          dfrtresp = dfrtjresp + dfrtmresp + dfrtcresp
          dsresp = dcrtresp + dfrtresp + dhresp

! ....... Calculate fraction of crop/forest carbon that is labeled
          crplblstg = 0.0
          forlblstg = 0.0
          if (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) .gt. 0.0) then
            crplblstg = carbostg(CRPSYS,LABELD) / &
                     (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD))
          endif
          if (carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD) .gt. 0.0) then
            forlblstg = carbostg(FORSYS,LABELD) / &
                     (carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD))
          endif

! ....... Calculate the daily delta 13C/14C values for ouput
          if (labtyp .eq. 1  .or.  labtyp .eq. 2) then
            !  Delta 14C output  (1)   Delta 13C output (2)
            if ((newoiresp(LABELD) + newoiresp(UNLABL)) .gt. 0.) then
              ddeloi = delLCout(labtyp,newoiresp(LABELD), newoiresp(UNLABL), &
                                ddeloi)
            endif
            if ((newoeresp(LABELD) + newoeresp(UNLABL)) .gt. 0.0) then
              ddeloe = delLCout(labtyp,newoeresp(LABELD), newoeresp(UNLABL), &
                                ddeloe)
            endif
            if ((newslitrsp(LABELD) + newslitrsp(UNLABL)) .gt. 0.) then
              ddsrfclit = delLCout(labtyp,newslitrsp(LABELD), &
                                   newslitrsp(UNLABL), ddsrfclit)
            endif
            if ((newmnrlrsp(LABELD) + newmnrlrsp(UNLABL)) .gt. 0.) then
              ddsmnrl = delLCout(labtyp,newmnrlrsp(LABELD), &
                                 newmnrlrsp(UNLABL), ddsmnrl)
            endif
            if ((newhetresp(LABELD) + newhetresp(UNLABL)) .gt. 0.) then
              ddhetresp = delLCout(labtyp,newhetresp(LABELD), &
                                   newhetresp(UNLABL), ddhetresp)
            endif
            if ((sum(newhetresp) + dcrtresp + dfrtresp) .gt. 0.0) then
              ddsoilresp = delLCout(labtyp,(dcrtresp * crplblstg) + &
                                    (dfrtresp * forlblstg) + &
                                    newhetresp(LABELD), &
                                    (dcrtresp * (1.0 - crplblstg)) + &
                                    (dfrtresp * (1.0 - forlblstg)) + &
                                    newhetresp(UNLABL), ddsoilresp)
            endif
            if (dcmresp .gt. 0.0) then
              ddcmresp = delLCout(labtyp,(dcmresp * crplblstg), &
                                  (dcmresp * (1.0 - crplblstg)), ddcmresp)
            endif
            if (dfmresp .gt. 0.0) then
              ddfmresp = delLCout(labtyp,(dfmresp * forlblstg), &
                                  (dfmresp * (1.0 - forlblstg)), ddfmresp)
            endif
            if (dcgresp .gt. 0.0) then
              ddcgresp = delLCout(labtyp,(dcgresp * crplblstg), &
                                 (dcgresp * (1.0 - crplblstg)), ddcgresp)
            endif
            if (dfgresp .gt. 0.0) then
              ddfgresp = delLCout(labtyp,(dfgresp * forlblstg), &
                                 (dfgresp * (1.0 - forlblstg)), ddfgresp)
            endif
            if ((carbostg(CRPSYS,LABELD)+carbostg(CRPSYS,UNLABL)) .gt. 0.0) then
              ddccarbostg = delLCout(labtyp,carbostg(CRPSYS,LABELD), &
                                   carbostg(CRPSYS,UNLABL), ddccarbostg)
            endif
            if ((carbostg(FORSYS,LABELD)+carbostg(FORSYS,UNLABL)) .gt. 0.0) then
              ddfcarbostg = delLCout(labtyp,carbostg(FORSYS,LABELD), &
                                   carbostg(FORSYS,UNLABL), ddfcarbostg)
            endif
            call wrtdels(time, curday, ddeloi, ddeloe, ddsrfclit, &
                         ddsmnrl, ddhetresp, ddsoilresp, ddcmresp, &
                         ddfmresp, ddcgresp, ddfgresp, ddccarbostg, &
                         ddfcarbostg, doiresp, doeresp, dslitrsp, &
                         dmnrlrsp, dhresp, dcrtjresp, dcrtmresp, &
                         dfrtjresp, dfrtmresp, dfrtcresp, dsresp, &
                         dmresp, dgresp)
          endif
        endif

! ..... Write to dc_sip.out
        if (time .ge. strplt) then
          npptot = sum(mcprd) + sum(mfprd)
          nee = npptot - dhresp
          tlai = (rleavc * 2.5) * btolai
          systemc = aglivc + bglivcj + bglivcm + rleavc + frootcj + &
                    frootcm + fbrchc + rlwodc + crootc + stdedc + &
                    wood1c + wood2c + wood3c + strucc(1) + metabc(1) + &
                    strucc(2) + metabc(2) + som1c(1) + som1c(2) + &
                    som2c(1) + som2c(2) + som3c
! ....... Since each layer of the soil profile for this exercise has
! ....... the same texture use the field capacity and wilting point of
! ....... the top Century soil layer to compute the water holding
! ....... capacity for the various soil widths
          wc_2cm = (afiel(1) - awilt(1)) * 2
          wc_3cm = (afiel(1) - awilt(1)) * 3
          wc_5cm = (afiel(1) - awilt(1)) * 5
          wc_10cm = (afiel(1) - awilt(1)) * 10
          wc_15cm = (afiel(1) - awilt(1)) * 15
          wc_30cm = (afiel(1) - awilt(1)) * 30
          call wrtdcsip(time, curday, trandly, evapdly, intrcpt, &
                        sublim, outflow, runoffdly, ppt(curday), saccum, &
                        melt, snow, snlq, petdly, stemp, wc_2cm, &
                        wc_3cm, wc_5cm, wc_10cm, wc_15cm, wc_30cm, &
                        dhresp, mcprd(ABOVE), mcprd(BELOWJ), mcprd(BELOWM), &
                        mfprd(LEAF), mfprd(FROOTJ), mfprd(FROOTM), mfprd(FBRCH), &
                        mfprd(LWOOD), mfprd(CROOT), npptot, nee, aglivc, &
                        bglivcj, bglivcm, rleavc, frootcj, frootcm, &
                        fbrchc, rlwodc, crootc, tlai, stdedc, wood1c, &
                        wood2c, wood3c, strucc(ABOVE), metabc(ABOVE), &
                        strucc(2), metabc(2), som1c(1), som1c(2), &
                        som2c(1), som2c(2), som3c, systemc)

        ! CH4 oxidation is converted from g/ha to g/m^2
        !call wrtmethane(time, curday, aglivc, bglivcj, bglivcm, &
        !                prev_mcprd1, prev_mcprd2, prev_mcprd3, &
        !                Com, ppt(curday) - watr2sat - irradd, &
        !                irradd, watr2sat, avgst_10cm, TI, Cr, Eh, &
        !                Feh, CH4_prod, CH4_Ep, CH4_Ebl, CH4_oxid/10000.0)
          if(wrtmeth.gt.0) then
            write(buffr,'(i6,i4)') cyear, curday
            buffr(1:6) = ADJUSTL(buffr(1:6))
            write(wrtmeth,'(a,13f12.4,f12.6,2f12.4,4f12.6)') buffr(:10), &
               aglivc, bglivcj, bglivcm, prev_mcprd, &
               Com, ppt(curday) - watr2sat - irradd, irradd, watr2sat, &
               avgst_10cm, TI, Cr, Eh, Feh, CH4_prod, CH4_Ep, &
               CH4_Ebl, CH4_oxid/10000.0
          endif

        endif
200   continue
! ... END DAILY LOOP

! ... If we are at the end of the year reset the curday back to 1
      if (month .eq. 12) then
        curday = 1
      endif

! ... Annual co2 accumulators (10/92)
      ast1c2 = ast1c2 + st1c2(UNLABL) + st1c2(LABELD)
      ast2c2 = ast2c2 + st2c2(UNLABL) + st2c2(LABELD)
      amt1c2 = amt1c2 + mt1c2(UNLABL) + mt1c2(LABELD)
      amt2c2 = amt2c2 + mt2c2(UNLABL) + mt2c2(LABELD)
      as11c2 = as11c2 + s11c2(UNLABL) + s11c2(LABELD)
      as12c2 = as12c2 + s12c2(UNLABL) + s12c2(LABELD)
      as21c2 = as21c2 + s21c2(UNLABL) + s21c2(LABELD)
      as22c2 = as22c2 + s22c2(UNLABL) + s22c2(LABELD)
      as3c2  = as3c2  + s3c2(UNLABL)  + s3c2(LABELD)

! ... Monthly output averaged for *.bin file
      stemp   = stempmth/dysimo(month)
      anerb   = anerb/dysimo(month)
      agdefacm(month) = agdefacsum/dysimo(month)
      bgdefacm(month) = bgdefacsum/dysimo(month)
      aminrl  = aminrl/(ntspm*dysimo(month))
      agdefac = agdefacm(month)
      bgdefac = bgdefacm(month)

      htran(month) = tran;         hpttr(month) = pttr
      !On Nov 15, 2012, at 1:55 PM, billp wrote:
      !I think that the index we use should be the total AET( evap plus tran) over pet
      !htran(month) = tran + evap;   hpttr(month) = pet

! ... Annual production accumulator
      cproda = cproda + cprodc + cprodf
      eproda = eproda + eprodc + eprodf

! ... Net Mineralization
      do iel = 1, nelem

! ..... Net mineralization for the mineralizing compartments
! ..... The structural component of litter and the wood compartments
! ..... are not mineralizers.  They should not be added into cmn or sumnrs.
        cmn = metmnr(SRFC,iel) + metmnr(SOIL,iel) + &
              s1mnr(SRFC,iel) + s1mnr(SOIL,iel) +   &
              s2mnr(SRFC,iel) + s2mnr(SOIL,iel) + s3mnr(iel)
        sumnrs(iel) = sumnrs(iel) + cmn

! ..... soilnm is net mineralization in the soil.
        soilnm(iel) = soilnm(iel) + s1mnr(SOIL,iel) + s2mnr(SOIL,iel) + &
                 s3mnr(iel) + metmnr(SOIL,iel) + strmnr(SOIL,iel) + w3mnr(iel)

! ..... Total net mineralization
        tnetmn(iel) = tnetmn(iel) + cmn + &
                      strmnr(SRFC,iel) + strmnr(SOIL,iel) + &
                      w1mnr(iel) + w2mnr(iel) + w3mnr(iel)
      enddo

! ... Add calculation for annet which is used in the N deposition equations
! ... in eachyr, cak - 06/25/02
! ... Compute annual actual evapotranspiration
      annet = annet + evap + tran

! ... Stream flow accumulators, cak - 04/08/03
      strmac = strmac + stream

! ... Calculate monthly respiration from decomposition for output
      if (month .eq. 1) then
        respmth = resp
        respsum = respmth
      else
        respmth = resp    - respsum
        respsum = respsum + respmth
      endif

! ... Calculate monthly autotrophic respiration for output
      if (month .eq. 1) then
        arspmth(1,:) = cautoresp
        arspmth(2,:) = fautoresp
        arspsum = arspmth
      else
        arspmth(1,:) = cautoresp(:) - arspsum(1,:)
        arspmth(2,:) = fautoresp(:) - arspsum(2,:)
        arspsum = arspsum + arspmth
      endif

! ... Compute output variables for printing or plotting.
      call savarp

      return
      end
