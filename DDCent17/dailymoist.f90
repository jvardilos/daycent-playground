
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine dailymoist(curday, agdefacsum, bgdefacsum, bgwfunc,  &
                            frlech, co2val, accum, evapdly, intrcpt,  &
                            melt, outflow, petdly, runoffdly, sublim, &
                            trandly, avgstemp, watr2sat, prev_bgprd,  &
                            Com, avgst_10cm, TI, Cr, Eh, Feh,     &
                            CH4_prod, CH4_Ep, CH4_Ebl, CH4_oxid)
      USE ISO_C_Binding
      use calflow;
      implicit none

      include 'cflows.inc'
      include 'const.inc'
      include 'doubles.inc'
      include 'dovars.inc'
      include 'fertil.inc' ! trace gas
      include 'jday.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'timvar.inc'
!      include 't0par.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'sitsoil.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

! ... FORMAL PARAMETERS
      integer :: curday  ! same as doy
      real    :: bgwfunc, agdefacsum, bgdefacsum
      real    :: frlech(MAXIEL)
      real    :: co2val
      real    :: accum, evapdly, intrcpt, melt, outflow, petdly
      real    :: runoffdly, sublim, trandly
      real    :: avgstemp
      real    :: watr2sat, prev_bgprd, Com, avgst_10cm, TI, Cr, Eh
      real    :: Feh, CH4_prod, CH4_Ep, CH4_Ebl
      double precision :: CH4_oxid

! ... This routine loops calls the water budget routine (h2oflux), the
! ... decomposition routine (decomp) and the trace gas routines at a
! ... daily timestep.
!     Note:  routine uses f90+ features and is no longer f77 compliant

! ... wfluxout[] - total net flux thru the bottom of layer each day (cm H2O)
! ...              (positive is downward, negative is upward)
! ... nitrate[]  - layers of the nitrogen pool (gN/m2)
! ... ammonium   - the ammonium pool, no layer structure (gN/m2)
! ... newminrl   - mineralization that has occurred in the current day (gN/m2)
! ... co2val     - CO2 effect on transpiration.   Added 8/14/98 -mdh

! ... LOCAL VARIABLES
      integer       :: ii
      integer       :: isdecid
      integer, save :: isagri = 0

      ! print status variables
      ! Fortan needs to test nflux.out, daily.out, summary.out exist
      integer, save :: unflux = 0, udaly = 0, usumry = 0 ! output logical units if open
      logical, save :: testout = .true. ! Check output files on first pass

      real          :: stream1
      real          :: amovdly(10)
      real          :: rprpet
      real          :: pottransp
      real          :: baseflow
      real          :: snlq1, snlq2
      real          :: wfluxout(MXSWLYR)
      real          :: co2_conc(MXSWLYR)
      real          :: newCO2
      real          :: critflow
      real          :: frlechd(MAXIEL)
      real          :: avgwfps
      real          :: tmns, tmxs
      real          :: croplai, treelai, totlai
      real          :: trandep
      real          :: CO2resp
      real          :: maxrwcf
      real          :: hwstress
      real          :: thermtemp, tmns_mlt, tmxs_mlt
      real          :: litrcrbn
      real          :: srfctemp
      real          :: soilresp2
      double precision, save   :: inipswc =0 ! daily initial profile soil water (for balance)
      double precision         :: deltaswc ! change in profile soil water (for balance)
      double precision         :: dN2lyr(MXSWLYR), dN2Olyr(MXSWLYR)
      double precision         :: newminrl, inorglch
      double precision, save   :: nflux_sum, nflux_sum2
      double precision         :: Nn2oflux, Dn2oflux, Dn2flux, NOflux
      double precision         :: nit_amt
      double precision, parameter   :: SCALE = 10000.0 ! convert gN/m^2 to gN/ha
      character (len=10), parameter :: subname = 'dailymoist'
      real :: dstr = 0   ! ,dmet,dsm1,dsm2,dsm3
      real :: efscltef;

! ... Fortran to C prototype

! ... FUNCTION DECLARATIONS
      real      ramp
      external  ramp


     ! write(*,'(a,i3,"/",i5,1p2e10.3)') 'ENTERING FUNCTION DAILYMOIST...',doy,cyear; call flushSTD(0)
! ... decodt is equivalent to 4 times a day (vs. 4 times a month)
! ... -mdh 8/31/94  (Formerly, decodt=dt/ntspm, monthly version)

      decodt = dt/(dysimo(month)*ntspm)

! ... ***********************************************************************

! ... Soil moisture

      if(testout) then
        inquire(file='nflux.out', exist=testout)
        if(testout) then
          unflux = 70
          open(unit=unflux,file='nflux.out', ACCESS='APPEND', IOSTAT=ii)
!          write(70,'(a10,1x,a8,4(1x,a12),1x,a16,5(1x,a12))') 'time',
!     &    'dayofyr','nit_N2O-N','dnit_N2O-N','dnit_N2-N', 'NO-N',
!     &    'CUM-N2O(gN/ha)','CUM-NO(gN/ha)'
          if (ii.gt.0) call abortrun('opening nflux.out in dailymoist')
        endif

        inquire(file='daily.out', exist=testout)
        if(testout) then
          udaly = 71
          open(unit=udaly,file='daily.out', ACCESS='APPEND', IOSTAT=ii)
!          write(80,'(a10,1x,a8,4(1x,a12),2(1x,a7),2(1x,a12))')
!     &    'time','dayofyr','PET(cm)','agdefac','bgdefac','srfctemp(C)',
!     &     'snow', 'snlq', 'thermunits', 'srad'
          if (ii.gt.0) call abortrun('opening daily.out in dailymoist')
        endif

        inquire(file='summary.out', exist=testout)
        if(testout) then
          usumry = 72
          open(unit=usumry,file='summary.out', ACCESS='APPEND', IOSTAT=ii)
!           write(90,'(a10,1x,a8,3(1x,a8),5(1x,a12))') 'time','dayofyr',
!     &    'tmax','tmin','ppt','N2Oflux','NOflux','CH4_oxid','NIT','CO2resp'
          if (ii.gt.0) call abortrun('opening summary.out in dailymoist')
        endif
        testout = .false. ! don't check next time
      endif

! ... This code was pulled from subroutine cycle. -mdh  8/31/94.

      cwstress = 0.0
      twstress = 0.0

      if (doy .eq. 1) then
        nflux_sum = 0.0
        nflux_sum2 = 0.0
      endif
      snlq1 = snlq
! ... Amount of water in the soil at the beginning of the day
      if(inipswc .eq. 0) inipswc = sum(swc(1:numlyrs+1))  ! sum(asmos(1:nlayer+1))
      newminrl = 0.0

      call fertopt(doy, month)

! ... Do not accumulate thermal units for a crop requiring a
! ... vernalization period until day length is increasing (Jan in
! ... the northern hemisphere, Jul in the southern hemisphere)
      if ((frtcindx .eq. 6) .and. (crpgrw .eq. 1) .and. &
           hrsinc .and. (.not. accumdd)) then
        accumdd = .true.
      endif

! ... Accumulate thermal units for the growing degree day implementation, cak - 04/17/03
      if (accumdd) then
!        thermunits = thermunits + max(0.0, avgtemp(doy) - basetemp)
! ..... Use the day length to calculate the thermal units, cak - 08/29/05
        if (daylength(doy) .lt. 12.0) then
          tmns_mlt = ((12.0 - daylength(doy)) * 3.0 + 12.0) / 24.0
        else
          tmns_mlt = ((12.0 - daylength(doy)) * 1.2 + 12.0) / 24.0
        endif
        tmns_mlt = min(0.95, tmns_mlt)
        tmns_mlt = max(0.05, tmns_mlt)
        tmxs_mlt = 1.0 - tmns_mlt
! ..... Set an upper limit on the calculation for accumulating growing
! ..... degree days such that when the maximum temperature for the day
! ..... is capped at basetemp(2) degrees C, cak - 05/21/2008
        thermtemp = tmxs_mlt * max(basetemp(2), tempmax(doy)) + &
                    tmns_mlt * tempmin(doy)
        thermunits = thermunits + max(0.0, thermtemp - basetemp(1))
! ..... Stop growth when the thermal units (growing degree days)
! ..... required to reach senesence or maturity for a non-grain
! ..... producting annual occur
! ..... Crop growth is stopped by a harvest event rather than the
! ..... growing degree day accumulator, cak - 05/19/2008
!        if ((thermunits .ge. ddbase) .and. &
!            ((frtcindx .eq. 3) .or. (frtcindx .eq. 4))) then
!          crpgrw = 0
!        endif
! ..... For a grain producing crop reaching ddbase starts the grain
! ..... filling period, cak 06/02/05
        if ((thermunits .ge. ddbase) .and. (frtcindx .ge. 5)) then
          if (.not. grnfill) then
            grnfill = .true.
            gwstress = 0.0
            grnfldys = 1
          endif

          !  Check if the grain filling is complete
          if (thermunits .lt. (ddbase + mnddhrv)) then
          ! Keep accumulating
          elseif (thermunits .ge. (ddbase + mxddhrv)) then
            ! Stop plant growth, set flag to trigger harvest event
            crpgrw = 0
            grnhrvt = .true.
            cgrwdys = 0
          else
            ! Use grain water stress term to trigger the harvest event
            hwstress = ramp(gwstress/grnfldys, 0.0, mnddhrv, 1.0, mxddhrv)
            if (hwstress .lt. (thermunits - ddbase)) then
              ! Keep accumulating
            else
              ! Stop plant growth, set flag to trigger harvest event
              crpgrw = 0
              grnhrvt = .true.
              cgrwdys = 0
            endif
          endif
          grnfldys = grnfldys + 1
        endif
      endif

      petann = petann + petdly

      amovdly = 0.0

! ... If the grass/crop plant functional type is annual calculate a
! ... dynamic value for claypg from 1 to claypg_const based on the
! ... number of days since planting and FRTC(3), claypg_const is the
! ... value of claypg as read for the current CROP option
      if ((cursys .eq. CRPSYS) .or. (cursys .eq. SAVSYS)) then
        if ((frtcindx .eq. 2) .or. (frtcindx .ge. 4)) then
          claypg = nint(ramp(real(plntcnt), 0.0, 1.0, frtc(3), &
                             real(claypg_const)))
        endif
      endif

! ... Calculate a dynamic value for nlaypg based on the crop and/or tree
! ... option used, cak - 01/29/03
      if (cursys .eq. CRPSYS) then
        nlaypg = claypg
      else if (cursys .eq. SAVSYS) then
! ..... For crops and grasses a leaf area of 1 = 100 grams of biomass
        croplai = aglivc * 2.5 * 0.01
        treelai = rleavc * 2.5 * btolai
        totlai = croplai + treelai
        if (totlai .gt. 0.0) then
          nlaypg = nint(line(treelai/totlai, 0.0, real(claypg), 1.0, real(tlaypg)))
        else
          nlaypg = min(claypg, tlaypg)
        endif
        if (nlaypg .lt. min(claypg, tlaypg)) then
          nlaypg = min(claypg, tlaypg)
        endif
        if (nlaypg .gt. max(claypg, tlaypg)) then
          nlaypg = max(claypg, tlaypg)
        endif
      else
! ..... This is a tree system
        nlaypg = tlaypg
      endif

! ... Calculate depth for transpiration, cak - 01/29/03
      trandep = sum(adep(1:nlaypg))

! ... Calculate biomass values to be used in soil surface temperature
! ... calculations, this code has been moved from the potprod
! ... subroutine, cak - 01/28/2010
      if (cursys .eq. FORSYS) then
        aglivb = rleavc * 2.5           ! Live biomass
        ! Second mod to remove effect of woodc -rm 1/91
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.0 ! Surface litter biomass
        stdead = 0.0                    ! Standing dead biomass
        woodb = (rlwodc + fbrchc) * 2.0 ! Wood biomass

      elseif (cursys .eq. SAVSYS) then
        aglivb = (rleavc + aglivc) * 2.5 ! Live biomass
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.0 ! Surface litter biomass
        stdead = stdedc * 2.5            ! Standing dead biomass
        woodb = (rlwodc + fbrchc) * 2.0  ! Wood biomass

      else
        aglivb = aglivc * 2.5            ! Live biomass
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.5 ! Surface litter biomass
        stdead = stdedc * 2.5            ! Standing dead biomass
        woodb = 0.0                      ! Wood biomass
      endif

      litrcrbn = som1c(SRFC) + som2c(SRFC) + metabc(SRFC) + strucc(SRFC)
      call watrflow(doy, month, nlayer, nlaypg, watertable, watrflag, &
                  avgtemp(doy), tempmin(doy), tempmax(doy), &
                  solrad(doy), rhumid(doy), windsp(doy), ppt(doy), &
                  aglivb, sfclit, stdead, rwcf, avh2o, asmos, snow, &
                  snlq, amovdly, petdly, evapdly, trandly, stream1, &
                  basef, pottransp, baseflow, accum, melt, intrcpt, &
                  outflow, tmelt, sublim, wfluxout, time, strplt, &
                  co2val, tmns, tmxs, runoffdly, trandep, soiltavewk, &
                  daylength(doy), woodb, elitst, pmxtmp, pmntmp, &
                  pmxbio, srfctemp, stsys, ststart, stamt, litrcrbn, &
                  watr2sat, srad(doy), sradadj(month), inipswc)

! ... Sum water added to keep the soil at saturation into the water input output variables
      !Ken check the new definitions of the precipitiation
      ppt(doy) = ppt(doy) + watr2sat
      annppt = annppt + watr2sat

! ... Accumulate the months srfctemp, returned via watrflow
      stempmth = stempmth + srfctemp

! ... Find the relative water content in the wettest soil layer to use in
! ... the calculation of water stress on potential growth, cak - 12/06/04
! ... Grass/crop water stress
      maxrwcf = maxval(rwcf(1:claypg))
      cwstress = cwstress + min(1.0, maxrwcf)
! ... If this is a grain crop and we are within the grain filling period
! ... calculate the water stress term for the harvest water stress calculation
      if (grnfill) then
        ! On Sep 5, 2013, at 12:12 PM, Parton,William wrote:
        ! I think that it would be good to use the same water stress equations for grain filling.
        gwstress = gwstress + 1.0/(1.0 + exp(cwscoef(2) * (cwscoef(1)-maxrwcf)))
!     &      1.0/(1.0+exp(12.*(0.533078-maxrwcf))) ! old stress parameters
      endif
! ... Tree water stress
      maxrwcf = maxval(rwcf(1:tlaypg))
      twstress = twstress + min(1.0, maxrwcf)

! ... If there is snow melting into the ground use the melt value returned
! ... from the watrflow subroutine to determine how much snow is melting
! ... into the ground, cak - 10/21/02
! ... ppt(doy) includes any irrigation
      if (melt .gt. 0.0) then
! ..... melt represents the amount of water draining into the soil when
! ..... there is snow on the ground, both precipitation and irrigation
! ..... amounts have been taken into account in the snow calculations,
! ..... cak - 12/13/02
        rprpet = melt / petdly
      else
        rprpet = (ppt(doy) + avh2o(3))/ petdly
      endif

! ... Accumulate daily stream flow, drainage from each layer, pet,
! ... evaporation, and transpiration by month
      stream(1) = stream(1) + stream1
      pet = pet + petdly
      evap = evap + evapdly + intrcpt + sublim
      tran = tran + trandly
      pttr = pttr + pottransp
      amov(1:nlayer) = amov(1:nlayer) + amovdly(1:nlayer)
! ... Accumulate runoff for month
      runoff = runoff + runoffdly
! ... Accumulate precipitation + irrigation for the month
      pptmonth = pptmonth + ppt(doy)

! ... Change CO2 respiration value passed to the trace gas submodel
! ... so that we are passing only soil respiration, cak - 10/22/2007
      newCO2 = sum(mt2c2) + sum(st2c2) + sum(s12c2) + sum(s22c2) + sum(s3c2)

      ! Methane production uses CO2 respiration from all crop/grass litter and SOM pools
      !  add surface to newCO2
      Com = newCO2 + sum(mt1c2) + sum(st1c2) + sum(s11c2) + sum(s21c2)

      dstr = sum(st2c2) ! pre decomposition structural
      sdco2sum = 0.0   ! sum of soil CO2 flow with cultivation
      ntdco2sm = 0.0  ! sum of soil C flow without cultivation

! ... *********************************************************************
! ... Decomposition Submodel
      !Moved the entire decomposition loop into the decomp routine   -klk April 2016
      ! Removed P and S chemistry from decomp to subroutine pschem.  -rm  6/91
      call decomp(decodt,decsys,amovdly,newminrl, bgwfunc, agdefacsum, bgdefacsum, &
           avgstemp, avgwfps, ppt(doy), rprpet,petdly, &
           time, cyear, month, doy, aceqcnt)

! ... *********************************************************************

! ... Update single precision occlud and secndy with double precision values cak - 03/20/02
      ! moved outside loop; float variables aren't used for decomposition KLK 11 Oct 15
      occlud = real(occlud_double)
      secndy = real(secndy_double)

! ..... Change CO2 respiration value passed to the trace gas submodel
      ! so that we are passing only soil respiration, cak - 10/22/2007
      soilresp2 = sum(mt2c2) + sum(st2c2) + sum(s12c2) + sum(s22c2) + sum(s3c2)
      ! Methane production uses CO2 respiration from all crop/grass litter and SOM pools
      !  add surface to newCO2
      Com = soilresp2 + sum(mt1c2) + sum(st1c2) + sum(s11c2) + sum(s21c2) - Com

      newCO2 = soilresp2 - newCO2 !   max(soilresp2 - newCO2, 0.000001)

      ! difference the structural accumulator to avoid having modify declig
      ! and tell the difference between soil litter and other pools.
      dstr = sum(st2c2) - dstr
! ... sdco2sum - soil decomposition CO2 flow sum. Needed to calculate the
!                effective decomposition effect.
! ... ntdco2sm - no-till soil decomposition CO2 flow sum
      sdco2sum = sdco2sum + dstr
      if(cltfac(4).gt.0) ntdco2sm = ntdco2sm + dstr/cltfac(4)
      !** Bill Parton requested to restrict the CO2 respiration effect
      !** on dentrification by capping ntdco2sm at 2.5 gC/m2/d. -mdh 11/12/2014
      !** ntdco2sm = min(2.5,ntdco2sm)

      ! ratio the total and no-till CO2 flows to get the weighted clteff
      if(ntdco2sm .gt. 0) then
        efscltef = sdco2sum / ntdco2sm
      else
        efscltef = 1.0
      endif

      newCO2 = max(sdco2sum, 0.000001)
      CO2resp = newCO2
      if (avgwfps .gt. 0.60) then
        newCO2 = newCO2 / bgwfunc
      endif

!      critflow = minlch/dysimo(month)
      critflow = minlch
      frlechd(1:nelem) = frlech(1:nelem)

! ... *********************************************************************

! ... Trace Gas Model

!      write(*,*) 'Time = ', time, ' doy = ', doy

! ... Are we running a decidious forest?
      if ((cursys .eq. FORSYS) .and. (decid .eq. 1)) then
        isdecid = 1
      else
        isdecid = 0
      endif
! ... Once cultivation, fertilization, or harvesting occurs in a system the
! ... agricultural effect needs to be "turned on".  Once this effect on
! ... methane oxidation has been invoked it stays in place.
      if (isagri .eq. 0) then
!        if (dofert .or. dohrvt) isagri = 1
        if (docult) isagri = 1
      endif

      call trace_gas_model(newminrl, ammonium, nitrate, texture, &
                           sand, silt, clay, afiel(1), bulkd, maxt, &
                           ppt(doy), snow, avgwfps, stormf, basef, &
                           frlechd, stream, inorglch, critflow, &
                           wfluxout, newCO2, co2_conc, efscltef, time, &
                           NOflux, Nn2oflux, Dn2oflux, Dn2flux, CH4_oxid, &
                           isdecid, isagri, aglivc, rleavc, btolai, &
                           crpstg(N), forstg(N), nit_amt, nreduce, &
                           cyear, doy, pHscalar(month), dN2lyr, dN2Olyr, &
                           esrsnk(N), minerl, nlayer, &
                           prev_bgprd, Com, avgst_10cm, TI, Cr, Eh, &     ! methane
                           Feh, CH4_prod, CH4_Ep, CH4_Ebl, watertable, &  ! methane
                           watrflag, bglivcj+bglivcm, tmxbio)             ! methane

! ... *********************************************************************
! ... Write to output files

      nflux_sum = nflux_sum+(Nn2oflux+Dn2oflux)*SCALE
      nflux_sum2 = nflux_sum2 + NOflux*SCALE

! ... Accumulate yearly trace gas output, cak - 09/23/02
      ! NOTE: The units on the following variables are gC/ha instead of gC/m2
      !       nit_amt_year, N2O_year, NO_year, N2_year,
      !       nit_amt_month, N2O_year, NO_year, N2_year
      N2O_year = N2O_year + (Nn2oflux+Dn2oflux)*SCALE
      NO_year = NO_year + NOflux*SCALE
      N2_year = N2_year + Dn2flux*SCALE
      nit_amt_year = nit_amt_year + nit_amt*SCALE
      CH4yrox = CH4yrox + CH4_oxid/SCALE ! CH4oxid_year  convert g/ha to g/m2
      CH4mnem = CH4mnem + CH4_Ep + CH4_Ebl
      CH4mnpr = CH4mnpr + CH4_prod

! ... Accumulate monthly trace gas output, cak - 05/14/04
      N2O_month = N2O_month + (Nn2oflux+Dn2oflux)*SCALE
      NO_month = NO_month + NOflux*SCALE
      N2_month = N2_month + Dn2flux*SCALE
      nit_amt_month = nit_amt_month + nit_amt*SCALE
      CH4mnox = CH4mnox + CH4_oxid/SCALE ! CH4oxid_month  convert g/ha to g/m2
      CH4yrem = CH4yrem + CH4_Ep + CH4_Ebl
      CH4yrpr = CH4yrpr + CH4_prod
! ... Growing season accumulator for N2O flux, cak - 06/06/2008
      n2oacc = n2oacc + (Nn2oflux+Dn2oflux)
      n2omth(month) = n2omth(month) + (Nn2oflux+Dn2oflux)

      if (time .ge. strplt) then

        if(unflux.gt.0) write(unflux, &
           '(f10.4,1x,i4,4(1x,f12.4),1x,f16.4,5(1x,f12.4))') time,doy, &
            Nn2oflux*SCALE, Dn2oflux*SCALE, Dn2flux*SCALE, NOflux*SCALE, &
            nflux_sum, nflux_sum2

        if(udaly .gt. 0) write(udaly, &
           '(f10.4,1x,i4,4(1x,f12.4),1x,2(f7.4),2(1x,f12.4))')time,doy, &
           petdly,agdefac,bgdefac,srfctemp,snow,snlq,thermunits,srad(doy)

        if(usumry .gt. 0) write(usumry, &
           '(f10.4,1x,i4,3(1x,f8.2),5(1x,f12.4))') time,doy, &
           tempmax(doy), tempmin(doy),ppt(doy), &
           (Nn2oflux+Dn2oflux)*SCALE, NOflux*SCALE, CH4_oxid, nit_amt*SCALE, &
          CO2resp*SCALE

        call wrtsoiln(time, doy, ammonium, nitrate)
        call wrtco2(time, doy, co2_conc)
        call wrtdn2lyr(time, doy, dN2lyr)
        call wrtdn2olyr(time, doy, dN2Olyr)
        call wrtwflux(time, doy, wfluxout)

        call wrtcflows(time, doy, som11tosom21, som12tosom22,    &
                       som12tosom3, som21tosom11, som21tosom22,     &
                       som22tosom12, som22tosom3, som3tosom12,      &
                       metc1tosom11, metc2tosom12, struc1tosom11,   &
                       struc1tosom21, struc2tosom12, struc2tosom22, &
                       wood1tosom11, wood1tosom21, wood2tosom11,    &
                       wood2tosom21, wood3tosom12, wood3tosom22)

      endif

! ....*********************************************************************

! ... Update state variables and accumulators and sum carbon isotopes
      call flowup(time)
      call sumcar

      ! Accelerate the equilibrium
      if(doaceq) call EqAclrat("AEQ") ! aceq

!      Moved bal_npool call to trace_gas_model where the nitrate profile is modified.
!      The minerl array should be rebalanced as soon as possible.

! ... *********************************************************************

! ... Report the water balnce at the end of the day.
      snlq2 = snlq
      deltaswc = inipswc                ! save the initial water content
      inipswc  =  sum(swc(1:numlyrs+1)) ! Get the current
      deltaswc = deltaswc - inipswc     ! now do the difference

      if (time .ge. strplt) then
        call watrbal(doy, time, ppt(doy), accum, melt, deltaswc, &
                     evapdly, trandly, sublim, intrcpt, &
                     outflow, snlq1, snlq2, snow, runoffdly)
      endif

! ... *********************************************************************

! ... If the minimum temperature for the day is lower than the tmpkill
! ... parameter for the crop and the crop has accumulated at least 1/2 of
! ... the base thermal units a killing frost has occurred
      if ((tempmin(doy) .le. tmpkill) .and. &
          (thermunits .ge. (ddbase/2.0) .and. (frtcindx .ge. 3))) then
        plntkill = .true.
      endif

      return
      contains
        include 'line.f'


        subroutine fertopt(doy, month)
          integer :: doy, month

! ... LOCAL VARIABLES
          real, save :: baseminrl(MAXIEL) =0.
          character (len=7), parameter :: subname = 'fertopt'
          integer :: iel, clyr
          real    :: fertiel

          ! Fertilization option
          ! moved into a subroutine so we can add Yao and Ram's fertilizer options
          ! without rewriting dailymoist.

          ! This code has been relocated from simsom so that fertilization can occur
          ! on a specific day, cak - 04/17/03

    !     Add any delayed release nitrogen to the ammonium pool


          ! Add an optional multiplier on feramt for N, cak - 04/05/04
          ! collapsed the following fertilizer code to modify N based on Ninput then do the
          ! sums rather than pasting 3 copies of the accumulator sums  KLK 4/2014
          if (dofert .and. doy .eq. fertday) then
           ! if(crfteff .eq. -1) then
           !    crfn = feramt ! controlled release fertilizer
           !    crfmewf = 0;
           ! endif

            do iel = 1, nelem
              fertiel = feramt(iel) ! cache addition so we can modify N with Nscalar
              ! capture the soil mineral when fertilization starts (net fert == 0)
              if(fertnet(iel) .eq. 0) baseminrl(iel) = sum(minerl(1:claypg, iel), &
                         MASK= minerl(1:claypg, iel) .gt. 0.0)

              if (iel .eq. N) then
                if (Ninput .eq. 1 .or. Ninput .eq. 3) fertiel = fertiel*Nscalar(month)
                clyr = SRFC
                call update_npool(clyr, fertiel, frac_nh4_fert, frac_no3_fert, &
                                  ammonium, nitrate, subname)
              endif
              esrsnk(iel)      = esrsnk(iel) - fertiel
              minerl(SRFC,iel) = minerl(SRFC,iel) + fertiel
              fertot(1,iel)    = fertot(1,iel) + fertiel
              fertac(1,iel)    = fertac(1,iel) + fertiel
              fertmth(month,iel) = fertmth(month,iel) + fertiel
              hrvfert(iel)     = hrvfert(iel) + fertiel ! sum fertilizer between harvest

              fertnet(iel) = fertnet(iel) + fertiel ! net fertilizer ! bminerl(iel) = sum(minerl(:,iel), 1)
            end do
            fertcnt = 1
            nreduce = ninhib

          elseif(sum(fertnet) .gt. 0) then
            do iel = 1, nelem
              fertiel = sum(minerl(1:claypg, iel), DIM= 1, &
                         MASK= minerl(1:claypg, iel) .gt. 0.0)
              if(fertiel < baseminrl(iel)) then
                fertnet(iel) = 0.0;
                baseminrl(iel) = 0.0;
              else if(fertiel < fertnet(iel)) then
                fertnet(iel) = fertiel - baseminrl(iel)
              endif
            end do
          endif
          return
        end subroutine fertopt
      end subroutine dailymoist
