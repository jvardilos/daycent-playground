
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine eachyr

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'zztim.inc'

! ... Perform tasks that only need to be done once a year.

! ... Function declarations
      real     fracis, ramp
      external fracis, ramp

! ... Local variables
      real     lfncmax, lfncmin, lfncon

! ... Added for savanna model (plot3.inc) BO

! ... Reset annual accumulators to zero
      call annacc
      fnhrvc = 0.0;     fnhrve = 0.0 ! make the fruit/nut harvest annual

! ... Set climate scalars as indicated, cak - 10/18/05
      if (wthinput .gt. 0) call climscale(time)

! ... Wet-dry fixation of N

! ... Determine annual precipitation and annual PET
      ! set precipyr to be the total average precip;
      ! Changed prcann to the actual yearly precipitation  KLK 11/2103
      ! prcann = sum(prcurr);   petann = sum(pevap)

      precipyr  = sum(precip * precscalar)
      prcann   = 0.
      agdefacm = -1.
      bgdefacm = -1.

! ... Scale the OMAD inputs as indicated, cak - 04/06/04
      if (OMADinput .gt. 0) call omadscale(time, OMADstart)

! ... Scale the N inputs as indicated, cak - 04/06/04
      if (Ninput .gt. 0) call nscale(time, Nstart)

! ... N fixation in atmosphere
      baseNdep = max(epnfa(INTCPT)+epnfa(SLOPE)*MIN(precipyr,80.0), 0.0)

!      wdfxs = epnfs(INTCPT)+epnfs(SLOPE)*MIN(precipyr,100.0)
! ... Now using annual ET in wdfxs calculation, cak - 02/21/02
! ... Use annual ET unless it is the first timestep
! ... No longer using the intercept in the calculation.
      if (annet .eq. 0.0) then
        wdfxs = epnfs(SLOPE)*MIN(precipyr,100.0)
      else
        wdfxs = epnfs(2) * (annet - epnfs(1))
      endif
! ... Reset annual accumulator for evapotranspiration
      annet = 0
      if (wdfxs .lt. 0.)  then
        wdfxs = 0.0
      endif

! ... This output varible represents only the non-symbiotic soil
! ... N-fixation, the atmospheric N deposition is added to this output
! ... variable in the simsom subroutine, cak - 04/05/04
!      wdfx = wdfxa+wdfxs
      wdfx = wdfxs

! ... Atmospheric S deposition
      satmt = max(0.0, satmos(1) + satmos(2)*precipyr)

! ... Determine what fraction of the carbon in new plant tissue is labeled
      if (labtyp .eq. 0) then
        cisofr = 0.0
        cisotf = 0.0
      elseif (labtyp .eq. 1) then
        cisofr = fracis(time,labyr)
        cisotf = cisofr
!      elseif (labtyp .eq. 2) then
! ..... cropin has set cisofr
! ..... treein has set cisotf
      endif

! ... Initialize co2 effects
      call co2eff(time)

! ... Implement pH shift as indicated, cak - 08/02/02
      if (phsys .gt. 0) call phshift(time)

! ... Added effect of co2 for forest; done here because not calcualted
! ... dynamically based on biomass like grassland/crop
! ... Direct CO2 effects only C/E ratio of leaves.
      ccefor(IMIN:IMAX, 1:FPARTS-1, 1:nelem) =
     &      cerfor(IMIN:IMAX, 1:FPARTS-1, 1:nelem)
      ccefor(IMIN:IMAX,1,1:nelem) =
     &      ccefor(IMIN:IMAX,1,1:nelem)+co2cce(FORSYS,:,1:nelem)
!      ccefor(IMIN,LEAF,1:nelem) = cerfor(IMIN,LEAF,1:nelem) *
!     &                        co2cce(FORSYS,IMIN,1:nelem)
!      ccefor(IMAX,LEAF,1:nelem) = cerfor(IMAX,LEAF,1:nelem) *
!     &                        co2cce(FORSYS,IMAX,1:nelem)
!
!      ccefor(IMIN,2:FPARTS-1,1:nelem) = cerfor(IMIN,2:FPARTS-1,1:nelem)
!      ccefor(IMAX,2:FPARTS-1,1:nelem) = cerfor(IMAX,2:FPARTS-1,1:nelem)

! ... Calculate leaf death rate multiplier for continuous forests 11/20/92
! ... Initialize LDRMLT to 1.0
      ldrmlt = 1.0

! ... Change leaf death rate multiplier if you have floating C/E ratios.
      if (ccefor(IMIN,LEAF,N) .ne. ccefor(IMAX,LEAF,N)) then
        if (rleavc .gt. 0) then
          lfncon = rleave(N) / rleavc
          lfncmin = 1 / ccefor(IMIN,LEAF,N)
          lfncmax = 1 / ccefor(IMAX,LEAF,N)
          ldrmlt = 1 + (maxldr - 1) *
     &             (lfncon - lfncmin) / (lfncmax - lfncmin)
        endif
      endif

      if (cursys .ne. FORSYS) then
        ! Determine the lignin fraction of plant residue added this year.
        call cmplig(cursys,fligni,wdlig)
      endif

! ... Compute SITPOT as a function of the annual precipitation, cak - 05/02/03
! ... sitpot_m is the SITPOT parameter value as read from the tree.100 file
! ... for the current tree, cak - 11/21/01
      sitpot = ramp(precipyr, 20.0, 1500.0, 90.0, 3250.0)
      sitpot = sitpot * sitpot_m

! ... Add code to reset the frstday or plntday values as necessary to allow a
! ... grass/crop that grows over the December/January boundary using the
! ... growing degree day implementation that has not yet started to grow to
! ... start growth once the growth criterion have been met, cak - 06/19/03
      if (frstschd .and. frtcindx .eq. 3 .and. frstday .gt. 1) then
        frstday = 1
      endif
      if (plntschd .and. frtcindx .ge. 4 .and. plntday .gt. 1) then
        plntday = 1
      endif

      return
      end
