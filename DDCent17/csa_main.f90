
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


!                           DISCLAIMER
!
!        Neither the Great Plains System Research Unit - USDA (GPSR) nor
!     Colorado State University (CSU) nor any of their employees make
!     any warranty or assumes any legal liability or responsibility for
!     the accuracy, completeness, or usefulness of any information,
!     apparatus, product, or process disclosed, or represents that its
!     use would not infringe privately owned rights.  Reference to any
!     special commercial products, process, or service by tradename,
!     trademark, manufacturer, or otherwise, does not necessarily
!     constitute or imply endorsement, recommendation, or favoring by
!     the GPSR or CSU.  The views and opinions of the authors do not
!     necessarily state or reflect those of GPSR or CSU and shall not
!     be used for advertising or product endorsement.

      program main

! ... Century Soil Organic Matter Model
!     Simulation of carbon, nitrogen, phosphorous, and sulfur cycling
!     As of Dec. 1991, uses a 1 month time step
!     Project - Soil Fertility in the Great Plains
!     Modeler - Bill Parton
!     Programmers - Vicki Kirchner, Becky McKeown, Laura Harding,
!                   Melannie Hartman

!       Modifications
!  9/2010  K. Killian
!     Added code to report run date and time for debugging

! ... State variables and flows are grams/m2.

      use calflow;

      implicit none
      include 'cflows.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'jday.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

! ...               (unit 1) = plot/print file used by modaid (unformatted)
! ... <site>.100    (unit 7) = parameter values and initial values for
! ...                          state variables; see subroutine sitein.
! ... fix.100       (unit 8) = fixed parameter values values for
! ...                          state variables; see subroutine fixin.
! ...               (unit 9) = a file of weather data read in subroutines
! ...                          wthini, weathr
! ... c14data      (unit 10) = a data file specifying years in which labeled
! ...                          carbon is added via plant growth and what
! ...                          fraction of the growth is labeled.
! ... nscale.dat   (unit 20) = a data file specifying years in which N input
! ...                          scalars are used and the scalar values
! ... omadscale.dat(unit 30) = a data file specifying years in which organic
! ...                          matter input scalars are used and the scalar
! ...                          values
! ... phscale.dat  (unit 40) = a data file specifying years in which pH
! ...                          scalars are used and the scalar values
! ... precscale.dat(unit 50) = a data file specifying years in which
! ...                          precipitation scalars are used and the scalar
! ...                          values, precipitation scalar are multipliers
! ... tmaxscale.dat(unit 55) = a data file specifying years in which
! ...                          maximum temperature scalars are used and the
! ...                          scalar values, maximum temperature scalar are
! ...                          addends
! ... tminscale.dat(unit 60) = a data file specifying years in which
! ...                          minimum temperature scalars are used and the
! ...                          scalar values, minimum temperature scalar are
! ...                          addends
! ... nflux.out    (unit 70) = N2/N2O fluxes computed by Trace Gas Model
! ... daily.out    (unit 80) = pet, defac, stemp, and snowpack water content
! ...                          computed by Trace Gas Model
! ... summary.out  (unit 90) = tmax, tmin, prec, N2O flux, NO flux, CH4, and
! ...                          gross nitrification computed by Trace Gas Model

! ... If you're getting floating point errors mentioned after you exit
! ... Century, uncomment the following lines, recompile, run Century
! ... in dbx with the 'catch FPE' option to find the offending code.
! ... You can also run Century outside of dbx, in which case you will
! ... get messages on your screen giving you an fpe error code (see
! ... the Floating Point Programmer's Guide, p20) and a not-very-
! ... useful-to-3rd-or-4th-generation-language-programmers location.
! ... The error handler 'mysigal' is a Fortran callable C routine
! ... written by Martin Fowler; it can be replaced by any user written
! ... handler or any of several library handlers, the most useful
! ... probably being SIGFPE_ABORT.  The calls to ieee_handler won't
! ... compile using poa's binaries.

!      external mysigal
!      ieeer=ieee_handler('set','invalid',SIGFPE_ABORT)
!      ieeer=ieee_handler('set','division',mysigal)
!      ieeer=ieee_handler('set','overflow',mysigal)
!      ieeer=ieee_handler('set','underflow',SIGFPE_ABORT)

! ... You probably won't want to uncomment the following line; inexact
! ... floating point signals occur all over the place.

!      ieeer=ieee_handler('set','inexact',mysigal)

! ... Fortran to C prototype
      INTERFACE

        SUBROUTINE closefiles()
          !MS$ATTRIBUTES ALIAS:'_closefiles' :: closefiles
        END SUBROUTINE closefiles

        SUBROUTINE wrttgmonth(time, N2O_month, NO_month, N2_month, &
                              CH4mnox, nit_amt_month, pptmonth)
          !MS$ATTRIBUTES ALIAS:'_wrttgmonth' :: wrttgmonth
          REAL             time
          REAL             N2O_month
          REAL             NO_month
          REAL             N2_month
          REAL             CH4mnox
          DOUBLE PRECISION nit_amt_month
          REAL             pptmonth
        END SUBROUTINE wrttgmonth

        SUBROUTINE wrtyearsum(time, N2O_year, NO_year, N2_year, &
                              CH4yrox, nit_amt_year, annppt)
          !MS$ATTRIBUTES ALIAS:'_wrtyearsum' :: wrtyearsum
          REAL             time
          REAL             N2O_year
          REAL             NO_year
          REAL             N2_year
          REAL             CH4yrox
          DOUBLE PRECISION nit_amt_year
          REAL             annppt
        END SUBROUTINE wrtyearsum

        SUBROUTINE wrtyrcflows(time, asom11tosom21, asom12tosom22, &
                               asom12tosom3, asom21tosom11, &
                               asom21tosom22, asom22tosom12, &
                               asom22tosom3, asom3tosom12, &
                               ametc1tosom11, ametc2tosom12, &
                               astruc1tosom11, astruc1tosom21, &
                               astruc2tosom12, astruc2tosom22, &
                               awood1tosom11, awood1tosom21, &
                               awood2tosom11, awood2tosom21, &
                               awood3tosom12, awood3tosom22)
          !MS$ATTRIBUTES ALIAS:'_wrtyrcflows' :: wrtyrcflows
          REAL time
          REAL asom11tosom21
          REAL asom12tosom22
          REAL asom12tosom3
          REAL asom21tosom11
          REAL asom21tosom22
          REAL asom22tosom12
          REAL asom22tosom3
          REAL asom3tosom12
          REAL ametc1tosom11
          REAL ametc2tosom12
          REAL astruc1tosom11
          REAL astruc1tosom21
          REAL astruc2tosom12
          REAL astruc2tosom22
          REAL awood1tosom11
          REAL awood1tosom21
          REAL awood2tosom11
          REAL awood2tosom21
          REAL awood3tosom12
          REAL awood3tosom22
        END SUBROUTINE wrtyrcflows

      END INTERFACE

! ... Local variables
      integer :: stim(8) = 0
      integer :: flag754, nancnt ! function and returned valu for error checking
      logical :: drtybin
      logical :: wrbin = .false.
!      character bfr*80
      character (len=12) runtim, rtime
      real    :: tmp, getevi

      ! start the real time timer.
      call date_and_time(TIME=runtim, VALUES=stim)

! ... do a "brute force" initialization of all common block variables, cak - 06/04/02
      call default()

! ... Obtain startup information from user, do initializations based on
! ... answers to Modaid questions
      call detiv(wrbin)
      call message('   Model is running... '//runtim(:2)//':'//runtim(3:4)//':'//runtim(5:))

!      call readblk() ! ... Read first event block;  extra readblk call removed 2018/07/26

      ! Adjust crop.100 and fix.100 parameters for daily production -mdh 1/95
      ! removed adjustpar call; This didn't adjust parameters any more;
      ! moved PET message to message to wthini

! ... Write out starting values
      if (wrbin) call wrtbin(time) ! don't write if not open
      call listvar(21,time,' ',' ')

      cyear = strtyr
      doy   = 1
      month = 1
      do while (cyear <= tend)

        ! Perform annual tasks
        if (month .eq. 1) then
          ! if we have finished the block, read the next block
          if (cyear .eq. blktnd) call readblk()

          call eachyr  ! Perform annual tasks  (must be called AFTER readblk KLK)
        endif

  ! ... The main driver for the model; call decomposition, growth, etc.
        call simsom()
        drtybin = .true.

  ! ... Update monthly production output, cak - 10/01/03
  !     Removed from main since the code didn't work!!!!  *****   KLK 3/3/2015
  !      - Last month's growth was always zero since accumulators are cleared
  !        on harvest/last
  !      - generated out of bounds
  !         . before frst command growth frstmth = 0
  !         . every january month = 1
  !         agcmth(month)  = 0; agcmth(month)  = max(0.0, agcacc  - sum(agcmth))
  !         eliminates out of bounds as long as there is only 1 growing season
  !         but the final month is still incorrect.


  ! ... Write yearly output
  ! ... Add output for the N2 flux for the year and convert fluxes to
  ! ... g/m^2, cak - 01/16/03
        if ((time .ge. strplt) .and. (month .eq. 12)) then
          call wrtyearsum(time, N2O_year/10000, NO_year/10000, &
                          N2_year/10000, CH4yrox, nit_amt_year/10000, annppt)
          call wrtyrcflows(time, asom11tosom21, asom12tosom22, &
                           asom12tosom3, asom21tosom11, asom21tosom22, &
                           asom22tosom12, asom22tosom3, asom3tosom12, &
                           ametc1tosom11, ametc2tosom12, astruc1tosom11, &
                           astruc1tosom21, astruc2tosom12, &
                           astruc2tosom22, awood1tosom11, awood1tosom21, &
                           awood2tosom11, awood2tosom21, awood3tosom12, &
                           awood3tosom22)
        endif

  ! ... Write monthly trace gas output, cak - 05/14/42
        if (time .ge. strplt) then
          call wrttgmonth(time, N2O_month/10000, NO_month/10000, &
                          N2_month/10000, CH4mnox, nit_amt_month/10000, pptmonth)
        endif

        ! eliminate all the finite precision correction for time.
        ! integer year and month won't drift
        time = cyear + month*dt

  ! ... Write out values
        if ((tplt - time) .lt. dt * 0.5) then
          if (wrbin) call wrtbin(time) ! don't write if not open
          call listvar(21,time,' ',' ')
          drtybin = .false.
          tplt = time + dtpl
        endif

        ! Update clock
        month = month+1
        if(month .eq. 13) then
          month = 1
          cyear = cyear +1
          doy = 1
        endif
      end do

      !log IEEE 754 NaN/Inf plotx values which generally indicate a bad run
      nancnt = flag754()

!...Write out final values if the run doesn't end on a print time
      if (drtybin) then
         if (wrbin) call wrtbin(time) ! don't write if not open
         call listvar(21,time,' ',' ')
      end if
      call wrtsite(0,0,' ') ! output a site file if one was scheduled.


! ... Close data files


      close(unit=9) ! Close the weather file

      if (labtyp .gt. 0) close(unit=10) ! Close the c14data file if necessary
      if (Ninput .gt. 0) close(unit=20) ! Close the nscale.dat file if necessary
      if (OMADinput .gt. 0) close(unit=30) ! Close the omadscale.dat file if necessary

! ... Close the phscale.dat file if necessary
      if (phsys .gt. 0) then
        close(unit=40)
      endif
! ... Close the precscale.dat file if necessary
      if (wthinput .eq. 4 .or. wthinput .eq. 5) then
        close(unit=50)
      endif
! ... Close the tmaxscale.dat file if necessary
      if (wthinput .eq. 2 .or. wthinput .eq. 3 .or. &
          wthinput .eq. 5) then
        close(unit=55)
      endif
! ... Close the tminscale.dat file if necessary
      if (wthinput .eq. 1 .or. wthinput .eq. 3 .or. &
          wthinput .eq. 5) then
        close(unit=60)
      endif

      tmp = getevi('close') ! Close the EVI file
      close(unit=15) ! Close the schedule file
      close(unit=70) ! Close N2/N2O flux file, nflux.out (unit=70)
      close(unit=71) ! Close daily.out
      close(unit=72) ! Close summary.out
      close(unit=73) ! Close methane.out
! ... Close *.out and *.csv files opened by the initsw subroutine
      call closefiles()

      ! removed endfile(unit=1)
      !  This is an anachronism at best especially when the file is being closed
      if (wrbin) close(unit=1) ! Close binary file

      runtim = rtime(stim)
      if(nancnt > 0) then
        call message('Execution finished '//runtim)
      else
        call message('Execution success. '//runtim)
      endif
      STOP

      end


      character*12 function rtime(ctim)
        integer, intent(inout) :: ctim(8)
        integer, dimension(8)  :: rtim
        integer i, timzn, stim2
        !integer, parameter, dimension(12) :: dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        ! convert an elapsed time into an output string.
        ! DayCent runs ought to be short.
        !   The string is sized to handle only 2 digit hour intervals.
        !   The first day of the month rollover will also only handle 1 day (24 - 48 hours)
        !   The latter can be extended by uncommenting the DOM and leap year calculation

        ! updated for daylight savings, date and changes KLK 7 Nov, 16

        ! integer, parameter, dimension(12) :: dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        !    'VALUE(1)':       year
        !    'VALUE(2)':       month of year
        !    'VALUE(3)':       day of the month
        !    'VALUE(4)':       time offset with respect to UTC in minutes
        !    'VALUE(5)':       hour of the day
        !    'VALUE(6)':       minutes of the hour
        !    'VALUE(7)':       seconds of the minute
        !    'VALUE(8)':       milliseconds of the second


        stim2 = ctim(2)
        rtim = ctim                     ! save the start in the difference
        call date_and_time(VALUES=ctim) ! get a new clock time
        rtim = ctim - rtim              ! this is now the elapsed clock time

        ! Get rid of the zulu year change that occur on new years unless you are GMT
        timzn = mod(rtim(4), 1440) ! time zone is a fraction of a day

        ! Simplifications to correct run time.
        if(rtim(1) .eq. 1  .or.  rtim(2) .eq. 1) then;     ! hours
          rtim(1:3) = 0; rtim(5) = rtim(5)+24; ! first of month date change; add 24 hours

          !rtim(1:3) = [0, 0, rtim(3)+dom(stim2)]; ! hours
          !if(stim2 .eq. 2  .and. &
          !   (mod(stim2,400).eq.0 .or. mod(stim2,4).eq.0 .and. mod(stim2,100).ne.0)) rtim(3) = rtim(3)+1
        endif;

        ! assume the run is going to take less that 1 day
        if(rtim(8) .lt. 0) then;          ! milliseconds
          rtim(8)= rtim(8)+1000; rtim(7)= rtim(7)-1;
        endif;

        if(rtim(7) .lt. 0) then;          ! seconds
          rtim(7)= rtim(7)+60; rtim(6)= rtim(6)-1;
        endif;

        !minutes
        rtim(6) = rtim(6) - timzn;   ! convert to  UTC to eliminate time change artifacts
        if(rtim(6) .lt. 0) then;
          rtim(6)= rtim(6)+60; rtim(5)= rtim(5)-1;
        endif;

!       construct the output format HH:MM:SS.sss
        write(rtime,'(i0.2,":",i2.2,":",i2.2,".",i3.3)') rtim(3)*24+rtim(5), &
              rtim(6),rtim(7),rtim(8)

        do while (index(' 0:',rtime(1:1)) .gt. 0);  rtime= rtime(2:);  end do
!        do while (rtime(1:3) .eq. '00:');  rtime= rtime(4:);  end do
        do i = 1,len_trim(rtime)
          if(rtime(i:i) .eq. ' ') rtime(i:i) = '0'
        enddo

        if(index(rtime,':') .eq. 0) rtime= trim(rtime)//' s'

        return
      end


!      integer function myhandler(sig, code, context)
!        integer sig, code(5), context
!        call abort()
!      end

