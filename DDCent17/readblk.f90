
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

! ... READBLK.F

      subroutine readblk()
! ... Reads the next block of events from the schedule file

      implicit none

! Modifications
!   19May2015: KLK
!    Added Century 4 block read error messages
!   Apr14: KLK
!    Added orchard modifications
!              added a tree plant event
!              allowed FRST or TFST to signal current system
!              Signal an error on FRST/TFST without a defined crop/tree
!              corrected a potential problem from bad block numbers, blknum
!              moved the parsing of the initial system to detiv
!   9/2002  K Killian
!   Dec12: KLK
!    included changes for Century orchard models and water table for Methane model
!   May11: K. Killian
!    Moved some Century upgrades
!      - restructured readblk to move the event processor to readevt
!      - some F90
!      - made sure all error messages are routed through message and can be trapped.
!        Moved the event reading code to a separate routine so that it
!        can be used to read global events

      include 'chrvar.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'zztim.inc'

! ... Local variables
      integer           :: blknum, ctemp, ttemp, ii
      integer           :: prevsys
      real              :: pltmo
      integer           :: lch, istat
      character (len=90):: buffr
      character (len=4) :: blklbl

      integer   ReadInt
      double precision   ReadReal
!      character (len=40) READSTRING


! ... Set the starting year of the block
! ... Using ANINT which rounds to the nearest whole number.  -rm
      strtyr = anint(time)

      ! Read the block number, ending year, and # of years set up
      lch = 0
      do while (lch .le. 0)
        blknum = ReadInt(15,buffr,lch,istat)
        blklbl = buffr; ! save the first 4 characters, the number, for error reporting
        ! Check for conversion errors starting the block
        if(istat .lt. 0) then ! End of file
          write(buffr,'(a,i6)') 'EOF reading block for:',cyear
          call abortrun(buffr)
        else if(istat.eq.0) then ! no numeric conversion
          ! Skip the 'Year Month Option' header
          if(index(buffr(lch:),'Option') .ne. 0 .or. &
                 index(buffr(lch:),'option') .ne. 0 .or. &
                 index(buffr(lch:),'OPTION') .ne. 0      ) then
            lch = 0
            cycle
          else ! Abort if we found about anything else
            call abortrun('reading Block #:'//trim(buffr))
          end if
        end if
      end do
      blknum = max(1, blknum)

      !read(15,*) blktnd
      lch = 0
      blktnd = ReadInt(15,buffr,lch,istat)
      if(istat.le.0) call abortrun('reading block '//blklbl//' Last year')

      !read(15,*) rptyrs
      lch = 0
      rptyrs = ReadInt(15,buffr,lch,istat)
      if(istat.le.0) call abortrun('reading block '//blklbl//' Repeat length')
      ! set a maximum BLOCK repeat length
      ! if(rptyrs .gt. 41297762) call abortrun('Repeat length exceeds maximum:'//buffr(:40))

      !read(15,*) strplt
      lch = 0
      strplt = ReadReal(15,buffr,lch,istat)
      if(strplt .lt. cyear) strplt = cyear ! output should start on or after the current year
      if(istat.le.0) call abortrun('reading block '//blklbl//' Output start year')

      ! read(15,*) pltmo
      lch = 0
      pltmo  = ReadReal(15,buffr,lch,istat)
      if(istat.le.0) call abortrun('reading block '//blklbl//' Output month')

      ! read(15,*) dtpl    ! removed fixed format read '(f6.3)' klk 19Oct2011
      lch = 0
      dtpl = ReadReal(15,buffr,lch,istat)
      if(istat.le.0) call abortrun('bad block '//blklbl//' Output interval'//trim(buffr))

! ... Add 1 more to blktnd for January of the ending year
      blktnd = blktnd + 1.

! ... Check dtpl;
      if (dtpl .lt. .01) then
        dtpl = dt  ! Reset to monthly if lower than lowest limit
      else

        ! Rounding; units are years, but it should represent an integral number
        ! of months as accurately as possible.
        dtpl = anint(12.*dtpl)/12.
      endif

! ... Set tplt, the next time to spit out output
!      pltmo = pltmo - 1
      tplt = strplt + real(pltmo)/12.
! ... Round tplt to 2 digits of precision, cak - 03/22/02
      tplt =  anint(tplt * 100.0) / 100.0
!      if (tplt-int(tplt) .eq. 0.0) then
!        tplt = tplt + dtpl
!      endif

! ... Determine the weather type and call wthini to initialize the block's weather
!      lch = 0; wthr = READSTRING(15,buffr,lch,ststrt,stend, .FALSE., 0)
!      if(lch .eq. -1  .and.  buffr .eq. '?--?') call abortrun('reading weather:'//trim(buffr)
      read(15,*) wthr  ! removed fixed format '(a1)' klk 19Oct2011
      if (wthr .eq. 'F' .or. wthr .eq. 'Y') then
        ! read(15,'(a)') wthnam ! lengthen input from '(a20)' klk 19Oct2011
        lch = 0; call READCLIN(15,wthnam,lch,.TRUE.) ! read with standerd dayCent parsing
        ! remove quotes which could have been used force paths
        wthnam = ADJUSTL(wthnam) ! left justify
        if    (wthnam(1:1) .eq. '"') then
          lch = index(wthnam(2:),'"')
        elseif(wthnam(1:1) .eq. "'") then
          lch = index(wthnam(2:),"'")
        else
          lch = 0
        endif
        if(lch.ne.0) wthnam = ADJUSTL(wthnam(2:lch))
        call wthini(cyear,doy) ! open the weather file for RAMS/Century linkage, ... mdh 12/96
      else if (wthr .eq. 'C') then ! continue weather file
      else if (wthr .eq. 'M'  .or.  wthr .eq. 'S') then
        call abortrun('mean and stochastic weather not supported in this version of DayCent')
      else
        call abortrun('unknown weather option "'//trim(wthr)//'"')
      endif

!...Read the events in the block
      call readevt(blknum)

! ... Reset evtptr, the array index into cmdary, timary, typary
      evtptr = 1

! ... Set up cursys, the current system(s) if this is the first block read
!     cursys: CRPSYS = crop/grass; FORSYS = forest; SAVSYS = savanna
!     Reset system for each block to handle changing systems.

! ... Store the value of the previous system
      prevsys = cursys

      ! temporary variables; system is active if anything is growing at the block start
      ctemp = 0;   if(crpgrw.eq.1) ctemp =1
      ttemp = 0;   if(forgrw.eq.1) ttemp =1
      do ii = 1, ttlind
        select case (cmdary(ii))
          case ('CROP')
            if (curcrp .eq. ' ') call cropin(typary(ii)) ! ensure a crop is loaded
            ctemp = CRPSYS
          case ('FRST','PLTM','LAST','SENM')
            ctemp = CRPSYS
          case ('TREE')
            if (curtre .eq. ' ') call treein(typary(ii)) ! ensure a tree is loaded
            ttemp = FORSYS
          case ('TFST','TPLT','TREM','THRV')
            ttemp = FORSYS
        end SELECT
      end do
      cursys = ctemp + ttemp

! ... If no crop or forest events were given in the block, use the previous system
      if (cursys .le. 0  .or.  cursys > SAVSYS) then
        cursys = prevsys
        if(cursys .le. 0  .or.  cursys > SAVSYS) call abortrun('undefined ecosystem, need a crop or tree defined')
        call message('undetermined system in block '//blklbl//': using previous system')
      endif

! ... Check if decsys needs to be updated.
      if (decsys .eq. CRPSYS .and. ttemp .eq. FORSYS) then
        decsys = FORSYS
      endif

      return
      end subroutine readblk



      subroutine readevt(blknum)

!      use flowbal  !aceq

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'jday.inc' ! prob
      include 'timvar.inc'

      integer   blknum

!...Functions
      integer          ReadInt, isevt, reldate
      double precision ReadReal
      character        READSTRING*40

!...Local variables
      integer   i, j, ststrt, stend, lch, istat, mxstp, parserr, nevt
      integer   yrnm, stpnm, lstyr, lstsp, argtyp, lchlyr, repeat, prevevnt
      integer   day
      integer   harvd, cultd, pltd, fertd, tremd, tlasd, lstmnth
      real      realarg
      character (len=5)   :: evntcmd
      character (len=28)  :: evntyp ! event option
      character (len=90)  :: buffr
      character (len=132) :: mssg
!      character tagint*3
      logical   echprt
      logical   deftstp, mnthin
      logical            :: clrafrt

      integer, parameter :: schunit =15
      integer, parameter :: TSTYR = NDAY

!...Reads the next block of events from the schedule file

!=======================================================================
!
! The following block reads the events from the schedule file.
!   Uses subroutines that allow free form number and string input.
!   It also simplifies the input format by relaxing the year input requirement
!   (it will assume the same year as the previous command)
!
!  History
!     Jul12: K. Killian
!              - added a persistant auto fertilizer, AFRT, event
!              The aufert value specified by AFRT persists until reset by
!              another AFRT event or cleared by a FERT event. There is also
!              an option to clear the on the last day of the block.
!              The AFRT option is the real aufert value. If the value ends
!              with a C, another AFRT will be scheduled at the end of the block.
!              - Corrected a bug that mis parsed plank lines in year - step combinations
!     May11: K. Killian
!              finished making schedule debugging a user input
!                first line 'echo' prints out events read
!                first line 'DEBUG' prints out events read and scheduled
!              Corrected a data check that prevented negative global dates.
!              comment changes and reordered code for clarity
!              changed end of block flag back from -999 to the X event
!             -for readability added some local variables for event inputs
!             -Modifications to handle daily input
!              Default is to read the block in months. The switch to DOY will
!              be done if the step input is > 12 and <= NDAYS or month/day.
!              the mode can be forced by inputting DAY or MONTH at block start.
!              month to DOY conversion occurs AFTER input is complete.
!              new variables
!              day       date value
!              deftstp   DEFault Time STeP switching
!              mnthin    month input
!              maxstp    maximum time step, number of month or weeks
!     Dec10: K. Killian
!            - extended relative dates to model repeating events at an interval
!              relative to the start of the schedule. Example a forest clear
!              at long intervals. This saves having to write VERY long blocks
!              that can overflow the event buffer.
!            - the word 'debug' at the start of the block triggers the debug
!              event parsing list for that block
!     Jun09: K. Killian
!              added relative dates for global events, block 0
!              suppress date order check for block 0 since they are sorted later
!     Apr09: K. Killian
!              A zero month now generates an error condition
!              for readability added some local variables for event inputs
!              removed the floating argument array
!     Mar08: K. Killian
!              New error message for possible missing numeric arguments
!     09/15/07: Kendrick Killian
!              added check for an over length event
!     07/10/03: Kendrick Killian
!              Changed variables for readability and isevt redefinition
!     04/17/03: Kendrick Killian
!              changed the error reporting variable from blknum to parserr
!              since this routine kills the run itself
!     07/15/03: Kendrick Killian
!              changed the error reporting and allowed multiple global
!              events in a single record to share years.
!     07/10/03: Kendrick Killian
!              added the tillage event.
!     05/2000 K. Killian
!              Moved the event reading code to a separate routine so that it
!              can be used to read global events
!     01/14/99: Kendrick Killian
!              added a drainage event.
!     04/03/95: Kendrick Killian
!              added error checking on Century events and dates
!    Written : 1/25/95 : Kendrick Killian
!
!  New External References
!    character function READSTRING (UNIT,BUFFR,lch,ststrt,stend,quotd)
!          double precision function      ReadReal   (UNIT, BUFFR, LCH, ISTAT)
!    integer function   ReadInt    (UNIT, BUFFR, LCH, ISTAT)
!
!  Additional comments
!    1) input commands are converted to upper case
!    2) input fields can be separated by commas, white space or endline
!    3) UNIX like comments can be inserted in the input stream
!       - # is the comment character
!       - comments extend to the end of the line
!    4) Century events can be entered without the year.  The year is
!        pulled from the previous event.
!    5) Century commands with NO date arguments generate a fatal error
!
!========================== BEGIN MODIFIED CODE ========================

!...Read and save the events for this block
      argtyp = 0
      lchlyr = 0
      lstsp  = 0
      lstyr  = 0
      parserr = 0
      repeat = 1
      ttlind = 0
      clrafrt = .FALSE.
      echprt  = .FALSE.
      dbgpsch = .FALSE.

!     time step conversion controls
      deftstp = .TRUE.
      mnthin  = .TRUE.
      mxstp  = MONTHS

!     set buffer pointer (lch) beyond the buffer to force an input line read
      lch = len(buffr)+1
!      number of events we can store
      nevt = lenary
      if(glevpt .gt. 0 .and. glevpt .lt. lenary) nevt = glevpt

      istat = 0
      do while(istat .GE. 0 .and. ttlind .le. nevt )
        evntyp= ""
        day   = 0
        yrnm = ReadInt(schunit,buffr,lch,istat)

!       respond to an empty line
        if(istat .eq. -1  .and.  lch .le. 0) then

!         End of file
          if(lch .eq. -1) then
            write(mssg,*) 'EOF reading block in readevt',blknum
            call abortrun(mssg)
          end if

!         EXIT for global events
          if(blknum .eq. 0) exit

!         SKIP a blank line
          istat=1
          cycle
        endif

!       NO CONVERSION: check for control commands or relative global dates

        ! look for global, block 0, dates relative to run start/end
        if (istat .eq. 0  .and.  lch .gt. 0  .and.  blknum .eq. 0) then
          yrnm = reldate(buffr, lch, istat)
          if(istat .gt. 1) then ! istat > 1 then we have a repeated command
            repeat = istat
            istat = 1
          endif
        endif

        if (istat .eq. 0) then
          i = index('echo ECHO debugDEBUG',buffr(lch:lch+4))
          if(i .ne. 0) then
            write(mssg,*) " debugging block ",blknum
            call message(mssg)
            echprt  = .TRUE.
            if(i .gt. 10) dbgpsch = .TRUE.
            lch = lch + index(buffr(lch:),' ') +1
            cycle
          endif

          if(deftstp) then
            if(index('monthMONTHMonth',buffr(lch:lch+4)) .ne. 0) then
              mnthin  = .TRUE.
              deftstp = .FALSE.
              mxstp  = MONTHS
              lch = 0
              cycle
            elseif(index('dayDAYDay',buffr(lch:lch+2)) .ne. 0) then
              mnthin  = .FALSE.
              deftstp = .FALSE.
              mxstp  = TSTYR
              lch = 0
              cycle
            endif
          endif

          if(argtyp .eq. 3) then
            mssg = 'missing event time or '//trim(buffr)//' argument'
          else
            call message('missing date: '//buffr(lch:lch+10))
          end if
          call FILRERR(schunit,0,mssg,buffr)
        end if

!        if (yrnm .eq. -999) exit  ! check for end of block

!...Fill in the date parameters
        ! reset lch so ReadInt doesn't return a blank line as a 0  KLK 11Jul2012
        if(lch .gt. len(buffr)) lch =0
        stpnm =  ReadInt(schunit,buffr,lch,istat)
        if(buffr(lch:lch) .eq. '/') then
          mnthin=.FALSE.
          lch = lch+1
          i = 0
          day = ReadInt(schunit,buffr,lch,i) ! read date portion

          if(i .eq. 0) then
            day = 0 ! do a month only conversion
          else if(i .lt. 0) then
            ! error out if date read fails
            call FILRERR(schunit,0,'bad month/date DOY',buffr)
          endif
        endif

!...check for the explicit timestep argument
        if (istat.eq.0) then
          stpnm = yrnm
!...global block input
          if (blknum .eq. 0) then
!           a global event record may have multiple events from the same year
!           check a previous event for a missing year
!           NOTE: This previous year code can be defeated if the command on a
!                 second line is farther to the right. This may or may not be
!                 what is intended.
            if (ttlind .gt. 0  .and.  lch .gt. lchlyr) then
              yrnm = lstyr
            else
              call abortrun('No year specified for global event:'//trim(buffr))
            endif
          else if (ttlind .gt. 0) then
            yrnm = lstyr
!...bad input
          else
            call message('Warning: assume first event occurs in year 1:'//trim(buffr))
              yrnm = 1
          endif
        else if (istat .lt. 0) then
          call FILRERR(schunit,istat,'I/O error reading timestep',buffr)
        endif
        lchlyr = lch

!...Read the command
        evntcmd = READSTRING(schunit,buffr,lch,ststrt,stend,.true.,1)
        if(lch.le.0) call FILRERR(schunit,lch,' Missing event',buffr)

        ! convert month/date to DOY
        if(day .gt. 0) then
          stpnm = ifrstdy(stpnm) + day -1
        elseif(day .lt. 0) then
          stpnm = ilstdy(stpnm) + day +1
        endif

        ! Look for dated exit command, 'X'
        if(evntcmd.eq.'X' .or. evntcmd.eq.'x') exit

        ! use parserr to mark an error  (parserr < 0  bad input)
        if(stpnm .gt. MONTHS .and. stpnm .le. TSTYR .and. deftstp) then
          mnthin = .FALSE.
          mxstp  = TSTYR
        endif

!...check date input consistancy
!     use parserr to mark an error  (parserr < 0  bad input)
        if (blknum .ne. 0) then
          ! for a block check  0 < sequence year <= repeat length
          if (yrnm .le.0) then
            write(mssg,'(a,i5)') 'ERROR: block year less than 0 ',yrnm ! timary(1,ttlind)
            call message(mssg)
            call message(buffr)
            yrnm = 1
            parserr = -1
          elseif (yrnm .gt. rptyrs) then
            write(mssg,'(a,i5)') 'ERROR: year exceeds rotation length ',yrnm
            call message(mssg)
            call message(buffr)
            yrnm = rptyrs
            parserr = -1
          endif
        endif
        if(stpnm.le.0 .or. stpnm.gt.mxstp) then
          write(mssg,'(a,i5)') 'ERROR: illegal time step ',stpnm
          call message(mssg)
          call message(buffr)
          stpnm = 0
            parserr = -1
        endif

!       check event order if we have more than 1 event
!         skip the check on global events since we will sort them later
        if (ttlind.gt.0  .and.  blknum .ne. 0) then
          if (lstyr.gt.yrnm  .or.                                        &
     &        (yrnm  .eq. lstyr  .and.  lstsp.gt.stpnm)) then
            call message('ERROR: dates out of sequence')
            call message(buffr)
            write(mssg,'(3(i5,a),i3)') lstyr,'/',lstsp," comes after ",yrnm,'/',stpnm
            call message(mssg)
            parserr = -1
          endif
        endif

!       argtyp the event argument type  istat = event number
        argtyp = isevt(evntcmd)

!       generate an error if this is not an event
        if(argtyp .eq. 0) then
          call FILRERR(schunit,0,'unknown CENTURY event: '//evntcmd,buffr)

!       Events with NO arguments PLTM FRST LAST SENM TFST TLST TFLR ',
        else if(argtyp .eq. 1) then
          evntyp= " "

        else ! Fill in additional information   (as required)
             ! All options are treated as strings KLK

          ! NOTE: set the case option to match the nextopt read case in chkdata
          evntyp= READSTRING(schunit,buffr,lch,ststrt,stend,.true.,0) ! changed to case sensitive 1 -> 0
          if(evntyp .eq. ',') lch = 0
          if(lch .le. 0) call FILRERR(schunit,lch,' Missing '//evntcmd//' option', buffr)

!         *** SPECIAL CASE events ***
          if(argtyp .eq. 3) then !check number argument convert (EROD, DRAN)
            i=0
            realarg = ReadReal(0,evntyp,i,istat)
!            write(*,*) "real type ",i,lch,event,realarg," '",evntyp,"'"
            if(istat .le. 0 .or. (i.lt.len(evntyp) .and. evntcmd .ne. 'DRAN')) &
            call FILRERR(schunit, istat,'bad '//evntcmd//' real data "'//trim(evntyp)//'" ',buffr)

!         *** SPECIAL CASE events ***
          ! Look for a block autofert clear flag
          else if(evntcmd.eq.'FERT') then
            clrafrt = .FALSE.! FERT events clear AFRT

          else if(evntcmd .EQ. 'ACEQ') then
              if(blknum .eq. 0  .and.  SCAN(evntyp, 'ARar') .gt. 0) then
                write(mssg,*) yrnm-strtyr+1
                mssg = adjustl(mssg)
                i = index(evntyp, '(') + 1
                j = MAX(len_trim(evntyp), index(evntyp, ')') - 1)
                ! remove enclosing parenthesis;
                ! they are not required but are accepted both to protect ','s and
                ! because the argument is similar to other arguments in parenthesis
                call EqAclrat("SET"//evntyp(i:j)//trim(mssg)//'I') ! parse with global auto schedule

              else
                call EqAclrat("SET"//trim(evntyp))                  ! parse without auto schedule
              endif

          else if(evntcmd.eq.'AFRT') then
            clrafrt = .FALSE. ! This clears the block
            i=0
            realarg = ReadReal(0,evntyp,i,istat)
            if (istat .gt. 0   .and.   i .le. len(evntyp)) then
              if(blknum.gt.0 .and. evntyp(i:i).eq.'C') clrafrt = .TRUE.
              evntyp(i:) = ' '
            endif

!          else if(evntcmd.eq.'TILL') then ! Combine TILL events
!            ! all TILL events need to be coded so they can be combined by tillin
!            !      to a single event leaving only one entry per time step
!
!!           convert block year to base 32 (1-1023 years maximum)
!            mssg(1:3) = tagint(yrnm)
!            mssg(1:2) = mssg(2:3)
!!           convert block "month". Using two digits coding day numbers
!            mssg(3:5) = tagint(stpnm)
!            mssg(3:5) = mssg(4:5)//" "
!!            call tillin(mssg(1:4),typary(ttlind), istat)
!!            typary(ttlind) = mssg
!            call tillin(mssg(1:4),evntyp, istat)
!            evntyp = mssg
!            ttlind = ttlind +istat
          endif

        endif

!     record the last entry read. These need to stay in input units
        lstyr  = yrnm
        lstsp  = stpnm

!     record the final location of the event
65    continue
      prevevnt = ttlind
      do while (prevevnt >=1)
        ! check previous events on the same date
!write(*,*) "check ", yrnm,stpnm, trim(evntcmd),"  ",trim(evntyp)," against", &
!prevevnt,timary(1:,prevevnt), trim(cmdary(prevevnt))," ", trim(typary(prevevnt))
        if(timary(1,prevevnt) == yrnm    .and. timary(2,prevevnt) == stpnm) then
           ! Is the event itself a duplicate todey?
           if(cmdary(prevevnt)   == evntcmd .and. typary(prevevnt)   == evntyp) exit; ! found a duplicate; don't need to search further

           prevevnt = prevevnt -1; ! decrement to a new event

        else
          ! prevevnt = prevevnt -1; ! decrement to a new event
          prevevnt = 0; ! ordered list (at least it's supposed to be) so terminate the search
        endif
      end do

      if (prevevnt > 0) then ! was there a duplicate?
        !print a warning message and
        write(buffr,'(i5," event ",i6,i4, 2(2x,a))') blknum, yrnm, stpnm, trim(evntcmd), evntyp
        call message("skipping duplicate event block"//trim(buffr))

      else
        ttlind = ttlind + 1 ! new event location

        ! Can this be stored without overflowing the array
        if(ttlind .gt. nevt) then
          write(buffr,'(i5," number ",i5,":  event \",i6,i4, 2(2x,a))') blknum, &
                       ttlind,yrnm, stpnm, trim(evntcmd), evntyp
          buffr = adjustl(buffr)
          call abortrun("event block "//trim(buffr)//"\ exceeded total event count")
        endif

       ! store the event description in the command arrays
       !  This ould be cool as a structure but I'm not inclined to fight for approval
        timary(1,ttlind) = yrnm
        timary(2,ttlind) = stpnm
        cmdary(ttlind) = evntcmd
        typary(ttlind) = evntyp
      endif


! ***************************************************
        if(echprt) then
          mssg = 'Read #'
          write(mssg(14:),'(i5,"  ->",i5,i4)') ttlind, yrnm,stpnm
          if(.not. mnthin) then
            call DOYdate(stpnm, .FALSE., i, day)
            write(mssg(len_trim(mssg)+2:),'("(",i2,"/",i2,")")') i,day
          endif
          mssg(len_trim(mssg)+3:) = evntcmd//' "'//evntyp
          mssg(len_trim(mssg)+1:) = '"'
          call message(mssg)
        end if
! ***************************************************
        if(repeat > 1) then
          yrnm = reldate(buffr, lch, repeat)
          goto 65
        endif
      enddo

!     exit the routine


!     Log the schedule mode
! ---- handle month to DOY corrections ------------
       if(mnthin  .and.  ttlind .gt. 0) then
         write(mssg,*) ttlind,"Monthly Events in Block ",blknum
         call message(mssg)
         i = 1
         lstmnth = 0
         harvd = 0
         cultd = 0
         pltd  = 32
         fertd = 0
         tremd = 0
         tlasd = 0
         do while (i.le. ttlind) ! use this so we can modify the limits
           ! if the month changed, clear the month and cult date flags
           if(lstmnth .ne. timary(2,i)) then
             harvd = 0
             cultd = 0
             pltd  = 32
             fertd = 0
             lstmnth = timary(2,i)
           endif
           ! assign a day to this this event
           timary(2,i) = evntday(timary(2,i), cmdary(i), typary(i))

           ! save the day of this months cult and harvest events
           if    (cmdary(i) .eq. 'HARV') then
             harvd = timary(2,i)
           elseif(cmdary(i) .eq. 'CULT') then
             cultd = timary(2,i)
           elseif(cmdary(i) .eq. 'FERT') then
             fertd = timary(2,i)
           elseif(cmdary(i) .eq. 'PLTM' .or. cmdary(i) .eq. 'TPLT') then
             pltd  = min(timary(2,i),pltd)
           elseif(cmdary(i) .eq. 'TREM') then
             tremd  = timary(2,i)
           elseif(cmdary(i) .eq. 'TLST') then
             tlasd  = timary(2,i)
           end if

           ! if we CULT and HARV this month; move the CULT after the harvest
           if(cultd .ne. 0  .and.  harvd .ne. 0) then
             ! search backwards until we find the cult index
             lch = i
             do while (cmdary(lch) .ne. 'CULT')
               lch = lch -1
             end do
             timary(2,lch) = harvd +1 ! reschedule the cult
             harvd = 0  ! clear the rescheduled harvest date
           end if

           ! if we FERT and plant this month; move the FERT to the plant
           if(pltd .ne. 32  .and.  fertd .ne. 0) then
             ! search backwards until we find the cult index
             lch = i
             do while (cmdary(i) .ne. 'FERT')
               lch = lch -1
             end do
             timary(2,lch) = pltd +1 ! reschedule the fert
             fertd = 0  ! clear the rescheduled fert date
           end if

           ! if we TLSt and TREM this month; move the TREM after the end of growth
           if(tremd .ne. 0  .and.  tlasd .ne. 0) then
             ! search backwards until we find the removal index
             lch = i
             do while (cmdary(lch) .ne. 'TREM')
               lch = lch -1
             end do
             timary(2,lch) = tlasd +1 ! reschedule the removal until after last
             tremd = 0  ! clear the rescheduled removal day
           end if

           ! code to replicate events.
           ! NOT needed with current DailyDayCent distributed IRRI and GRAZ events.
           ! if(timary(2,i) .gt. 1000) call repeatevnt()
           i = i+1
         end do ! if(i .lt. ttlind) goto 200 ! enddo
         i = evntday(-999, 'X', '') ! clear the event day assignments between blocks

         ! if(echprt) call prntevt('Translated: ', 1,ttlind)

         call sortevt(1, ttlind)

         if(echprt) then
          ! write(mssg,*) ttlind," events scheduled:"; call message(mssg);
            call message('')
            call prntevt('Scheduled: ', 1,ttlind)
         endif

       end if

       if(clrafrt) then
         ttlind = ttlind + 1
         if(ttlind .gt. nevt) then
           write(buffr,'(i5)') blknum
           buffr = adjustl(buffr)
           call abortrun("adding AFRT event, block"//trim(buffr)//", exceeded total event count")
         endif
         timary(1,ttlind) = rptyrs
         timary(2,ttlind) = 365
         cmdary(ttlind) = 'AFRT'
         typary(ttlind) = "0.0"
       end if

!...   return if the read was successful or non fatal message
      if(parserr .ge. 0) return
!
!...abort on illegal event input
      call abortrun('parsing events')
!============================ End Modified Code ========================
      contains
        integer function evntday(month, event, typaray)
          integer month
          character*(*) event, typaray
          ! modify to allow certain events to have multiple monthly entries.
          ! 'CULT', 'HARV', 'TREM', 'FERT', 'OMAD'
          ! can be shifted by a day to allow for a second event in the month
          !  ignore  'GRAZ' and 'irri' since Daycent spreads them out and there
          !  is no way to avoid the collision
          !
          !       convert monthly values to a day based on the monthly event order
          !  Events 'GRAZ', 'TREM', 'FIRE', 'HARV' and 'THRV' occur after plant
          !    growth is complete so set them to the last day of the month.
          !  Irrigation, IRRI, is a special case. Return the first week *1000
          !    plus the last week so read block can schedule a modified event
          !    in any week of the month.
          !  All unlisted events occur in the middle of the month.

          include 'const.inc'
          include 'jday.inc'
           ! idysimo(12)=(31,28,31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
           ! ilstdy(12)= (31,59,90,120,151,181,212,243,273,304,334,365)
           ! ifrstdy(12)=( 1,32,60, 91,121,152,182,213,244,274,305,335)

          integer m !, winmonth
          integer, save :: lmonth = -999
          integer, save :: lcult, lfert, lharv, ltrem

          ! clear the month memory at the start/end of a block
          if(month .eq. -999 .and. event .eq. 'X') then
            lmonth = month
            evntday = 0
            return
          endif

          m = modulo(month-1,12)+1 ! Make sure this a valid month

          if(lmonth .ne. month) then
            lcult = ifrstdy(m) +2
            lfert = ifrstdy(m) +2
            lharv = ilstdy(m)  -1
            ltrem = ilstdy(m)  -1
            lmonth = month
          endif

         select case (event)
           ! definitions occur as early as possible
           case ('TREE','CROP','OTRF','THRV')
             evntday = ifrstdy(m)

           case ('IRRI','GRAZ')
           ! Special Cases - redefine evntday
           ! DayCent spreads these out (I think)
!             evntday = (ilstdy(m)-ifrstdy(m))*1000 + evntday
             evntday = ifrstdy(m)+1

           ! cultivation occured before growth
           case ('CULT')
             evntday = lcult
             lcult   = lcult -1

           ! growth starts in the first week of the month
           case ('FRST','PLTM','TFST','TPLT')
             evntday = ifrstdy(m) +3

           case ('FERT','AFRT')
             evntday = lfert
             lfert   = lfert +1

           ! harvest and tree removal occur at the end of the month
           case ('HARV')
             evntday = lharv
             lcult   = lharv +1

           case ('TREM')
             evntday = ltrem
             lharv   = ltrem +1

           ! SENM occurs before harvest which probably isn't done. Delay until the last day
           case ('SENM')
             evntday = ilstdy(m)

           ! initial FLOD occurs just after harvest
           case ('FLOD')
             if(typaray(1:1) .eq. '0'  .or.  typaray(1:1) .eq. 'w') then
               evntday = ilstdy(m)  -10
             else
               evntday = ifrstdy(m) + 5
             endif

           ! Century 4 removal and harvest events occur after growth.
           ! Cindy  putthem on the next to last day
           case ('LAST','TLST','FIRE')
             evntday = ilstdy(m) -1

           case default
             ! by default assume the event will occur in the middle
             ! intentionally leave FLOD for mid month to avoid puddle at plant and harvest
             evntday = ifrstdy(m) +14
         end select
        return
        end function evntday

        subroutine repeatevnt()
          integer :: j
!           This is probably a cleaner implementation since the current
!           smeared applications is hidden and I don't think can be overridden.

          ! treat events (IRRI) that translate to multiple weekly events

          ! we need to insert lch entries to the middle of the stack
          lch = timary(2,i)/1000
          timary(2,i) = timary(2,i) - lch*1000

          ! make room for a block of events in the stack
          timary(:,i+1+lch:ttlind+lch) = timary(:,i+1:ttlind)
          cmdary(i+1+lch:ttlind+lch) = cmdary(i+1:ttlind)
          typary(i+1+lch:ttlind+lch) = typary(i+1:ttlind)
          do j=ttlind, i+1, -1
            timary(:,j+lch) = timary(:,j)
            cmdary(j+lch) = cmdary(j)
            typary(j+lch) = typary(j)
          enddo
          ttlind = ttlind+lch

          ! store the repeat count as a non-printing character in the command
          ! IRR and GRA are unique so we don't need the character
          typary(i)= typary(i)(:len(typary(i))-1)//char(lch) ! (char(lch +ichar("0")))

          ! NOW duplicate the irrigation blocks
          do j= 1, lch
            timary(1,i+j) = timary(1,i)
            timary(2,i+j) = timary(2,i) +j
            cmdary(i+j) = cmdary(i)
            typary(i+j) = typary(i)
          enddo

          i = i+ lch  !    don't reprocess the added events
        end subroutine repeatevnt
      end subroutine readevt


      subroutine readglobalevt()

      implicit none
      include 'chrvar.inc'

        call readevt(0)

        !======================== global events ========================
        if (ttlind .gt. 0) then
          ! sort events since events relative to start dates can change the order
          call sortevt(1, ttlind)
          glevpt = lenary+1-ttlind
          timary(:,glevpt:lenary) = timary(:,1:ttlind)
          cmdary(glevpt:lenary)   = cmdary(1:ttlind)
          typary(glevpt:lenary)   = typary(1:ttlind)

          timary(:,1:ttlind) = 0
          cmdary(1:ttlind)   = " "
          typary(1:ttlind)   = " ";   !        ttlind = 1
        else
          glevpt = ttlind ! no global events
        endif
        !======================== End global events ====================
        return
      end subroutine readglobalevt


      integer function isevt(buffr)
      implicit none
        character*(*) buffr

! This routine encapsulates the event search function, returning a unique
! numeric value corresponding to the argument type required by the event.
!    0  Not an event
!    1  no argument
!    2  string argument
!    3  real argument
!
! Readblk has been modified so not to duplicate this code.
! If the string argument is a previously parsed event token, (4 characters)
! then case conversion is NOT done.
! If it is a longer string then the character tools are called to parse and
! uppercase the tokens. (this may be extraneous because of other code)
!
! programming note: it was easier to assume that the default is no argument
!     and to put those events at the end of the event list.
!
! ADDING EVENTS
!  Events are added to the parser by adding the event nemonic to the correct
!  list, depending on its argument type. Then increase the event count e?cnt
!  to reflect the number of events in the revised list. These values are in
!  the declaration portion of this routine.
!  The code in this routine or in readblk will not require any changes. There
!  will need to be a " case ('EVENT')" block added to schedl handle the new
!  event
!
! History
!   May11: K. Killian
!     removed the unused numarg from the argument list
!     Used F90 declaration conventions to convert the event list into 3 separate
!       list. making it easier to add an event to appropriate list
!   Apr11: K. Killian
!     Added a the accelerated equilibrium command ACEQ for testing
!   Mar09: K. Killian
!     Added a fixed area savanna with fertilizer preference event, OTRF
!   08/07: K. Killian
!     Added fruit harvest commands THRV and TFLR
!   07/07: K. Killian
!     Recoded so the function returns the argument type and parameter has the
!     unused event number. Changes some variables for readability
!   written 6/2003  K Killian

        integer   lch,ststrt,stend
        character READSTRING*40, evt*4

        integer   charend,  numevt

        character (len=62), parameter :: wc =                            &
     &             '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'//              &
     &             'abcdefghijklmnopqrstuvwxyz'

        ! type 1 events, those requiring NO argument
        integer, parameter :: e1cnt = 7
        character (len=e1cnt*5), parameter :: c1el = 'PLTM FRST '//      &
     &    'LAST SENM TFST TLST TFLR '

        ! type 2 events, those requiring a .100 option argument
        integer, parameter :: e2cnt = 15
        character (len=e2cnt*5), parameter :: c2el = 'CROP CULT '//      &
     &    'FERT FIRE GRAZ HARV IRRI IRIG OMAD TREE TREM AFRT OTRF TPLT ACEQ ' ! TILL

        ! type 3 events, those requiring a numeric argument
        ! NOTE:: numeric arguments checked to be numbers regardless of storage.
        integer, parameter :: e3cnt = 4
        character (len=e3cnt*5), parameter :: c3el = 'DRAN EROD FLOD THRV '

        isevt = 1
        charend = len_trim(buffr)

        if     (charend .lt.4) then
!          return 0 if it is shorter than an event
           isevt = 0
           return ! return the unknown event
        else if(charend .eq.4) then
!          scan the event list if we have a likely candidate
           numevt = INDEX(c2el//c3el//c1el,buffr(:4))
!           write(*,*) "single ",charend,"  '",buffr(:4),"'  ",numevt
        else
           lch = 1
           numevt = 0

!          look through the long string for any likely candidates
10         evt = READSTRING(0, buffr, lch, ststrt, stend, .TRUE., 1)
!           write(*,*) "LONG ",lch,charend,"  '",buffr(ststrt:stend),"'"
           if(stend-ststrt+1 .eq. 4) numevt= INDEX(c2el//c3el//c1el,evt)
           if (numevt .eq. 0 .and. lch .gt. 0) goto 10
        end if

!       convert the character number to an event number
!          a 0 event number, no event, generates an event type 0
        numevt = (numevt+4)/5
        if(numevt .eq. 0) then
          isevt = 0
        else if(numevt .le. e2cnt) then
          isevt = 2
        else if(numevt .le. e2cnt+e3cnt) then
          isevt = 3
        end if

        return
      end function isevt


      integer function reldate(buffr, lch, istat)
        character buffr*(*)
        integer   lch, istat

!=======================================================================
!
! This routine handles relative year input. This format provides a format
!   to index dates based on the schedule start and end dates. This format
!   was created to handle relative dates in global events. Made a stand alone
!   routine to allow the format to be applied to other date inputs
!    ++n    n years from start of the run    ++0 is the first year
!    --n    n years from the end of the run, --1 is the last year
!    %n     every n years from the start of the run
!
!   Changes lch and istat. Caller should handle usual error conditions.
!
!  History
!     Dec10:  K. Killian
!       added a repeating (modular) offset from run start
!     June09: K. Killian
!       created this routine to handle relative dates in global events.

        include 'timvar.inc'

        integer   ReadInt
        integer   relyr, repcnt, lchin, istin
        save      relyr, repcnt

!       set a return value
        reldate = 0
        lchin   = lch
        istin   = istat

!       repeat section for a modular repeating section
        if(repcnt .gt. 0  .and.  istat .gt. 1) then
          repcnt = repcnt +1
          reldate = strtyr -1 + repcnt* relyr
          istat = istat -1
          if(istat .eq. 1) repcnt = 0
          return
        end if

!       evaluate whether the error was caused by a relative date flag
        if    (buffr(lch:lch+1) .eq. '++') then
          reldate = 1
          lch = lch+2
        elseif(buffr(lch:lch+1) .eq. '--') then
          reldate = 2
          lch = lch+2
        elseif(buffr(lch:lch) .eq. '%') then
          reldate = 3
          lch = lch+1
        else
          return ! returns istat letting caller handle error
        endif

!       read the relative year
        istat=0
        relyr = ReadInt(0,buffr,lch,istat)

!       did the number parse still fail?  caller handles error
        if(istat .le. 0) return

!       offset relative to start
        if    (reldate .eq. 1  .and.  relyr .ge. 0) then
          reldate = strtyr + relyr
          istat = 1

        elseif(relyr .le. 0) then
!         resignal the original condition if the offset is bad
          istat = istin
          lch   = lchin
          reldate = 0

!       offset relative to end
        elseif(reldate .eq. 2) then
          reldate = tend - relyr
          istat = 1

!       interval offset relative to start
        elseif(reldate .eq. 3) then
          istat  = (tend - strtyr) / relyr
          repcnt = 1
          reldate = strtyr -1 + repcnt* relyr
        else
        end if

        return
      end function reldate


!     convert the week DOY to month
      integer function doy2mnth(DOY, leapyr)
        integer DOY
        logical leapyr
        integer sdc ! spring day count ! invers: Mod[spday[test] + 59, 365] + 1
        if(leapyr) then
          sdc = modulo(DOY - 61, 366)
        else
          sdc = modulo(DOY - 60, 365)
        endif
        doy2mnth=modulo((sdc/153)*5 + int(modulo(sdc,153)/30.5)+2,12) +1

        return
      end function doy2mnth


      Subroutine DOYdate(DOY, leapyr, mn, dy)
        integer, intent (In)  :: DOY
        logical, intent (In)  :: leapyr
        integer, intent (Out) :: mn, dy

      ! Function converts DOY and leap year to a month and day
      !
      ! The method is, starting in March, break the year up into 5 month
      ! sequences of 31, 30, 31, 30, 31 days. Calculate the sequence, sequence
      ! month, and the remainder day. Combine sequence and month and reorder
      ! to conventional format.

      ! Coded to use integer arithmetic instead of slower floating point
      ! and modular arithmetic.
      !  dy = modulo(DOY - 61, 366) = (DOY -61 .lt. 0)? DOY -61 : DOY + 305
      !  mn = dy/30.5              = (2*dy)/61
      !  dy = int(dy - mn*30.5) +1 - (2*dy - mn*61)/2 +1
      !  mn = modulo(cc*5 + mn +2, 12) +1 = cc*5 +mn +(mn<=9)? 3: -9;

        integer :: cc

        ! convert the DOY to a spring day (leap year sensitive)
        if(leapyr) then  !  sdc = modulo(DOY - 61, 366) leap year
          dy = DOY - 61
        else
          dy = DOY - 60
        endif
       if(dy .lt. 0) dy = DOY + 305 ! =  DOY -60 +365 = DOY -61 +366

        cc = dy/153                      ! 153 day cycle count
        dy = dy - cc * 153               ! remainder: days into cycle

        mn = (2*dy)/61                   ! number of average (30.5 day) months
        dy = (2*dy - mn*61)/2 +1         ! remainder  days into month, date
        mn = cc*5 + mn +3; if (mn>12) mn=mn-12 ! total months converted to convention

        return
      end subroutine DOYdate


      character*5 function datestr(DOY, leapyr)
        integer :: DOY
        logical :: leapyr

      ! Function converts DOY and leap year to a month/day string
      ! This is the same routine as DOYdate with the output written as a string

        integer :: cc, mn, dy

        ! convert the DOY to a spring day (leap year sensitive)
        if(leapyr) then  !  sdc = modulo(DOY - 61, 366) leap year
          dy = DOY - 61
        else
          dy = DOY - 60
        endif
       if(dy .lt. 0) dy = DOY + 305 ! =  DOY -60 +365 = DOY -61 +366

        cc = dy/153                      ! 153 day cycle count
        dy = dy - cc * 153               ! remainder: days into cycle

        mn = (2*dy)/61                   ! number of average (30.5 day) months
        dy = (2*dy - mn*61)/2 +1         ! remainder  days into month, date
        mn = cc*5 + mn +3; if (mn>12) mn=mn-12 ! total months converted to convention

        call DOYdate(DOY, leapyr, mn, dy)
        write(datestr,'(i2,"/",I2)') mn,dy

        if(datestr(2:2) .eq. " ") datestr(:2) = ADJUSTR(datestr(:2))
        if(datestr(4:4) .eq. " ") datestr(4:) = ADJUSTL(datestr(4:))
        return
      end function datestr



      subroutine sortevt(fevt, levt)

!      sort the events.
!      Required since month to week conversion may not leave the list sorted
!      This is a simple insertion sort.
!      Not normally efficient, since it is an N2 sort, but this is the best
!      case for a nearly sorted array.
!
!      modified Feb 2015 KLK
!      changed indexing on inner loop to handle an especially pathological dual
!       crop monthly schedule and limit access to input data range.
!
      implicit none
      integer   fevt, levt

      include 'chrvar.inc'
      include 'const.inc'

      ! local variables
      character (len=5)   :: ecmd
      character (len=28)  :: etyp
      integer   i, j, etim(2)
      logical :: test

       !write(*,'(a,2i5)') "sortevt ",fevt, levt;
       if(levt-fevt .eq. 1) return ! trivial array; just return

       ! write(*,'(2(a,i5),i4,a,a4,3a,g12.5)') 'start', 1,'   ->', timary(:,1), &
       !         ' "',cmdary(1),'" "', trim(typary(1)),'"';

       do i = fevt+1, levt

         test = timary(1,i-1).gt.timary(1,i) .or.                            &
                (timary(1,i-1) .eq. timary(1,i)  .and.  timary(2,i-1) .gt. timary(2,i))

         ! write(*,'(3(a,i5),i4,a,a4,3a,L12,L3)') 'sort ', i-1," < ",i,'   ->', &
         !           timary(:,i),' "',cmdary(i),'" "', trim(typary(i)),'"', test;

         if(test) then
           ! store the early event for movement down the stack
           etim = timary(:,i)
           ecmd = cmdary(i)
           etyp = typary(i)

           j = i
           do while(j .ge. fevt+1)
             if(timary(1,j-1).gt.etim(1) .or. &
                (timary(1,j-1).eq.etim(1) .and. timary(2,j-1).gt.etim(2))) then
               ! Move this event up the stack to make room for misplaced event
               timary(:,j) = timary(:,j-1)
               cmdary(j)   = cmdary(j-1)
               typary(j)   = typary(j-1)
               j = j-1
             else
               exit
             end if
           end do

           ! place the event we are moving
           timary(:,j) = etim
           cmdary(j)   = ecmd
           typary(j)   = etyp
         endif

         ! call prntevt('       Sort: ', fevt,levt)
       enddo
       return
      end subroutine sortevt

      subroutine prntevt(tag, fevt, levt)
        implicit none
        integer   fevt, levt
        character (len=*) :: tag

        include 'chrvar.inc'

        ! local variables
        integer   l
        character (len=132) :: mssg

        do l=fevt,levt
          write(mssg(:15),'(3i5)') l,timary(:,l)
          mssg(16:) = " "//cmdary(l)//' "'//typary(l)
          mssg(len_trim(mssg)+1:) = '"'
          call message(tag//trim(mssg))
        end do
      end subroutine prntevt

!spday[d_]   := Mod[d - 60, 365]      (* DOY to spring day *)
!spdayl[d_]  := Mod[d - 61, 366]      (* leap DOY to spring day *)
!sd2doy[sd_] := Mod[sd + 59, 365] + 1 (* spring day to DOY *)
!sd2dly[sd_] := Mod[sd + 60, 366] + 1 (* spring day to leap DOY *)
!sm[m_]      := Mod[m + 9, 12]        (* month to spring month *)

!(* day to month conversion  first use doy->spday and then convert \
!spday -> month
!   conversion has 2 parts
!    1   is the 153 day repeat number
!    2   is the day within the repeat sequence *)
!
!dy2m[sd_] := Mod[IntegerPart[sd/153]*5 + IntegerPart[Mod[sd, 153]/30.5 +1] +1, 12] +1
!dy2m[spday[test]]
!dy2m[spdayl[testl]]
!
!(*month to FIRST spring day    note: no leap year correction*)
!mfsd[sm_] := 153*Floor[sm/5.] + Floor[Mod[sm, 5]*30.5 - 0.5] + 1
!mfsd[sm[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}]]
!mfld[sm_] :=  Mod[153*Floor[sm/5.] + Floor[(Mod[sm, 5] + 1)*30.5 - 0.5], 365]
!sd2doy[mfld[sm[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}]]]
!sd2dly[mfld[sm[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}]]]
