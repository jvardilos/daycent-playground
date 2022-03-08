
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... SCHEDL.F

      subroutine schedl(curday)

!     Determine the current scheduling options from the schedule file stack
!
!   The structure was modified to use the event parser for both block and
!   global events. The routine is basically a do while loop. The do while control
!   logic is to complex to sit in a statement so it occupies the first section of
!   the routine. The control logic evaluates the active block, global events and
!   returns if no events are current.
!   the second section is the event parser. It is "called" by setting the event
!   and jumping (goto) to the case statement. It returns to the top of the event
!   of the do-while loop.
!
!
!    Modifications
!     Sep14: K. Killian
!              Made the case parser a local subroutine. Greatly helps readability!
!              Modified the FLOD event to also set the DRAIN parameter to 0.3.
!              If the DRAN event is also set on the same day it will use that
!              value instead of the default.
!     Oct10: K. Killian
!              Modified the DRAN event so that valueF will set the default
!              flooded drain state. A drain on the same day will set the
!              DRAIN for that FLOD.
!              Changed the default flooded drain value to 0.0
!     Sep14: K. Killian
!              Modified the FLOD event to also set the DRAIN parameter to 0.3.
!              If the DRAN event is also set on the same day it will use that
!              value instead of the default.
!     Aug12: K. Killian
!              Added the FLOD event for flooded rice . This sets
!              FLOD    watertable    watrflag
!              0         0 (false)     0
!              1         1 (true)      0  rain accumulates without drainage
!              2         1 (true)      0  automatically add water to saturation
!              A DRAN event on the same date will set flooded  it will use that
!              value instead of the default.
!     Jul12: K. Killian
!              added a persistant auto fertilizer, AFRT, event. The AFRT option
!              is the real aufert value.
!              The value is saved in the local variable savafrt which is used
!              to set aufert if the value is non-zero. savafrt persists until
!              reset by another AFRT event or cleared by a FERT event.
!     Apr11: K. Killian
!             F95 converted event evaluator from i ..elseif  to select .. case
!             F95 added save attributes to variable declaration
!             simplified the control loop portion; comment changes
!     Dec10: K. Killian
!             modified global relative event input to simplify periodic events
!             %n     every n years from the start of the run
!             Intended for long repeat events. Disturbances to natural ecosystems
!             due to catastrophic events are likely candidates.
!     Mar09: K. Killian
!            - Dropped fltary from the common. For those events the floating
!              value is converted from the character string. This is a little
!              less efficient but the conversion load is small and it saves
!              managing this extra, virtually empty, array
!     Jul08:  K. Killian
!              Modifications to the global event index code
!     5/2005  K. Killian
!              Provide an undo function to the DRAN event which restores DRAIN
!                DRAN -1 swaps the current value and undodran(1)
!                DRAN -2 loads the site file value
!              undo's are limited through extends since previous value is not passed
!    11/2004 : Kendrick Killian
!              initialize the forest CO2 effect after reading in the tree
!    03/2004 : Kendrick Killian
!              calculate block year crtyr from timvar integer variable cyear
!    07/2003 : Kendrick Killian
!              changed the array logic so global events superceed block events
!              corrected a bug treating mixed global/block events
!    09/2002 : Kendrick Killian
!              corrected a bug with blocks having less than 2 events
!              simplified the crtyr calculation
!    06/2000 : Kendrick Killian
!              added a method of running a "global" events.
!    01/14/99 : Kendrick Killian
!              added a drainage event.
!
!NOTE Several century changes have been copied for reference and commented

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'zztim.inc'

! ... Argument declarations
      integer curday

! ... Local variables
      integer       :: crtyr, evpt
      real          :: daydran
      real,save     :: satdran = 0.01 !default saturation, FLOD, DRAIN value
      logical       :: fndthrv
      logical       :: newsenm  !local to flag SENM set in this call
      logical       :: inflod, indran ! local flags that this step changed DRAIN
      ! calculations can change aufert. Save values to restore the next use
      real, save    :: savedfert = 0.0 ! restore normal FERT events
      real, save    :: savafrt = 0.0   ! restore persistant AFRT events
      integer, save :: nxtaceq(2) = (0, 0) ! date of next equilibrium acceleration
      integer, parameter :: lco = 28;
      character(len=132) :: mssg
      character(len=lco), save :: curcult = ' ', curfert = ' '
      character(len=lco), save :: curfire = ' ', curgraz = ' '
      character(len=lco), save :: curharv = ' ', curirri = ' '
      character(len=lco), save :: curomad = ' ', curtrm  = ' ', curafrt = ' '

!...Check forest doflst to reset forgrw; done here so that forest grows
!     through the last month of growth
      if (doflst) then
        forgrw = 0
        dofngw = .false.
!        doflst = .false.
      endif

! ... Reset do variables to false
      doaceq  = .false.
      docult  = .false.
      doerod  = .false.
      dofert  = .false.
      doflst  = .false.
      dofone  = .false.
      dofrst  = .false.
      dograz  = .false.
      dohrvt  = .false.
      dolast  = .false.
      doomad  = .false.
      doplnt  = .false.
      if(dosene) senmschd = .false. ! clear pending senescence if handled
      dosene  = .false.              ! clear triggered SENM
      dotrem  = .false.
      fndthrv = .false. ! reset tree harvest flag found flag; leaves dothrv intact
      dofire(CRPSYS) = .false.
      dofire(FORSYS) = .false.
      dofire(SAVSYS) = .false.
      aufert = savafrt ! 0.0
      harmth = 0
      evpt   = 0
      newsenm= .false.              ! No events read yet so this must be false

      daydran= -1.                  ! No DRAN event this step
      inflod = .false.
      indran = .false.
      if(irrcnt  .eq. 0) domirri = .false.
      if(senecnt .ne. 0) senecnt = mod(senecnt + 1,31) ! if active, update the days since senescence  cak - 04/17/03

!...Calculate block year, crtyr, from the new integer year
!   this replaces a number of floating point conversions, patches, and bad code.
      crtyr = mod((cyear - strtyr), rptyrs)+1
      ! write(*,*) "schedl",crtyr, curday

!...return if there are no events in the block (pathological but ...)  KLK
      if (ttlind .lt. 1) return

      ! auto schedule an equilibrium event; equivalent to an ACEQ event.
      ! schedule if the count is active, there is a day for the event, and the date matches
      if (aceqcnt > 0 .and. nxtaceq(2) > 0  .and.  nxtaceq(1) == cyear .and. nxtaceq(2) == doy) then
        call EqAccel()
      endif

10    continue

!...Determine if evtptr needs reset to 1 (corrected for trivial blocks KLK)
!     if the event pointer is beyond the last one   OR
!        block date is beyond the last block event     then reset the pointer
      if (evtptr .gt. ttlind .or. crtyr .gt. timary(1,ttlind) .or.       &
            (crtyr .eq. timary(1,ttlind) .and.                           &
             curday .gt. timary(2,ttlind))) then
       evtptr = 1
       if(evpt .eq. 1) return ! special case:  return if we just did 1
      endif

!...Look if an event in timary matches the current time
!     First check the block events
      if(crtyr.eq.timary(1,evtptr) .and. curday.eq.timary(2,evtptr))then
        evpt     = evtptr
        if(dbgpsch) mssg = "      "
        call parsevnt(evpt,cmdary(evpt),typary(evpt))
        evtptr   = evtptr +1
        goto 10        ! return the event loop control

!     check for a global event
      else if (glevpt .gt. 0  .and.  glevpt .le. lenary) then
        ! Does the next global event match the current time
        if(timary(1,glevpt).eq.cyear.and.timary(2,glevpt).eq.curday)then
          evpt   = glevpt
          if(dbgpsch) mssg = "global"
          call parsevnt(evpt,cmdary(evpt),typary(evpt))
          glevpt = glevpt + 1
          goto 10        ! return the event loop control
        end if
      end if

      ! dothrv is persistant; this sets dothrv if a new event was seen.
      ! the event loop can clear it
      if(fndthrv) dothrv = fndthrv

      if (inflod .or. indran) call setdrain() ! FLOD and DRAN can set drain parameter
      return

! ===========================================================================
      contains
         subroutine parsevnt(evpt,curevt,evtype)
           integer       :: evpt
           character*(*) :: curevt, evtype

           integer :: istat, lch

           !... functions
           double precision :: ReadReal
           integer          :: ReadInt

           ! ***************************************************
           if(dbgpsch) then
             mssg = "Execute event "
             write(mssg(15:),'(i4,i6,4i5)') evpt,cyear,month,curday,timary(:,evpt)
             mssg(47:) = curevt//' "'//trim(evtype)//'"'
             ! mssg(len_trim(mssg)+2:) = curevt//' "'//trim(evtype)//'"'
             call message(mssg)
           end if

           ! ***************************************************
           select case (curevt)
             case ('CROP')
               if (curcrp .ne. evtype) then
                 call cropin(evtype)
                 call co2eff(time)
               endif

             case ('PLTM')
               doplnt = .true.
               plntday = timary(2,evpt)
               plntschd = .true.
               saveplntday = plntday
               if(.not. newsenm) senmschd = .false. ! new crop, clear previous SENM
               senecnt = 0 ! stop the no growth period after senescence

             case ('HARV')
               dohrvt = .true.
               if (curharv .ne. evtype) then
                 call harvin(evtype,curharv)
               endif
               hrvtday = timary(2,evpt)
               harvschd = .true.

             case ('FRST')
               dofrst = .true.
!               crpgrw = 1
               frstday = timary(2,evpt)
               frstschd = .true.
               savefrstday = frstday
               if(.not. newsenm) senmschd = .false. ! new crop, clear previous SENM
               senecnt = 0 ! stop the no growth period after senescence

             case ('LAST')
               dolast = .true.
               lastday = timary(2,evpt)

             case ('SENM')
               !dosene = .true. this is now only true when the event is triggered
               seneday = timary(2,evpt)
               senmschd = .true.
               newsenm= .TRUE. !this is today's SENM;

             case ('FERT')
               dofert = .true.
               aufert = savedfert
               savafrt = 0.0  !  ANY FERT clears the persistant autofert
               if (curfert .ne. evtype) then
                 call fertin(evtype,curfert)
                 savedfert = aufert
               endif
               fertday = timary(2,evpt)

             case ('AFRT')
               if (curafrt .ne. evtype  .or.  savafrt .eq. 0.0) then
                 lch   = 1
                 istat = 0
                 savafrt = max(ReadReal(0,evtype,lch,istat), 0.)
                 aufert  = savafrt
               endif

             case ('CULT')
               docult = .true.
               if (curcult .ne. evtype) then
                 call cultin(evtype,curcult)
               endif
               cultday = timary(2,evpt)

             case ('OMAD')
               doomad = .true.
               if (curomad .ne. evtype) call omadin(evtype,curomad, labtyp)
               omadday = timary(2,evpt)

             case ('IRRI')
               ! monthly irrigation; retained for compatibility
               domirri = .true.
               if (curirri .ne. evtype) then
                 call irrgin(evtype, curirri, domirri, irintvl)
               endif
               irrcnt = irintvl ! set the irrigation count
               irriday = timary(2,evpt)

             case ('IRIG')
               ! general dayCent appropriate irrigation
               ! - single irrigation, 0, is amount specified on day scheduled
               ! - auto irrigation ,1-4, is repeated for airrdy days or until
               !   another irrigation event, plant harvest, senescence, last, tlst
               domirri = .false.
               if (curirri .ne. evtype) then
                 call irrgin(evtype, curirri, domirri, irintvl)
               endif
               irrcnt = irintvl
               irriday = timary(2,evpt)

             case ('GRAZ')
               dograz = .true.
               if (curgraz .ne. evtype) then
                 call grazin(evtype,curgraz)
               endif
               grazday = timary(2,evpt)

             case ('EROD')
               doerod = .true.
               istat  = 0
               lch    = 1
               psloss = ReadReal(0,evtype,lch,istat) ! fltary(evpt, 1)
               erodday = timary(2,evpt)

             case ('DRAN') ! Modify the drain variable
               istat = 0
               lch   = 1
               daydran = ReadReal(0,evtype,lch,istat)
                 ! set water table flags
                 if(lch.le.len(evtype)) then
                 if(evtype(lch:lch) .eq. 's') then
                   satdran = daydran
                 endif
               endif
               indran  = .true.

             case ('FIRE')
               dofire(cursys) = .true.
               if (curfire .ne. evtype) then
                 call firein(evtype,curfire)
               endif
               fireday = timary(2,evpt)

             case ('OTRF')
               istat = 0
               lch = 1
               if(evtype(1:1) .eq. '(') lch =2
               otfrac = ReadReal(0,evtype,lch,istat)
               if(otfrac .ne. -1 .and. (otfrac < 0.0  .or.  otfrac > 1.0)) then
                 call message('ignoring undefined orchard fraction" '//trim(evtype))
                 otfrac = -1 ! 0<= otfrac <=1
               endif
               otffrc = ReadReal(0,evtype,lch,istat)
               if(istat .le. 0 .or. otffrc .lt. 0. .or. otffrc .gt. 1.) then
                  call message('assuming area weighted orchard fertilizer fraction: '//trim(evtype))
                  otffrc = otfrac ! assume a uniform fertilizer application
               endif
               if(cursys .ne. SAVSYS) then
                  call message('OTRF is the tree allocation event for mixed (SAVANNA) orchards: '//trim(evtype))
               endif

             case ('TREE')
               if (curtre .ne. evtype) then
                 call treein(evtype)
                 call co2eff(time)

                 frstage = 0
!                update the new leaf C/E ratio, ccefor, with the current CO2 level
!                Did the calculation here since schedl had all the values
                 ccefor(IMIN,LEAF,1:nelem)= cerfor(IMIN, LEAF, 1:nelem) * &
                                            co2cce(FORSYS, IMIN, 1:nelem)

                 ! PROBABLY UNNEEDED
                 ! as writtent pltlig uses tree and should be updated
                 ! BUT using the tree in the plant lignen is PROBABLY INCORRECT KLK 7 APR 2017
                 call cmplig(cursys,fligni,wdlig) ! Recalculate PLANT lignin
               endif

             case ('TREM')
               dotrem = .true.
               if (curtrm .ne. evtype) call tremin(evtype,curtrm)
               tremday = timary(2,evpt)

             case ('TFST')
               dofone = .true.
               forgrw = 1
               foneday = timary(2,evpt)
               dothrv = .false. ! cancel any remaining harvest flags

             case ('TPLT')
               dofone = .true.
               forgrw = 1
               dothrv = .false. ! cancel any remaining harvest flags
               call trplntin(evtype) !  read the initial conditions

             case ('TFLR')
               dofngw = .true.
               dothrv = .false. ! cancel any remaining harvest flags

             case ('THRV')
               istat = 0
               lch   = 1
               tfhrvfr  = ReadReal(0,evtype,lch,istat)
               fndthrv = .true.
!              check for the harvest fraction
               if(tfhrvfr .gt. 1.0  .or.  tfhrvfr .lt. 0.) call abortrun   &
               ("fruit harvest fraction out of bounds"//trim(evtype))

             case ('TLST')
               doflst = .true.
               flstday = timary(2,evpt)

             case ('FLOD')
               istat = 0
               lch   = 1
               watrflag = ReadInt(0,evtype,lch,istat)
               inflod = .true.
               if(watrflag .gt. 0) then ! set water table
                 watertable = 1
!                 watrflag = min(max(watrflag-1,0),1) ! this is really a logical (bool)
                 watrflag = min(watrflag-1,1) ! this is really a logical (bool)
               else           ! input flag < 0 clears water table; this was Cindy's DRAN event
                 watertable = 0
                 watrflag   = 0
               endif

             case ('ACEQ')
               call EqAccel()

           end select

         end subroutine parsevnt

         subroutine EqAccel()
           integer :: EqAcldat
           integer :: newyr

           doaceq = .true.;
           newyr = EqAcldat()
           if(newyr > 0) then
             nxtaceq = (/ cyear + newyr,  doy /)
           else
             nxtaceq = 0
           endif
           return
         end subroutine EqAccel


         subroutine setdrain !(drain, daydran, inflod, indran, watertable,
           ! logical    :: inflod, indran ! local flags that this step changed DRAIN
           ! integer    :: watertable
           ! float      :: drain, daydran

           if(inflod) then                  ! DRAN is different under FLOD
             if(watertable .gt. 0) then     ! set anaerobic if water table set
               undodran(3) = drain          ! save the current value
               if (indran) then             ! FLOD with DRAN to set drain value
                 drain = daydran            ! Set DRAIN to the users input
               else
                 drain = satdran            ! set DRAIN to default anaerobic value
               endif
             else                           ! drain the flood
               drain   = undodran(3)        ! Set DRAIN pre-FLOD value
             endif

           else                             ! dran without FLOD
             ! This code was in the DRAN block but
             if(daydran .eq. -2) then       ! reset DRAIN to the stored site value
                undodran(1) = drain
                drain = undodran(2)

             else if (daydran .eq. -1) then ! UNDO DRAIN
                ! swap the undodran(1) and DRAIN
                daydran = undodran(1)
                undodran(1) = drain
                drain = daydran

             else                           !       set DRAIN
                undodran(1) = drain
                drain = min(1., max(0., daydran))
!               check for the drain fraction
                if(drain .lt. 0.  .or.  drain .gt. 1.0 ) call abortrun    &
          &      ("Error: illegal drain fraction"//typary(evpt))
             endif
           endif
         end subroutine setdrain

      end
