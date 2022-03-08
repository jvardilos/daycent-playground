
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... IRRGIN.F

      subroutine irrgin(tomatch,curirri,domirri, irintvl)

      implicit none
      include 'parcp.inc'

! ... Argument declarations
      character (len=28)   :: tomatch, curirri
      integer, intent(out) :: irintvl
      logical, intent(in)  :: domirri

! ... Read in the new irrigation type
!
!       auirri
!            = 0 Fixed Irrigation; applied regardless of available water
!            = 1 irrigate top 30 cm to field capacity
!            = 2 irrigate with a specified amount of water applied
!            = 3 irrigate top 30 cm to field capacity plus PET
!            = 4 irrigate rooting zone to field capacity
!                ALL auto irrigation applied when available water fraction < fawhc
!
! Allowed immediate fertilizer input where the values are coded in match
!  All fields are preceeded by a real value. No field separator is required
!  but they can be separated by commas for readability
!  Defined irrigation fields:
!        A auirri   irrigation type
!        F fawhc    fraction available water to initiate irrigation
!        C water    added in cm type 0 irramt   type 2 irraut
!        L aintvl   irrigation interval, number of days to maintain auto irrigation
!
! modifications
!     4/2015 K. Killian
!       made irraut an optional variable to maintain compatibility. There is no reason
!       to maintain 2 irrigation amounts

!...Functions
      double precision parsImedFld
      real             chkdata

! ... Local variables
      integer            :: unit, is, me, ms, nv, istat, sordr
      real               :: tmpval, stk(1)
      character (len=16) :: varfnd
      logical            :: warn
      character (len=13),parameter :: keylist = "AFCL?afcl";
      character (len=4), parameter :: evnt = "IRIG";

!...stack variables
      ! Number of lines to read for each type
      integer, parameter           :: REQLNS = 3, OPTLN = 1, SDEPTH = 6
      real, save                   :: stack(REQLNS+OPTLN, SDEPTH)
      integer, save                :: Tstack(SDEPTH) = 0
      character (len=28), save     :: stackname(SDEPTH) = '####'
      integer                      :: stakfind, m
      character (len=6), parameter :: callr = 'irrgin'

        ms = 1
        me = len_trim(tomatch)

        if(tomatch .ne. curirri) warn = .true.

        ! write(*,*) "irrgin=", tomatch, curirri, domirri, irintvl

        if(tomatch(:1) .eq. '(') then
          if(tomatch(me:me) .ne. ')') call  &  ! OOOOPS  the end of the data was lost
                     abortrun('unterminated immediate irrigation '//tomatch)

          unit = 5;
          is = 1; ! not used for immediate input but define for consistancy
          m  = 2
          me = len_trim(tomatch(:me-1))
          if(me<=m .and. warn) call message("Warning: ignoring empty immediate irrigation: "//trim(tomatch))

          ! clear the irrigation variables
          auirri  = 0
          fawhc   = 0.8 ! changed to ensure auto irrigation does something
          irramt  = 0
          irintvl = 1;

          sordr=0 ! set the sort order to 0
          do while (m.le.me) ! the last character is the terminator

            !treat tab as a soft delimiter
            if(tomatch(m:m) .eq. '	') then; m=m+1; cycle; endif

            ms = m ! save the field start for error reporting
            tmpval = parsImedFld(nv, tomatch, sordr,istat,m, ms,me,keylist,evnt)

            select case (nv)
              case (1) ! A  irrigation type
                auirri = int(tmpval)
              case (2) ! F  fraction available water
                fawhc = tmpval
              case (3) ! C  Cm     water added
                irramt = tmpval
              case (4) ! L  length auto irrigation interval
                irintvl = int(tmpval)
              case DEFAULT
                call abortrun('Unknown '//evnt//' key '//tomatch(ms:)//' in event '//tomatch)
                !call message('Warning: Unknown fertilizer key '//tomatch(ms:)//' in event '//tomatch)
                !m = m+1; cycle; ! exit
            end select

            if(istat.eq.0) then          ! the value wasn't parseable
              call abortrun('Missing immediate irrigation value before: '//tomatch(m:m)//" in "//trim(tomatch))
              !call message('Warning: Missing immediate fertilizer value: '//tomatch(m:m)//" in "//trim(tomatch))
              !m = m+1; cycle;
            endif
          end do
        else

          ! Check the stack for the data
          is = stakfind(tomatch,SDEPTH,Tstack,stackname)

          if(is.lt.0) then
            is = abs(is)
            stackname(is) = tomatch ! record the name in the name array

            unit = 11
            call oldopen(unit,'irri.100')   ! open the input file
            call findopt(unit,'irri.100',tomatch,REQLNS)  ! find the option in the irri file
          else
            unit = 0
          endif

          m = 1
          auirri = int(chkdata(unit,callr,'auirri',m,stack(1,is)))

          fawhc  = chkdata(unit,callr,'fawhc',m,stack(1,is))

          varfnd = ''
          irramt = chkdata(unit,callr,varfnd,m,stack(1,is))
          if(unit .ne. 0) then
            if(index(varfnd,'irraut') .ne. 0) then
              ! varfnd = 'findname'
              ! value = chkdata(unit,'cropin',varfnd,m,stack(1,is))
              m=m-1 ! We are going to reset the stack
              tmpval = chkdata(unit,callr,'irramt',m,stack(1,is)) ! read without the stack
              if(auirri .eq. 0  .and.  tmpval .gt. 0) then
                irramt = tmpval
                stack(m-1, is) = irramt
              endif
            elseif(index(varfnd,'irramt') .eq. 0) then
              ! This is just an obtuse way of triggering the standard error message
              irramt = chkdata(unit,callr,'reread',0,stk) ! reuse value
              irramt = chkdata(unit,callr,'irramt',0,stk)
            endif
          endif

          ! Added the auto irrigation length, KLK - 04/14/15
          varfnd = ' '
          nv = irintvl
          irintvl = chkdata(unit,callr,varfnd,m,stack(1,is))
          if(unit .ne. 0) then
            if(varfnd .ne. 'aintvl') then
             irintvl = nv
             stack(m-1, is) = nv
            endif
          endif

        endif

        if(unit .ne. 0) close(unit)
        curirri = tomatch

        ! check the irrigation interval against other inputs.
        if(domirri) then   ! reset things for the monthly 4 part application
          irintvl = 28;
          if(auirri .eq. 0) irramt = irramt/4. ! get the weekly application rate
        elseif(auirri .eq. 0) then
          if(irramt .eq.0) then ! This is a clear event.
            irintvl = 0  ! don't need multiple zero applications stop irrigation period
          else
            irintvl = 1  ! ignore irrigation interval input for
          endif
        endif

! write(*,*) "return irrgin:",unit," auirri=",auirri,"   fawhc=",fawhc,"   irramt=",irramt, "   irintvl",irintvl
        return

      end subroutine irrgin
