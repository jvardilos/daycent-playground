
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine fertin(tomatch,curfert)

      implicit none

! ... Argument declarations
      character (len = 28) :: tomatch, curfert

!...Read in the new fert type
!
! Allowed immediate fertilizer input where the values are coded in tomatch
!  All fields are preceeded by a real value. No field separator is required
!  but they can be separated by commas for readability
!  Defined fertilzer fields:
!        N nitrogen
!        P phosphorus
!        S Sulfer
!        F Frac_NH4
!        I nitrogen inhibitor   integer part is ninhtm, fractional part is ninhib
!        A auto fert (stand alone)
!
! modifications
!     7/2016 K. Killian
!       Improved immediate parsing error trapping and reporting
!       Identify a decimal point used as a delimiter as an error
!     6/2015 K. Killian
!       corrected a bug in fert.100 parsing
!     4/2015 K. Killian
!       minor recoding to use select, do while loop, and the findopt subroutine
!       corrected resulting unlabeled input indexing
!     7/2013 K. Killian
!       corrected a bug that allowed an option without a number to be ignored.
!        An immediate input without a number generates a warning
!     2/2013 K. Killian
!       chkdata was modified to allow '_' in frac_nh4. frac_nh4 input is now enabled.
!       changed default NH4 fraction
!     6/2007 K. Killian
!       adapted to respond to the fertilizer immediate format
!         use the Century input management routines stakfind and chkdata
!         instead the older read and ckdata
!       Modified input so that frac_no3_fert is defined as 1- frac_nh4_fert.
!         this simplifies immediate input and simplified the possibility of
!         the N fraction normalization error. frac_no3_fert is now an optional
!         input to maintain file compatibility
!
!
      include 'const.inc'
      include 'fertil.inc'

!...Functions
      double precision ReadReal, parsImedFld
      real             chkdata

!...Local variables
      integer            :: sordr, nv, ms, me, is, m, unit, istat
      real               :: scaleopt
      double precision   :: tmpval
!      character (len=16) :: varfnd
      character (len=13),parameter :: keylist = "NPSFIA?npsfia";
      character (len=4), parameter :: evnt = "FERT";

!...stack variables
      ! Number of lines to read for each type
      integer, parameter           :: REQLNS = 7, OPTLN = 1, SDEPTH = 6
      real, save                   :: stack(REQLNS+OPTLN, SDEPTH)
      integer, save                :: Tstack(SDEPTH) = 0
      character (len=28), save     :: stackname(SDEPTH) = '####'
      integer                      :: stakfind
      character (len=6), parameter :: callr = 'fertin'

      ! write(*,*)"fertinST  tomatch '",trim(tomatch),"'  curfert '",curfert,"'"
!... Is this an immediate mode??
!    Note: Immediate mode data is in memory so there is nothing to be gained
!    by putting it on the stack
      tmpval = 0
      scaleopt = -1
      ms = 1
      me = len_trim(tomatch)

      if(index(tomatch,'*') .gt. 0) then
        if(tomatch(:1) .eq. '('  .and.  tomatch(me:) .eq. ')') then
          ms = 2
          me = len_trim(tomatch(:me -1))
        endif

        ! parse the initial part as a number
        istat = 0
        m = ms
        scaleopt = ReadReal(0, tomatch, m, istat)
        if(istat .le. 0) then ! numeric parse worked
          call message ('Warning: missing weighting factor in fertilizer option: '//trim(tomatch))
          scaleopt = 1
        endif

        do while(index(' 	*',tomatch(m:m)) > 0  .and.  m < me)
          m = m+1
        end do

        if(m >= me) then
          call abortrun("Missing fertilizer option in "//trim(tomatch))
        endif

        ms = m
      endif

      if(tomatch(:ms) .eq. '(') then
        if(tomatch(me:me) .ne. ')') call &  ! OOOOPS  the end of the data was lost
                  abortrun('unterminated immediate fertilizer '//tomatch(ms:me))

        m  = 2
        me = len_trim(tomatch(:me-1))
        if(me<=m) call message("Warning: empty immediate fertilizer: "//trim(tomatch))

! clear the fertilizer amounts
        feramt(N) = 0.
        feramt(P) = 0.
        feramt(S) = 0.
        !On Mar 1, 2013, at 8:53 AM, delgro wrote:
        !yes - what Ernie suggests should be our default assumption
        ! On Thu, 28 Feb 2013 19:58:11 -0700, Ernie Marx <erniemarx@gmail.com> wrote:
        ! What do y'all think about changing the default to 80% NH4 and 20% NO3?
        frac_nh4_fert = 0.8d0
        ninhtm = 0
        ninhib = 1

        sordr=0 ! set the sort order to 0
        do while (m .le. me) ! the last character is the terminator
          !treat tab as a soft delimiter
          if(tomatch(m:m) .eq. '	') then; m=m+1; cycle; endif

          ms = m ! save the field start for error reporting
          tmpval = parsImedFld(nv, tomatch, sordr,istat,m, ms,me,keylist,evnt)

          select case (nv)
            case (1:3)
              feramt(nv) = tmpval
              ! check for a labeled form separated by commas (6N,2P)
            case (4) ! H
              frac_nh4_fert = tmpval
            case (5) ! I
              ninhtm = int(tmpval)
              ninhib = tmpval - ninhtm
            case (6) ! A
              aufert = tmpval
            case DEFAULT
              call abortrun('Unknown '//evnt//' key '//tomatch(ms:)//' in event '//tomatch)
              !call message('Warning: Unknown fertilizer key '//tomatch(ms:)//' in event '//tomatch)
              !m = m+1; cycle; ! exit
          end select

          if(istat.eq.0) then          ! the value wasn't parseable
            call abortrun('Missing immediate "//evnt//" value before: '//tomatch(m:m)//" in "//trim(tomatch))
            !call message('Warning: Missing immediate fertilizer value: '//tomatch(m:m)//" in "//trim(tomatch))
            !m = m+1; cycle;
          endif
        end do
      else

        unit = 0
        ! Check the stack for the data
        is = stakfind(tomatch(ms:me),SDEPTH,Tstack,stackname)

        if(is.lt.0) then
          unit = 11
          is = abs(is)
          call oldopen(unit,'fert.100') ! open the input file
          call findopt(unit,'fert.100',tomatch(ms:me),REQLNS) ! find the option
        endif

        ! record the name in the name array
        stackname(is) = tomatch(ms:me)

        m = 1
        feramt(N) = chkdata(unit,callr,'feramt',m,stack(1,is))
        feramt(P) = chkdata(unit,callr,'feramt',m,stack(1,is))
        feramt(S) = chkdata(unit,callr,'feramt',m,stack(1,is))
        aufert    = chkdata(unit,callr,'aufert',m,stack(1,is))

        ninhib    = chkdata(unit,callr,'ninhib',m,stack(1,is))
        ninhtm    = nint(chkdata(unit,callr,'ninhtm',m,stack(1,is)))
        frac_nh4_fert = chkdata(unit,callr,'frac_nh4',m,stack(1,is))

        if (feramt(N) .gt. 0.0  .and.                                    &
            (frac_nh4_fert .gt. 1. .or. frac_nh4_fert .lt. 0.)) call     &
             abortrun('frac_nh4_fert out of bounds in fert.100!')

        ! don't read FRAC_NO3.
        ! This is the last value we don't even have to pretend to read it
        ! varfnd = ' '
        ! value = chkdata(unit,callr,varfnd,m,stack(1,is))
        ! if(varfnd .ne. 'FRAC_NO3') then; ; endif

        if (unit.ne.0) close(unit)

        ! scale the fertilizer rates if a factor was input
        if(scaleopt .ge. 0.)feramt = feramt * scaleopt

      endif

      frac_no3_fert = 0.5d0 - frac_nh4_fert + 0.5d0 ! make sure NO3 is normalized

!     set the auto fertilization
!      if (aufert .lt. 0) then; dosafrt = .TRUE.; aufert = abs(aufert); endif


      !write(*,*) "fert levels ",tomatch(ms:me),": ",feramt,aufert,ninhib,ninhtm,frac_nh4_fert,frac_no3_fert
      curfert = tomatch
      return

      end subroutine fertin


      double precision function parsImedFld(nv, tomatch, sordr,istat,m, ms,me,keylist,evnt)
        integer, intent(in)               :: ms, me
        integer, intent(inout)            :: sordr,nv, m, istat
        character (len=*), intent(in)     :: tomatch,keylist,evnt

!...parse immediate Field
!   Parses the value, provides delimiter screening and initial error checking
!   returns the value as the return and the index integer based on the keylist order
!   in one of the formal parameters.
!   (probably ought to make this a module and return a structure)
!
!   This is only part of the search loop, but since the rest of the loop
!   assigns the read value to the proper variable It is hard to encapsulate.
!
!   nv       the argument key number
!   tomatch  input string that is being parsed
!   sorder   sort order, this is the assignment order for non keyed input
!   istat    conversion status, returned to allow calling routine to accept a key
!            without a value
!   m        Pointer to current character in tomatch
!   ms       start of active field in tomatch
!   me       end of parsable data (not including close ')'
!   keylist  A list of variable key letters in the form UPPER?lower
!   evnt     Event type being parsed; used for error reporting
!
!   Called by
!   fertin, irrigin, treein
!

        double precision ReadReal

          sordr=sordr+1 ! increment sort order regardless of delimiter

          istat = 0
          parsImedFld = ReadReal(0, tomatch, m, istat)
          ! write(*,*) "tmpval =",parsImedFld, " del= '",tomatch(m:m),"'  istat = ",istat,"  m=",m,' read:',tomatch(ms:m-1)

          ! ERROR bad numeric conversion; this shouldn't happen
          if(istat.lt.0) call abortrun('parsing immediate '//evnt//': '//trim(tomatch))

          ! process a valid number

          ! find the variable for this number
           nv=mod(index(keylist,tomatch(m:m)),(len(keylist)+1)/2)
          if(nv .gt. 0) then
            m=m+1 ! key found; accept the key character
            ! allow 'key,', construction
            ! NOTE; nested if because fortran ifs don't short circuit.
            if(m .lt. me) then; if(tomatch(m:m) .eq. ',') m=m+1; endif
          else if(tomatch(m:m) .eq. '.' .and. index('0123456789.',tomatch(m-1:m-1))>0) then ! period as a delimiter is ambiguous
            ! this is probably unintended; error the parse and let the user try again
            call abortrun('ambiguous immediate '//evnt//' ['//trim(tomatch(ms:me))// &
             '] in:  '//trim(tomatch))
          else
            nv=sordr
          endif
          ! write(*,*) "Immediate ",evnt,parsImedFld," '",tomatch(m:m),"' ",m,istat,'nv=',nv,"  '",tomatch(m:me),"'"

         return
      end function parsImedFld
