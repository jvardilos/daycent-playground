
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... OMADIN.F

      subroutine omadin(tomatch,curomad, labtyp)

! ... Read in the new omad type
!
!     modifications
!     On Dec 11, 2014, at 4:06 PM, Keough,Cynthia <Cindy.Keough@colostate.edu> wrote:
!     OMADTYP input parameter:
!          1 organic matter added to the structural and metabolic pools, astlbl = C13/C14 concentration
!          2 organic matter added to som2c pool, astlbl = C13/C14  concentration
!          3 organic matter added to the structural and metabolic pools, astlbl = C13/C14 fraction
!          4 organic matter added to som2c pool, astlbl = C13/C14 fraction
!     30Nov12 KLK rewrote comments
!             changed omadtyp to Match Cindy K range definition
!               1 litter
!                   organic matter added to the structural and metabolic pools
!               2 decomposed material
!                   organic matter added to som2c pool
!             the fractional part is the mixing fraction between the types
!     6/2007 K. Killian
!       used the new stack management routine stakfind
!     02/2006  K. Killian, C Keough   30Nov12 KLK rewrote this comment
!       Add code to read the new flag, OMADTYP, which indicates the fraction
!       of the added material that is that is decomposed, or composted.
!       Decomposed material, 1, will be added directly to som2c. Traditional
!       omad, 0, is added to structural and metabolic litter pools.
!     03/2004  K. Killian
!       Added an "unit conversion" variation of the immediate mode. With this
!       option the lead number MULTIPLIES ASTGC in the the omad.100 option.
!       Use it to convert between the input practice and Century's g/m2 C.
!       -- Format   weight*option.
!          Imbedded spaces are allowed but you have a limited field.
!     4/2003  K. Killian
!       minor changes to ensure compatability with g77
!     1/23/2002 KLK
!       Added an immediate mode input. Parenthesis are required to signal the
!       immediate mode. The immediate mode format is:
!          (newASTGC, OPTION)
!        where
!           newASTGC  is the the amount to apply
!           OPTION    is the omad.100 option
!       newASTGC replaces ASTGC read from OMAD.100. However, a valid option
!       MUST be specified so the remaining ratios must be read from OMAD.100.
!       NOTE: The numeric parser can read the input without a delimiter, comma
!       or space but the delimiter is required if the lead characters of OPTION
!       is a D, C (or a number) that could be included in scientific notation.
!      (Examples:
!          (25M)      and   (25,M) are equivalent
!          (100E2P) is parsed as 1.0E4 with OMAD type 'P' not as (100,E2P).
!

! ... Argument declarations
      implicit none
      character*28        :: tomatch, curomad
      integer, intent(in) :: labtyp

      include 'const.inc'
      include 'parcp.inc'

!...function definitions
      double precision ReadReal
      real             chkdata
      character (len=16)  :: rawval

! ... Local variables
      integer   is, m, unit, istat
      real      newastgc
      logical   weight, lblconc
      character (len=28) ::  match
      character (len=16) ::  varfnd
      character (len=6), parameter :: callr = 'omadin'

!...stack variables
      ! Number of lines to read for each type
      integer, parameter           :: REQLNS = 6, OPTLN = 1, SDEPTH = 6
      real, save                   :: stack(REQLNS+OPTLN, SDEPTH)
      integer, save                :: Tstack(SDEPTH) = 0
      character (len=28), save     :: stackname(SDEPTH) = '######'
      integer                      :: stakfind

      newastgc = -1
      weight   = .FALSE.
      match    = tomatch
      lblconc  = .FALSE.

!... Is this an immediate mode??
!    For OMAD, this is a hybrid; immediate weighting of a defined profile
!    parse the weighting factor then continue with the rest of the routine
      m=index(match,'(')
      if(m .ne. 0) then
        is= index(match,')')
        if(is .eq. 0) then
!         OOOOPS  looks like the end of the data was lost
          call abortrun('unterminated immediate organic addition `'//tomatch)
        endif

        ! read the application rate (astgc weighting value) from the command line
        m=m+1
        istat=0
        newastgc = ReadReal(0, match, m, istat)
        if(istat .le. 0) then
          ! ERROR bad numeric conversion!! Internal error exceeded buffer
          call abortrun('parsing OMAD weighting factor: '//tomatch)
        endif

        if(m .ge. is-1) then
          call abortrun("Missing OMAD option in "//trim(tomatch))
        else if(tomatch(m:m) .eq. '*') then
          m=m+1
          weight = .TRUE.
        end if

        ! remove the immediate mode characters from the string
        match = adjustl(match(m:is-1))

   !write(*,*) "Immediate OMAD ",tomatch,istat, newastgc," * '", trim(match),"'"
      end if



      unit = 0
!...Check the stack for the data
      is = stakfind(match,SDEPTH,Tstack,stackname)

      if(is.lt.0) then
        unit = 11
        is = abs(is)
!... open the input file and find the removal option
        call oldopen(unit,'omad.100')
!       find the option in the file
        call findopt(unit,'omad.100',match,REQLNS)
      endif

!     record only the base option in the stack name array
      ! astgc modifiers take little time and tomatch must have changed
      stackname(is) = match

      m = 1
!       Did they specify the optional omadtyp
        varfnd = ' '
        omadtyp = chkdata(unit,callr,varfnd,m,stack(1,is))

!       if this is a file read then process the input
        if (unit .ne. 0) then
!         found omadtyp.   Check the range
          if (index(varfnd,"omadtyp") .ne. 0) then
            if    (omadtyp.ge.1. .and. omadtyp.le.2.) then
              omadtyp = omadtyp - 1.0
              lblconc = .TRUE.
            elseif(omadtyp.ge.3. .and. omadtyp.le.4.) then
              omadtyp = omadtyp - 3.0
              lblconc = .FALSE.
            else
              call abortrun("composted fraction, omadtyp, out of range in "//trim(tomatch))
            endif

!         did NOT find the optional input; use the default
          else
!           mark this line for the next read
            omadtyp = chkdata(unit,callr,'reread',m,stack(1,is))

!           the default is a litter type concentration
            omadtyp = 0.
          endif
          stack(m-1,is)  = omadtyp
        endif

        astgc = chkdata(unit,callr,'astgc',m,stack(1,is))
        if(weight) then
!           write (*,*) " Weighting OMAD C addition to ",newastgc,"*",astgc," = ",newastgc * astgc
           astgc = newastgc * astgc
        else if(newastgc .ne. -1) then
           astgc = newastgc
!           write (*,*) " Changing OMAD C addition to ",astgc
        end if


        astlbl = chkdata(unit,callr,'astlbl',m,stack(1,is))
        ! Convert concentrations to fractions
        if(lblconc) then
          if (labtyp .eq. 1) then !  fraction of C labeled 14C
            ! fraclabl = ((astlbl / 1000.0) + 1.0) / 1000.0
            astlbl = FRAC_C14 * (1.0 + astlbl / 1000.0)
          elseif (labtyp .eq. 2) then ! 13C fraction if 13C labeling
            astlbl = PEEDEE * (1.0 + astlbl / 1000.0)
          endif
          astlbl = astlbl / (1.0 + astlbl)
          stack(m-1,is)  = astlbl
        endif

        astlig = chkdata(unit,callr,'astlig',m,stack(1,is))
        if ((astlig .lt. 0) .or. (astlig .gt. 1.0)) then
          call abortrun('ASTLIG out of range [0.0, 1.0]; read:'//trim(rawval()))
        endif

        astrec(N) = 1.0 / chkdata(unit,callr,'astrec',m,stack(1,is))
        astrec(P) = 1.0 / chkdata(unit,callr,'astrec',m,stack(1,is))
        astrec(S) = 1.0 / chkdata(unit,callr,'astrec',m,stack(1,is))

        if (unit.ne.0) close(unit)

!       Modify the curomad so that it remembers to input value KLK
        curomad = tomatch  !  = match
!        write(*,*) "stack",is, stackname(is), stack(:,is)

      return
      end subroutine omadin
