!               Copyright 1998 Colorado State University
!                       All Rights Reserved

! Check input data
!   This reads the default .100 files and checks the input record agaist the targer
!     This is an upgrade of the ckdata routine. Made it a module to allow more
!     descriptive methods to access the state variables.
!

 module chkdatstat

   implicit none
  ! private                      ! hide everything except the function calls
  ! public flowdat, nflows, stackhead, ierr

    ! State variables (should these be module variables)
    integer, parameter       :: lo = 16         ! number of characters alloted to .100 types
    character (len=20), save :: file            ! open .100 file name
    character (len=4),  save :: dbgfile = ''    ! debug flag (file name)
    character (len=lo), save :: par             ! saved input variable name
    character (len=lo), save :: opt             ! saved .100 option
    character (len=80), save :: string          ! stored input record
    integer, save            :: lenval          ! length of value read
    integer, save            :: istat           ! conversion flag on value
    integer, save            :: lpnam, lpar     ! length of saved variable name and string
    real, save               :: value           ! last parsed value
    logical, save            :: reread=.FALSE.  ! next call should return current record
    logical, save            :: leof=.FALSE.    ! last read hit EOF

    character (len=26), parameter :: upcase  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character (len=37), parameter :: lcword  = 'abcdefghijklmnopqrstuvwxyz0123456789_'
    character (len=3),  parameter :: del     = '(,)'
    character (len=5),  parameter :: ws      = ' 	'//achar(10)//achar(12)//achar(13) ! [# \t\f\r] perl \s = [space, tab, linefeed, formfeed, return]

!    integer, save        :: lpar2, lnam2    ! length of saved variable name
!    character, save      :: par2*lo         ! saved input variable name

 end module chkdatstat


      subroutine findopt(unit,routin,expect,m)

      use chkdatstat;

!      locates an option (name) in the file.
!          unit     the i/o unit to read
!          routin   is the input file name for file open and error reporting
!          expect   is the option name to be found MUST BE LEFT JUSTIFIED
!          m        is the number of values in the input option
!
!      Modifications
!        7/2016  K. Killian
!          corrected handling 'expect' string length and it no longer modifies itself
!
!        3/2015  K. Killian
!          copied to stand alone routine for clarity

      integer,          intent(IN) :: unit
      character(len=*), intent(IN) :: routin
      character(len=*), intent(IN) :: expect
      integer,   intent(INOUT)     :: m


!...Local variables
      integer            :: j,  sl, lxp
      character          :: mssg*180
      character(len=10)  :: crpchar
      integer (kind=4)   :: crophash
      logical            :: cimatch

!...Functions
      character(len=lo)  :: nextopt, ucexpct, ucopt

        opt = ' ';
        if (m .eq. 0) call abortrun('chkdata - called with zero skip length')

        leof = .FALSE. ! this should handle an EOF

        ! expect = ADJUSTL(expect)
        lxp = len_trim(expect)

        cimatch = .false.        ! case insensitive match found
        ucexpct = expect(:lxp)   ! copy the search string for the case insensitive search
        call setcase(ucexpct, 1) ! upper case the option

        ! write(*,*) unit,routin,'  :  expect "',expect,'"',lxp,'"',expect(:lxp),'"', m
        do while (opt .ne. '?--?')
          sl = 0; opt = nextopt(unit,sl,string)
          if (expect(:lxp) .eq. opt) return
          ! try an upper case match if one hasn't been found
          if (.not. cimatch) then
            ucopt = opt;
            call setcase(ucopt, 1);
            cimatch = (ucexpct .eq. ucopt) ! found an upper case match
          endif
        end do

        ! CASE INSENSITIVE match found.
        ! rewind and locate that option
        if(cimatch) then
          opt = ' ';
          rewind abs(unit)
          do while (opt .ne. '?--?')
            sl = 0; opt = nextopt(unit,sl,string)
            ucopt = opt; call setcase(ucopt, 1);
            if (ucexpct .eq. ucopt) then
              call message("Warning: event type '"//trim(opt)//"' not found in "//trim(file)//": Using '"//expect(:lxp))
              return
            endif
          end do
        endif

        !  SEARCH FAILED
        ! log was what was found and terminate the run
        call message("Option '"//expect(:lxp)//"' not found in file "//trim(file))
        call message('     the following options were located:')
        rewind abs(unit)
        sl = 0
!       loop through the file printing the option names
225     par = nextopt(unit,sl,string)
        if (par .eq. '?--?') call abortrun('Unexpected EOF searching '// &
                   trim(file)//" for event option '"//expect(:lxp)//"'")

          sl = m
          mssg = "'"//trim(par)//"'"
          j = 17
          mssg(j:) = ':'
          if(file(:4) .eq. "crop") then
            mssg(j:) = '('//trim(crpchar(crophash(par)))//'): '
            j = 31
          endif
          call message(mssg(:j)//string(1:30))
        goto 225

      end subroutine findopt


      function rawval()
        use chkdatstat;

        implicit none
        integer :: ii
        character (len=lo)  :: rawval

        ii = VERIFY(string, " 	") ! search for first nonwhite, space/tab, character
        !ii=1;
        !do while ((string(ii:ii) .eq. ' ' .or. string(ii:ii) .eq. '	') .and. ii .lt. len(string))
        !  ii = ii+1;
        !end do

        rawval = string(ii:lenval)
        return

      end function rawval


      real function chkdata(unit,routin,expect,m,stk)

      use chkdatstat;

      integer,   intent(IN)    :: unit
      character, intent(IN)    :: routin*(*)
      character                :: expect*(*)
      integer,   intent(INOUT) :: m
      real,      intent(INOUT) :: stk(*)

!...reads the input .100 files and returns the variable values
!   the subroutine has three modes.  The input variables take on different uses
!   depending on the mode
!       1) unit =0  returns the values from a stack array.  This option is
!                   included to simplify the calling routines.
!           unit     is zero
!           routin   unused
!           expect   unused
!           m        is the index in the stack array; incremented before return
!              NOTE: generates a fatal error if m is zero
!           stk      is the stack array
!
!       2) unit >0  returns the requested variable from the input unit
!                   the value is also stored in the stack attay stk if m >0
!           unit     the i/o input unit to read
!           routin   is the input file name for file open and error reporting
!           expect   is the variable expected to be returned
!                        if expect is a blank string chkdata returns the variable found
!           m        is the index in the stack array; incremented before return
!              NOTE: if m = 0 the stack is not set and m not incremented
!           stk      is the stack array
!
!     uses ReadReal to read the input record and return the data's value
!
!       modifications
!         2/2013 K. Killian
!         - modified to allow '_' in variable names for frac_nh4.
!        02/2013  K. Killian
!         -  explicitly limit routin interactions with lstroutin to 4 characters
!        10/2010  K. Killian
!         -  clear par before loading new values to prevent character carryover
!        10/2009  K. Killian
!         -  prevent stack counter from incrementing during 'findname'
!        04/2009  K. Killian
!         - modifications to 'findname' option to ensure that the names can
!           be read at the end of the option
!         - identify EOF during an option read
!        12/2008  K. Killian
!         - removed the option to return a previously read line. This was used
!           as a reread that didn't work for arrays.
!           It also allows array element to be labeled by the name without
!           including the unchecked index
!         - Added a 'findname' option to recover variable name/value without
!           effecting the stack or file position (automatic reread)
!         - insure that rereads are to same file
!        05/2007  K. Killian
!           moved the skip line function to nextopt where the reading is done
!        09/2006  K. Killian
!           Corrected a bug when routin is actually the file name
!        08/2004  K. Killian
!           Made a change that should simplify handling optional values.
!             Saved the value and name from the last call. Then if the next
!             call rerequested the same data then just return it.
!           Changed the call to ReadReal so it skips blanks.
!           Changed the error messages.
!        03/2004  K. Killian
!           Make sure that "expect" is left justified. The OMAD unit conversion
!           option would leave a non justified option name.
!         4/2003 K. Killian
!           minor string restructuring for g77 compiler
!        10/2002 K. Killian
!           moved the option reading code to a subroutine so it can handle
!           variable length data options
!         5/2002 K. Killian
!           initialized istat = 0; non zero initial values in istat can
!           induce reconvert modes in readReal.
!         5/2000 K. Killian
!           changed the error checking to handle IO errors during the read
!         3/2000 K. Killian
!           corrected a bug in handling a ReadReal conversion error
!        12/1999 K. Killian
!           changed ReadReal to the new case. Corrected a bug in the error
!           printout. Does checking in lower case.
!        K. Killian  5/22/95


!...Functions
      double precision   :: ReadReal

!...Local variables
      integer i, j, lch, skipspac, lexp
      character          :: mssg*180
      logical            :: findnam, trncexp
      character          :: option*28


      chkdata = tiny(value)      ! initialize to an innocuous value
      option  = ' '
      !write(*,*) "chkdata(",unit,m,stk(1),"  '",routin,"'   '",expect,"'"

      if (expect .eq. 'echo') then ! set the echo flag
        dbgfile = routin(:4)
        call message("")
        call message("Echo "//routin(:4)//" option: "//trim(opt))
        return
      elseif (dbgfile .ne. routin(:4)) then ! clear the debug on subroutine change
        dbgfile = ""
      endif

!
      ! find a name from the list
      ! Not sure this has ever been used
      findnam = .FALSE.
      if (expect .eq. 'findname') then
        findnam = .TRUE.
        expect  = " "
      endif

      if (expect(:2) .eq. 'O:') then
        option = expect
        expect = ' '
      endif


      !...return a value from the stack variable
      !     Put this option before the reread so it just returns the stack value
      if (unit .eq. 0) then
        if (m .eq. 0) call abortrun('(internal) chkdata - arguments unit and m both zero')
        chkdata = stk(m)
        if(.not. findnam) m = m+1
        reread = .false.  ! clear the reread flag
        return
      endif

      ! find the end of the expected variable name
      trncexp = .false.
      lexp = index(expect,'(')-1; if(lexp <= 0) lexp = len_trim(expect);
      ! is there an alias; It will searched for so this can also be a truncated version
      i = index(expect,'|');
      if(i > 0) then
        trncexp = .true.
        lexp = i-1
        if(lexp < -1) call abortrun('malformed .100 variable name:'//expect)
      endif

!     handle the reread condition
      if(expect .eq. 'reread') then
         reread = .TRUE.
         chkdata = value
         goto 150

!     If we are reading optional variables and hit the of the file.
      elseif(leof .and. expect .eq. ' ') then
        reread  = .false.                 ! clear the reread parameter
        par     = ' '; lpnam=1; lpar=1    ! clear the saved parameter
        if(expect .eq. ' ') expect= 'EOF' ! return EOF file
        chkdata= -tiny(chkdata)           ! return EOF file

!     look for the return command
!     WATCH THIS: The file names really MUST be correct!!
      else if(reread) then ! .and. lstroutin .eq. routin(:4)) then
        ! clear the reread flag unless this is a name test read
        if(.NOT. findnam) then
          reread = .FALSE.
        endif
        chkdata = value

        if(m .ne. 0) then
          stk(m) = value
          m = m+1
        end if

        ! set/check the expect value
!        if     (istat.eq.0) then
!          mssg = "parsing data for "//trim(file)//";"//trim(expect)//" in line '"//trim(string)//"'"
!          call abortrun(mssg)
!
!        elseif (expect .eq. ' ') then
        if     (expect .eq. ' ') then
          expect = par ! (:lpar)

        elseif (expect(1:1) .eq. '-') then
          if(index(par,expect(2:lexp)).eq.0) then
            chkdata= tiny(chkdata) ! return fail
          endif

        elseif (index(par(:lpar),expect(:lexp)).eq.0) then
          call message ('bad record: '//trim(string))
          mssg = 'variable '//trim(expect)//' not found in '//trim(file)
          call abortrun(mssg)
        endif

      else

!...Read data from the input file
!       set lch = 0 so ReadReal will skip comments
101     skipspac = 0
        lch = 0
        istat = 0
        value   = ReadReal(unit, string, lch, istat)
        ! write(*,*) "value ",value," istat ",istat,"  string:",string
        chkdata = value
        lenval  = len_trim(string(:lch-1))

        if(istat.eq.0  .and. string .eq. 'debug') then
          dbgfile = routin(:4)
          call message("")
          goto 101 ! debug request for this option; redo the read
        else if(istat.le.0) then
           par = string(max(1,lch):)//' '
           i   = 2
           if(lch .le. 0) leof = .TRUE.  ! store the EOF condition
        else
           leof = .FALSE. ! this should handle an EOF
!VERIFY(STRING, SET, back)  the position of the first character in STRING not in SET.
!SCAN(STRING, SET, back)    the position of the first character from SET in the string STRING.
!                           If BACK is true, you will get the rightmost such character.

           !...extract the variable name
           par = ' '
           lpnam=0
           i = 1
           do while (i.le.len(par) .and. lch.le.len(string))
             j = index(lcword//del,string(lch:lch)) + index(upcase,string(lch:lch))
             if(string(lch:lch) .eq. '(') then
               skipspac = i
               if(lpnam .eq.0) lpnam = i-1
             else if(string(lch:lch) .eq. ')') then
               skipspac = 0
             endif

             ! terminate the parameter on a comment in parenthesis (malformed subscript)
             if(skipspac .ne. 0 .and. string(lch:lch).eq.' ') then
               i = i-1

             ! allowed characters in the name
             else if(j .ne. 0) then
               if(j<=26) then
                 par(i:i) = lcword(j:j)     ! lower case
               else
                 par(i:i) = string(lch:lch) ! copy the character
               endif

             ! first character is junk (probably a single quote)
             else if (i .eq. 1) then
               i=0
             ! remove terminal junk (single quotes, tabs, comments, ...)
             else
               par(i:) = ' '
               if(lpnam.eq.0) lpnam = i-1
               exit ! terminate the loop
             endif
             lch = lch+1
             i = i+1
           end do
           lpar = len_trim(par)
        endif

        ! debug print
        if (unit.ne.0  .and.  dbgfile .ne. '') then
          if(expect .eq. ' '  .and.  option .ne. ' ') then
            write(mssg,'(i3,2x,a,":",a,t26,"found: ",f10.4,2x,a,t62,2a)') m,routin, &
             option,chkdata,       ' "'//par(:lpar)//'"','-> "',trim(string)//'"'
          else
            write(mssg,'(i3,2x,a,":",a,t26,"found: ",f10.4,2x,a,t62,2a)') m,routin, &
              trim(expect),chkdata,' "'//par(:lpar)//'"','-> "',trim(string)//'"'
          endif
          call message(mssg)
        endif

        ! Test name for expected name
        mssg = ' '
        if (expect .eq. ' ') then
          expect = par
          if (findnam) then
            reread = .TRUE.
          endif

        else if(istat.eq.0) then
          mssg = "parsing data for "//trim(file)//";"//trim(expect)//" in line '"//trim(string)//"'"

              !(trncexp  .and.  index('|'//expect//'|', '|'//par(:lpnam)//'|')     ! truncated length compare
        else if(trncexp  .and.  index('|'//expect, '|'//par(:lpnam)) .eq. 0.) then ! search for an option
          mssg = "variables "//trim(file)//";"//trim(expect)// &
                 " not found in line '"//trim(string)//"'"

        else if(.not.trncexp  .and.  expect(:lexp) .ne. par(:lpnam)) then   ! full name compare
          mssg = "variables "//trim(file)//";"//trim(expect)// &
                 " not found in line '"//trim(string)//"'"
        endif
        if(mssg .ne. ' ') call abortrun(mssg)

!       prevent the findame command from incrementing the stack
        if (m.ne.0  .and.  .not. findnam) then
          stk(m) = chkdata
          m = m+1
        end if

      endif

150   continue

!      if(prntflg .and. m.ne.0) then
!        write(mssg,*) chkdata
!        call message(expect//': '//par//"   "//trim(mssg))
!      endif

      return
      contains


      end function chkdata



      function nextopt(unit,m,bffr)
      use chkdatstat;
      implicit none

      character (len=lo)  :: nextopt
      character*(*) bffr
      integer unit, m

!...steps through the input .100 files and returns the next option name
!
!    Option names must now be determined by the structure of the line. The
!    following rules are used to determine option names:
!    a) may be quoted
!    b) must be <= lo, (16) characters, excluding quote characters.
!    c) unquoted strings starting with a character other than numbers
!       (digits 0-9 plus -, ., and +,) are options
!    d) unquoted strings consisting of a number followed by an alpha string
!       without a white delimiter are options (ex:  15A)
!
!    Using these rules the following option line would not be recognized:
!       2957    Modern yield Corn cultivar       (not)
!       2957  # Modern yield Corn cultivar       (valid)
!       2957                                     (valid)
!       '2957'  Modern yield Corn cultivar       (valid)
!       2957A   Modern yield Corn cultivar       (valid)
!       A2957   Modern yield Corn cultivar       (valid)
!       MYCC    Modern yield Corn cultivar       (valid)
!
!     uses ReadReal to parse the numbers
!
!    Arguments
!       unit  logical unit to read
!       m     lines to skip between options. 0 skip for start of file
!             The routine will work with 0 but it just evaluates every data line
!       bffr  The logical record where the option is found
!
!    Modifications
!       written  10/2002 K. Killian
!       05/2007  K. Killian
!          moved the skip line read to nextopt where the reading is done
!       7/2003 K. Killian
!         converted the option length to lo characters
!         Used READSTRING instead of trying to duplicate the quoting and
!         comment capabilities.
!       1/2003 K. Killian
!         writes an error message if an option is partially quoted
!         comment changes .
!         expand TABS in options thus effectively terminating option names
!         Imbedded spaces (or tabs) were tolerated in Century 4 and you might
!         still be able to quote an option with spaces. However, TABS display
!         as 4 - 8 spaces thus visually filling out the field.
!      10/2002 K. Killian
!         comment changes and simplified parsing rules
!


! ... Function types
      double precision ReadReal
      logical          isblnk
      character        READSTRING*40

! ... Local variables
      real      vread
      integer   i, lstr, lch, instat, stbeg, stend
      character TAB*1, FRMAT*9
!        TAB is a tab character  "CHAR(9)"
      parameter (TAB = '	')

      data FRMAT/'(   (/)) '/

      if(m.gt.0) then
!...prepare the format statement
        write(FRMAT(2:4),'(I3)') m-1
!        FRMAT(2:4) = '  '//char(ichar('0')+mod(m-1,10))
!        if (m-1.ge.10)  FRMAT(3:3) = char(ichar('0')+mod((m-1)/10,10))
!        if (m-1.ge.100) FRMAT(2:2) = char(ichar('0')+(m-1)/100)
         read(abs(unit), FRMAT, end=800)
      endif


!20      read(abs(unit), '(a)', iostat=i, end=200) bffr
20      lch = 0
        ! NOTE: set the case option to match evntyp read case in readblk
        nextopt=  READSTRING(abs(unit),bffr,lch,stbeg,stend,.TRUE.,0)
        if(nextopt .eq. '?--?') return

!       Find the last non-commented, non-space character in the buffer
        lstr = len_trim(bffr)

!       We have a label If the first character is not consistant with a number.
!          NOTE: assume a leading white space COULD be part of a numeric entry
        if(index('0123456789-.+ '//TAB,bffr(1:1)) .le. 0) goto 500

! ...   It could be a number. Rather than rewrite the rules, use ReadReal to
!       parse the number and find the next "word" in the input line.
        instat = 0
        lch = 1
        vread = ReadReal(0, bffr, lch, instat)

!       If it was not a complete number (instat <=0) OR if the entire string
!       is a number (a number without a variable name) then it must be a label
        if(instat .le. 0  .or.  lch .ge. lstr) goto 500

!       If we succesfully parsed a leading number and we don't have a white
!       delimiter then we will call it a label.
!        NOTE: The foolproof test would be to test the variable but ...
        if(.not.isblnk(bffr(lch-1:lch-1))) goto 500

!       If we made it here this MUST be a data line so read another line
        goto 20

! ...   terminate any label at a TAB character   "CHAR(9)"
500     i = index(nextopt,CHAR(9))
        if(i .ne. 0) nextopt = nextopt(:i-1)
        return

! ...   return the EOF indicator
800     nextopt = '?--?'
        return
      end



      subroutine oldopen(unit,filnam)

      use chkdatstat;

      character*(*) filnam
      integer unit
!
!      The routine opens an existing file.  The local directory is tested then
!       the remote directory is tested if the library path is specified
!
!        unit is the logical unit
!        filnam is the file that should be read.
!
!      Modifications
!        3/2015  K. Killian
!          clears the chkdata state variables to prevent reread across file calls
!        4/2003  K. Killian
!          minor changes to ensure compatability with g77

!...Local variables
      integer   :: status, lnam

! this common contains the library path
      character filpath*256
      integer plen
      common /libpath/plen,filpath
      save /libpath/
!
! check the local directory
!
         lnam = len_trim(filnam)
         open(unit=unit, file=filnam, status='OLD', iostat=status)
!
! prepend the path and check that directory
!    look for null string end:  file=filpath(:index(filpath,char(0))-1)//filnam
!
         if (status.ne.0 .and. plen .gt. 0) then
            open(unit=unit, status='OLD', iostat=status,file=filpath(:plen)//filnam(:lnam))
         endif

         if (status.ne.0) call abortrun('opening file: '//filnam(:lnam))

         ! New file; clear the pending read state variables
         file      = filnam                    ! save the file root for debug
         dbgfile   = ""                        ! clear debug flag input file name
         reread    = .false.                   ! not going to reread the stored line
         par       = ' '; lpar = 1; value = 0; ! clear saved variable and value
         opt       = ' '                       ! saved .100 option
!         lstroutin = ' '                      ! clear file  stored data
         string    = ' ';                      ! clear stored input record (probably not needed)
         leof      = .false.                   ! not read EOF

         return

      end



      integer function stakfind(tomatch,SDEPTH,Tstack,stackname)

       integer SDEPTH
       integer Tstack(SDEPTH)
       character*(*) tomatch, stackname(SDEPTH)

!    A routine to manage the option stack for the input routines
!    It sorts the stack keeping the last option at the front of
!    the stack to help prevent trashing
!
!    Modifications:
!     6/2007   KLK
!      initial code

       integer i, j, m

       stakfind = 0
       j = 0
!...Check the stack for the data
       do 15 m = 1, SDEPTH
         stakfind = Tstack(m)

!        record the last unused slot
         if(stakfind .gt. SDEPTH  .or.  stakfind .le. 0) then
           stakfind = m
           goto 20

!        Look for a matching entry
         elseif (tomatch .eq. stackname(stakfind)) then
!          roll the match to the front of the list
           do 10 i = m, 2, -1
10           Tstack(i) = Tstack(i-1)
           Tstack(j+1) = stakfind

!          return the index into the data arrays
           return
         endif
15     continue

!      the stack was full; push the stack dropping the oldest entry
       stakfind = Tstack(SDEPTH)

!      push the stack dropping the oldest entry
20     continue
       do 25 i = SDEPTH, 2, -1
25       Tstack(i) = Tstack(i-1)

!      put the new entry in front
       Tstack(1) = stakfind

!      Signal that the option is not in the data arrays
       stakfind = -stakfind

      return
      end
