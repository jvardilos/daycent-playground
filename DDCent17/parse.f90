! Modification:
!   02/13/18 : Kendrick Killian
!     rewrite isblnk;
!     convert lcase to a more general setcase and use it in READSTRING
!   02/13/13 : Kendrick Killian
!     corrected a potential infinite loop condition parsing quotes in comments
!   10/05 K. Killian
!     use single byte integer instead of the byte extension for gfortran.
!     NOTE: Sun f77 still requires byte, convert in SUN Makefile
!   06/05 K. Killian
!     removed out of bounds tests on strings
!   04/05 K. Killian
!     Significant simplification of the uncomment routine
!     Changes to error messages
!   08/03 K. Killian
!     corrected a bug caused string pointer not skip to trailing white space
!  4/2003  K. Killian
!     minor changes to ensure compatability with g77
!  5/2000  K Killian
!    made sure all error messages are routed through message and can be trapped
! 12/1999  K Killian
!  new routines    len_trim and ladjust versions of the f90 intrinsics
!                  FILRERR file read error routine
!  READSTRING      added a quote character capability and the ability to
!                  specify the case returned. Currently used to disable
!                  token parsing inside of parenthesis
!  ReadReal        changed INTEG to RNUM logical to simplify logic (?)
!                  other minor logic changes
!  READCLIN        added the ability to specify the returned case
!
!     Copyright 1992 - 2004    Kendrick Killian     All Rights Reserved

!****   FUNCTION ISBLNK
      logical function isblnk(line)
      character*(*) line

        !verify finds characters not in the set
        isblnk = (VERIFY(line, ' '//'	') .eq. 0); ! use the modern intrinsic

      end function isblnk


!**** CHARACTER FUNCTION SETCASE
      subroutine setcase(strng, ccase)
        ! change string case to upper >0 or lower <0 using an array manipulation
        ! this allows the compiler to use processor array instructions

        character (len=*), intent(INOUT) :: strng
        integer,           intent(IN)     :: ccase

        if (ccase == 0) return ! get out if no change requested
        if (ccase .LT. 0) then    ! lower case the token
          call lowercase(strng, len(strng))
        else
          call uppercase(strng, len(strng))
        endif

        return
       contains

        subroutine lowercase(caray, maxch)
          integer, intent(IN) :: maxch
          character (len=1), dimension(maxch), intent(INOUT) :: caray
         ! NOTE: the difference between ascii upper and lower case is bit 5
         !       set bit 5, + x20, to upper cases, clear, - x20, lower case
          caray = merge(achar(ibset(iachar(caray),5)),caray, (caray>='A' .and. caray<='Z')) ! & ! lower case
                ! ibits(dim(iachar(caray),ucs)-1,0,8)<26)  !  bit manipulation mask
                ! Commented the obtuse but interesting bit manipulation for the mask
                !    ibits(dim(iachar(caray),[ul]cs)-1,0,8)<26) ! [upper/lower]
                ! dim()-1 indexes character number to a 0 based upper/lower case count,
                ! ibits transforms index<0 to large positive
                ! then test the mask
                ! not sure whether the mit manipulations are faster than the > part
                ! of the obvious character range test
        end subroutine lowercase

        subroutine uppercase(caray, maxch)
          integer, intent(IN) :: maxch
          character (len=1), dimension(maxch), intent(INOUT) :: caray
          caray = merge(achar(ibclr(iachar(caray),5)),caray, (caray>='a' .and. caray<='z')) ! & ! upper case
                ! ibits(dim(iachar(caray),lcs)-1,0,8)<26)  !  obtuse bit manipulation mask
        end subroutine uppercase
      end subroutine setcase



! 3456789012345678901234567890123456789012345678901234567890123456789012 4567890
      character (40) function READSTRING(UNIT,buffr,lch,ststrt,stend,      &
                                         quotd, ccase)
!
!=======================================================================
!
! Description
!   Parses character data separated by white spaces or hard delimiters.
!   see comments below.  Null record results in zero parsed field
!
! Arguments
!
! - Input
!   buffr  character*(*), input string
!   UNIT    integer,       IO unit used to read input lines
!   lch     integer,       is the index of the next character to parse
!                          updated by the routine as parsing occurs
!   quotd   logical,       true if parenthesis or single quotes
!                          dissable delimeter parsing
!   ccase   integer,       -1 for lower case, 0 no change, 1 upper case
!
! - Output
!   ststrt  integer,       index of string start
!   stend   integer,       index of string end
!
! History
!   08/12/03 : Kendrick Killian
!              corrected a bug caused lch to point to trailing white space
!   11/ 4/99 : Kendrick Killian
!              added a quote character to disable delimiter searching
!   03/31/95 : Kendrick Killian
!              Forced a soft terminator at the buffer end.  Corrects the
!              bug which causes right justified strings to be dropped
!   01/25/95 : Kendrick Killian
!              Modified routine to locate and return a single field
!    3/15/91 : Kendrick Killian
!              Written
!
! Error Conditions
!   None
!
! External References
!      subroutine READCLIN   reads an input line
!
! Additional comments
!   1) a data field is established/bounded by a pair of upstream and
!      downstream delimiters.
!   2) commas are considered as hard delimiters:
!      - hard delimiters delimits field with no exception.
!      - consecutive hard delimiters constitute a null field.
!   3) White spaces (blanks and tabs) are soft delimeters:
!      - soft delimiters that preceed or follow an established data
!        field are ignored.
!      - consecutive soft delimeters do not constitute a null field.
!   4) a record is assumed to be bounded by soft delimiters.
!   5) quotd .true. dissables token parsing inside parenthesis. This
!      is useful for parsing subscripts or similar structures
!
!     Copyright 1992 - 2004    Kendrick Killian     All Rights Reserved
!========================== BEGIN EXECUTABLE CODE ======================
!
      character    buffr*(*), curquot*1
      integer      UNIT, lch, ststrt, stend
      logical      quotd
      integer      qdepth, pdepth, ccase

      integer I2, I, LENLIN

! Set Initial Values
      READSTRING = " "
      curquot    = " "
      I = lch - 1
      qdepth = 0
      pdepth = 0

! Set status (I2) flag to a hard delimiter condition
!    I2 condition codes are:
!        -1  Field terminated with a white space
!         0  Last Field terminated with a hard delimiter
!        >0  tracking through a token; location of last hard delimeter

      I2    = 0
      stend = 0

! determine record specified length
      LENLIN = LEN(buffr)

! read a record if required

! write(*,*) "READSTRING: ",4,trim(buffr),lch,ststrt,stend,quotd, ccase
20    I= I+1
      if (I.GT.LENLIN .OR. I.EQ.0) then
        if (I2.GT.0) then
          ststrt = I2
          stend = LENLIN
          GOTO 90
        else if (UNIT.LE.0) then
          lch = -1
          return
        else
          ! write(*,*) "READCLIN:",UNIT
          call READCLIN(UNIT,buffr,lch,.true.)
!          write(*,*) "read:",lch,buffr
          if (lch.lt.0) then
            ! call FILRERR(UNIT,lch,"read an End-Of-File ",buffr)
            ststrt = 1
            stend = 4
            buffr = '?--?'
            READSTRING = buffr
!           write(*,*) "EOF RETURN: ",ststrt,stend,READSTRING
            return
          endif
          I = lch
        endif
      endif


! Check to see if we have quoted the string.
        if (quotd .and. ISQUOTE(buffr(I:I))) then
          qdepth= qdepth+1
           if (I2 .LE. 0) I2 = I
           if(buffr(I:I) .EQ. curquot) then
              curquot = ""
              qdepth = 0
           else
              curquot = buffr(I:I)
              qdepth = 1
           endif

! Check to see if we have a parentheses which is a special quote character
        else if (quotd .and. (buffr(I:I).EQ.'(' .or. buffr(I:I).EQ.')')) then
          if (I2 .LE. 0) I2 = I
          if (buffr(I:I).eq.'(') then
            pdepth= 1+pdepth
          else if (pdepth.le.1) then
            pdepth= 0
          else if (pdepth.gt.0) then
            pdepth= pdepth-1
          endif

! if we have quoted the string then get the next character
        else if ( qdepth + pdepth .gt. 0) then
!          if (I2 .LE. 0) I2 = I
          GOTO 20

! Check to see if we have encountered a hard delimeter
        else if (ISHARD(buffr(I:I))) then
!
! Check if the previous field was terminated with a hard terminator
!  If so this is a null field
          if (I2 .EQ. 0 ) then
            ststrt = lch
            stend  = lch

! Check to see if we are stepping over a data field
! If so we need to terminate the field
          else if (I2 .GT. 0 ) then
            ststrt = I2
            stend = I-1
          endif

!  Found a field
          GOTO 90

! Have we encountered a white character with an intervening field
        else if ((I2 .GT. 0) .AND. ISWHITE(buffr(I:I))) then
          ststrt = I2
          stend = I-1

!         skip over trailing white space  (2 ifs avoid string bound overflow)
30        if(I .GE. LENLIN) goto 90
                 if(.not.ISWHITE(buffr(I+1:I+1))) goto 90
          I = I+1
          GOTO 30


! If the character is not white and we are starting a new field then
! we need to store its location
! (NOTE: don't worry about hard terminators since they won't get here

        else if ((I2 .LE. 0) .AND. .NOT.ISWHITE(buffr(I:I))) then
          I2 = I
        endif

! If we are here we are either stepping over white spaces after
! terminating a field or we are stepping over a valid field
!  In either case    KEEP GOING

      GOTO 20
!
! Clean up and end processing
!
   90 CONTINUE
      lch = I+1

! Remove quote characters from the returned string
      if(ISQUOTE(buffr(ststrt:ststrt)) .and. stend-1.gt.ststrt) then
        ststrt=ststrt+1
        if(ISQUOTE(buffr(stend:stend))) stend = stend-1
      endif

      READSTRING = buffr(ststrt:stend)
      call setcase(READSTRING, ccase)  !  change the case the output token

!     write(*,*) 'READSTRING RETURN: ',READSTRING
      return
        contains

!  Logical statement functions to defines delimeters
        logical function ISWHITE(CHR);  character CHR;
           ISWHITE = (CHR.EQ.' ' .OR. CHR.EQ.'	');
        end function
        logical function ISHARD(CHR);   character CHR;
           ISHARD  = (CHR.EQ.',');
        end function
        logical function ISQUOTE(CHR);  character CHR;
           ISQUOTE = (CHR.EQ."'" .or. CHR.EQ.'"');
        end function

      end function READSTRING


!2345678901234567890123456789012345678901234567890123456789012345678901234567890
      double precision function ReadReal(UNIT, buffr, lch, ISTAT)
!
!=======================================================================
!
! This is a general purpose routine for parsing a numeric character
! string and converting it to a number.
!
! Arguments
! - Input
!   UNIT,  INTEGER,       unit to read string from if lch = 0 or
!                         lch > len(buffr), a new buffr will be read and the
!                         conversion will start with the new line.
!                         UNIT <= 0 no read
!   buffr, character*(*), character string to be converted.
!   lch,   INTEGER,       is the index of the first character to parse
!                         this is updated as parsing occurs
!
!                         NOTE: lch > len(buffr) or lch = 0 Read a new buffer.
!                               1 <= lch <= len(buffr) the routine will return
!                               the conversion even if the buffer is empty
!
!                               lch < 0 and Cfp or Cip
!                               will return the previous value providing access
!                               to the other output form. This allows output of
!                               the real value in an incomplete integer read
!                               or the full precision of an integer.
!                               lch will NOT be modified!
!                               ISTAT will be set to 1
!
! - Output
!   lch,   INTEGER,       the index of the next token
!   ISTAT, integer,       conversion status
!                         -1, Conversion error
!                          0, no numerical conversion (alpha field)
!                          1, valid integer (and real)
!                          2, valid floating point and integer approximation
!                          3, valid floating point and integer overflow
!
! External References
!      subroutine READCLIN   reads an input line
!
! History
!   Modified :  6/10/15 : Kendrick Killian
!             made Cex, Cdec, ND default integers. They don't need to be
!             big but it silences compile warnings.
!   Modified : 10/19/11 : Kendrick Killian
!             made the repeat request depend on lch instead of ISTAT.
!             ISTAT is OUT only and will be cleared on any call.
!   Modified : 08/19/03 : Kendrick Killian
!             changed the behavior when lch = 0
!             It now reads the new line with skipblnk set to true.
!             This is in contrast to lch>len(buffr) which will return
!             a blank line error with comment or null lines.
!   Modified : 05/03/00 : Kendrick Killian
!             changed the behavior when UNIT > 0 to return the READCLIN
!             status codes through LCH
!                 LCH =  0, Empty line
!                 LCH = -1, End of File
!   Modified : 02/12/95 : Kendrick Killian
!             Substantial recoding to reduce the possibility of misreporting
!             unusual strings (strings starting with .,-,+,e)
!   Modified : 01/25/95 : Kendrick Killian
!             Modified arguments
!   Modified : 01/23/95 : Kendrick Killian
!             Removed FORTRAN77 non-compliant internal unformatted READ
!   Written  : 03/27/91 : Kendrick Killian
!
!     Copyright 1992 - 2014    Kendrick Killian     All Rights Reserved
!========================== BEGIN EXECUTABLE CODE ======================

      character, intent(INOUT) ::  buffr*(*)
      integer,   intent(IN)    :: UNIT
      integer,   intent(INOUT) :: lch
      integer,   intent(OUT)   :: ISTAT
      integer  ReadInt

!  local variables
      logical :: NOSPAC, RDEC, RNUM, NEG, NEGE, skipblnk
      integer :: Cex, Cdec, ND
      integer :: I, J, LOCE, LENLIN
      integer, save          :: Cip =0
      DOUBLE PRECISION, save :: Cfp =0.

!  Overflow/underflow value for a 32 bit integer
      integer ovrflw, Ndig
      PARAMETER (ovrflw = 2147483647, Ndig = 9)

!  Statement functions
      character CHR*1
      logical  ISDIGIT,ISWHITE,ISHARD
!     logical  ISUPC,ISLOWC

!     ISUPC(CHR) = CHR.GE.'A' .AND. CHR.LE.'Z'
!     ISLOWC(CHR) = CHR.GE.'a' .AND. CHR.LE.'z'
      ISDIGIT(CHR) = CHR.GE.'0' .AND. CHR.LE.'9'
      ISWHITE(CHR)  = CHR.EQ.' ' .OR. CHR.EQ.'	'
      ISHARD(CHR)  = CHR.EQ.',' .OR. CHR.EQ.';'

      data skipblnk/.FALSE./

      if (lch .lt. 0  .and.  Cfp .ne. 0) then
        ReadReal  = Cfp
        Cfp = 0.
        ISTAT = 1
        return
      endif
      RDEC = .TRUE.
      GOTO 10

      ENTRY ReadInt(UNIT, buffr, lch, ISTAT)
      if (lch .lt. 0  .and.  Cip .ne. 0) then
        ReadInt  = Cip
        ISTAT = 1
        Cip = 0
        return
      endif
      RDEC = .FALSE. ! return type

!     set initial values value
10    ISTAT = 0
      ReadReal  = 0.
      ReadInt = 0
      LENLIN = LEN(buffr)
!     If we are parsing a buffer, UNIT=0, and request a read lch=0
!      then set the lch to character one and avoid a trivial error return
      if (UNIT .EQ. 0 .and. lch .EQ. 0) lch= 1
      skipblnk = .FALSE.
      if(lch .EQ. 0) then
        skipblnk = .TRUE.
        lch = LENLIN +1
      end if

      I     = lch-1
      J     = 0
      ND    = 0
      Cex   = 0
      Cip   = 0
      Cfp   = 0.
      Cdec  = 0
      LOCE  = 0
      NEG    = .FALSE.
      NEGE   = .FALSE.
      RNUM   = .FALSE.
      NOSPAC = .TRUE.

! check for numeric entries.  Increment I until a non numeric character

20    I = I+1

!     Test for end of line
      if (I .GT. LENLIN) then

!       parsed something so jump out of the loop
        if (J.NE.0) GOTO 30

!       signal an error in case we don't read
        lch= -1
        if(UNIT.GT.0) call READCLIN(UNIT,buffr,lch,skipblnk)
        if (lch.LE.0) then
           ISTAT = -1
           return
        endif
        I = lch-1
        GOTO 20
      else
!       write (*,'(a,i3,3a,l2,f18.0,i10,3i4,3L3,2i4)')
!     1 ' character i=',I ,' "',buffr(I:I),'"',RNUM,
!     1    Cfp, Cip, ND, Cdec, Cex, NEG, NEGE, NOSPAC, J, LOCE
!   Accept a digit
        if (ISDIGIT(buffr(I:I)) .AND. NOSPAC) then
          J = ICHAR(buffr(I:I)) - ICHAR('0')
          if (LOCE.GT.0) then
!   record an exponent digit
            Cex = Cex * 10 + J
          elseif (ND .GE. Ndig) then
!   Suppress excess precision  (record decimal point shift on integer)
            if (.not.RNUM) Cdec = Cdec+1
          elseif (ND+J .EQ. 0)  then
!   Suppress leading zero's (record decimal shift if is a decimal
            ISTAT = 1
            if (RNUM) Cdec = Cdec-1
          else
!   record a significant digit
            ISTAT = 1
            Cip = Cip * 10 + J
            ND = ND + 1
            if (RNUM) Cdec = Cdec-1
          endif
          J=I
          NOSPAC = .TRUE.
!  Accept the first decimal point
        else if ((buffr(I:I).EQ.'.') .AND. .not.RNUM .AND. NOSPAC) then
          RNUM = .TRUE.
          if (J.EQ.0) J = I
!  Accept a plus sign
        else if ((buffr(I:I).EQ.'+') .AND. (J.LE.LOCE)) then
          if (J.EQ.0) J = I
!  Accept a minus sign
        else if ((buffr(I:I).EQ.'-') .AND. (J.LE.LOCE)) then
          if (J.EQ.0) J = I
          if (LOCE.EQ.0) then
            NEG  = .TRUE.
          else
            NEGE = .TRUE.
          endif
!   Step over spaces
        else if (ISWHITE(buffr(I:I))) then
          NOSPAC = ((LOCE.ge.J).or.(J.eq.0))
!   Check for an E or D in an exponent field
        else if (LOCE.EQ.0 .AND. ISTAT.NE.0 .AND.                        &
          index("EeDd",buffr(I:I)).NE.0) then
          RNUM = .TRUE.
          NOSPAC = .TRUE.
          LOCE = I
        else
!        write(*,*) ' decode ',buffr(lch:I),lch,I,J,ND,LOCE,NOSPAC,ISTAT
!        Be sure and step over a hard delimiter
          if (ISHARD(buffr(I:I))) I = I + 1
          GOTO 30
        endif
      endif
! step to the next character
      GOTO 20

! was anything accepted
30        lch = I
          if (ISTAT .EQ. 0) then
!           Found something (.+-) but they were not part of a number
            if(J.NE.0) lch = J
!           return and report a character conversion
            return
          endif

!  update the buffer pointer when we terminated on an E
          if (LOCE .GT. J) lch = LOCE


!   if you have reached here start the conversion to a real number
!         start by combining the mantissa and the exponents
!      write (*,*) ' parsed',Cfp,Cip,Cdec,Cex,ND,NEG,NEGE,ISTAT
      if (NEGE) Cex = -Cex
      Cex = Cex+Cdec

      if (NEG) Cip = -Cip
!     nasty this DBLE(Cip) tends to generate an inexact if the IEEE error
!     trapping is enabled
      Cfp =  DBLE(Cip) * 10.d0**Cex

      if (Cfp.gt.ovrflw) then
        Cip = ovrflw
        ISTAT = 3
      else if (RNUM) then
        ISTAT = 2
        Cip = DBLE(Cfp)
      endif

      if (RDEC) then
        ReadReal = Cfp
!        write (*,*) ' --- Real Return --',ReadReal, lch, ISTAT
      else
        ReadInt  = Cip
!        write (*,*) ' --- Integer Return --',ReadInt, lch, ISTAT
      endif

      return
      end


      subroutine READCLIN(UNIT,buffr,lch,skipblnk)
!
!=======================================================================
!
! Description
!   Reads a character buffer from the specified unit.
!   The routine skips comment lines and sets the case of alpha characters
!
! - Input
!   UNIT,   integer,       IO unit to read input lines from
!   skipblnk logical,      F forces the routine to return on a blank line
!
! - Output
!   buffr, character*(*), input buffer
!   lch      integer,       is the index of the first character to parse
!                           lch = 0   blank line
!                           lch = -1  End of file
!
! History
!   Modifications :
!     02/13/13 : Kendrick Killian
!              corrected a potential infinite loop in parsing quoted strings
!     02/12/03 : Kendrick Killian
!              added the ability to read from standard input
!    9/10/2002 : K. Killian
!              suppressed any nonprintable characters at the end of the buffer
!     06/20/00 : Kendrick Killian
!              removed an obscure bug resulting from comments starting in
!              the overflow variable, pad.
!     05/03/00 : Kendrick Killian
!              allow the routine to report an empty line with lch = 0
!              coding changes to use len_trim
!     10/08/98 : Kendrick Killian
!              added a choice of output case and return on blank line option
!     09/30/98 : Kendrick Killian
!              added a blank line and lch = 0 return on end of file
!              stops on a severe read error
!     03/14/95 : Kendrick Killian
!              added a 5 character over read to check for input overflow
!    Rewritten : 1/25/95 : Kendrick Killian
!               Recoded to remove special field processing
!
! Error Conditions
!   None
!
! External References
!   NONE
!
! Additional comments
!   1) input commands are converted to upper case
!   2) UNIX like comments can be inserted in the input stream
!      - # is the comment character
!      - comments extend to the end of the line
!
! Modifications
!   7/23/2018 by K. Killian
!    Left justify buffr; matches '*' input (this probably should be an input option)
!   9/10/2002 by K. Killian
!    suppressed any nonprintable characters at the end of the read buffer
!!=======================================================================
!
      character (len=*) :: buffr
      integer lch, UNIT
      logical skipblnk

!     local variables
      character PAD*(4), curqote*1
      integer i, lenp
      integer, save :: lenb = 0
      logical READY, trimcmnt
      character (len=40) :: ioerrmsg

      character(*), parameter :: prtyp= "abcdefghijklmnopqrstuvwxyz`0123456789-=[]\;',./" // &
                                        'ABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*()_+{}|:"<>?'

        if(UNIT.LE.0) then
           lch = -1
           return
        end if

        INQUIRE(UNIT=UNIT,OPENED=READY)
        if(READY) GOTO 10
           call FILRERR(UNIT,-1,"Unopened file",buffr)
           return

10      buffr = ' '
        PAD   = ' '
        READ(UNIT,'(a,a)',iostat=lch, IOMSG=ioerrmsg) buffr,PAD
        buffr = ADJUSTL(buffr) ! match the behavior of '*' input
!        write(*,'(i5,5a)') lch,"  '",buffr,"'   '",PAD,"'"

        ! Abort the run on a severe file error: Section 9.4.1.4i states
        ! Execution of an input/output statement containing the IOSTAT= specifier
        ! causes the variable specified in the IOSTAT= specifier to become defined
        !   (1)   With a zero value if neither an error condition, an end-of-file
        !         condition, nor an end-of-record condition occurs,
        !   (2)   With a processor-dependent positive integer value if an error
        !         condition occurs,
        !   (3)   With a processor-dependent negative integer value if an
        !         end-of-file condition occurs and no error condition occurs, or
        !   (4)   With a processor-dependent negative integer value different
        !         from the end-of-file value if an end-of-record condition occurs
        !         and no error condition or end-of-file condition
        if(lch .GT. 0) call FILRERR(UNIT,lch, "fatal read error "//trim(ioerrmsg), buffr)


!        reject any non printable characters at the end of the buffer.
         ! do with a scan of keyboard characters avoid ASCII dependence
         ! SCAN pad for the last position with a character from SET
         ! hopefully this to eliminates any DOS/UNIX problems.
        lenp = SCAN(PAD, prtyp, back=.true.)

        lenb = SCAN(buffr, prtyp, back=.true.)
        buffr(lenb+1:) = " "

!       IF we have a truly blank line
        if(lenb+lenp .EQ. 0) then

!         RETURN blank buffr on an "end of file" mark
!           ONLY do this for blank lines (no values in buffr/PAD) because
!           Certain misformed files can generate an EOF on the last valid line
!           with g77. Luckily the next read gives an EOF with a blank line.
          if (lch .LT. 0) return

!         are we skipping blank lines?
          if(skipblnk) GOTO 10

!         then return it!
          return
        endif

!  Check for comments and buffer overflow
        call uncomment (buffr,curqote,lenb,trimcmnt)
        if(curqote .ne. ' ') &
           call FILRERR(UNIT,0,"Input truncated, unterminated quotation ",buffr(:lenb))
        if(trimcmnt) lenp = 0

!       remove comment data from the pad
!        write(*,*) 'uncomment pad:"',PAD,'"'
        if(lenp>0) call uncomment (PAD,curqote,lenp,trimcmnt)

!       generate an error if there is a printable character in the pad.
        if (lenp>0) then
          do i = 1,lenp
            if (index(prtyp, PAD(i:i)) .ne. 0) then
!             This is a fatal error, so put most of the error message in buffr
              call FILRERR(UNIT,0,"Input truncated at '"//PAD//"'", buffr//PAD)
            endif
          end do
        endif
!        write(*,*) "uncommented:",buffr(:lenb)

!  check for a blank line
        lch = 0
 50     lch = lch +1
!  NOTE: the second string is a TAB  This does not comply with FORTRAN 77
!        standards.  If this is a problem: ASCII tab is CHAR(9)
        if(lch .GT. lenb) then
          if(skipblnk) GOTO 10
          lch = 0
          return
        end if
        if(buffr(lch:lch).EQ.' ' .OR. buffr(lch:lch).EQ.'  ') GOTO 50
      return
      end


      subroutine uncomment(buffr,curqote,lenb,trimcmnt)
      character (*)  ::  buffr
      character (1)  ::  curqote
      integer   lenb
      logical   trimcmnt

!     This is a routine to remove characters after a comment character.
!     Comment characters  escaped or quoted in single or double quotes are
!     treated as part of the string.
!
!     returns a .TRUE. value and the MODIFIED string if characters were trimmed

      integer    j, sq, cq
      character (len=1), parameter  ::  Comment = '#'


      ! initialize the return values character
      curqote = ' '
      trimcmnt  = .FALSE.
      lenb = len_trim(buffr)

      sq = 0
20    sq = sq+1
      ! write(*,*) "comments? ",sq,"  ",trim(buffr(sq:))
      j = INDEX(buffr(sq:),Comment)      ! Check for comment character
        if(j .eq. 0) return ! return if we didn't find one
        j = j +sq-1 ! this is its true location
        ! write(*,*) "found comment? ",j,"  ",trim(buffr(j:))

        ! is the comment character escaped
        !   NOTE: use the ASCII code for '\', 92, since the character itself
        !         tends to drive compilers NUTS
        if(j .eq. 1) then
           ! Do Nothing   The first character can't be quoted
        elseif(buffr(j-1:j-1) .eq. char(92)) then
          ! If so first remove the escape character and restart search
          buffr(j-1:) = buffr(j:)
          sq = j
          goto 20

        else
          ! Search to see if the comment character is quoted
          do while (sq .lt. j)

            ! write(*,*) j, sq, (buffr(sq:sq).EQ."'" .or. buffr(sq:sq).EQ.'"'), buffr(sq:)
            if (buffr(sq:sq).EQ."'" .or. buffr(sq:sq).EQ.'"') then
              curqote = buffr(sq:sq)

              ! search for the close quote character
              cq = INDEX(buffr(sq+1:), curqote)

              ! if there is not another quote then return curqote uncleared
              !  this may be an error unless you are doing multi buffer strings
              if (cq.eq.0) return

              ! close the quote
              curqote = ' '
              sq = cq+ sq

!             if the comment character was quoted, search for another comment
              if (sq .gt. j) goto 20
              ! the quote was closed BEFORE the comment. Is there another?
            endif
            sq = sq+1
          end do
        endif

      ! Remove commented characters
       ! write(*,*) "remove comments ",trim(buffr(j:))
      buffr(j:) = ' '
      trimcmnt = .TRUE.
      lenb = j
      end


      subroutine c2f(Cstr, result)
! This is a routine to convert a C string to a FORTRAN character variable
!     logical*1 Cstr(1)  !     byte Cstr(1)
      integer, parameter :: byte  =  selected_int_kind(1)
      integer (byte)     :: Cstr(1)
      character (*)      :: result
      integer n

       result = ' '
       do 10 n = 1, LEN(result)
         if (Cstr(n) .eq. 0) return
10       result(n:n) = CHAR(Cstr(n))
      end


      subroutine FILRERR(UNIT,lch,MSGE,buffr)
      integer UNIT, lch
      character (*) :: MSGE, buffr
      character (len=256) :: mssg

!****   parameter
      integer, parameter :: defunit = 128

        if(lch .EQ. -1) then
          write (mssg,*) 'I/O failure ',MSGE," unit ",UNIT
        else
          if(UNIT.eq.defunit) then
            mssg= trim(MSGE)//' reading: standard I/O'
          else
            inquire(UNIT=UNIT,NAME=mssg)
            mssg= trim(MSGE)//' reading: '//trim(mssg)
          endif

          if (lch .eq. 0) then
!             write(*,*) mssg
             call message (mssg)
             mssg= 'Line: '//trim(buffr)
          endif
        endif

!        write(*,*) mssg ; stop
        call abortrun(mssg)
      end
