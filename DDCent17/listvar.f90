 module listnames

   implicit none

     include 'table.inc'

   contains
 ! methods to access table names;  they can/will change with the storage method
     INTEGER FUNCTION search(word,tbl)
       INTEGER,intent(out) :: tbl
       CHARACTER*(*),intent(in) :: word
       character (len=12) :: erindx

       ! locate the variable name
       ! This search requires a space at the end of every variable name
       ! search = index(' '//table//' ',' '//trim(word)//' ')
       ! search = search/tblen +1       ! raw variable number
       search = 1;
       do while (search <= MAXVAR)
         if(table(search) .eq. word) exit
         search = search+1 !search the string array
       end do
! write(*,*) 'search:',word,search
       if(search > MAXVAR) then
         search = -1;  tbl = -1;
       else
         !convert the variable number to a table-variable pair
         tbl = 1
         do while (tbl < MAXINC  .and.  search > tvals(tbl))
           search = search - tvals(tbl)
           tbl = tbl+1
         end do
 ! write(*,*) "search word:",word,"  location ",index(table,trim(word))/tblen +1, " (", &
 ! table((search-1)*tblen+1:search*tblen),")",tbl,search
         if(search > tvals(tbl)) then ! error check the results
           write(erindx, "(2i5)") tbl,search;
           call abortrun("incorrect index "//erindx//"generated in search for:"//word)
         endif
       endif

       RETURN
     END FUNCTION search

     character*16 FUNCTION tblvnam(tbl,var)
       integer, intent(in) :: tbl,var
       integer :: tbloc
       character (len=12) :: erln

       if(tbl > MAXINC  .or.  var>tvals(tbl)) then
         write(erln, "(2i5)") tbl,var;
         call abortrun ("list variable location "//erln//" outside bounds");
       endif

       ! tbloc = (sum(tvals(1:tbl-1))+var)*tblen; tblvnam = table(tbloc-tblen+1:tbloc);
       tbloc = (sum(tvals(1:tbl-1))+var); tblvnam = table(tbloc);
! write(*,*) "tblvnam word:",tbl,var, sum(tvals(1:tbl-1))+var, tbloc, tblvnam

     END FUNCTION tblvnam
  end module listnames

      subroutine listvar(lis,time,biname,ascnam)
      use listnames
      implicit none

      INTEGER   lis
      REAL      time
      CHARACTER biname*(*), ascnam*(*)
!
!  This routine opens the lis file and writes the records
!     lis     The unit number
!             lis  = 0 input specifies the input command mode
!     time    the date for the output record;
!             time returns a negative value on I/O error
!     biname  the binary/evt file name used to label the lis file;
!             In command mode biname contains the command

!              inputime       direct a number to start/end times in order
!              startim        set stime if not set
!              notim          set etime if not set
!              cleartim       clear stime and etime

!
!     ascnam  the name of the output file to open
!             if the file name is ' ' or '-', no file is opened;
!     Modified:
!        7/2007  K. Killian
!           corrected a formatting bug that printed negative numbers as 0
!        4/2006  K. Killian
!           increased output list and corrected full buffer warning
!        4/2003  K. Killian
!           minor changes to ensure compatability with g77
!        KLK  4/16/2002
!           Modified the output format to provide scientific output for
!           overflow and underflow values

      logical opn
      INTEGER i
      REAL    outvar      !      output variable
      CHARACTER outfrmt*12
      character*10 :: crpchar

!****   Object common block
         INTEGER, PARAMETER :: maxvr=84
         REAL,    PARAMETER :: notim = 0.00137
         INTEGER  numchsn, found(2,maxvr), rntim, crpv(2)
         REAL      stime, etime
         CHARACTER  del*1
         COMMON /varlist/ stime,etime,rntim,numchsn,found,crpv,del
      CHARACTER  lnout(maxvr)*13, ctim*8
      CHARACTER (len=5), parameter::  chfrmt='(85a)'

      include 'outval.inc'

!**** Command mode
      if(lis .eq. 0) then
         if(biname(1:9).eq.'delimiter') then
           if(len(biname) .ge. 10) then
             del = biname(10:10)
             if(del.eq.'t' .or. del.eq.'T') then
!               NOTE:   the string is a TAB
                del = CHAR(9)
                call message('setting the delimiter to tab')
             else if(del.eq.'c') then
                del = ','
                call message('comma separated values')
             else if(del .eq. ' ') then
                call message('space delimited lis file')
             else if(index('0123456789+-.EDed*',del) .ne. 0) then
                call message("ignoring numeric delimiter '"//del//"'")
                del = ' '
                time = -1
             else
                call message('lis file delimiter is "'//del//'"')
             endif
           endif

         else if(biname.eq.'inputime') then
!         This is a number so assign it to a date
            if (stime .eq. notim) then
              stime = time
            elseif (etime .eq. notim) then
              etime = time
            else
              time = -9
            endif

         else if(biname.eq.'startim') then
           if (stime .eq. notim  .or.  stime .eq. -9999) stime = time

         else if(biname.eq.'endtim') then
           if (etime .eq. notim  .or.  etime .eq. -9999) etime = time

         else if(biname.eq.'cleartim') then
!           call listval
           stime = notim
           etime = notim
           rntim = 0
!           time = 0

         else
           time = -2
           call message (' unknown listvar command')
         end if
         return
      end if

!**** return if there are no output variables specified
      if(numchsn .le. 0) return

      inquire(unit=lis,OPENED=opn)
      if(.not. opn) then

!****   Check for a valid file name. If none is given mark the text file as
!        inactive; this keeps us from opening or writing to an unnamed file
        if(ascnam.eq.' ' .or. ascnam.eq.'-') then
           call message(' no ascii file ')
           time = -1.
           numchsn = -1
           return
        endif

!****	Open ascii file
        open(unit=lis,file=ascnam, status='UNKNOWN', iostat=i)
        if (i.gt.0) call abortrun('opening output file '//trim(ascnam))
        write(lis,'(a)') trim(biname)

!****	Create the variable header record
        do 40 i = 1, numchsn
          lnout(i) = tblvnam(found(1,i), found(2,i))
          if(del .eq. ' ') lnout(i)= adjustr(lnout(i))  ! right justify the column name
40      continue

!****   write the header record
        if(del .eq. ' ') then
          write(lis,chfrmt) '    time',(del//trim(lnout(i)),i=1,numchsn)
        else
          write(lis,chfrmt) 'time',(del//trim(lnout(i)),i=1,numchsn)
        endif
        write(lis,'("")')
        return
      end if

        if (rntim .eq. 0) then
!****   Convert input run times back to code ull and
!       Add very small value to etime b/c of real numbers in binary
           if (stime .eq. -9999) stime = notim
           if (etime .eq. -9999) then
              etime = notim
           else if (etime .ne. notim) then
              etime = etime * (1.000001)
           endif
           rntim = 1
        endif

!       return if we are outside the print window
!       print if(stime=notim or time>stime) .and. (etime=notim or time<etime)
        if ((stime.ne.notim .and. time.lt.stime) .or. &
            (etime.ne.notim .and. time.gt.etime)) return

        write(ctim,'(f8.2)') time
        if(del .ne. ' ') ctim = ADJUSTL(ctim)
!       Modified the output format KLK  4/16/2002
        do 80 i = 1, numchsn
          if (found(1,i) .eq. 1) then
            outvar = vals1(found(2,i))
          elseif (found(1,i) .eq. 2) then
            outvar = vals2(found(2,i))
          elseif (found(1,i) .eq. 3) then
            outvar = vals3(found(2,i))
          elseif (found(1,i) .eq. 4) then
            outvar = vals4(found(2,i))
          elseif (found(1,i) .eq. 5) then
            outvar = vals5(found(2,i))
          endif

          lnout(i) = ' '
          if(found(1,i) .eq. crpv(1) .and. found(2,i) .eq. crpv(2)) then
            lnout(i) = crpchar(outvar)
            lnout(i) = ADJUSTR(lnout(i))
          else
!           zero the output for very small output values.
            if(abs(outvar) .lt. 1.0e-9) then
               lnout(i)(len(lnout(i))-4:) = '0.0  '
            else
!              standard format
               outfrmt = '(4x,   f9.3)'

!              scientific output for values that won't display in f9.3
               if(abs(outvar).ge.100000. .or. abs(outvar).lt.0.1) then
                 outfrmt = '(2x,1PE11.4)'
!                reduce precision for negative values (the delta 13 values)
                 if(outvar .lt. 0.) outfrmt(11:) = '3)'
               endif

               write(lnout(i),outfrmt) outvar
            endif
          endif
80      continue

!
!****   write a data record
        if(del .ne. ' ') then
          write(lis,chfrmt) trim(ctim),(del//trim(adjustl(lnout(i))),i=1,numchsn)
        else
          write(lis,chfrmt) trim(ctim),(del//lnout(i),i=1,numchsn)
        endif
       RETURN
      END subroutine listvar


!**** SUBROUTINE GETLIST
      SUBROUTINE getlist(unitnum,bffr,varcnt,errpt,istat)
        use listnames
        CHARACTER bffr*(*)
        INTEGER unitnum, varcnt, errpt, istat
!
!  Description
!    locates variable names and adds them to the ouput tables
!
!  Arguments
!    Input
!    unitnum  INTEGER,       IO unit to read input lines from  0 for internal
!    buffr    CHARACTER*(*)  input string,  unit 0
!    varcnt   INTEGER        number of variables found
!    errpt    INTEGER        Error reporting conditions
!                         -1)  no reporting,   0) list,  1) verbose,  2) Fatal
!    istat    INTEGER        unknown word, error, count
!                          0)  normal operation,  >0 number of unrecognized words
!
!       Modifications
!   4/2005  K. Killian
!...A new output variable reports bad variables, errors, vs duplicate entries to
!   allow better list processing
!

!****   Object common block
         INTEGER, PARAMETER :: maxvr=84
         REAL,    PARAMETER :: notim = 0.00137
         INTEGER  numchsn, found(2,maxvr), rntim, crpv(2)
         REAL      stime, etime
         CHARACTER  del*1
         COMMON /varlist/ stime,etime,rntim,numchsn,found,crpv,del
         save   /varlist/

!****   Local variables
        INTEGER nv, nt, i, lch, ststrt, stend, fnum
        CHARACTER word*15, READSTRING*40
        logical, save :: wrn=.TRUE.

        varcnt= numchsn
        fnum= 0
        nv= 0
        word= ' '
        istat = 0

!****   Obtain list of variables from input file
        if (unitnum .eq. 0) then
          lch = 1
        else
          lch = len(bffr)+1
        endif
10      word = READSTRING(unitnum,bffr,lch,ststrt,stend,.true.,-1)

        if (lch .gt. 0) then

!         keep track of how many tokens we have parsed
          fnum= fnum+1

!         find the word in the output variable list
          nv= search(word,nt)
          if(word .eq. 'crpval') then
            crpv(1) = nt
            crpv(2) = nv
          endif

! Check for a non variable error for this word
          if (nv.eq.-1) then

!           increment the error count
            istat = istat+1

!            if (errpt.gt.0 .or. (errpt.eq.0 .and. fnum.ge.2))               &
            if (errpt.gt.0 .or.                                              &
               (errpt.eq.0 .and. (lch.lt.len(bffr) .or. fnum.ge.2)))         &
               call message('ERROR: Variable "'//word(:stend-ststrt+1)       &
                            //'" not found')
            if (errpt .ge. 2) stop 'Abnormal Termination'
            goto 10
          else

! make sure that the new variable is a not already on the list
            do i = 1, numchsn
              if (nt.eq.found(1,i) .and. nv.eq.found(2,i)) then
                 if(errpt.gt.0)                                              &
                  call message('Warning: duplicate output '//word)
                 goto 10
              end if
            end do

! add the variable to the list
            if (numchsn .lt. maxvr) then
              numchsn = numchsn+1
              found(1,numchsn) = nt
              found(2,numchsn) = nv
              goto 10
            endif

            if (wrn) then
              call message('Warning: Too many output variables. last'//      &
                   ' accepted: '//trim(tblvnam(found(1,maxvr), found(2,maxvr))))
              wrn = .FALSE.
            endif

            goto 10
          endif

        endif
        varcnt = numchsn-varcnt
        RETURN
      END SUBROUTINE getlist



!****   INTEGER FUNCTION flag754
      integer function flag754()

!****   A routine to search for IEEE exception values.
!       the presence of these values indicates that the simulation went bad.
!       This is generally the result of bad input values. The routine is
!       in listvar since it is similar to search.
!       Returns the count of IEEE 754 exceptions

!       IEEE 754 specifies (NaN .ne. NaN) is true
!       IEEE 754 specifies a/Inf = 0 and a+Inf = Inf

        use listnames
        include 'outval.inc'

!****   Local variables
        integer i
        character mssg*24 ! increase length for conversion to C string



        flag754 = 0
        mssg = "IEEE754 EXCEPTIONS:"

          do i = 1, tvals(1)
            if (test754(vals1(i))) then
              flag754 = flag754+1
              if(flag754 .eq. 1) call message(mssg)
              write(mssg,'(2a,f4.0)') tblvnam(1,i),"=",vals1(i)
              call message(mssg)
            endif
          end do

          do i = 1, tvals(2)
            if (test754(vals2(i))) then
              flag754 = flag754+1
              if(flag754 .eq. 1) call message(mssg)
              write(mssg,'(2a,f4.0)') tblvnam(2,i),"=",vals2(i)
              call message(mssg)
            endif
          end do

          do i = 1, tvals(3)
            if (test754(vals3(i))) then
              flag754 = flag754+1
              if(flag754 .eq. 1) call message(mssg)
              write(mssg,'(2a,f4.0)') tblvnam(3,i),"=",vals3(i)
              call message(mssg)
            endif
          end do

        return

      CONTAINS
        logical function test754(v)
        real    v
!        test754(v) = (ieee_is_nan(v) .or. IEEE_IS_FINITE(v)) .eq. 0.0)
!        test754(v) = (v .ne. v .or. 1.0e37/(1.+abs(v)) .eq. 0.0)
          test754 = (ISNAN(v)  .or.  abs(v) .gt. huge(v))
        end function test754
      end function flag754


      function crophash(curcrp)
      integer (kind=4) :: crophash
      character (len=*) :: curcrp

!...Determine the 'numerical hash' of the curcrp, for use as an output
!     The new hash code is the following
!     It consists of 2 fields.
      ! Field 1 (21 bits) codes up to the first 4 alphanumeric
      ! (0-9 A-Z case insensitive) characters from the name.
      ! Other characters are skipped.
      !
      ! Field 2 (8 bits) makes it possible to distinguish between names that
      ! are longer than 4 characters by summing of all the alphanumeric
      ! characters in name. The sum is formed using the following rules
      !    numbers    0-9      1 + number
      !    letters    A-Z     11 + sequence number of letter (case ignoring)
      !    other characters   skipped
      ! The maximum value should be 28 characters * 37 = 1036 which slightly
      ! exceeds the 10 bits allocated for storage. To keep

      ! storage
      ! the two fields are stored in an IEEE 32 bit floating point word
      !   S EEEEEEEE FFFFFFFFFFFFFFFFFFFFFFF
      !   0 1      8 9                     31
      ! field 2 is stored in the upper 10 bits S - 10
      ! field 1 is stored in the mantissa field 11-31
      ! to create a well formed floating point value, field 2 starts with
      ! a bias of 10 to ensure that there are bits in the exponent field.
      ! To ensure that the sum does not exceed the bits allocated, once
      ! sum of the stored characters exceeds 31, the redundant sum of the
      ! stored characters is subtracted from the total sum.

      integer          :: ii, jj, cs, crpch

      crophash = 10
      cs       = 0
      crpch  = 0
      do ii = 1, len(curcrp)
        if (curcrp(ii:ii) .ne. ' ') then
          jj = iachar(curcrp(ii:ii))
          if (curcrp(ii:ii).ge.'0' .and. curcrp(ii:ii).le.'9') then
            jj = 1 + jj - iachar('0')
          elseif(curcrp(ii:ii).ge.'A' .and. curcrp(ii:ii).le.'Z') then
            jj = 11+ jj - iachar('A')
          elseif(curcrp(ii:ii).ge.'a' .and. curcrp(ii:ii).le.'z') then
            jj = 11+ jj - iachar('a')
          else
            cycle
          endif

          if(ii .le. 4) then
             crpch = crpch*37+jj ! keep the first 4 characters
             cs = cs + jj
          endif

          crophash = crophash + jj ! keep the sum as a hash function
      !write(*,'(a, i3,2i10,Z20)') curcrp(ii:ii),jj,crophash,crpch,crpch

          !if (curcrp(ii:ii) .ge. '0' .and. curcrp(ii:ii) .le. '9') then
          !  crpval = crpval + ((ichar(curcrp(ii:ii)) - ichar('0'))/10.)
          !else
          !  crpval = crpval + (ichar(curcrp(ii:ii)) - ichar('A')) + 1
          !endif
        endif
      end do
      if(cs.gt.31) crophash = crophash - cs
      crophash=ior(ishft(crophash,21), crpch)
      !write(*,'(a,i20,g30.9)') curcrp, crophash,transfer(crophash,0.0)
      return
      end


      function crpchar(crpval)
      character*10 :: crpchar
!      integer (kind=4) :: crpval ! type cast the normal 32 bit float
      real (kind=4) :: crpval ! type cast the normal 32 bit float

      ! convert a crop hash value to a nice printable string
      ! quote string to aid importing into database, excel ...
      ! does not use write so you can write crpchar without rentrance problems
      integer          :: ii, jj, crpch, cs
      integer  (kind=4),  save :: lastcv = 0 ! retain last hash
      character (len=10), save :: lastch     ! last translation
      integer, parameter :: ach0 = iachar('0')

      crpchar = ' '
      ! does this match the last value?
      if(crpval .eq. lastcv) then
        crpchar = lastch  ! return the stored value
      else

        crpch = ibits(transfer(crpval,0),0,21) ! extract the char pattern

        crpchar = "'" ! lead quote
        cs = 0   ! sum of stored characters
        ii = 5   ! decode from right to left
        do while (crpch .gt. 0) ! loop over nonzero values
          ! divide with remainder
          jj = crpch
          crpch = crpch / 37
          jj = jj-crpch*37
          cs = cs + jj  ! the retained digit sum
          ! convert integer io character
          if    (jj .eq. 0) then   ! this shouldn't happen
             crpchar(ii:ii) = '*'
          elseif(jj .le. 10) then  ! digit
             crpchar(ii:ii) = achar(jj -1 + ach0)
          else                     ! letter
             crpchar(ii:ii) = achar(jj -11 +iachar('A'))
          endif
          ii = ii - 1
        end do

        ! process the character sum
!        jj = ishft(transfer(crpval,0),-21) ! shift the sum into an integer
        jj = ibits(transfer(crpval,0),21,10) ! extract the sum into an integer
        jj = jj -10            ! remove bias
        if(cs .gt. 31) jj = jj + cs ! add retained digit sum
        ! convert number to string (Don't think write is reentrant)
        ii = 10
        do while(jj>0)
          cs = jj
          jj = jj /10
          crpchar(ii:ii) = achar(cs -jj*10 +ach0)
          ii = ii-1
        end do
        crpchar(6:) = ADJUSTL(crpchar(6:)//"'") ! shift digits next to name
        crpchar(2:) = ADJUSTL(crpchar(2:))      ! shift digits next to name
!      write(*,'(2a,Z10,i12,g20.8)')'crpchar:',crpchar,crpval, &
!       crpval,transfer(crpval,0.0)
        lastcv = crpval
        lastch = crpchar
      endif

      return
      end function crpchar


      BLOCK DATA listval

!****   Object common block
!     found    INTEGER(2,*)   output tables
!     numchsn  INTEGER        number of entries in the output tables
         INTEGER, PARAMETER :: maxvr=84
         REAL,    PARAMETER :: notim = 0.00137
         INTEGER  numchsn, found(2,maxvr), rntim, crpv(2)
         REAL      stime, etime
         CHARACTER del*1
         COMMON /varlist/ stime,etime,rntim,numchsn,found,crpv,del

         data del, rntim, numchsn, found/' ', 0, 0, maxvr*0, maxvr*0/
         data stime, etime/notim, notim/
      END

!elemental logical function ieee_is_nan(x)
!real,intent(in):: x
!ieee_is_nan = isnan(x)
!end function
!elemental logical function ieee_is_finite(x)
!real,intent(in):: x
!ieee_is_finite = .not. (isnan(x) .or. abs(x) > huge(x))
!end function


