   program build
    implicit none

!****   Build include files for use with the list100.f program and a user's model.
!       Inputs:
!            file name;    up to MAXINC include files to mine for output variables
!            -I directory; adds a directory to the include search path
!       Outputs:
!            outval.inc
!                    used by the model to write output variables to
!                    binary file and by list100 to retrieve those
!                    variables from the binary file
!            table.inc
!                    used by listvar, parameters describing the variable label array
!            table.f
!                    used by listvar, Block Data containing the variable labels
!   Modified 07Jul2018 K. Killian
!    Added an include file search feature similar to that used by the compiler/linker
!    a series of -I directory sets a search path for the input include files
!    table.inc and outval.inc will be opened in the first (-I) directory or the default directory
!    MAJOR CHANGE: Build OVERWRITES existing output files instead of aborting.
!                  Files checked out of version control loose modification dates.
!                  Shouldn't interrupt a clean build to backup files in version control.
!   Modified 31Jul2011 K. Killian
!    changed to handle longer lines of F90 input
!    maximum variable length upgraded from 6 to 10 characters for Daycent
!    some simplification of the continue code
!   Modified 30Oct2010 K. Killian
!    converted to f90
!    used f90 string functions to build table.f output lines
!    convert integers to strings without using WRITE
!    eliminated partial line output to table.f
!    modified parser to accept dimensioning in integer/real statements
!   Modified 7/8/2000 K. Killian
!    added the ability to read arguments
!    added the save line to outval.inc
!    program ignores comment lines in the include files
! commands useful for debugging
!  mv outval.INC outval.inc ; mv table.F table.f ; mv table.INC table.inc
!  mv outval.inc outval.INC ; mv table.f table.F ; mv table.inc table.INC

    integer, parameter:: MAXINC = 5, MAXV = 500, maxpth = 8
    integer, parameter :: tblen=16

    integer, dimension(MAXINC,MAXV,3) ::   indx
    character (len=16) :: label(MAXINC,MAXV)
    character (len=16) :: vartyp(MAXINC)


    !* Local variables
    integer   ostat
    integer count, fnum, i, j, incmax, tvals(MAXINC), lcmnd, maxvar
    character*255 fname(MAXINC), tname
    character (len=256), dimension(maxpth) :: incpath
    character (len=256) :: inlin   ! input line to parse
    character (len=256) :: outval = 'outval.inc'  ! output file name
    character (len=256) :: table  = 'table.inc'   ! output file name
    character (len=25)  :: cmnd; ! executable command
    ! character     itochar*16
    logical error, goahead
    integer :: npath  =0
    integer :: incnum =0

    !* Initialization
    fnum = 11
    tvals     = 1
    fname     = 'dummy'
    inlin = '';
    incpath = '';


 ! parse the command line arguments
call GET_COMMAND_ARGUMENT(0,cmnd,lcmnd,i) ! get cmnd the executable name
call commandline()                        ! parse the arguments

call findinc(incpath,npath,fname,incnum)  ! do the search for the include files

    !* Check for existance of old files
    error = .false.

    call outfnam (outval)

!****   Open the outval.inc include file
    open(unit=fnum, file=outval, status='unknown', iostat=ostat)
     if(ostat .gt. 0) call abortrun('Output file '//trim(outval)//' could not be opened')

    write(fnum,'(2("!",15x,a/)//)') 'Copyright 1993 Colorado State University', &
                                    '          All Rights Reserved'

    vartyp = 'real'
    count  = 1
    incmax = 1
    !* Determine the number of variables and the labels in each .inc
    do i = 1, incnum
      call parser(i,count,fname(i),MAXINC,MAXV,indx, label, vartyp(i))
      call putcount(fnum,fname(i),i)
      tvals(i) = count
      incmax   = max(incmax, count)
    enddo
    do i = incnum+1, MAXINC
       call putcount(fnum,fname(i),i)
    end do
    write(fnum,'(/)')


    !* Write out the data type of all include files to outval.inc
    !write(fnum,*) '     real vals1, vals2, vals3, vals4, vals5'
    if(incnum .lt. MAXINC) vartyp(incnum+1:MAXINC) = 'real'
    do i = 1, MAXINC
      call putype(vartyp(i),fnum,i,tvals(i))
    enddo

!***add a save statement to outval.inc
23  count = 12
    tname = '      save'
    do i = 1, incnum
       j = index(fname(i), '.')-1
       if(j<0) j = len_trim(fname(i))
       if(count+j+4 .gt. 72) then
          write(fnum,'(a)') tname
          goto 23
       endif
       tname(count:) = ' /'//fname(i)(1:j)//'/,'
       count = count+j+4
    end do
    tname(count-1:) = ' '
    write(fnum,'(a)') trim(tname)

!***Close the new include file
    close(unit=fnum)

!***Open table.f90
    call outfnam (table)
    count = 0;
    maxvar = sum(tvals(1:incnum));
    open(unit=fnum,file=table,status='unknown', iostat=ostat)
     if(ostat .gt. 0) call abortrun('Output file '//trim(table)//' could not be opened')


!***Write common and type declarations for table.inc
    write(fnum,'(5x,"integer, parameter :: MAXINC =",i3)') incnum
    write(fnum,'(5x,"integer, parameter :: tblen =",i3)') tblen
    write(fnum,'(5x,"integer, parameter :: tvals(MAXINC) = (/",i4)',  &
          ADVANCE="NO") tvals(1)
    do i=2,incnum
      write(fnum,'(",",i4)',ADVANCE="NO") tvals(i)
    end do
    write(fnum,'("/)"/)')
    write(fnum,'(5x,"integer, parameter :: MAXVAR =",i5)') maxvar

!***Write the labels for each include file
    write(*,*) "table.f ",incnum

!c  common/tables/tvals(  3), table(  3,281)
    write(fnum,'(5x,"character (len=",i2,"), dimension(MAXVAR), parameter :: table = (/ &")') tblen

    do i = 1, incnum
      write(fnum,'(a)') "!"
      call wrtable(i,fnum,tvals(i),MAXINC,MAXV,indx, count, maxvar, label)
    end do

!***Close table
    close(unit=fnum)

     contains

       ! SUBROUTINE PUTCOUNT
       subroutine outfnam(file)
         character (len=*) :: file
         error = .false.

         inquire(file=file,exist=goahead)
         if (.not. goahead  .and.  incpath(1) .ne. ' ') then
           file = trim(incpath(1))//'/'//file
           inquire(file=file,exist=goahead)
         endif
         if (goahead) then
           write(*,*) 'Overwriting '//trim(file)
         else
           write(*,*) 'Creating    '//trim(file)
         endif
       end subroutine outfnam

       ! SUBROUTINE PUTCOUNT
       subroutine putcount(fnum,label,num)
       integer fnum, num
       character*(*) :: label
       integer       :: cend

       cend = index(label, '.')
       if(cend .eq. 0) then
         cend = len_trim(label)
       else
         cend = cend -1
       endif

       write(fnum,"(6x,'common /',a,'/ vals',i1)") label(:cend), num

       return
       end subroutine putcount

       ! SUBROUTINE putype
       subroutine putype(vartyp,fnum,num,count)
       integer fnum, num, count
       character*(*) :: vartyp
       character*3   :: cntstr

       write(cntstr,'(i3)') count
       cntstr = adjustL(cntstr)

       write(fnum,"(6x,a,' :: vals',i1,'(',a,')')") trim(vartyp), num, trim(cntstr)
       return
       end subroutine putype

       ! SUBROUTINE PUTCOUNT
       subroutine findinc(incpath,npath,fname,incnum)
         character (len=*), dimension(*) :: incpath
         character (len=*) :: fname(MAXINC)
         integer           :: npath, incnum

         character (len=256) :: inlin   ! input line to parse
         integer           :: lend, lenf
         logical           :: goahead !,error

file:    do i = 1, incnum
           lenf = len_trim(fname(i))
           inquire(file=fname(i), exist=goahead)
           if(goahead) cycle file

path:      do j = 1, npath
             lend = len_trim(incpath(j))

             if(lend > 0  .and. lend < 255) then ! make sure we stay inside the defined string
               if(incpath(j)(lend:lend) .ne. '/') then
                 lend = lend +1
                 incpath(j)(lend:) = '/'
               endif
             endif

             if(lend+lenf .gt. 255) cycle ! skip the inquire if the file neme will truncate

             if(lend > 0) inlin = incpath(j)(:lend)//fname(i)
             inquire(file=inlin, exist=goahead)
             if(goahead) then
               if(lend > 0) fname(i) = inlin ! copy back if path string was defined
               cycle file;
             endif
           end do path
         end do file

       end subroutine findinc

       subroutine commandline ()
         integer     :: nargs
         character (len=8)   aswitch ! argument switch
         integer     :: windx, clen
         logical     :: dpars = .false.

         windx  = 0
         aswitch = '  '

        ! ... Get command line arguments
         call retarg(nargs,'arg count')
         if (nargs .eq. 0) call abortrun("input include files")

        ! ... Process command line arguments
        ! ... Add a new argument to indicate reading an extended <site>.100 file
        ! ... which was created from a site.nc file from a Gridded Century run
        ! ... cak - 10/05/01
         do
           if(inlin .eq. ' ') call retarg(0, inlin)
           if(dpars)  call message('retarg: '//trim(inlin)//'  aswitch: '//trim(aswitch)) ! echo the input if requested

           ! SPECIAL CASES
           ! look for debug argument command
           if(inlin .eq. 'dpars'  .or.  inlin .eq. '--dpars') then
             call message("debugging parameter parsing")
             write(*,*)"debugging parameter parsing"
             dpars = .true.              ! set the debug logical
             inlin = ''; cycle;         ! clear the input line and reread

           else if(inlin .eq. cmnd(:lcmnd)) then ! skip the command (make dependancy)
             inlin = " ";    ! clear the input line
             cycle;          ! reread
           else if(inlin(:2) .eq. '-c') then ! Immediately read input from command file
             if(len_trim(inlin) .gt. 2) then
               inlin = inlin(3:)           ! Appended command file name
             else
               call retarg(0, inlin)       ! load the command file name
             endif
             inquire(file=inlin,exist=goahead) ! is there a file by that name?
             ! abort if there is nothing to open
             if (.not. goahead) call abortrun('command file not found: '//trim(inlin))
             call retarg(-21, inlin)     ! have the read package open the command file
             inlin = " ";    ! clear the input line
             cycle;          ! reread
           end if

           ! determine if this the next switch
           if(aswitch .eq. ' ') then
             if(inlin .eq. '?--?') then
               if(dpars)  call message("EOF")
               exit
             elseif(inlin(:2) .eq. '--') then
               aswitch = inlin
               inlin  = ' '
             elseif(inlin(:1) .eq. '-') then
               aswitch = inlin(:2)
               inlin   = inlin(3:)
             endif
           endif

!           ! These switches don't require an argument
!           select case (aswitch)
!             case ('-h', '-?')
!               call manual()
!
!           end select

           if(inlin .eq. ' ') cycle ! remaining switches need to look at the arguments

           ! ---- BARE WORD,  No Pending flag; treat as an include file
           if(aswitch .eq. ' ') then
             if(incnum<MAXINC) then
               aswitch = '-f'
             else
               call message('Skipping unknown bare word argument: '//trim(inlin))
               inlin = ''
               cycle
             endif
           endif

           ! At this point everything should be defined
           clen = len_trim(inlin)
           if(dpars)  call message('** case aswitch: '//aswitch//'  inlin: '//inlin(:clen))

           select case (aswitch)

             case ('-f')
               if(incnum<MAXINC) then
                  incnum = incnum+1
               else
                  call abortrun("too many include files")
               endif
               fname(incnum) = inlin(:clen)

             case ('-I') ! add to search path
               if(npath<maxpth) then
                  npath = npath+1
               else
                  call abortrun("too many directories in search path")
               endif
               incpath(npath) = inlin(:clen)

!             case ('-g')
!               ext_grid = .true.
!               schnam = inlin(1:clen)
!               if (index(schnam,'.sch').eq.0) schnam(clen+1:clen+4) = '.sch'
!               inquire(file=schnam,exist=goahead)
!               if (.not. goahead) call abortrun('schedule file not found; '//trim(schnam))

             case default
               call message('Skipping Unknown argument: '//aswitch)

           end select
           aswitch = ' '
           if(inlin(:1) .ne. '-') inlin   = ' '
         end do

         !... perform any close operations
         call retarg(-9,'')

       end subroutine commandline

   end program build



!****   SUBROUTINE PARSER
        subroutine parser(incnum,count,filename,MAXINC,MAXV,indx, label, vartyp)
        integer :: incnum, count
        character*(*) filename
        character (len=*) :: vartyp

        integer   MAXINC, MAXV
        ! integer, dimension(:,:,:) :: indx
        integer, dimension(MAXINC,MAXV,3) :: indx
        character (len=*) :: label(MAXINC,MAXV)

        !****   Local variables
        character (len=16) :: wrdrt
        integer   aptr, i, ii, j, k, istat, vcnt, lintyp, ptr
        character inrec*132, word*25
        logical   cont, cont90

!****   Function declarations
        integer ReadInt

!****   Open the output include file
        write(*,*) "reading:    ",trim(filename)
        open(unit=10,file=filename,status='OLD', iostat=istat)
        if(istat .gt. 0) call abortrun('Include file "'//filename//'" could not be opened')

!****   Initialize ptr into label and indx arrays
        aptr = 0

!****   Initialize the count of variables in this file
        count = 0
        lintyp = 0
        cont=.false.
        cont90=.false.

!****   Find line that begins common block
10      read(10,'(a)',end=60) inrec

        ! remove anything that is a f90 comment
        i=index(inrec, '!')
        if(i.gt.0) inrec(i:) = ''
        if (len_trim(inrec).eq.0) goto 10

        cont=.false.
        if(inrec(1:1) .eq. 'c' .or. inrec(1:1) .eq. 'C' .or.            &
           inrec(1:1) .eq. '1') goto 10

        i=1
        if    (cont90) then
          cont=.true.
          cont90 = .false.
          if(inrec(:6) .eq. '     &') inrec(:6) = ' ' ! redundant but accepted
        elseif(inrec(6:6) .ne. ' ' .and. inrec(:5) .eq. ' ') then
          cont=.true.
          inrec(6:6) = ' '
        end if

        inrec = adjustl(inrec) ! remove leading spaces

        ! if(inrec(len_trim(inrec):) .eq. '&') then
        k= index(inrec, '&')
        if(k .gt. 0) then ! f90 continue mark
          inrec(k:) = ' '
          cont90=.true.
        else              ! continue is one line at a time
          cont90 = .false.
        endif

!       evaluate the command if we are not in continue
        if(.not. cont) then

          if((inrec(:6) .eq. 'COMMON' .or. inrec(:6) .eq. 'common')) then
!****       On first line, move past 'common/label' characters
            i = index(inrec, '/',BACK=.true.)
            lintyp = 1
            ii = index(inrec,'/')
            filename = inrec(ii+1:i-1);  ! write(*,*) "common block ", filename
            i = i +1

          else if(inrec(:4) .eq. 'SAVE' .or. inrec(:4) .eq. 'save') then
            ! skip a save statement
            goto 10

          else
            i = max(index(inrec,'DIMENSION'), index(inrec,'dimension'),     &
                    index(inrec,'REAL'), index(inrec,'real'),               &
                    index(inrec,'INTEGER'), index(inrec,'integer'),         &
                    index(inrec,'LOGICAL'), index(inrec,'logical'),         &
                    index(inrec,'DOUBLE PRECISION'), index(inrec,'double precision'))

            if(i.ne.0) then
              ii = index(inrec(i:), ' ') +i-1
              if ((inrec(i:i) .eq. 'D' .or. inrec(i:i) .eq. 'd') .and.  &
                   (inrec(i+1:i+1) .eq. 'O' .or. inrec(i+1:i+1) .eq. 'o')) then
                do while (inrec(ii:ii) .eq. ' ')
                  ii = ii+1
                end do
                ii = index(inrec(ii:), ' ') +ii-1
              endif
              vartyp = inrec(i:ii-1)
              i = ii
            endif
            lintyp = 2
          endif

          if(i.eq.0) then
            lintyp = 0
            goto 10
          endif
          inrec = adjustl(inrec(i:))
          i=1
        endif

!       remove leading spaces
        if(inrec(1:1) .eq. ',') inrec = adjustl(inrec(2:))
        i=1

        lenlin = len_trim(inrec)
        if (lenlin.eq.0) goto 10 ! skip a blank line

!****   Parse each line of the common block
! 40      if (len_trim(inrec).gt.0) then

!****   Get each word from the current line
50      if (i .le. lenlin) then
          call nextword(inrec,i,word)

          istat = len_trim(word)
          if (istat.gt.0) then
            wrdrt = word

            k=index(word,'(')
            if(k.ne.0) then
              wrdrt = wrdrt(:k-1)
            endif

            ptr = 1
            if(aptr.ne.0) then
              do while (label(incnum,ptr) .ne. wrdrt .and. ptr .le. aptr)
                ptr = ptr + 1
              end do
            endif

!****       Store the label
            if(ptr .gt. aptr) then
              label(incnum,ptr) = wrdrt
              indx(incnum,ptr,1) = 1
              indx(incnum,ptr,2) = 0
              indx(incnum,ptr,3) = 0
              count = count + 1
              aptr = ptr
! write(*,*) count, incnum,ptr,label(incnum,ptr) ,indx(incnum,ptr,:)
            endif

            if(k.ne.0) then
!****         Calculate the total number of values in variable
              k = k+1
              vcnt = 1
              do j=1,3
                istat = 0
                indx(incnum,ptr,j) = ReadInt(0,word,k,istat)
                if(indx(incnum,ptr,j) .gt. 0) vcnt =                    &
                    vcnt * indx(incnum,ptr,j)
              enddo
              count = count + vcnt -1
            endif

          endif
!****     Check if there could be more words in line
          goto 50
        endif
!****   Read the next line
        goto 10

60      continue
        close(unit=10)

        return
        end subroutine parser


!****   SUBROUTINE NEXTWORD
        subroutine nextword(inrec,pos,neword)
        character inrec*(*), neword*(*)
        integer pos

!****   Local variables
        integer newpos, lenlin
        logical inparen

!****   Initialization
        neword = ' '
        newpos = 1
        inparen = .false.
        lenlin = len_trim(inrec)

!****   Skip over any leading white space
10      if (pos .le. lenlin .and. inrec(pos:pos) .eq. ' ') then
           pos = pos + 1
           goto 10
        endif

!****   Build neword from line
20      if (pos .gt. lenlin) then
          return
        elseif (inrec(pos:pos) .eq. '(') then
          inparen = .true.
        elseif (inrec(pos:pos) .eq. ')') then
          inparen = .false.
        elseif ((inrec(pos:pos) .eq. ',' .and. .not. inparen) .or.        &
                (inrec(pos:pos) .eq. ' ')) then
          pos = pos + 1
          return
        endif
        neword(newpos:newpos) = inrec(pos:pos)
        newpos = newpos + 1
        pos = pos + 1
        goto 20

        end subroutine nextword



!**** SUBROUTINE WRTABLE
      subroutine wrtable(incnum,fnum,tval,MAXINC,MAXV,indx, tblcnt, maxvar, label)
        integer incnum, fnum, tval, tblcnt, maxvar
        integer   MAXINC, MAXV, indx(MAXINC,MAXV,3)
        character (len=*) label(MAXINC,MAXV)

       integer, parameter :: tblen=16

!****   Local variables
        integer count, i, j, k, m


!****   Write labels
        i = 0
        count =0
        do while (count .lt. tval)
          i = i+1
!****     Write out the data line
          if (indx(incnum,i,3) .ne. 0) then
            do m = 1, indx(incnum,i,3)
              do k = 1, indx(incnum,i,2)
                do j = 1, indx(incnum,i,1)
                  count = count+1
                  call wrtheadr(fnum,tblcnt,label(incnum,i),j,k,m)
               end do
              end do
            end do
          else if (indx(incnum,i,2) .ne. 0) then
            do k = 1, indx(incnum,i,2)
              do j = 1, indx(incnum,i,1)
                count = count+1
                  call wrtheadr(fnum,tblcnt,label(incnum,i),j,k)
              end do
            end do
          elseif (indx(incnum,i,1) .ne. 1) then
            do j = 1, indx(incnum,i,1)
              count = count+1
                  call wrtheadr(fnum,tblcnt,label(incnum,i),j)
            end do
          else
            count = count+1
            call wrtheadr(fnum,tblcnt,label(incnum,i))
          endif
        end do
        return
      CONTAINS

        subroutine wrtheadr(fnum,tblcnt,label, i1, i2, i3, i4)
        integer, optional :: i1, i2, i3, i4
        integer fnum, tblcnt
        character (len=*) :: label

!****   Local variables
        integer lhead, ndig
        character     indx*16
        character (len=tblen) :: vrnam

        vrnam = trim(label)
        lhead = len_trim(vrnam)

        if(present(i1)) then
          indx = itochar(i1,ndig)
          vrnam(lhead+1:) = "("// indx(:ndig)//")"
          lhead = lhead +2 +ndig
        endif
        if(present(i2)) then
          indx = itochar(i2,ndig)
          vrnam(lhead:) = ","// indx(:ndig)//")"
          lhead = lhead +1 +ndig
        endif
        if(present(i3)) then
          indx = itochar(i3,ndig)
          vrnam(lhead:) = ","// indx(:ndig)//")"
          lhead = lhead +1 +ndig
        endif
        if(present(i4)) then
          indx = itochar(i4,ndig)
          vrnam(lhead:) = ","// indx(:ndig)//")"
          lhead = lhead +1 +ndig
        endif

        tblcnt = tblcnt +1
        if(tblcnt .eq. maxvar) then
          write(fnum,'(5x,a,i5)') "'"//vrnam//"'/)     !", tblcnt
        else
          write(fnum,'(5x,a,i5)') "'"//vrnam//"',  &  !", tblcnt
        endif

        return
        end subroutine wrtheadr


        character *16 function itochar(intnum, ndig)
        integer n, intnum

        if(intnum .eq. 0) then
          itochar(1:1) = '0'
          ndig=1
        else
          ndig=11
          itochar = ''
          n = abs(intnum)

          do while (n .gt. 0)
            ndig = ndig-1
            itochar(ndig:ndig) = achar(IACHAR("0") + (n -int(n/10) *10))
            n=n/10
          end do

          if (intnum < 0) then
            ndig = ndig-1
            itochar(ndig:ndig) = '-'
          endif

          itochar = adjustl(itochar)
          ndig = 11 - ndig
        endif

        return
        end function itochar

      end subroutine wrtable
