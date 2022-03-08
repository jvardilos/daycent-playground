      program sitarchive
!
! compiled with: gfortran  -g -fno-underscoring -fcaller-saves   ../../list100/filext.f  messagefor.f90  retarg.f90   parse.f90  chkdata.f90  Obj/initlyrs.o Obj/watreqn.o Obj/initsw.o Obj/mssg.o Obj/initsite_tg.o  Obj/initsrad.o Obj/swpotentl.o  writsit.f90 sitein.f90  sitarchive.f -o DDCsitarchive
!    gfortran -g -fno-underscoring -fcaller-saves ../list100/filext.f message.f90 retarg.f90 parse.f90 chkdata.f90 Obj/initlyrs.o Obj/watreqn.o Obj/initsw.o Obj/mssg.o Obj/initsite_tg.o Obj/initsrad.o Obj/swpotentl.o writsit.f90 sitein.f90 sitarchive.f -o ../bin/DDCsitarchive
!
! This routine reads and writes records to a century binary site archive.
!
! 9 Mar 2013  K. Killian
!    corrected a bug that triggered a block copy instead of moving a single
!    record from one archive to a larger record number in a second archive.
!    copy/scan now log completions
! 27 July K. Killian
!    Added a copy option to move blocks of records to a different archive.
!    recoded the argument processing so that it just parses the arguments
!    The difference is that the file names are now more order dependant
!    file testing is done in the main routine
! 27 July K. Killian
!    added a record check scan mode
!  Modified   Jan 2008    K. Killian
!    recoded the argument processing to correct a problem with out of order file names
!  Modified   6/13/2000    K. Killian
!    Now using the machine independent argument routine.
!    Included binary file open/close in binrw so it can be used for testing.
!  Modified   2/23/2000    K. Killian
!    expanded to include the enhanced site parameters and a character field.
!    To distinguish the new archive format the extension has been changed to
!    .esa (enhanced site archive). New archives will be created in the enhanced
!    format. The routine will still read and write the old 240 line format.
!  Modified   8/17/1999    K. Killian
!    corrected a bug that allowed it to read a nonexistant file.
!  Modified   1/15/1999    K. Killian
!    Added a routine to completely generate the site file syntax.
!  Coded      1/14/1997    K. Killian
!    this routine combines the original routines dearch.f archive.f
!    It provides a common interface for dealing with the site archive.
!    This is similar to the interface written for the binary archive.
!
!  NOTE: has NO index maximum checking
!

!   nsitrec # site records
!   vrsn    1) normal site    2) enhanced (default)
        integer vrsn, mxrec, nsitrec(2)
        parameter (mxrec=290)

        integer nargs, i, lch, lwrt, lsit, sit, ivopt, ostat, flgvrb
        integer :: rindx = 0, windx = 0
        logical :: filexist, sitparread, soilread
        character (len=132)  srcnam, wrtnam, buffr
        character (len=80)   label

        integer, parameter :: MXSWLYR = 21

!      Functions
!       real ReadReal
!
        character*1 squote
!        parameter (squote= char(39))
        parameter (squote= '''')

!****       file unit numbers  sit is must be compatable with sitein
        parameter (sit=7)

        data nsitrec, label/240,mxrec, '-'/


           ! Set some variables default
           call defsitpar()

!****   Determine command line arguments
        call retarg(nargs,'arg count')
        if (nargs .lt. 2) then
          write(*,'(a)')' Site archive utility',                         &
     &      ' usage:  sitarchive [-r/w] <Record#> <archive.esa/.sarc> '  &
     &      //' <site.100>',                                             &
     &      '         move site date between ascii and binary '//        &
     &      'files as determined by file name.',                         &
     &      '         Binary records must have the '//                   &
     &      'action and record number specified',                        &
     &      '      -r/w #  read or write archive record'
          stop
        else
!****   get the command line and open the specified files
          call getcmdline(rindx,srcnam,lsit,                       &
     &                      windx,wrtnam,lwrt,vrsn, flgvrb)
        endif

        ! if(wrtnam .eq. 'STDOUT') call message('set quiet /dev/null')

!****   check to see if the site file exists
        filexist = .FALSE.
        if(srcnam.ne.'NONE') inquire(file=srcnam,exist=filexist)
        if(.not. filexist) call abortrun('Unknown site file: '//srcnam(:lsit))

!****  Execute a record scan
        if(rindx .lt. 0   .and.  windx .lt. rindx) then
          call blockrecord(srcnam, wrtnam, abs(rindx),abs(windx), label)
          stop;
        endif



!****  Read the data from the archive
        if (rindx .ne. 0) then
           call binsitrw(srcnam(:lsit), rindx, 0, label, 17) !read binary record
        else
          write(*,'("Reading         ",a)') srcnam(:lsit)

!****     Open the <site.100> file
          open(unit=sit, file=srcnam, status='OLD', iostat=ostat)
          if(ostat .gt. 0) call abortrun(' opening Site file '//trim(srcnam))
          read(sit,'(a)',end=90) buffr
          label = srcnam
          if(label(1:4) .eq. 'new_') then
             label = label(5:)
             i=lsit-2
          else
             i=lsit+2
          endif
          lch = index(buffr,'data from ')
          if(lch.ne.0) label(i:) = 'C'//buffr(lch+4:)
!
!**** read the first site 100 records
          rewind(sit)
          ivopt=0
          call sitein(ivopt,flgvrb, sitparread, soilread)
        endif


!****  Write the data
        filexist = .FALSE.
        if (windx .eq. 0) then
          if(wrtnam.ne.'STDOUT') then
            inquire(file=wrtnam,exist=filexist)
            if (filexist .eqv. .true.) then
              ! verify the overwrite
              write(*,'("Replacing       ",a)') wrtnam(:lwrt)
            else
              write(*,'("Creating        ",a)') wrtnam(:lwrt)
            endif
          endif

          call sitprint(trim(wrtnam),'"'//trim(label)//'"  '//trim(srcnam),  &
     &       rindx) ! write site.100 file
        else
           call binsitrw(wrtnam(:lwrt), 0, windx, label, 17) !output binary record
        endif

        stop

!****     error message for site file read
90        write(*,*) 'Fatal error reading site file ',srcnam(:lsit)
!****   Close files
        close(unit=sit)
      end


!****   SUBROUTINE GETCMDLINE
      subroutine getcmdline(rindx,srcnam,lsit,                     &
     &                      windx,wrtnam,lwrt,vrsn, flgvrb)
        integer rindx, lsit, windx, lwrt, vrsn, flgvrb
        character*(*) srcnam, wrtnam

!       Local variables
        integer istat, lch, larg, indx
        character (len=1) :: rwd
        character*128 argv
        integer ReadInt

        nrd   = 0
        lsit  = 6
        srcnam = 'NONE'
        wrtnam = 'STDOUT'
        rindx = 0
        windx = 0
        indx  = -1
        rwd   = ' '
        vrsn  = 0
        argv = ' '
        flgvrb = 0;

        call retarg(0, argv)

        do while (argv .ne. '?--?')
          larg = len_trim(argv)
           ! write(*,*) "  argv: ",argv(:larg);

          ! test the argument for a bare number
          lch = 0
          istat = 0
          indx = ReadInt(0,argv,lch,istat)
          if (istat.eq.3  .or.  indx.lt.0) then
             call abortrun('bad record index '//trim(argv))
          endif

          ! are we looking for a number?
          if(rwd .ne. ' '  .and.  istat .gt. 0  .and.  indx .ge. 0) then
            select case (rwd)
              case('r')
                rindx = indx
              case('w')
                windx = indx
              case('c', 's')
                rindx = -indx
                windx = -indx
                if(istat.eq.2 .and. argv(lch:lch) .eq. '.') then
                  lch = lch +1
                  windx = -ReadInt(0,argv,lch,istat)
                endif
                if(rwd .eq. 's') then
                  wrtnam = "ScAnIn"
                  lwrt = 6
                endif

                if(abs(rindx) .gt. abs(windx) .or. windx.eq.0 .or. rindx.eq.0) then
                  call abortrun("improper index range "//trim(argv))
                endif
              case default
                call abortrun("unknown option "//argv(:2))
            end select
            rwd = ' '
          else if(argv(:8) .eq. "--dbgsit") then
            flgvrb = 1

          else if(argv(:1) .eq. "-") then
            rwd = argv(2:2)
            ! is the record number here?
            if(larg .gt. 2) then
              argv = argv(3:)
              cycle
            else
              argv = " "
            endif

          else  ! These should be file names
            ! write(*,*) 'file name ',rwd,istat, trim(argv)

            if(srcnam .eq. 'NONE') then
              srcnam = argv
              lsit = larg
            else
              wrtnam = argv
              lwrt = larg
            endif

          endif

          call retarg(0, argv)
        end do

        ! A assume file names
        istat = 0
        call filext(srcnam, '.100.sarc.esa.dsa',0,larg,istat)
        nrd = 0
        if(wrtnam .ne. 'STDOUT'  .and.  wrtnam .ne. 'ScAnIn') call   &
              filext(wrtnam, '.100.sarc.esa.dsa',0,larg,nrd)

        ! check to see if the order is swapped
        if((istat .eq. 1  .and.  rindx .ge. 1  .and.  nrd .gt. 1) .or.    &
           (nrd   .eq. 1  .and.  windx .ge. 1  .and.  istat .gt. 1)) then
             argv = srcnam            ! swap names
             srcnam = wrtnam
             wrtnam = argv
             indx = nrd               ! swap types
             nrd = istat
             istat = indx
             indx = lsit              ! swap length
             lsit = lwrt
             lwrt = indx
        endif

        ! should have an archive if there is a read index
        if(rindx.ge.1 .and. istat.le.1) then
             call abortrun("missing index for "//trim(srcnam))
        endif

        ! other than SCAN there should be an archive if with a write index
        !  SCAN is different since windx is the end limit
        if(wrtnam .ne. 'ScAnIn'  .and.  windx .gt. 0  .and.  nrd .le. 1) then
          call abortrun("missing index for "//wrtnam(:larg))
        endif

        if(srcnam .eq. 'NONE') call abortrun("No site file to read")

        return
      end



!****   SUBROUTINE blockrecord
      subroutine blockrecord(srcnam, wrtnam, rindx, windx, label)
        character*(*) srcnam, wrtnam, label
        integer rindx, windx
        do i = rindx, windx
           call binsitrw(trim(srcnam), i, 0, label, -17) !read binary record
           if (wrtnam .eq. 'ScAnIn') cycle
           call binsitrw(trim(wrtnam), 0, i, label, -17) !read binary record
        end do

        if (wrtnam .eq. 'ScAnIn') then
          write(*,*) "Scanned records ", rindx, " - ",windx
        else
          write(*,*) "copied records ", rindx, " - ",windx, " from ",    &
          trim(srcnam)," to ",trim(wrtnam)
        end if

      end
