!cat listvar.f90 table.f > lv.f ;gfortran list100.f90 ../DailyDayCentEVI/parse.f90 ../DailyDayCentEVI/retarg.f90 ../DailyDayCentEVI/message.f lv.f  table.f -o ~/bin/DDClist100; rm lv.f
!cat listvar.f   table.f > lv.f; gfortran -static list100.f90 -O3 -g -fno-underscoring -fcaller-saves -Wall 32bit/message.o 32bit/retarg.o 32bit/parse.o lv.f -o /data/paustian/Century/bin/DDClist100; rm lv.f

      program list100
!23456789012345678901234567890123456789012345678901234567890123456789012
!
! Modified 6/13/00 by K. Killian
!  corrected bug parsing input files
!  Improved help message and output error message construction.
!  changed goto to cycle.
! Modified 6/13/00 by K. Killian
!    Now using the machine independent argument routine. Also remembers
!    binary file name from lis file to use if no other binay file specified.
! Mar 7, 99   K. Killian   Command line Interface modifications
! Modified the command line interface to expand input options
!   NOTES:
!     1) The first two numeric inputs will be interpreted as a dates
!        This means all numeric file names will be a problem
!     2) The binary and lis file types are determined by the extensions.
!     3) Extensions are added to file names input without them.  The files are
!        tested against binary, list, and variable files in that order.  I am
!        assuming that the extensions will be there since the century system
!        adds them by default.
!     4) if an existing list file is specified it WILL be deleted and recreated
!     5) If no list file is specified then we will use the binary file name
!     6) The variable file has no default extension.  However it must exist.
!     7) A -r extension RECREATES an existing list file; reading the dates
!        and variables from the specified file.  Variables can be added to
!        the list and the dated can be changed by command line inputs
!     8) An existing list file can be used for variable list.  This happens
!        automaticaly when you overwrite a list file.
!
! 4/25/03   K. Killian   Minor modifications for g77
! 2/20/97   K. Killian   Interface modifications
! Modified to use an internal character buffer instead of using the
!     non standard fortran that suppresses record ends.
!     presently the buffer length is 256 bytes or 20 output variables
!     A check was added to terminate the input list at a full buffer
! Added changes to make the file names more tolerent
!     allows file extensions to be specified on .bin and ASCII output files
!     allows ASCII file extensions other than .lis
! Modifications to the interactive input mode include:
!     allows over writting of the list file
!     reasks for mistyped output variables
!
!   compiled with: f77 -u -g list100.f chartool.o table.o -o  ~/bin/list100
!   Dos compiled : fl /AL /Ox /FPc /c list100.for table.for
!   Dos compiled : fl list100.obj table.obj

      integer in, narg, numchsn, lis
      real :: time=0
      character*150 ascname, biname

      include 'outval.inc'
      include 'table.inc'

!****	Initialization
        numchsn = 0
        ascname = '-'
        biname  = '-'
!****	Assign file unit numbers
        in = 1
        lis = 2


!****   Determine command line arguments or interactive
      call retarg(narg,'arg count')
      if (narg .eq. 0) then
        call interactive(in,biname,ascname)
      else
        call getcmdline(narg,in,biname,ascname)
      endif

!****	Open ascii file and write the header records
          call listvar(lis,time,biname,ascname)

!****	make sure there is something to output
!++++++++++++++++++++++++++++++++
      if (time .lt. 0) then
        write(6,'(/a,a)') 'Fatal Error: No output variables;',           &
     &                    ' no output file created'
        STOP
      endif


!****	Read binary file, writing to ascii
50    continue
       SELECT CASE (MAXINC)
          CASE (1)  !if (MAXINC .eq. 1) then
            read(unit=in,end=90) time,vals1
          CASE (2)  !else if (MAXINC .eq. 2) then
            read(unit=in,end=90) time,vals1,vals2
          CASE (3)  !else if (MAXINC .eq. 3) then
            read(unit=in,end=90) time,vals1,vals2,vals3
          CASE (4)  !else if (MAXINC .eq. 4) then
            read(unit=in,end=90) time,vals1,vals2,vals3,vals4
          CASE (5)  !else if (MAXINC .eq. 5) then
            read(unit=in,end=90) time,vals1,vals2,vals3,vals4,vals5
        END SELECT  !endif

        call listvar(lis,time,' ',' ')
      goto 50

!****	Close files
90    close(unit=in)
      close(unit=lis)

      end


!****   SUBROUTINE INTERACTIVE
      subroutine interactive(in,biname,ascname)
        integer in
        character*(*) biname, ascname
        character     YN*1,   bffr*150
        integer clen, len_trim, lch

!****   Local variables
        integer itemp, istat
        logical existascii, binok

!****   Print title
        write(6,'(/11x,a/7x,a/15x,a/a)') 'CENTURY List100',              &
                        'Binary to Ascii Utility','12/99',               &
                        '   Enter name of binary input file:'

!****   Obtain name of binary input file
10      read(5,'(a)') bffr
        clen = len_trim(bffr)
        if(clen.eq.0) stop

!****   Check that binary file exists
        call chkbin(in,bffr,biname,binok)
        if (.not. binok) then
          write (6,*) ' Unknown file:'//biname
          goto 10
        endif

!****   Obtain name of ascii output file
20      write(6,'(/a)') '   Enter name of ASCII output file:'
        read(5,'(a)') ascname

!****   If ascii file already exists, ask for it again
        clen = len_trim(ascname)
        if(clen.eq.0) goto 20
! Check to see if an extension was specified  (.??? format)
! add the .lis extension in none was specified
        if(ascname(clen-3:clen-3).ne.'.') ascname(clen+1:clen+4)= '.lis'
        inquire(file=ascname,exist=existascii)
        if (existascii) then
          write(6,'(/3a/a)') '  ASCII file ',ascname(:clen),' exists.',  &
                            '   Overwrite  (y/n)'
          read(5,'(a)') YN
          if (YN.ne.'Y' .and. YN.ne.'y') goto 20
        endif

!****   Obtain time interval
        write(6,'(/a/a)') '   Enter starting time,',                     &
                         '   <return> for time file ends:'
        read(5,'(i5)') itemp
        if (itemp .gt. 0) call listvar(0,real(itemp),'startim',' ')

        write(6,'(a/a)') '   Enter ending time,',                        &
                         '   <return> for time file ends:'
        read(5,'(i5)') itemp
        if (itemp .gt. 0) call listvar(0,real(itemp),'endtim',' ')

!****   Obtain list of variables
        write(6,'(/a/a)') '   Enter variables','   <return> to quit:'
100     CALL READCLIN(5,bffr,lch,.false.)
          if(lch.le.0) return
!      getlist(unitnum,bffr,varcnt,errpt,istat) vs (unitnum,bffr,varcnt,errpt)
          call getlist(0,bffr,itemp,1,istat)
        goto 100
        end


!**** SUBROUTINE GETCMDLINE
      subroutine getcmdline(narg,in,biname,ascname)
! Turn this routine into an event loop to process the inputs
        integer narg, in
        character*(*) biname, ascname
        character (len=256) bffr

!****  Local variables
        character ARGV*80
        logical islist
        logical :: replis = .false., cmndf = .false., varead = .false.
        logical :: binok  = .false., nolstfl=.TRUE.
        logical :: echo  = .false.

        integer :: i = 1
        integer :: lch, lenA, istat, varcnt
        real rt
        double precision ReadReal
        real, parameter :: notim = 0.0078125

!****  functions
        integer :: isvar

        bffr = ""

        do while (i .le. narg)

!****  Get the argument
          call retarg(i, ARGV); !        write(*,*) 'ARGV:',trim(ARGV)
          i=i+1

          if(cmndf) then
            call retarg(-2, ARGV)
            if(ARGV(:4) .eq. "?--?") call abortrun ("opening command file:"//trim(ARGV(5:)))
            i=1
            call retarg(narg,'arg count')
            cmndf = .false.
            cycle
          endif

          !if(ARGV(:4) .eq. "?--?") exit ! leave if this this the last argument

          lenA = len_trim(ARGV)
          if(echo) write(*,*) ' argument ',i,narg,' length ',lenA,'  ',ARGV(:lenA)

          if(ARGV .eq. '?') call manual     !****  bare '?' help request
          ! parse flags
          if(ARGV(1:1) .eq. '-') then
            SELECT CASE (ARGV(:2))
              CASE ('-h','-?') ! help
                call manual
              CASE ('-r')      !****  Check for a redo list file flag
                replis = .true.
                cycle
              CASE ('-d')      !****  Is this a delimiter flag
                rt = 1
                call listvar(0,rt,'delimiter'//ARGV(3:3),' ')
                if(rt.lt.0) call message ('ignoring an unknown delimeter')
                cycle
              CASE ('-c')      !****  command file flag
!          -c is special since the arguments MAY be in the file
                cmndf = .TRUE.
                cycle
              CASE ('-v')      !****  accept/ignore the variable flag
                ARGV = ' '
                cycle
              CASE ('-e')      !****  debug echo
                ARGV = ' '
                echo = .true.
                cycle
              CASE DEFAULT
                write(*,*) 'skipping unknown switch '//ARGV
                cycle
            END SELECT  !endif
          endif

          ! write(6,*) ' argument ',i,' length ',lenA,'  ',ARGV(:lenA)
!****  if the dates are not set check to see if this is a number
          lch = 1
          istat = 0
          rt = ReadReal(0, ARGV, lch, istat)
          if(istat.gt.0 .and. lch.gt.lenA) then
!         This is a number so assign it to a date
            call listvar(0,rt,'inputime',' ')
            if (rt .eq. -9) then
               write (*,*) 'Fatal Error: Unknown numeric argument ',ARGV(:lenA)
               STOP
            endif
            cycle
          endif

!****  Check if this is a binary file  (should normally be the first argument)
          if(.not.binok) then
            call chkbin(in,ARGV,biname,binok)
            if(binok) cycle
          endif

!****  Check if this is a variable
          call getlist(0,ARGV,varcnt,0,istat)
          ! write(6,*) ' call getlist ',ARGV(:lenA), "   ",varcnt
          if (varcnt.gt.0) cycle

!****  Check if this is a list file
          if (nolstfl .and. islist(ARGV,ascname,biname,replis,bffr)) then
             nolstfl = .FALSE.
             cycle
          endif

!****  If we haven't read a variable file, is this is a valid one
          if (.not. varead) then
             if (isvar(ARGV) .gt. 0) then
                varead = .TRUE.
                cycle
             end if
          end if

!****  if you are here then I don't know what to do with this argument
          if (ARGV(lenA-3:lenA-3).eq.'.') then
            ARGV = 'Warning: Missing file     '//ARGV(:lenA)
          else
            ARGV = 'Warning: Unknown Argument '//ARGV(:lenA)
          endif
          write (6,*) ARGV(:lenA+26)
        end do

        ! echo the file names
        if(echo) then
          write(*,*) "biname  ",trim(biname)
          write(*,*) "ascname ",trim(ascname)
        endif

        ! Make sure there is a binary file specified
        if(biname .eq. '-') then
          write(*,*) 'Fatal Error:  No binary file specified'
          stop
        end if

        ! If we have a binary file and no ASCII file then use the binary root
        if(binok  .and.  ascname .eq. '-') then
            ascname = biname(1+index(biname,'/',.TRUE.):)
            ascname(len_trim(ascname)-2:) = 'lis'
        endif

        if(replis) then
           lch = 1
           rt = ReadReal(0, bffr, lch, istat)
           call listvar(0,rt,'startim',' ')
           rt = ReadReal(0, bffr, lch, istat)
           call listvar(0,rt,'endtim',' ')

           call getlist(0,bffr(lch:),varcnt,-1,istat)
        endif
        return
      end


!**** integer function isvar
      integer function isvar(ARGV)
        character*(*) ARGV

        integer  istat
        character bffr*120
        logical filok

        isvar = 0
        inquire(file=ARGV,exist=filok)
! Is it an existing file?  If not then it is an error or a list root
        if (.not.filok  .or.  (index(ARGV,'.lis').ne.0)) return

! At this point we have definately found some sort of file
        open(unit=15,file=ARGV,status='OLD',iostat=istat)
        if(istat .eq. 0) then
          call getlist(15,bffr,isvar,1,istat)
          close(unit=15)
        endif
      end


!**** SUBROUTINE islist
      logical function islist(ARGV,ascname,biname,replis,bffr)
        character*(*) ARGV, ascname, biname, bffr
        logical replis

        integer alen
        logical lexist

        islist = .false.
        lexist = .false.
        ! write(*,*) 'islist(',trim(ARGV),"  ",trim(ascname),"  ",trim(biname),replis,"  ",trim(bffr)

! Is it an existing file?  If not then it is an error or a list root
        alen = len_trim(ARGV)
        if (ARGV(alen-3:alen).eq.'.lis') then
          ! accept a name that ends in list
          islist = .TRUE.
          ascname = ARGV
        elseif (index(biname,ARGV) .gt.0  .or.  ascname .eq. '-') then
          ! accept a name that matches a string in the binary name
          ! or if the list file name is not set
          islist = .TRUE.
          ascname = trim(ARGV)//'.lis'
        else
          ! accept an existing list file by that name?
          inquire(file=ARGV(:alen)//'.lis',exist=lexist)
          if(lexist) then
            ascname = ARGV(:alen)//'.lis'
            islist = .TRUE.
          end if
        end if

        ! return if this is not a list file
        if(.not.islist) return

        ! Is it an existing file?  If not print a message and return
        inquire(file=ascname,exist=lexist)
        if(.not. lexist) then
          write (6,*) 'Generating list file: ',trim(ascname)
          islist = .true.
          return
        endif


        ! are we going to try and replace an existing list file.
        if(replis) then

          ! If so then we need to open and read some information
          !    copy the critical data to an intermediate file
          write(*,*) ' opening list file: ',trim(ascname)
          open(unit=15,file=ascname,status='OLD',err=995)
          read(15,'(a)',err=995) bffr ! read the file name reported
          if (index(bffr,'.bin').ne.0) then
             ! remember the binary file name from this list if we have no other
             if (biname .eq. '-') biname = bffr
          elseif (index(bffr,'.sch')+index(bffr,'.evt') .eq. 0) then
             call FILRERR(15,1,'Not a list file ',bffr)
          endif

10        read(15,'(a)',err=998) bffr    ! read the variable list

          ! remove '   time             ' from the input line
          alen = index(bffr,'time')
          if(alen .eq. 0) goto 10       ! skip the blank line (without time)
          bffr = adjustl(bffr(alen+4:))

          close(unit=15,status='DELETE')
          !write(*,*) "isvar:",isvar, "  istat",istat," final bffr:'",trim(bffr),"'"
        endif

995     return
998   call FILRERR(15,-1,'unparsable list file ',bffr)
      end


!**** SUBROUTINE CHKBIN
      subroutine chkbin(in,argv,biname,binok)
        character*(*) argv, biname
        integer in
        logical binok

        integer bl
        real time

        bl = len_trim(argv)
! Check to see if the bin extension was specified
        if (bl.gt.4   .and.  argv(bl-3:bl).eq.'.bin') then
          inquire(file=argv,exist=binok)
          biname = argv
          if(.not. binok)  goto 100 ! missing file Generate an open error
        else
! add the .lis extension (structure modified for PC)
          inquire(file=argv(:bl)//'.bin',exist=binok)
          biname = argv(:bl)//'.bin' ! Accept the file name
        endif

        if(.not. binok) return

!****	Open and test the binary file  (Really necessary??)
        ! write(6,*) ' binary file ',trim(biname)
        open(unit=in,file=biname,form='UNFORMATTED',status='OLD',err=100)
        read  (in,err=110) time
        rewind (in)
        return

100     write(*,'(/2a)') 'Fatal Error opening binary file ',trim(biname)
        STOP
110     write(*,'(/2a)') 'Fatal Error reading binary file ',trim(biname)
        STOP
      end subroutine


      subroutine manual
       write(6,'(a)') 'Century output utility',                          &
        ' usage: list100 <.bin> <out.lis> [var|vfil] [date]',            &
        '  Users Manual',                                                &
        '    Arguments',                                                 &
        '      -h,  print help',                                         &
        '      -r   rebuilds an existing list file',                     &
        '      -d,  the next character is output delimiter (no space)',  &
        '      -dt  makes the file tab delimited.',                      &
        '      -c   "command_file" reads remaining arguments from a',    &
        '           disk file. Any argument, file names, switches and',  &
        '           variable names can be placed in the file.',          &
        '', &
        ' 1) list100 without arguments runs the interactive mode',       &
        '      described on page 6-2 of the Century manual.',            &
        ' 2) Output variables are specified by any combination of:',     &
        '       the command line or free format file,',                  &
        '       a list file.',                                           &
        '       a free format variable file,',                           &
        '       or by free format interactive input.',                   &
        '      Free format input reads space, comma, or line end',       &
        '      delimited names. Input is terminated by an empty line',   &
        '      OR end of file. Older one variable per line files are',   &
        '      still compatable.',                                       &
        ' 3) All input can be provided from the command line.',          &
        '    Commands are interpreted using these rules:',               &
        '    a) binary (.bin) and list (.lis) file extensions are optional.',&
        '    b) Default list file uses the binary file name as a root',  &
        '    c) File names must contain an alphabetic character.',       &
        '    d) The list file specified WILL BE OVER WRITTEN!!!',        &
        '    e) output variables are added without duplication',         &
        '    f) command files may include commands or just an output'//  &
        '       variable list ', &
        '    g) -r argument reads a list file for the binary file and',  &
        '        variables and then overwrites it',    &
        '       The binary file may be superseded by command line inputs',&
        '    h) The start date is the first numeric argument.', &
        '    i) The end date is the last numeric argument.',    &
        '    j) a date of -9999 indicates run start or end.',   &
        ' Usage:', &
        '   list100 -r <fil.lis>',                       &
        '   list100 <binary> variable list ',            &
        '   list100 <binary> <lisfil> -c commandfile '
        stop
       end
