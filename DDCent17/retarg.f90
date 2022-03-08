      subroutine retarg(imod,ARGV)
! *********************************************************************
!  retarg   Return argument
!           This is a machine independent wrapper to hide the ugly details
!           of the argument calls. It provides all the normal functions
!           needed to access the arguments.
!           It also provides a method of imitating a command line by
!           loading and parsing a character buffer. There is a current
!           limitation of 64 arguments in the buffer mode.
!
!      NOTE:   All the MAC and OSX codes work for the respective Absoft
!              compilers. The MAC designation is for the MPW implementation.
!              There is a unix compatability library available that uses the
!              UPPERCASED versions of the UNIX variables(case is critical)
!              (However, Absoft recomends using the native API's)
!              The OSX ABSOFT lines directly implement the Absoft structures.
!              These work only if MAIN is Absoft Fortran. Under OSX, It might
!              be cleaner to directly access the C getarg functions but ....
!              that's for later
!
!  Operating modes are based in the input values of imod and ARGV.
!  Reserved messages in ARGV are checked before the imod options.
!     All reserved messages return a value in imod.
!
!  ARGV = 'arg count'
!           returns the number of arguments in imod
!           If that count is negative there was a parsing error.
!           The successfully parsed arguments can be located as usual.
!  ARGV = 'current arg'
!           returns the number of the last returned argument in ARGV
!  ARGV = 'echo'
!           prints the full command line argument to terminal
!  imod > maximum number of arguments
!           returns the string '?--?' in ARGV
!  imod > 0
!           returns the requested argument in ARGV. Similar to getarg call
!  imod = 0
!           steps through the requested arguments returning them in ARGV.
!           Returns '?--?' when all arguments have been returned.
!  imod = -1
!           sets argument mode
!  imod = -2
!           opens an external 'argument' file. All requests are then
!           read from the file.
!  imod = -21
!           opens an external 'argument' file.
!           Open error causes an abnormal termination.
!  imod = -3
!           loads the internal buffer and sets to internal mode
!  imod = -4
!           reads from standard input (interactive mode)
!  imod = -5
!           returns argument 0, the executable path
!  imod = -9
!           closes files and any internal variables
!  imod < 0   (other values)
!           returns 'ERROR: Unknown Option'
!
! Modified
!  10/02/2016 K. Killian
!           use F 2006 intrinsics to access command line
!           added a command line echo call
!  3/12/2002 K. Killian
!           changed the Mac code to be compatable with OSX. The MAC/OSX code
!           is used for both the OSX and MPW environments. The only difference
!           is the common definition.
!  1/12/2000 K. Killian
!           included a close argument
!  9/24/2000 K. Killian
!           included the Mac argument code instead of the dummy iargc/getarg
!
! OSX f90   use UNIX getarg, iargc and compile with -lU77
! *********************************************************************
      !integer, intent(inout) :: imod
      integer                 :: imod
      character,intent(inout) :: ARGV*(*)
! UNIX      !      integer iargc      ! OSX  IARGC
! DOS       !      integer nargs

! MAC commons carrying  total number of arguments,  string pointer array
! ABSOFT__ARGS is initialized in MAIN; frt0globals is initialized in MPW
! MAC       !      COMMON /frt0globals/ globalargc, globalargv
! OSX ABSOFT!      COMMON /ABSOFT__ARGS/ globalargc, globalargv  !, environc
! NOTE: we are just ignoring the last element environ,  pointer to environment

! these are the Mac variable and pointer definitions
! MAC/OSX   !      INTEGER*4 globalargc, Cstrng_ptaray(0)  !, environ
! MAC/OSX   !      POINTER (globalargv,Cstrng_ptaray)      ! handle to array

!****   local variables
      logical whsp
      integer i, ststrt, stend, larg, istat
      character tmp*4

!****   function definitions
      integer len_trim
      character READSTRING*40

!****   parameter
      integer defunit, maxarg
      parameter (defunit=128, maxarg=64)

!****   local static variables
      integer,   save :: numarg, lch
      integer,   save :: crntarg = 0, sargc = -1
      integer*2, save :: locarg(maxarg)
      logical,   save :: bffrmode =.false., filmode =.false., intmode =.false.
!      character, save ::  argbffr*512 = '?--?'
      character (len=512), save :: argbffr = '?--?'

! write(*,*) 'retarg( imod=',imod,'ARGV= ',trim(ARGV),')'

        if(sargc.lt.0) then
!****   Determine command line arguments or interactive
! UNIX        get argument number!           sargc = iargc() ! OSX IARGC
! DOS         get argument number!           sargc = nargs() -1
! MAC/OSX     get argument number!           sargc = globalargc - 1
           sargc = COMMAND_ARGUMENT_COUNT()
           numarg = sargc
        endif

!****  If/elseif structure handles the input case information

!****   Return the argument count
        if(ARGV .eq. 'arg count') then
          imod = numarg
          numarg = abs(numarg)
          return

!****   Return the number of the currently accessed argument number
        else if(ARGV .eq. 'current arg') then
          imod = crntarg
          return

!****   Return the number of the currently accessed argument number
        else if(ARGV .eq. 'echo') then
          argbffr="command:"
          call GET_COMMAND(argbffr(10:), larg, istat)
          if(istat <= 0  .and.  larg >0) call message(argbffr(:min(larg+10, len(argbffr))))
          return

!****   Return the currently accessed argument number
!****   MODE SET: Open standard input for interactive input
        else if(ARGV .eq. 'stdin') then
           intmode = .TRUE.
           numarg = 10000
           crntarg = 0
           imod = defunit
           return

!****   CLEANUP MODE: clean up any internal junk
        else if(imod.eq.-9) then
!****      close the command file if necessary
           inquire(unit=4,opened=whsp)
           if(whsp) close(4)

!****   MODE SET: Retrieve system arguments
        else if(imod.eq.-1) then
           bffrmode = .FALSE.
           filmode  = .FALSE.
           intmode  = .FALSE.
           numarg = sargc
           crntarg = 0
           return

!****   MODE SET: Open a command disk file
        else if(imod.eq.-2  .or.  imod.eq.-21) then
           filmode = .TRUE.
           call comfileopen(imod, ARGV, numarg, crntarg)
           return

!****   MODE SET: get the arguments from an internal buffer
        else if(imod.eq.-3) then
           crntarg = 0
           argbffr = ARGV
           bffrmode = .TRUE.

!          Remove any leading spaces
10         if(argbffr(1:1).ne.' ') goto 15
              argbffr(1:) =argbffr(2:)
              goto 10
15         continue

!          count and store the location of any arguments
!          initialize the search parameters
           locarg(1)= 1
           numarg = 1
           whsp = .FALSE.
           locarg(maxarg)=len_trim(ARGV)

!          loop through the remaining characters in the string
           do 20 i=2,locarg(maxarg)
             if(argbffr(i:i).eq.' ') then
!               set the space
                whsp = .TRUE.

             else
!               non space character
!               if the last character whas white the record the word start
                if(whsp) then
!                  terminate loop if we have exceeded the argument count
!                  this leaves the remainder of the string as a single argument
                   if(numarg+1.eq.maxarg) goto 100
                   numarg = numarg+1
                   locarg(numarg) = i
                endif
                whsp = .FALSE.
             endif
20         continue
           locarg(numarg+1)=locarg(maxarg)


!****   get the executable
        else if(imod.eq.-5) then
               call GET_COMMAND_ARGUMENT(0, ARGV, larg, istat) ! Fortran 2003 intrinsic
               if(larg == 0   .or.  istat .ne. 0) then
                 ARGV = 'DDCentEVI'
                 call message('Warning: executable name unavailable')
               endif

!****   normal argument requests
        else if(imod.ge.0) then
           if(imod.gt.0) then
              crntarg = imod
           else
              crntarg = crntarg+1
           end if

           if(crntarg.lt.1 .or. crntarg.gt.numarg) then
!****          Unknown argument number
               ARGV = '?--?'
               crntarg = numarg+1

           else if(bffrmode) then
!****          return a string from the buffer
               ARGV = argbffr(locarg(crntarg):locarg(crntarg+1)-1)
               ! write(*,*) "bffrmode argument: ",ARGV

           else if(intmode) then
!****          return a string from STDIN
               tmp=READSTRING(defunit,argbffr,lch,ststrt,stend,.TRUE.,0)
               ARGV = argbffr(ststrt:stend)
               ! write(*,*) "intmode argument: ",ARGV

           else if(filmode) then
!****          return a string from the file
               tmp = READSTRING(4,argbffr,lch,ststrt,stend,.TRUE., 0)
               ARGV = argbffr(ststrt:stend)
               ! write(*,*) "filmode argument: ",ARGV

           else
!****          Get the argument from the system
! UNIX      get argument !               call getarg(crntarg, ARGV) !OSX GETARG
! DOS       get argument !               call getarg(crntarg, ARGV, istat)
! MAC/OSX fortran string !               CALL c2f(VAL(Cstrng_ptaray(crntarg+1)), ARGV)
               call GET_COMMAND_ARGUMENT(crntarg, ARGV, larg, istat) ! Fortran 2003 intrinsic
!              remove any enclosing quotes since DOS doesn't bother
               if(ARGV(:1).EQ.'''' .OR. ARGV(:1).EQ.'"') then
                 ARGV=ARGV(2:)
                 larg = larg-1
               endif
               if(ARGV(larg:larg).EQ.'''' .OR. ARGV(larg:larg).EQ.'"') ARGV(larg:larg)=' '
           endif

!****   Unknown control option
        else
          ARGV = 'ERROR: Unknown Option'
        end if

!****   Normal return
        return

!****   Buffer error
100     numarg = -maxarg

        return
      contains
        subroutine comfileopen(imod, ARGV, numarg, crntarg)
          character (len=*) :: ARGV
          integer :: imod, numarg, crntarg

          integer :: istat, lch
          character (len=8) :: tmp

            call message("opening command file: "//ARGV) ! log the file read
            open(unit=4,file=ARGV,status='old',iostat=istat)
            if(istat.ne.0) then
              write(tmp,'(i4)') istat  ! signal file is unreadable
              if(imod .eq. -21) then
                ! signal the file open error
                ! abort on open error
                call abortrun("opening file "//trim(ARGV)//"istat="//trim(tmp))
              else
                ARGV = "?--?"
              endif
            endif
            numarg = 0
            crntarg = 0

            lch=0
            numarg = -1 ! the loop will count the end so subtract one from the start
            do while (lch .ge. 0)
              numarg = numarg + 1
              tmp = READSTRING(4,argbffr,lch,ststrt,stend,.TRUE., 0)
            end do
            rewind(4)
          return
        end subroutine comfileopen
      end subroutine retarg


      subroutine parsfilnam(cmnd,path,location,pl)
! *********************************************************************
!  prsfilnm
!           This is a machine independent wrapper to hide the ugly details
!           of separating portions of a file path. It returns the location
!           and length of the requested file name part.
!           NOTE: termdir MODIFIES path
!
!  Operating modes are based on the input values of cmnd.
!
!  cmnd = 'file'
!           returns the file name from the path.
!  cmnd = 'dir'
!           returns the directory portion of the path.
!  cmnd = 'termdir'
!           makes sure a path string is terminated with the correct delimiter.
!  cmnd = 'root'
!           returns the file name without the extension. For Century these
!           extensions follow the DOS form of '.???'

      character, intent(In)    :: cmnd*(*)
      character, intent(InOut) :: path*(*)
      integer,   intent(Out) :: location(2), pl


!****   internal variables
        character delim*1
! MAC  !        parameter delim = ':'
! UNIX !
        parameter (delim = '/')
! DOS  !        parameter (delim = '\')
!  ASCII encoding because the slash is a special character on some systems
! DOS  !        parameter (delim = char(92))

!****   internal variables
        integer lendir, lenext, len_trim, i
        logical, parameter :: BACK = .true.

        pl = len_trim(path)
        lenext = INDEX(path, '.', BACK)
        lendir = index(path, delim, BACK)

        lenext = 0
        do 10 i = pl,1,-1
          if(path(i:i) .eq. delim) goto 20
          if(path(i:i) .eq. '.'  .and. lenext .eq. 0) lenext = i
10      continue
        i = 0
20      lendir = i

        if (cmnd .eq. "file") then
           location = [lendir+1, pl]

        else if(cmnd .eq. 'termdir') then
           if(path(pl:pl) .ne. delim) then
              pl = pl+1
              path(pl:) = delim
           end if

        else if (cmnd .eq. "dir") then
           if(lendir .eq. 0) then
              location = [1,1]
              pl = 0
           else
              location = [1, lendir]
              pl = lendir
           end if

        else if(cmnd .eq. 'root') then
!          remove an extension if one exists and cmnd = root
           if(lenext.gt.1) then
             lenext = lenext-1
           else
             lenext = pl
           endif
           location = [lendir+1,lenext]
           pl = lenext - lendir +1

        else
           location = [0,0]
           pl = 0
        end if

      return
      end

      character*256 function prsfilnm(cmnd,path,pl)
! *********************************************************************
!  prsfilnm
!           This is a machine independent wrapper to hide the ugly details
!           of separating portions of a file path.
!
!  Operating modes are based on the input values of cmnd.
!
!  cmnd = 'file'
!           returns the file name from the path.
!  cmnd = 'dir'
!           returns the directory portion of the path.
!  cmnd = 'termdir'
!           makes sure a directory is terminated with the correct delimiter.
!  cmnd = 'root'
!           returns the file name without the extension. For Century these
!           extensions follow the DOS form of '.???'

      character path*(*), cmnd*(*)


!****   internal variables
        character delim*1
! MAC  !        parameter delim = ':'
! UNIX !
        parameter (delim = '/')
! DOS  !        parameter (delim = '\')
!  ASCII encoding because the slash is a special character on some systems
! DOS  !        parameter (delim = char(92))

!****   internal variables
        integer pl, lendir, lenext, len_trim, i

        pl = len_trim(path)
        lenext = 0
        do 10 i = pl,1,-1
          if(path(i:i) .eq. delim) goto 20
          if(path(i:i) .eq. '.'  .and. lenext .eq. 0) lenext = i
10      continue
        i = 0
20      lendir = i

        if (cmnd .eq. "file") then
           prsfilnm = path(lendir+1:pl)
           pl = pl-lendir

        else if(cmnd .eq. 'termdir') then
           prsfilnm = path(:pl)
           if(path(pl:pl) .ne. delim) then
              pl = pl+1
              prsfilnm(pl:pl) = delim
           end if

        else if (cmnd .eq. "dir") then
           if(lendir .eq. 0) then
              prsfilnm = " "
              pl = 1
           else
              prsfilnm = path(:lendir)
              pl = lendir
           end if

        else if(cmnd .eq. 'root') then
!          remove an extension if one exists and cmnd = root
           if(lenext.gt.1) pl = lenext-1
           prsfilnm = path(lendir+1:pl)
           pl = pl - lendir

        else
           prsfilnm = "unknown option "//cmnd
!           pl = len_trim(prsfilnm)
           pl = 0
        end if

      return
      end
