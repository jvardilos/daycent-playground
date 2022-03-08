
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

!     Fortran routines to print to standard output. This provides an
!     internal capability to redirect output to a file. This version
!     routes the string and prints using C file handles. If we mix each
!     STDOUT calls from both languages
!         a it works but the order of the output is not certain
!     OR  b the output from one side gets lost.
!     Neither option is desirable.
!     FORTRAN output should be done a string then "printed" with these calls.

!     Modified to use Fortran 2008 iso C binding   KLK March 2012

      subroutine stdfil(filnam)
       use iso_c_binding, only: C_CHAR, C_NULL_CHAR
       !...Write a string to standard output without advance

       ! Argument declarations
       character :: filnam*(*)

       interface
         subroutine mssgfil(string) bind(C, name="mssgfil")
           use iso_c_binding, only: c_char
           character(kind=c_char) :: string(*)
         end subroutine mssgfil
       end interface

!     Redirect standard output to a file

      !open(6,status='UNKNOWN',file=string)
      call mssgfil(trim(filnam)//char(0))

      return
      end


      subroutine message(string)
       use iso_c_binding, only: C_CHAR, C_NULL_CHAR

       !...Write a string to standard output and terminate the record

       ! Argument declarations
       character :: string*(*)
       interface
         subroutine mssg(string) bind(C, name="mssg")
           use iso_c_binding, only: c_char
           character(kind=c_char) :: string(*)
         end subroutine mssg
       end interface

!     Done in a separate subroutine so that specific OS methods of
!     writing to standard output need only to make changes here.

!      if(string .eq. 'set quiet') then
!         open(6,status='UNKNOWN',file='century.log', iostat=ostat)
!         if(ostat .gt. 0) call abortrun('Opening log file century.log')
!         return
!      endif

         call mssg(trim(string)//char(0)) !write(*,'(a)") trim(string)

      return
      end


      subroutine abortrun(string)
       use iso_c_binding, only: C_CHAR, C_NULL_CHAR

       !...Write out the string on standard output and error out the run

       ! Argument declarations
       character :: string*(*)

       interface
         subroutine abortmssg(string) bind(C, name="abortmssg")
           use iso_c_binding, only: c_char
           character(kind=c_char) :: string(*)
         end subroutine abortmssg
       end interface

        !write(*,'(a)') trim(string), 'Abnormal Termination'
        call abortmssg(trim(string)//char(0))
        stop
      end


      subroutine writdel(lun,buffr)
       use iso_c_binding, only: C_CHAR, C_NULL_CHAR
       implicit none

       ! Argument declarations
       INTEGER   lun
       CHARACTER buffr*(*)
        interface
          subroutine mssg(buffr) bind(C, name="mssg")
            use iso_c_binding, only: c_char
            character(kind=c_char) :: buffr(*)
          end subroutine mssg
        end interface

       !  This replaces the space delimiter in a data string with the delimiter
       !   and writes it to open unit lun
       !   essentially it executed the grep command s/ +/del/g
       !     lun     The unit number
       !             lun  = -1 input specifies the delimiter
       !             lun  =  0 write to terminal
       !     buffr   the string buffer

       INTEGER i, j, l, m
       CHARACTER (len=1), save :: del=' '

        if(lun .eq. -1) then
          del = buffr(1:1)
          return
        endif


        l = len_trim(buffr)
        if(del .ne. ' ') then
          buffr = ADJUSTL(buffr)
          !grep loop
          i=index(buffr," ")
          do while(i .ne. 0  .and.  i .lt. l)
            j = i+1
            m = l
            do while(buffr(j:j) .eq. " "); j = j +1; end do
            buffr(i:i) = del
            i = i+1
            if(j > i) then
               l = l - (j -i)
               buffr(i:l) = buffr(j:m)
            endif
            i = index(buffr(i:)," ") +i -1
          end do
        endif

        if(lun .eq. 0) then
          call mssg(buffr(:l))
        else
          write(lun,*) buffr(:l)
        endif
      return
      end subroutine writdel


      subroutine flushSTD(UNIT)
      USE, INTRINSIC :: ISO_FORTRAN_ENV;

      integer :: UNIT

      ! a routine to flush a fortran unit.
      ! This uses the Fortran 2003 standard flush statement
      ! then it calls the POSIX standard fsync
      !   (on windows this may have to be FlushFileBuffers())
      !
      ! This is NOT used in standard operation but is Intended to
      ! encapsulate code used in a debugging situation.
      ! fnum returns POSIX file descriptor number corresponding to the open Fortran I/O unit

      ! ON A POSIX compliant system, linux, Mac OS  NOT MinGW WINDOWS
      ! uncomment the interface and fsync calls!

      ! Declare the interface for POSIX fsync function
  !    interface
  !      function fsync (fd) bind(c,name="fsync")
  !      use iso_c_binding, only: c_int
  !        integer(c_int), value :: fd
  !        integer(c_int) :: fsync
  !      end function fsync
  !    end interface

      integer :: iostat

        if(UNIT .eq. 0) then
          flush(OUTPUT_UNIT);
  !        iostat = fsync(fnum(OUTPUT_UNIT)) ! Flush and sync
          flush(error_unit);
  !        iostat = fsync(fnum(error_unit))  ! Flush and sync
        else
          flush(LU);
  !       iostat = fsync(fnum(UNIT))                 ! Flush and sync
        endif
      return
      end subroutine flushSTD
