
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine message(string)

! ... Argument declarations
      character*(*) string
      integer  ostat

!     Done in a separate subroutine so that specific OS methods of
!     writing to standard output need only to make changes here.

      if(string .eq. 'set quiet') then
         open(6,status='UNKNOWN',file='century.log', iostat=ostat)
         if(ostat .gt. 0) call abortrun('Opening log file century.log')
         return
      endif

         write(6,*) trim(string)

      return
      end


      subroutine abortrun(string)

      character*(*) string

!...Write out the string on standard output and error out the run
!     Done in a separate subroutine so that specific OS methods of
!     writing to standard output need only to make changes here.

      write(6,*) 'Error: '//trim(string)
      write(6,*) 'Abnormal Termination'
      stop
      end

      subroutine writdel(lun,buffr)
      implicit none

      INTEGER   lun
      CHARACTER buffr*(*)
!
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
          call message(buffr(:l))
        else
          write(lun,*) buffr(:l)
        endif
      return
      end subroutine writdel
