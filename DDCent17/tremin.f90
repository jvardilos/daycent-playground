
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... TREMIN.F

      subroutine tremin(tomatch,curtrm)

      implicit none

!...Read in the new tree removal type
!    Modifications:
!     4/2014   KLK
!      Use the chkdata input routine
!      Added removal for the fruit/nut pool;
!        included input stack for orchard management
!        stack is sorted so oldest entries are rolled off the stack

      include 'const.inc'
      include 'forrem.inc'

! ... Argument declarations
      character (len = 28) :: tomatch, curtrm

! ... Local variables
      integer   is, m, unit
      real      chkdata, value
      character varfnd*16

      ! stack variables
      ! Number of lines to read for each type
      integer, parameter      :: TREMLNS = 20,  TREMOPT = 5,  SDEPTH = 4
      real, save              :: stack(TREMLNS+TREMOPT, SDEPTH)
      integer                 :: stakfind
      integer, save           :: Tstack(SDEPTH) = 0
      character(len=28), save :: stackname(SDEPTH) = '####'

      unit = 0
!...Check the stack for the data
      is = stakfind(tomatch,SDEPTH,Tstack,stackname)

      if(is.lt.0) then
        unit = 11
        is = abs(is)
!... open the input file and find the removal option
        call oldopen(unit,'trem.100')
        call findopt(unit,'tremin',tomatch,TREMLNS)
      endif

!     record the name in the name array
      stackname(is) = tomatch

!        write(*,'(a,i5,a,7x,6(i4,2x,a6))') 'entry ',is,' is '//tomatch
!     &        ,(Tstack(i),stackname(min(SDEPTH,Tstack(i))),i=1,SDEPTH)
!        write(*,*) is, (stack(L,is), L=1,TREMLNS+TREMOPT)

      m = 1
        evntyp = int(chkdata(unit,'tremin','evntyp',m,stack(1,is)))

        remf(1) = chkdata(unit,'tremin','remf',m,stack(1,is)) ! LEAF
        remf(2) = chkdata(unit,'tremin','remf',m,stack(1,is)) ! fine branch
        remf(3) = chkdata(unit,'tremin','remf',m,stack(1,is)) ! large wood
        remf(4) = chkdata(unit,'tremin','remf',m,stack(1,is)) ! dead fine branch
        remf(5) = chkdata(unit,'tremin','remf',m,stack(1,is)) ! dead large wood

!       look for the new fruit/nut removal rate
        varfnd = ' '
        remf(6) = chkdata(unit,'tremin',varfnd,m,stack(1,is))
        if (unit.gt.0 .and. index(varfnd,'remf') .eq. 0) then
!         default to the leaf ratio
          remf(6) = remf(1)
          stack(m-1,is) = remf(6)
          value = chkdata(unit,'tremin','reread',m,stack(1,is))
        endif

!       the tree root die rates (not removal)
        fd(1) = chkdata(unit,'tremin','fd',m,stack(1,is))
        fd(2) = chkdata(unit,'tremin','fd',m,stack(1,is))

!       the return rates for the previously removed tree fractions
        retf(1,1) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(1,2) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(1,3) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(1,4) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(2,1) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(2,2) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(2,3) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(2,4) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(3,1) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(3,2) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(3,3) = chkdata(unit,'tremin','retf',m,stack(1,is))
        retf(3,4) = chkdata(unit,'tremin','retf',m,stack(1,is))

!       look for the fruit/nut return rates (require none or all 4)
        varfnd = ' '
        retf(4,1) = chkdata(unit,'tremin',varfnd,m,stack(1,is))
        if (unit.gt.0 .and. index(varfnd,'retf') .eq. 0) then
!         default to the leaf return ratios if the lines are missing
          retf(4,1) = retf(1,1)
          retf(4,2) = retf(1,2)
          retf(4,3) = retf(1,3)
          retf(4,4) = retf(1,4)

!         Set the values in the stack
          stack(m-1,is) = retf(4,1)
          stack(m,  is) = retf(4,2)
          stack(m+1,is) = retf(4,3)
          stack(m+2,is) = retf(4,4)
        else
          retf(4,2) = chkdata(unit,'tremin','retf',m,stack(1,is))
          retf(4,3) = chkdata(unit,'tremin','retf',m,stack(1,is))
          retf(4,4) = chkdata(unit,'tremin','retf',m,stack(1,is))
        endif


        if (unit.gt.0) close(unit)
        curtrm = tomatch

      return

      end
