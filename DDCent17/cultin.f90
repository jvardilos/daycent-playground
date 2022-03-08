
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... CULTIN.F

      subroutine cultin(tomatch,curcult)
      implicit none

! ... Argument declarations
      character (len = 28) tomatch, curcult

! ... Read in the new cult type
! modifications
!     Aug2011 KLK
!       checks to ensure input cultivation rates are between 0 and 1
!       used the new stack management routine stakfind

      include 'parcp.inc'

!...Local variables
      integer m, unit, is, i

!...stack variables
!     Number of lines to read for each type
      integer,parameter        :: CULTLNS = 11,  OPTLN = 0,  SDEPTH = 8
      real, save               :: stack(CULTLNS+OPTLN, SDEPTH)
      integer, save            :: Tstack(SDEPTH)=0
      character (len=28), save :: stackname(SDEPTH) = '######'
      integer  stakfind
      real     chkdata

      unit = 0
!...Check the stack for the data
      is = stakfind(tomatch,SDEPTH,Tstack,stackname)

      if(is.lt.0) then
        unit = 11
        is = abs(is)
!... open the input file and find the removal option
        call oldopen(unit,'cult.100')
!       find the option in the file
        call findopt(unit,'cult.100',tomatch,CULTLNS)
      endif

!     record the name in the name array
      stackname(is) = tomatch


      m = 1
        cultra(1) = chkdata(unit,'cultin','cultra',m,stack(1,is))
        cultra(2) = chkdata(unit,'cultin','cultra',m,stack(1,is))
        cultra(3) = chkdata(unit,'cultin','cultra',m,stack(1,is))
        cultra(4) = chkdata(unit,'cultin','cultra',m,stack(1,is))
        cultra(5) = chkdata(unit,'cultin','cultra',m,stack(1,is))
        cultra(6) = chkdata(unit,'cultin','cultra',m,stack(1,is))
        cultra(7) = chkdata(unit,'cultin','cultra',m,stack(1,is))
        do i=1,7
          if(cultra(i) .lt. 0.  .or. cultra(i) .gt. 1.) call     &
          abortrun("cult input "//tomatch//"cultra < 0 or >1")
        end do

        clteff(1) = chkdata(unit,'cultin','clteff',m,stack(1,is))
        clteff(2) = chkdata(unit,'cultin','clteff',m,stack(1,is))
        clteff(3) = chkdata(unit,'cultin','clteff',m,stack(1,is))
        clteff(4) = chkdata(unit,'cultin','clteff',m,stack(1,is))

        if (unit.ne.0) close(unit)
        curcult = tomatch

      return
      end
