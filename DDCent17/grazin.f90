
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... GRAZIN.F

      subroutine grazin(tomatch,curgraz)

      implicit none

! ... Argument declarations
      character*28 tomatch, curgraz

! ... Read in the new graze type
! modifications
!     6/2007 K. Killian
!       used the new stack management routine stakfind

      include 'const.inc'
      include 'parcp.inc'


! ... Local variables
      integer   is, m, unit
      ! real      value
      ! character varfnd*16

!...stack variables
!     Number of lines to read for each type
      integer, parameter :: GRAZLNS = 11,  OPTLN = 0,  SDEPTH=6
      real, save         :: stack(GRAZLNS+OPTLN, SDEPTH)
      integer, save      :: Tstack(SDEPTH) = 0
      character (len=28) :: stackname(SDEPTH) = '######'
      integer  stakfind
      real     chkdata

      unit = 0
!...Check the stack for the data
      is = stakfind(tomatch,SDEPTH,Tstack,stackname)

      if(is.lt.0) then
        unit = 11
        is = abs(is)
!... open the input file and find the removal option
        call oldopen(unit,'graz.100')
!       find the option in the file
        call findopt(unit,'graz.100',tomatch,GRAZLNS)
      endif

!     record the name in the name array
      stackname(is) = tomatch


      m = 1
        flgrem = chkdata(unit,'grazin','flgrem',m,stack(1,is))
        fdgrem = chkdata(unit,'grazin','fdgrem',m,stack(1,is))
        gfcret = chkdata(unit,'grazin','gfcret',m,stack(1,is))
        gret(N) = chkdata(unit,'grazin','gret',m,stack(1,is))
        gret(P) = chkdata(unit,'grazin','gret',m,stack(1,is))
        gret(S) = chkdata(unit,'grazin','gret',m,stack(1,is))
        grzeff = int(chkdata(unit,'grazin','grzeff',m,stack(1,is)))

        fecf(N) = chkdata(unit,'grazin','fecf',m,stack(1,is))
        fecf(P) = chkdata(unit,'grazin','fecf',m,stack(1,is))
        fecf(S) = chkdata(unit,'grazin','fecf',m,stack(1,is))
        feclig = chkdata(unit,'grazin','feclig',m,stack(1,is))

        if (unit.ne.0) close(unit)
        curgraz = tomatch

      return
      end
