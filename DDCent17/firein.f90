
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... FIREIN.F

! modifications
!   4Aug2011 K. Killian
!     introduced crop stack using the new stack management routine stakfind
!     Code the full name in tomatch instead of just the first 5 characters

      subroutine firein(tomatch,curfire)
      implicit none

! ... Argument declarations
      character (len=28) tomatch,curfire

      include 'const.inc'
      include 'parcp.inc'

! ... Read in the new fire type

!...Local variables
      integer   unit, is, m
      real      chkdata

!...stack variables
      ! Number of lines to read for each type
      integer   FIRELNS, FIREOPT, SDEPTH
      parameter (FIRELNS = 20,  FIREOPT = 0,  SDEPTH=6)
      real      stack(FIRELNS, SDEPTH)
      integer   Tstack(SDEPTH), stakfind
      character stackname(SDEPTH)*28
      data      stackname, Tstack/SDEPTH*'####', SDEPTH*0/
      SAVE      Tstack, stack, stackname

      unit = 0
      ! Check the stack for the data
      is = stakfind(tomatch,SDEPTH,Tstack,stackname)

      if(is.lt.0) then
        unit = 11
        is = abs(is)
!... open the input file and find the removal option
        call oldopen(unit,'fire.100')
        ! find the option in the file
        call findopt(unit,'fire.100',tomatch,FIRELNS)
      endif

      ! record the name in the name array
      stackname(is) = tomatch

      m = 1
        flfrem =    chkdata(unit,'firein','flfrem',m,stack(1,is))
        fdfrem(1) = chkdata(unit,'firein','fdfrem',m,stack(1,is))
        fdfrem(2) = chkdata(unit,'firein','fdfrem',m,stack(1,is))

        ! Next 14 parameters read from the fire.100 file have been added
        ! to handle changes for burning dead wood and returning carbon
        ! to the som3c pool as charcoal, cak - 01/02
        fdfrem(3) = chkdata(unit,'firein','fdfrem',m,stack(1,is))
        fdfrem(4) = chkdata(unit,'firein','fdfrem',m,stack(1,is))

        fret(1,1)   = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(1,N+1) = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(1,P+1) = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(1,S+1) = chkdata(unit,'firein','fret',m,stack(1,is))

        fret(2,1)   = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(2,N+1) = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(2,P+1) = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(2,S+1) = chkdata(unit,'firein','fret',m,stack(1,is))

        fret(3,1)   = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(3,N+1) = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(3,P+1) = chkdata(unit,'firein','fret',m,stack(1,is))
        fret(3,S+1) = chkdata(unit,'firein','fret',m,stack(1,is))

        frtsh =   chkdata(unit,'firein','frtsh',m,stack(1,is))
        fnue(ABOVE) = chkdata(unit,'firein','fnue',m,stack(1,is))
        fnue(BELOW) = chkdata(unit,'firein','fnue',m,stack(1,is))

        if (unit.ne.0) close(unit)
        curfire = tomatch

      return
      end
