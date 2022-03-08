
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... HARVIN.F

      subroutine harvin(tomatch,curharv)

      implicit none

! ... Argument declarations
      character (len=28) tomatch, curharv

! ... Read in the new harvest type
!     modifications
!     Dec/2011 K. Killian
!       used the new stack management routine stakfind

      include 'parcp.inc'
      include 'param.inc'

! ... Local variables
      integer   is, m, unit
      real      chkdata

!...stack variables
!     Number of lines to read for each type
      integer   HARVLNS, HARVOPT, SDEPTH
      parameter (HARVLNS = 6,  HARVOPT = 0,  SDEPTH=4)
      real      stack(HARVLNS, SDEPTH)
      integer   Tstack(SDEPTH), stakfind
      character stackname(SDEPTH)*28
      data      stackname, Tstack/SDEPTH*'####', SDEPTH*0/
      SAVE      Tstack, stack, stackname

      unit = 0
!...Check the stack for the data
      is = stakfind(tomatch,SDEPTH,Tstack,stackname)

      if(is.lt.0) then
        unit = 11
        is = abs(is)
!... open the input file and find the removal option
        call oldopen(unit,'harv.100')
!       find the option in the file
        call findopt(unit,'harv.100',tomatch,HARVLNS)
      endif

!     record the name in the name array
      stackname(is) = tomatch


      m = 1
        aglrem = chkdata(unit,'harvin','aglrem',m,stack(1,is))
        bglrem = chkdata(unit,'harvin','bglrem',m,stack(1,is))
        flghrv = int(chkdata(unit,'harvin','flghrv',m,stack(1,is)))
! ..... Add a check for himax > 0.0 if grain is to be harvested,  cak - 10/02/00
        if (flghrv .eq. 1) then
          if (himax .le. 0.0) call abortrun('harv.100'// &
             "grain HIMAX must be greater than 0.0.")
        endif
        rmvstr = chkdata(unit,'harvin','rmvstr',m,stack(1,is))
        remwsd = chkdata(unit,'harvin','remwsd',m,stack(1,is))
        hibg = chkdata(unit,'harvin','hibg',m,stack(1,is))

        if (unit.ne.0) close(unit)
        curharv = tomatch

      return
      end
