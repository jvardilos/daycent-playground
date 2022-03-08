
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine crop(time, bgwfunc, tfrac, tavedly, curday, avgstemp)
      use calflow;

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot2.inc'
      include 'timvar.inc'

c ... Argument declarations
      real             time
      real             bgwfunc
      real             tfrac
      real             tavedly
      real             avgstemp
      integer          curday

c ... Driver for calling all of crop code.

c ... If microcosm selected, skip the rest of the crop code
      if (micosm .eq. 1) then
        goto 999
      endif

c ... Update flows so direct absorption will be accounted for
c ... before plant uptake.
      call flowup(time)
      call sumcar

c ... Grow (growth checks crpgrw and exactly what should be done)
      ! Remove avgstemp from the growth parameter list CAK 31Aug2011
      call growth(tfrac, tavedly, month)

c ... Fall of standing dead
      call falstd(pltlig, tfrac)

c ... Death of roots
      call droot(pltlig, tfrac, avgstemp)

c ... Death of shoots
      call dshoot(bgwfunc, tfrac, curday)

c ... Update state variables and accumulators and sum carbon isotopes.
      call flowup(time)
      call sumcar

999   continue

      return
      end
