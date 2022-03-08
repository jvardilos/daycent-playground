
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine trees (bgwfunc, tavewk, tfrac, tavedly, avgstemp)
      use calflow;

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot2.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Argument declarations
      real             bgwfunc, tavedly, tavewk, tfrac
      real             avgstemp

c ... Simulate forest production for the month.

c ... Update flows so direct absorption will be accounted for
c ... before plant uptake
      call flowup(time)
      call sumcar

      ! Remove avgstemp from the growth parameter list CAK 31Aug2011
      call treegrow(tfrac, tavedly)

c ... Death of tree parts
      call wdeath(tavewk, bgwfunc, tfrac, avgstemp)

c ... Update state variables and accumulators and sum carbon isotopes.
      call flowup(time)
      call sumcar

      return
      end
