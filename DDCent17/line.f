
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      real function line(x, x1, y1, x2, y2)

        implicit none

        real x, x1, y1, x2, y2 ! Arguments

        !  This function is the generic equation of a line from two points.
        !  Given 2 known points and a new X, calculate Y.
        !  slope = (y2 - y1) / (x2 - x1)
        !  y = slope * (x - x2) + y2

        ! formatted to be included as a local function through CONTAINS
        !  co2eff.f dailymoist.f90 froota.f froota.f90
        !  grem.f growth.f prelim.f ramp.f treegrow.f90

        line = (y2 - y1) / (x2 - x1) * (x - x2) + y2
        return
      end
