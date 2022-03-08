
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

      real function carctanf(x,a,b,c,d)

        implicit none

        real x, a, b, c, d ! Arguments

        real, parameter :: PI=3.141592653589793

        !
        ! release 1.0  (first formal release of modaid) james m. vevea,  NREL
        ! reformatted f90; name changed to distinguish this from the C99 complex atan KLK
        ! james m. vevea
        ! natural resource ecology lab
        ! colorado state university
        ! fort collins, colorado  80523

        ! this is functionally equivalent to the catanf routine described in:
        ! some graphs and their functional forms
        ! technical report no. 153
        ! william parton and georgo innis  (1972)
        ! natural resource ecology lab.
        ! colorado state university
        ! fort collins, colorado  80523

        carctanf = b + (c / PI) * atan(PI * d * (x - a))

        return
      end function carctanf
