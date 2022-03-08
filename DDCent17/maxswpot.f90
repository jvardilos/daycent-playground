
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      real function maxswpot(numlayers)

      implicit none
      include 'parfx.inc'
      include 'plot1.inc'
      include 'site.inc'

! ... Argument declarations
      integer numlayers

! ... This function calculates the soil water potential of the wettest
! ... soil layer in the plant rooting zone.

! ... Arguments:
! ...   numlayers - number of soil layers in the plant rooting zone

! ... Local variable explanations:
! ...   b      - slope of retention curve; exponent value for power function
! ...   BAR2CM - conversion factor for bars to centimeters H2O
! ...   psis   - "saturation" matric potential of "ilyr" (cm H2O ?)
! ...   swptnl - soil water potential of the current layer (bars)
! ...   theta  - volumetric soil water content * 100
!
!  Refactored original calculation:
!    to simplify the minimum
!     a Moved the root absorption limit of 30 to top setting the maximium value
!     b removed swptnl = 80.0 when asmos = 0 since it is above the 30 limit
!     c Used max function
!     d grouped all sand clay calculations
!     e after that expressed the remaining loop in fortran matrix instructions

! ... Local variables
      real, parameter  :: BAR2CM = 1024
      real             :: psis, swptnl, thetas
      double precision :: b

      ! as long as sand and clay dor not depend on lyr then these calculations
      ! are layer invariant and can be removed from the loop or array calculations
      thetas = (-14.2 * sand) - (3.7 * clay) + 50.5         ! 36.3 <= thetas <= 50.5   ~45
      b = (-0.3 * sand) + (15.7 * clay) + 3.10              ! 2.8  <= b      <= 18.8   ~8.25
      psis = 10.0**((-1.58 * sand) - (0.63 * clay) + 2.17)  ! 2.54 <= psis   <= 3.17

      ! Calculate the minimum soil water potential for the profile
!      swptnl = minval((psis/((((asmos(1:numlayers)/adep(1:numlayers)) *100) / thetas)**b)) / BAR2CM, asmos(1:numlayers) >0)
      swptnl = maxval(asmos(1:numlayers)/adep(1:numlayers), asmos(1:numlayers) >0)
      swptnl = (psis/((100 * swptnl / thetas)**b)) / BAR2CM

! *** Original loop
!    base   - base value for power function
!    lyr    - current soil layer
!    thetas - volumetric soil water content at saturation for layer (% volume)
!       maxswpot = 30.0
!       do lyr = 1, numlayers
! c ..... Calculate the soil water potential for the current soil layer
!         if (asmos(lyr) .gt. 0) then
!           thetas = (-14.2 * sand) - (3.7 * clay) + 50.5
!           b = (-0.3 * sand) + (15.7 * clay) + 3.10
!           psis = 10.0**((-1.58 * sand) - (0.63 * clay) + 2.17)
!           base =  (((asmos(lyr) / adep(lyr)) * 100) / thetas)
!           swptnl = (psis/ BAR2CM) / (base**b)
!           maxswpot = min(maxswpot, swptnl)
!         endif
!       end do

      ! Place a limit the maximum water potential of the wettest layer because roots
      ! are not able to extract water below this level, cak - 03/20/2008
      maxswpot = min(swptnl, 30.0)

      return
      end function maxswpot
