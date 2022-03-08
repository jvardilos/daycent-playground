
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine tcalc(avgstemp, teff, tfunc)

      implicit none

! ... Argument declarations
      real :: avgstemp, teff(4), tfunc

! ... computes the effect of temperature on decomposition.
      ! Older versions of Century used exponential and density function.
      ! Created 10/95 - rm
      !
      ! The temperature effect is now being computed using an arctangent curve.
      ! CAK - 03/16/01
      !
      ! Called From:  calcdefac
      !
      ! Variables
      !   AVGSTEMP:  weighted average of the average soil temperature in the
      !              second and third soil layers
      !   TEFF(1):   "x" location of inflection point
      !   TEFF(2):   "y" location of inflection point
      !   TEFF(3):   step size (distance from the maximum to minimum point)
      !   TEFF(4):   slope of line at inflection point


! ... Local variables
      real, save::  normalizer = 0;

       ! normalizer is the value of the numerator at 30 deg C
       ! since teff is a fix parameter and shouldn't change store the normalization.
        if(normalizer .eq. 0) normalizer = carctanf(30.0, teff(1), teff(2), teff(3), teff(4))

        tfunc = max(0.01, carctanf(avgstemp, teff(1), teff(2), teff(3), teff(4)) / normalizer)

      return
      contains
       include 'catanf.f'
      end
