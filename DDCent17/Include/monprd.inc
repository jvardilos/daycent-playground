
!               Copyright 1996 Colorado State University
!                       All Rights Reserved


! ... Code specific to GRIDDED CENTURY
! ... Monthly production variables.

      common/monprd/nit_amt_year, nit_amt_month, mcprd, mfprd,           &
     &              mrspdyflux, cmrspdyflux, fmrspdyflux,                &
     &              grspdyflux, cgrspdyflux, fgrspdyflux,                &
     &              mrspTempEffect, mrspWaterEffect,                     &
     &              N2O_year, NO_year, N2_year,                          &
     &              N2O_month, NO_month, N2_month,                       &
     &              stempmth, annppt, pptmonth,                          &
     &              cwstress, gwstress, twstress

      double precision nit_amt_year, nit_amt_month  ! stored as gC/ha

      real mcprd(3)
      real fmrspdyflux(7), fgrspdyflux(7), mfprd(7) ! These should be dimensioned FPART
      real mrspdyflux(2), cmrspdyflux(3)
      real grspdyflux(2), cgrspdyflux(3)
      real mrspTempEffect(2,2), mrspWaterEffect(2)
      real N2O_year, NO_year, N2_year               ! stored as gC/ha
      real N2O_month, NO_month, N2_month            ! stored as gC/ha
      real stempmth, annppt, pptmonth
      real cwstress, gwstress, twstress


      save /monprd/
