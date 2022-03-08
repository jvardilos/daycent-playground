
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine mthacc(month, agdefacsum, bgdefacsum)

      implicit none
      include 'const.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

! ... Argument declarations
      real    :: agdefacsum, bgdefacsum
      integer :: month

! ... Reset monthly accumulators.

! ... Local variables

! ... Initialize monthly accumulators
      stream = 0 ! array

      pet = 0
      evap = 0
      tran = 0
      pttr = 0
      rain = 0
      agdefacsum = 0.0
      bgdefacsum = 0.0
      anerb      = 0.0
      irract = 0.0
      runoff = 0.0
      amov(1:nlayer) = 0 ! array

! ... Initialize monthly co2 accumlators (10/92)
      mt1c2 = 0.0 ! array
      mt2c2 = 0.0 ! array
      st1c2 = 0.0 ! array
      st2c2 = 0.0 ! array
      s11c2 = 0.0 ! array
      s12c2 = 0.0 ! array
      s21c2 = 0.0 ! array
      s22c2 = 0.0 ! array
      s3c2  = 0.0 ! array
      wd1c2 = 0.0 ! array
      wd2c2 = 0.0 ! array
      wd3c2 = 0.0 ! array

! ... Initialize monthly accumulator for volatilization of N during
! ... harvest, senescence, and return from grazing animal waste,
! ... cak 01/02
      volpl = 0.0

! ... Initialize monthly accumulator for symbiotic N fixation to track
! ... fixation for both grasses and trees as necessary, cak - 10/15/02
      nfix = 0.0

! ... Initialize monthly C production, cak - 11/20/03
      cprodc = 0.0
      cprodf = 0.0
      eprodc(1:MAXIEL) = 0.0
      eprodf(1:MAXIEL) = 0.0

! ... Initialize monthly accumulator for soil surface temperature,
! ... cak - 11/20/03
      stempmth = 0.0

! ... Initialize monthly accumulators for maintenance respiration
! ... cak - 05/14/04
      mrspflux(CRPSYS) = 0.0
      mrspflux(FORSYS) = 0.0
      mrspmth(CRPSYS) = 0.0
      mrspmth(FORSYS) = 0.0
        cmrspflux = 0.0 ! array
        fmrspflux = 0.0 ! array
      sumrsp = 0.0

! ... Initialize monthly accumulators for growth respiration
! ... cak - 01/16/2007
      grspflux(CRPSYS) = 0.0
      grspflux(FORSYS) = 0.0
      grspmth(CRPSYS) = 0.0
      grspmth(FORSYS) = 0.0
        cgrspflux = 0.0 ! array
        fgrspflux = 0.0 ! array

! ... Initialize monthly respiration from decomposition output variables
! ... cak - 02/21/2007
      respmth(UNLABL) = 0.0
      respmth(LABELD) = 0.0
      aminrl = 0

! ... Initialize monthly soil respiration output variables
! ... cak - 02/21/2007
      srspmth(CRPSYS) = 0.0
      srspmth(FORSYS) = 0.0

! ... Initialize monthly autotrophic respiration output variables
! ... cak - 03/27/2007
      arspmth = 0.0 ! array

! ... Initialize monthly trace gas accumulator output variables,
! ... cak - 05/14/04
      N2O_month = 0.0
      NO_month = 0.0
      N2_month = 0.0
      CH4mnox = 0.0
      CH4mnpr = 0.0
      CH4mnem = 0.0
      nit_amt_month = 0.0
      pptmonth = 0.0

      agcmth(month)  = 0
      bgcjmth(month) = 0
      bgcmmth(month) = 0
      fcmth(month)   = 0

      return
      end
