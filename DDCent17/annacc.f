
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine annacc

      implicit none
      include 'cflows.inc'
      include 'const.inc'
      include 'monprd.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'param.inc'

      ! Reset annual accumulators.
      ! NOTE: The annet annual accumulator is reset in eachyr as it is used
      !       in the calculation for non-symbiotic soil N before fixation is reset

      ! Local variables

      ! Initialize annual removal accumulators
      prcann = 0.0
      petann = 0.0
      nfixac = 0.0
!      cgrain = 0.0;      egrain = 0.0
      cgracc = 0.0
      snfxac = 0.0
      accrst = 0.0
      shrema = 0.0
      shrmai = 0.0
      sdrema = 0.0
      sdrmai = 0.0
      creta  = 0.0
      resp   = 0.0
      cautoresp = 0.0
      fautoresp = 0.0
        accrste = 0.0
        ereta   = 0.0
        shrmae  = 0.0
        sdrmae  = 0.0
        egracc  = 0.0
        ! Initialize mineralization accumulators
        tnetmn  = 0.0
        sumnrs  = 0.0
        soilnm  = 0.0

      ! Initialize litter inputs
      clitad = 0.0; elitad = 0.0 ! yearly litter additions

      ! Initialize annual C production
      cproda = 0.0
      eproda = 0.0

      ! Initialize cinputs
      cinput = 0.0

      ! Reset minimum total non-living C, an annual value
      totc = 1000000

      ! Initialize co2 accumulators (10/92)
      ast1c2 = 0.0
      ast2c2 = 0.0
      amt1c2 = 0.0
      amt2c2 = 0.0
      as11c2 = 0.0
      as12c2 = 0.0
      as21c2 = 0.0
      as22c2 = 0.0
      as3c2 = 0.0

      mrspann = 0.0  ! Initialize annual maintenance respiration accumulators
      grspann = 0.0  ! Initialize annual growth respiration accumulators
      srspann = 0.0  ! Initialize annual soil respiration accumulators

      strmac = 0.0   ! Initialize stream accumulators (04/03)

      ! Initialize annual accumulators for N volatilization (04/14/04)
      voleac = 0.0
      volgac = 0.0
      volpac = 0.0

      tgzrte = 0.0  ! Initialize E return from grazing accumulators (04/14/04)

      ! Initialize accumulators for N deposition and non-symbiotic soil N fixation
      wdfxas = 0.0
      wdfxaa = 0.0
      wdfxa = 0.0

      ! Reset accumulators for yearly trace gas output, cak - 09/23/02
      N2O_year = 0.0
      NO_year = 0.0
      N2_year = 0.0
      CH4yrox = 0.0
      CH4yrem = 0.0
      CH4yrpr = 0.0
      nit_amt_year = 0.0
      annppt = 0.0

      ! Modify irrtot output so it is annual rather than simulation long cak - 06/25/2008
      irrtot = 0.0

      ! Modify fertot output so it is annual accumulator for mineral fertilizer additions
      ! add annual accumulators for organic matter additions C, omadtot,
      ! and E omaetot(1..3), cak - 10/20/2008
      omadtot = 0.0
      fertot  = 0.0
      omaetot = 0.0

      ! Initialize annual accumulators that track decomposition carbon flows
      ametc1tosom11  = 0.0
      ametc2tosom12  = 0.0
      astruc1tosom11 = 0.0
      astruc1tosom21 = 0.0
      astruc2tosom12 = 0.0
      astruc2tosom22 = 0.0
      asom11tosom21  = 0.0
      asom12tosom22  = 0.0
      asom12tosom3   = 0.0
      asom21tosom11  = 0.0
      asom21tosom22  = 0.0
      asom22tosom12  = 0.0
      asom22tosom3   = 0.0
      asom3tosom12   = 0.0
      awood1tosom11  = 0.0
      awood1tosom21  = 0.0
      awood2tosom11  = 0.0
      awood2tosom21  = 0.0
      awood3tosom12  = 0.0
      awood3tosom22  = 0.0

      return
      end
