
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine inprac(system)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

      ! Argument declarations
      integer system

! ... Initialize annual production accumulators.
      ! changes for orchard  KLK Apr 2015
      ! - changed the termination conditions for the accumulators. Tying these to the
      !   crop and tree system ignored the savanna system
      ! - clear the fertilizer, omad and N2O long with the mineral uptake when
      !   neither system is growing
      ! - capture the production accumulators whenever either the tree or crop
      !   quits growing. This means an orchard will probably get 2 updates but
      !   since they are not cleared until neither is growing the end of growing
      !   season should be correct unless do something difficult like a winter cover grass
      ! ** These changes need to be evaluated for GDD termination and other edge conditions

      ! In the CROP system or SAVANNA, if it is the last month
      ! of the growing season reset the accumulators for the grasses.
      !if (dolast .and. (.not. crpgrw) .and. system .eq. CRPSYS) then
      if (dolast .and. crpgrw .eq. 0) then
        ! Aboveground carbon production
        agcisa = 0.0         ! array
        agcprd = agcacc
        agcacc = 0.0
        ptagc = 0.0
        ! Belowground carbon production
        bgcisja = 0.0        ! array
        bgcjprd = bgcjacc
        bgcjacc = 0.0
        bgcisma = 0.0        ! array
        bgcmprd = bgcmacc
        bgcmacc = 0.0
        ptbgc = 0.0
          ! N, P, and S uptake by plants
          eupaga = 0.0 ! array
          eupbga = 0.0 ! array

        ! capture the production accumulators
        fertprd = fertac ! array
        omadprd = omadac
        omadpre = omadae  ! array
        n2oprd = n2oacc
      endif

      ! In the FOREST system or SAVANNA, if it is the last month
      ! of the growing season reset the accumulators for the trees.
      ! if (doflst .and. (.not. forgrw) .and. system .eq. FORSYS) then
      if (doflst .and. forgrw .eq. 0) then
        ! Total forest carbon
        fcprd = fcacc
        fcacc = 0
        ! Leaf carbon production
        alvcis = 0.0        ! array
        rlvprd = rlvacc
        rlvacc = 0.0
        ! Fine root carbon production
        afrcisj = 0.0       ! array
        frtjprd = frtjacc
        frtjacc = 0.0
        afrcism = 0.0       ! array
        frtmprd = frtmacc
        frtmacc = 0.0
        ! Fine branch carbon production
        afbcis = 0.0        ! array
        fbrprd = fbracc
        fbracc = 0.0
        ! Large wood carbon production
        alwcis = 0.0        ! array
        rlwprd = rlwacc
        rlwacc = 0.0
        ! Coarse root carbon production
        acrcis = 0.0        ! array
        crtprd = crtacc
        crtacc = 0.0
        ! FRUIT/NUT carbon production
        afncis = 0.0        ! array
        frnprd = frnacc
        frnacc = 0.0
        ! N, P, and S uptake by plants
            eupprt(1:FPARTS-1, :) = 0.0    ! array

        ! capture the production accumulators
        fertprd = fertac ! array
        omadprd = omadac
        omadpre = omadae  ! array
        n2oprd = n2oacc
      endif

!      if (((.not. forgrw) .and. system .eq. FORSYS)  .or.  &
!          ((.not. crpgrw) .and. system .eq. CRPSYS)  .or.  &
!          ((.not. forgrw) .and. (.not. crpgrw) .and. system .eq. SAVSYS)) then
      ! assumes the other growth will will be zero in forsys and crpsys

      ! Reset the mineral uptake, fertilizer, omad and N2O flux if no growth is occurring
      ! do it here instead of exclusively for crops   KLK
      if (crpgrw .eq. 0  .and.  forgrw .eq. 0) then
        ! Add code to reset the growing season accumulators for fertilizer
        ! and organic matter addition and N2O flux, cak - 06/06/2008
        ! Fertilizer addition
            fertac = 0.0     ! array
            fertmth = 0.0    ! array
        ! Organic matter addition
        omadac = 0.0
          omadae = 0.0      ! array
            omadmth = 0.0   ! array
            omadmte = 0.0   ! array
        ! N2O flux
        n2oacc = 0.0
          n2omth = 0.0      ! array

          ! N, P, and S uptake by plants
          eupprd = eupacc   ! array
          eupacc = 0.0      ! array
      endif

      return
      end
