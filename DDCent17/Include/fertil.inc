
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      common/fertil/nitrate, ammonium, frac_nh4_fert,frac_no3_fert,      &
     &              aufert, feramt, ninhib, ninhtm, nreduce, Nscalar,    &
     &              OMADscalar, otffrc, safcdr, hrvirr, hrvafert,        &
     &               hrvfert, hrvomadc, hrvomade, afrtup, fertnet

      double precision :: nitrate(21)
      double precision :: ammonium, frac_nh4_fert, frac_no3_fert
      real    :: aufert, feramt(3), fertnet(3)
      real    :: hrvafert(3), hrvfert(3), hrvomadc, hrvomade(3)
      real    :: hrvirr
      real    :: ninhib, Nscalar(12), nreduce, OMADscalar(12)
      real    :: otffrc, safcdr
      real    :: afrtup(3)
      integer :: ninhtm
!      integer :: texture

      save/fertil/

! included Melannie's npool common.      5Mar14
!      double precision nitrate(21), ammonium, frac_nh4_fert, frac_no3_fert
!      integer          texture

!  otffrc  orchard tree fertilizer fraction, 0 < 1.0
!          this is a kludge that allocates the proportion of ammendments
!          to the tree component
!  safcdr  savanna auto fertilization crop demand ratio
!  afrtup  Automatic Fertilization uptake
! hrvafert(3)   ! harvest auto fert
! hrvfert(3)    ! harvest regular fert
! hrvomadc      ! clear harvest omad c
! hrvomade(3)   ! clear harvest omad minerals
