
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... IRRIGT.F

      real function irrigt(pptdly, petdly)

      implicit none
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'

! ... Argument declarations
      real, intent(in) :: pptdly, petdly

! ... Simulate irrigation
!
!       auirri
!            = 0 Fixed Irrigation; applied regardless of available water
!            = 1 irrigate top 30 cm to field capacity
!            = 2 irrigate with a specified amount of water applied
!            = 3 irrigate top 30 cm to field capacity plus PET
!            = 4 irrigate rooting zone to field capacity
!                ALL auto irrigation applied when available water fraction < fawhc

! ... Local variables
! ... Also, since this is a function, irrtot should not be modified
! ... in this function (suggestion).
      real      :: toth2o
      real,save :: lstnlaypg = 0, rtwcap=0, rtwilt=0, rtwfld = 0

        ! Check tavedly so that irrigation water is not added as snow
        ! Allow irrigation if the temperature is below freezing, cak - 09/16/02
        irrigt = 0.0 ! initialize the irrigation

        ! Add amount given by user
        if (auirri .eq. 0) then
          irrigt = irramt

        ! Add amount automatically
        else if (auirri .eq. 2) then
          if(avh2o(1) .le. fawhc*awhc) irrigt = irramt

        ! Add amount automatically to field capacity
        else if(auirri .eq. 1  .or.  auirri .eq. 3) then
          if(avh2o(1) .le. fawhc*awhc) then
            ! Eliminated the total soil water loop.   KLK 21 Apr 2015
            ! To within round-off, the same water deficit is derivable from the values
            ! used to trigger irrigation; it is the difference between the available
            ! water (setasmos) and the water holding capacity, awhc, (prelim).
            irrigt = awhc - avh2o(1) - pptdly ! current water deficit
            if(auirri .eq. 3) irrigt = irrigt + petdly
            irrigt = max(irrigt,0.0)
          endif

        ! Irrigate automatically to field capacity in rooting zone
        else if (auirri .eq. 4) then
          !rzwilt =  sum(awilt(1:nlaypg) * adep(1:nlaypg))
          !toth2o =  sum(asmos(1:nlaypg)) - rzwilt
          !adepsum = sum(afiel(1:nlaypg) * adep(1:nlaypg)) - rzwilt

          ! recalculate the water capacity and wilting point water any time the
          ! number of layers change
          if(lstnlaypg .ne. nlaypg) then
             rtwfld = sum(afiel(1:nlaypg) * adep(1:nlaypg))
             rtwilt = sum(awilt(1:nlaypg) * adep(1:nlaypg))
             rtwcap = rtwfld - rtwilt
             lstnlaypg = nlaypg
          endif

          toth2o =  sum(asmos(1:nlaypg))
          if ((toth2o - rtwilt) .le. rtwcap * fawhc) irrigt = max(rtwfld - toth2o - pptdly, 0.)
        endif

      irrtot = irrtot + irrigt
      hrvirr = hrvirr + irrigt

      return
      end
