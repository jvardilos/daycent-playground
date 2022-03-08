
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine potfor(tavedly, petdly, tfrac, tavemth, srad, daylength)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dynam.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'site.inc'

! ... Argument declarations
      real    daylength, tavedly, petdly, tfrac
      real    tavemth
      double precision srad

      ! Compute monthly potential production for forest

      ! Outputs
      ! pforc gC/m^2/time
      ! h2ogef(2)

      ! RESPPT is an array of respiration values for the forest system
      ! production components.

      ! Added savanna model, pcropc now pforc
      !                    tgprod now tfprod (BO)

      ! The change we made to the equation that is calculating growth drought
      ! stress was too severe, back it off some.  CAK - 08/29/2011


! ... Function declarations
      real     gpdf, frespr, laprod, pprdwc, lacalc
      external gpdf, frespr, laprod, pprdwc, lacalc

! ... Local variables
      real     frlive, lai, potprd
      real, parameter :: w2lang = (3600 / 41840.00)

! ... Estimate potential production based on temp & h2o

      if (tavedly .gt. 0.0) then

! ..... Calculate temperature effect on growth.
        potprd = gpdf(tavedly, ppdf(1,2), ppdf(2,2), ppdf(3,2), ppdf(4,2))

! ..... Added to match version 3.0 -lh 4/93
        potprd = potprd * .8

! ..... Calculate moisture effect on growth
        if (petdly .ge. .01) then
! ....... Calculate potential growth based on the relative water content
! ....... of the wettest soil layer, cak - 12/06/04
          !this should be h2ogef(2) not h2ogef(1)   KLK 9Jun2011
          !converted to use tree specific water stress   KLK 31Jan2011
          h2ogef(2) = 1./(1.+exp(twscoef(2)*(twscoef(1)-twstress)))
        else
          h2ogef(2) = 0.01
        endif

      ! For large wood, calculate the percentage which is live (sapwood)
        frlive = sapk / (sapk + rlwodc)

      ! Calculate LAI needed for its effect on production -rm 5/91
      ! Theoretical LAI maximum based on large wood biomass, cak - 07/24/02
      ! Include fine branch in the woody component for LAI calculation, cak - 10/20/2006
        lai = lacalc(fbrchc, rlwodc, maxlai, klai)

      ! Use solar radiation value calculated by Peter Thornton's subroutines
        ! for calculating potential production, cak - 06/18/2009
        ! Convert the solar radiation value from W/m^2 to langleys/day
        ! removed intermediate value    KLK 9Jun11
        pforc = srad * daylength * w2lang * & ! solrad
                prdx(2) * tfrac * potprd * h2ogef(2) * laprod(lai,laitop) * co2cpr(FORSYS)

! ..... Compute carbon allocation fractions for each tree part,  mdh 5/11/01
! ..... Call moved from treegrow subroutine, cak - 07/01/02
        if (pforc .gt. 0.01) then
          call treeDynC(pforc, tree_a2drat, tavemth, tree_cfrac)
        else
          pforc = 0.0
        endif

      else
        pforc = 0.
      endif

      return
      end
