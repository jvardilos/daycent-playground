!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine cropDynC(rtsh, fracrc)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'potent.inc'

! ... Argument declarations
      real fracrc, rtsh

! ... Compute carbon allocation fractions for aboveground and belowground
! ... plant parts.
! ...   agprod - estimated aboveground plant production
! ...   rtsh   - root/shoot ratio

! ... Function declarations
      real     froota, rtimp
      external froota, rtimp

! ... Local variables
      integer  iel, lyr
      real     agprod, availm(MAXIEL), demand, eavail(MAXIEL), maxNfix
      real     rimpct, totale, toler, effprc

      toler = 1.0E-30

! ... Estimate the fraction of carbon going to the roots
      if (frtcindx .eq. 0) then
! ..... Use Great Plains equation for root to shoot ratio
!        Limit the great plains equation          KLK 8/15/2002
!        Make sure precipitation is above zero in the Great Plains equation.
!        Since these are input parameters this may not be obvious. First ensure
!        precip is positive. Then check that it is larger than a positive zero.
!        rtsh = (bgppa + grwprc*bgppb)/(agppa + grwprc*agppb)
!        fracrc = 1.0 / (1.0/rtsh + 1.0)
        effprc = max(grwprc,0.1)
        if(agppb.ne.0.) effprc = max(effprc, -agppa/agppb*1.8)
        rtsh = (bgppa + effprc*bgppb)/(agppa + effprc*agppb)
        ! Make sure the rtsh does not exceed an unreasonable value
        rtsh = min(rtsh,5.0)
        fracrc = min(rtsh / (1.0 + rtsh), 0.99)

      elseif (frtcindx .eq. 1 .or. frtcindx .eq. 3) then
! ..... A perennial plant (grass)
        fracrc = (cfrtcw(1)+cfrtcw(2)+cfrtcn(1)+cfrtcn(2))/4.0
      elseif (frtcindx .eq. 2 .or. frtcindx .ge. 4) then
! ..... An annual plant (crop)
        fracrc = (frtc(1)+frtc(2))/2.0
      endif

! ... Estimate total production
!      tgprod = agprod / (1 - fracrc)
! ... Estimate aboveground production, cak - 08/22/03
      agprod = tgprod * (1 - fracrc)

! ... Determine nutrients available to plants for growth.
      do 20 iel = 1, nelem
        availm(iel) = 0.0
! ..... Nutrients available to grasses/crops are in the top claypg layers,
! ..... cak 01/29/03
        do 30 lyr = 1, claypg
          if (minerl(lyr,iel) .gt. toler) then
            availm(iel) = availm(iel) + minerl(lyr, iel)
          endif
30      continue
20    continue

! ... Calculate impact of root biomass on available nutrients
      rimpct = rtimp(riint, rictrl, bglivcj+bglivcm)

! ... Calculate soil available nutrients, based on a maximum fraction
! ... (favail) and the impact of root biomass (rimpct), adding storage.
      ! converted to use crop favail  KLK 30Jan13
      do iel = 1, nelem
        eavail(iel) = (availm(iel) * favail(iel,CRPSYS) * rimpct) +
     &                 crpstg(iel)
      end do

! ... Compute the minimum and maximum C/E ratios
      call fltce(nelem, aglivc, co2cce)

! ... Estimate the demand
      do 50  iel = 1, nelem
! ..... Initialize fixation to 0
        maxNfix = 0.0
! ..... N FIXATION
        if (iel .eq. N) then
          maxNfix = snfxmx(CRPSYS) * (tgprod / 2.5)
        endif
! ..... DEMAND based on the maximum E/C ratio.
        demand = 0.0
        demand = demand + ((agprod / 2.5) *
     &                    (1.0 / cercrp(IMIN,ABOVE,iel)))
        demand = demand + ((tgprod - agprod) / 2.5) *
     &                    (1.0 / cercrp(IMIN,BELOW,iel))
        totale = eavail(iel) + maxNfix

! ..... New calculation -mdh 5/10/01
        crop_a2drat(iel) = min(1.0, totale / demand)
        crop_a2drat(iel) = max(0.0, crop_a2drat(iel))
50    continue

! ... New way of calculating fracrc, the fraction of root carbon. -mdh 3/21/00
! ... crop_a2drat(iel) = ratio of available mineral to mineral demand

      fracrc = froota(crop_a2drat,h2ogef(1),CRPSYS)
      ! not sure rtsh is used but this will stop a divide by zero
      if(fracrc .le. 0.99) then
        rtsh = fracrc/(1 - fracrc)
      else
        rtsh = 99
      endif

      return
      end
