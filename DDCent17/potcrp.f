
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine potcrp (cancvr, tavedly, petdly, tfrac, srad,
     &                   curday, daylength)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'evivars.inc'

! ... Argument declarations
      integer curday
      real    cancvr
      real    daylength, tavedly, petdly, tfrac
      double precision srad

      ! Compute monthly production potential based upon montly precip and
      ! restrict potential production based upon method specified by grzeff.

      ! The change we made to the equation that is calculating growth drought
      ! stress was too severe, back it off some.  CAK - 08/29/2011

! ... Function declarations
      real     gpdf, pprdwc
      external gpdf, pprdwc

! ... Local variables
      real     agprod, aisc, bgp, bgprod, bioc, biof,
     &         bioprd, fracrc, potprd, ratlc, rtsh, sdlng,
     &         shdmod, subcan, temp1, temp2, temp3
      real, parameter :: w2lang = (3600 / 41840.00)

! ... Compute shading modifier for savanna
      if (cursys .eq. SAVSYS) then
        if (cancvr .le. 0.001) then
          aisc = 0.
        else
          aisc = 5 * exp(-.0035 * (rleavc*2.5)/cancvr)
        endif
        subcan = aisc/(aisc + 1.)
        shdmod = (1.0-cancvr) + (cancvr*subcan)
      else
        shdmod = 1.0
      endif

! ... Calculate moisture effect on growth
      if (petdly .ge. .01) then
! ..... Calculate potential growth based on the relative water content
! ..... of the wettest soil layer, cak - 12/06/04
        h2ogef(1) = 1.0/(1.0 + exp(cwscoef(2) *
     &                             (cwscoef(1)-cwstress)))
      else
        h2ogef(1) = 0.01
      endif

! ... Estimate plant production:
      if (tavedly .gt. 0.0) then

! ..... Calculate temperature effect on growth
        potprd = gpdf(tavedly, ppdf(1,1), ppdf(2,1), ppdf(3,1),
     &                ppdf(4,1))

! ..... Calculate biof
        if (bioflg .eq. 1) then

! ....... Calculate maximum potential effect of standing dead on plant growth
! ....... (the effect of physical obstruction of litter and standing dead)
          bioc = stdedc + .1*strucc(SRFC)
          if (bioc .le. 0.) then
            bioc = .01
          endif

          if (bioc .gt. pmxbio) then
            bioc = pmxbio
          endif
          bioprd = 1. - (bioc/(biok5+bioc))

! ....... Calculate the effect of the ratio of live biomass to dead biomass
! ....... on the reduction of potential growth rate.  The intercept of this
! ....... equation ( highest negative effect of dead plant biomass ) is equal
! ....... to bioprd when the ratio is zero.
          temp1 = (1. - bioprd)
          temp2 = temp1*0.75
          temp3 = temp1*0.25
          ratlc = aglivc/bioc
          if (ratlc .le. 1.0) then
            biof = bioprd+(temp2*ratlc)
          elseif (ratlc .le. 2.0) then
            biof = (bioprd+temp2)+temp3*(ratlc-1.)
          else
            biof = 1.0
          endif
        else
          biof = 1.0
        endif

! ..... Restriction on seedling growth
! ..... sdlng is the fraction that prdx is reduced
        if (aglivc .gt. fulcan) then
          seedl = 0
        endif

!       do not latch full canopy.  Can cause hayed/grazed crops to over produce
!        if (seedl .eq. 1) then
           sdlng = min(1.0, pltmrf + aglivc*(1-pltmrf) /fulcan)
!        else
!           sdlng = 1.0
!        endif

! ..... Compute total production, cak - 08/22/03
! ..... Added EVI option, rm - 10/09
! ..... This equation is from:
! .....     Potter et al., "Terrestrial Carbon Sinks for the United States
! .....        Predicted from MODIS Satellite Data and Ecosystem Modeling",
! .....        Earth Interactions, Volume 11 (2007), Paper No. 13.
! ..... eMax is the light utilization efficiency term (g C/MJ PAR)
        if (eviFlag) then
! ....... Compute C assimilation into plant
! .......   srad is Watts/m2 convert to  MJ/day/m2  daylength * 60 * 60 / 10^6
          tgprod = srad * daylength * potprd * h2ogef(1) * .0036 *
     &          eviday * eMax
! ....... Convert from C to biomass
          tgprod = tgprod * 2.5
        else
! ....... Use the solar radiation value as calculated by Peter Thornton's
! ....... subroutines in the calculation of potential production,
! ....... cak - 06/18/2009
! ....... Convert the solar radiation value from W/m^2 to langleys/day
          tgprod = srad * daylength * potprd * h2ogef(1) * w2lang *
     &         prdx(1) * biof * shdmod * sdlng * co2cpr(CRPSYS) * tfrac
        endif

! ..... Dynamic carbon allocation routines for crop/grsss, used to compute
! ..... root/shoot ratio, cak - 07/01/02
! ..... Do not call the dynamic carbon allocation routine when there is no
! ..... production, cak - 09/09/02
        if (tgprod .gt. 0.0) then
          call cropDynC(rtsh, fracrc)
        else
          tgprod = 0.0
          agp = 0.0
          pcropc = 0.0
          goto 40
        endif

!c ..... Change root/shoot ratio if burning occurs
!        if (firecnt .ge. 1) then
!          rtsh = rtsh + frtsh
!          firecnt = firecnt + 1
!          if (firecnt .gt. 5) then
!            firecnt = 0
!          endif
!        endif

! ..... Use the fraction of carbon allocated to the roots rather than
! ..... root to shoot ratio to determine amount of aboveground and
! ..... belowground production, cak - 08/22/03
!        bgprod = agprod * rtsh
        bgprod = tgprod * fracrc
        agprod = tgprod - bgprod

        agp = agprod
        bgp = bgprod
        tgprod = agp + bgp
      else
! ..... No production this month
        tgprod = 0.0
        agp = 0.0
        pcropc = 0.0
        goto 40
      endif

! ... Determine if grazing occurs
      if (grazcnt .ge. 1) then
        call grazrst(agp, bgp, flgrem, gremb, grzeff, rtsh, tgprod)
        grazcnt = grazcnt + 1
        if (grazcnt .gt. 31) then
          grazcnt = 0
        endif
      endif

! ... Update accumulators & compute potential C production
      ptagc = ptagc + agp/2.5
      ptbgc = ptbgc + bgp/2.5
      pcropc = tgprod / 2.5

40    continue

      return
      end

!Definitions
! One langley is one thermochemical calorie per square centimetre.
! In SI units, one langley is 41840.00 J/m² (or joules per square metre).
!                             = 11.622 watt-hours/m2.
! watt = J/s
