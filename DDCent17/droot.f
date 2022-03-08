
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine droot(pltlig, tfrac, avgstemp)
      use calflow;

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      pltlig(3)
      real      tfrac
      real      avgstemp

c ... Simulate death of roots for the month.

c ... Function declarations
      real      gpdf, maxswpot
      external  gpdf, maxswpot

c ... Local variables
      integer   iel
      real      fr14, recres(MAXIEL), rdeath, rtdh
      real      srfclittr, soillittr, tempeff, tmpajr, watreff
      real      accum(ISOS), cturn, temp
      ! set the Century 4 default perennial senescence.
      ! this could be replaced by a variable
      real, parameter :: rtsen = 0.30

c ... Death of roots

c ... Add code to age fine roots, juvenile fine roots age to the mature
c ... fine root pool.  Modify this subroutine so that the death rate of
c ... roots is a function of soil water potential and soil temperature,
c ... cak - 06/28/2007
c ... See:  A Model of Production and Turnover of Roots in Shortgrass Prairie
c ...       Parton, Singh, and Coleman, 1978
c ...       Journal of Applied Ecology

c ... Cap the temperature effect on fine roots at -2 and +28 degrees C
      if (avgstemp .gt. 28.0) then
        temp = 28.0
      else if (avgstemp .lt. -2.0) then
        temp = -2.0
      else
        temp = avgstemp
      endif

c ... Soil temperature effect on root death rate
      tempeff = (temp - 10.0)**2 / 4.0 * 0.00175 + 0.1
      tempeff = min(tempeff, 0.5)
c ... Soil water potential effect on root death rate
      watreff = maxswpot(claypg)
      watreff = carctanf(watreff, 35.0, 0.5, 1.0, 0.05)
c ... Root death is driven by the maximum of the soil temperature
c ... effect and the soil water potential effect on root death rate,
c ... cak - 06/28/2007
      rtdh = max(tempeff, watreff)

      rdeath = 0. ! make sure rdeath is defined.
      accum = 0.0
!      ! transfer carbon from the Upon root death for perennial plants
!c ... carbohydrate storage pool to the C source/sink, CAK - 05/14/2014
!      liveCtotal = bglivcj + bglivcm
!      liveCremoved = 0.0
      if (bglivcj .gt. 0.0) then
c ..... Death of juvenile fine roots
        if(dosene) then
          if(frtcindx .eq. 2  .or.  frtcindx .ge. 4) then ! is this an annual
            rdeath = fsdeth(2) * bglivcj ! annual remove the same fraction as shoots
          else
            ! remove the fixed fraction as long as it is smaller than the shoot removal
            ! if rtsen becomes a user input this probably should be removed!!!
            rdeath = min(fsdeth(2), rtsen) * bglivcj
          endif
        else
          rdeath = min(rdrj * tfrac * rtdh, 0.95) * bglivcj
        endif
        recres(:nelem) = bglivej(:nelem)/bglivcj
        fr14 = bglcisj(LABELD)/bglivcj
c ..... A fraction of the dead roots are transferred to the surface
c ..... litter layer, the remainder goes to the soil litter layer
c ..... cak - 05/14/2007
        srfclittr = rdeath * rdsrfc
        soillittr = rdeath - srfclittr
        call partit(srfclittr,recres,SRFC,bglcisj,bglivej,
     &              pltlig(BELOWJ),fr14)
        call partit(soillittr,recres,SOIL,bglcisj,bglivej,
     &              pltlig(BELOWJ),fr14)
      endif


c ... Soil temperature effect on aging of juvenile roots
      ! handle juvenile roots in one step so we keep dead roots from maturing
      ! this uses the same fr14 and recres values as death calculation
      ! given flow, it is easier to not mature dead roots than it is to kill
      ! newly matured roots.    KLK Nov 2012
      if (bglivcj .gt. rdeath) then
        tmpajr = gpdf(temp, 37.0, 0.0, 3.0, 3.0)
        cturn = cmxturn * tmpajr * tfrac * (bglivcj - rdeath)
        if (cturn .gt. 0) then
          call csched(cturn, fr14, 1.0,
     &                bglcisj(UNLABL), bglcism(UNLABL),
     &                bglcisj(LABELD), bglcism(LABELD),
     &                1.0, accum)
          do iel = 1, nelem
            call flow(bglivej(iel), bglivem(iel),time,cturn*recres(iel))
          end do
        endif
      endif

      if (bglivcm .gt. 0.0) then
c ..... Death of mature fine roots
        if(dosene) then
          if(frtcindx .eq. 2  .or.  frtcindx .ge. 4) then
            rdeath = fsdeth(2) * bglivcm ! annual remove the same fraction as shoots
          else
            rdeath = rtsen * bglivcm
          endif
        else
          rdeath = min(rdrm * tfrac * rtdh, 0.95) * bglivcm
        endif

        do 15 iel = 1, nelem
          recres(iel) = bglivem(iel)/bglivcm
15      continue
        fr14 = bglcism(LABELD)/bglivcm
c ..... A fraction of the dead roots are transferred to the surface
c ..... litter layer, the remainder goes to the soil litter layer
c ..... cak - 05/14/2007
        srfclittr = rdeath * rdsrfc
        soillittr = rdeath - srfclittr
        call partit(srfclittr,recres,SRFC,bglcism,bglivem,
     &              pltlig(BELOWM),fr14)
        call partit(soillittr,recres,SOIL,bglcism,bglivem,
     &              pltlig(BELOWM),fr14)
      endif

      return
      contains
       include 'catanf.f'
      end
