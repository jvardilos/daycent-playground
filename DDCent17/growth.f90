
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine growth(tfrac, tavedly, month)
      use calflow;

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'zztim.inc'
!      include 'timvar.inc'
      include 'evivars.inc'

! ... Argument declarations
      real             tfrac, tavedly
      integer          month

! ... Simulate production for the month.
! ...   tfrac   - fraction of month over which current production event
! ...             occurs (0-1)
! ...   tavedly - mean air temperature over production period (deg C)

      ! Remove avgstemp from the growth parameter list.
      ! plant root respiration rates are controlled by the aboveground photosynthesis
      ! so use air temperature to effect on root maintenance respiration CAK 31Aug2011

! ... Function declarations
      real      fsfunc, maxswpot, rtimp
      external  fsfunc, maxswpot, rtimp

! ... Local variables
      integer          iel, lyr
      real             accum(ISOS)
      real             agfrac, amt, availm(MAXIEL)
      real             calcup, cfrac(CPARTS-1)
      real             euf(CPARTS), fsol
      real             gnfrac, rimpct
      real             tm, uptake(4,MAXIEL)
      real             fraclblstg
      real             soilm(MAXIEL)
      real             cropNfix
      real             namt
      real             cprodcdy, eprodcdy(MAXIEL)
      real             agnpp, mrspReduce, rootadj
      real             totprod
      real ::          eup(CPARTS, MAXIEL)
      logical          latgrow ! late season perennial growth storage
      character        subname*10
      double precision frac_nh4, frac_no3
      real, parameter :: toler = 1.0E-30

      if (tfrac .lt. 0.0 .or. tfrac .gt. 1.0) then
        write(subname,*) tfrac
        call abortrun('growth tfrac   ('//subname//') out of bounds')
      endif
      if (tavedly .le. -999.0) then
        write(subname,*) tavedly
        call abortrun('growth tavedly ('//subname//') out of bounds')
      endif

      subname = 'growth    '

      uptake = 0.0
      accum = 0.0

! ... labeled material fraction in the crop/grass carbohydrate storage pool
      if (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) .gt. 0.0) then
        fraclblstg = carbostg(CRPSYS,LABELD) / &
                    (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD))
      else
        ! Cindy Keough May 8, 2014,
        ! When storage pool goes to zero set the labeled fraction for to that of
        ! new material entering from the current crop.
        fraclblstg = cisofr
      endif

! ... Temperature effect on maintenance respiration for aboveground
! ... components
      mrspTempEffect(CRPSYS,SRFC) = 0.1 * exp(0.07 * tavedly)
! ... Bound maintenance respiration temperature effect between 0.0 and 1.0, cak - 09/16/02
      mrspTempEffect(CRPSYS,SRFC) = max(0.0, min(1.0, mrspTempEffect(CRPSYS,SRFC)))

! ... Temperature effect on maintenance respiration for belowground
! ... components
      mrspTempEffect(CRPSYS,SOIL) = 0.1 * exp(0.07 * tavedly)
! ... Bound maintenance respiration temperature effect between 0.0 and 1.0,
! ... cak - 09/16/02
      mrspTempEffect(CRPSYS,SOIL) = max(0.0, min(1.0, mrspTempEffect(CRPSYS,SOIL)))

! ... Add a soil water term to the root maintenance respiration
! ... equation, cak - 06/27/2007
! ... Calculate the soil water potential of the wettest soil layer
! ... in the crop rooting zone
      mrspWaterEffect(CRPSYS) = maxswpot(claypg)
      if (mrspWaterEffect(CRPSYS) .le. 76.0) then
        mrspWaterEffect(CRPSYS) = (80.0 - mrspWaterEffect(CRPSYS)) / 80.0
      else if (mrspWaterEffect(CRPSYS) .gt. 76.0) then
        mrspWaterEffect(CRPSYS) = 0.05
      endif
      mrspWaterEffect(CRPSYS) = min(1.0, mrspWaterEffect(CRPSYS))
      mrspWaterEffect(CRPSYS) = max(0.0, mrspWaterEffect(CRPSYS))

! ... Reduce the amount of maintenance respiration based on
! ... predicted annual aboveground NPP and root biomass.
! ... Calculate a predicted NPP value (gmC/m2/yr) based on average
! ... annual precipitation, use 50% of this predicited value as an
! ... approximation of aboveground production
      agnpp = (-40.0 + 3.0 * precipyr) * 0.5
! ... Adjust this value based on current root biomass
      if ((bglivcm+bglivcj) .gt. 150.0) then
        rootadj = 1.0
      elseif ((bglivcm+bglivcj) .lt. 0.0) then
        rootadj = 0.05
      else
        rootadj = line((bglivcm+bglivcj), 0.0, 0.05, 150.0, 1.0)
      endif
      agnpp = agnpp * rootadj
! ... Use two line functions to linearly decrease maintenance
! ... respiration as the amount of carbohydrate stored in the
! ... carbohydrate storage pool gets smaller based on the predicted
! ...  annual aboveground NPP calculated above, cak - 01/08/2010
      if (agnpp > 0.000000001) then
        amt = sum(carbostg(CRPSYS,:))
        if (amt .le. cmrspnpp(1) * agnpp) then
          mrspReduce = cmrspnpp(2)
        else if (amt .lt. cmrspnpp(3) * agnpp) then
          mrspReduce = line(amt, cmrspnpp(1) * agnpp, cmrspnpp(2), cmrspnpp(3) * agnpp, cmrspnpp(4))
        elseif (amt .gt. cmrspnpp(5) * agnpp) then
          mrspReduce = cmrspnpp(6)
        else
          mrspReduce = line(amt, cmrspnpp(3) * agnpp, cmrspnpp(4), cmrspnpp(5) * agnpp, cmrspnpp(6))
        endif
        mrspReduce = max(0.0, mrspReduce)

! ...   Added maintenance respiration (mrspflux) calculation. -mdh 2/99
        cmrspdyflux(ABOVE) = ckmrspmx(ABOVE) * mrspTempEffect(CRPSYS,SRFC) * aglivc * &
                             tfrac * mrspReduce
        cmrspdyflux(BELOWJ) = ckmrspmx(BELOWJ) * mrspTempEffect(CRPSYS,SOIL) * &
                              mrspWaterEffect(CRPSYS) * bglivcj * tfrac * mrspReduce
        cmrspdyflux(BELOWM) = ckmrspmx(BELOWM) * mrspTempEffect(CRPSYS,SOIL) * &
                              mrspWaterEffect(CRPSYS) * bglivcm * tfrac * mrspReduce

        ! Maintenance respiration reduces carbohydrate storage pool, mdh - 9/4/01
        mrspdyflux(CRPSYS) = sum(cmrspdyflux)
        call csched(mrspdyflux(CRPSYS),fraclblstg,1.0, &
                    carbostg(CRPSYS,UNLABL),csrsnk(UNLABL), &
                    carbostg(CRPSYS,LABELD),csrsnk(LABELD), &
                    1.0,cautoresp)

      else
        mrspReduce = 0.0
        cmrspdyflux = 0.0
        mrspdyflux(CRPSYS) = 0
      endif

! ... Determine actual production values, restricting the C/E ratios
      if (crpgrw .eq. 1 .and. pcropc .gt. 0.0 .and. .not. (senecnt .gt. 0)) then

! ..... Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, bglivcj+bglivcm)

! ..... Calculate carbon fraction in each part
        cfrac(ABOVE) = agp / tgprod
        cfrac(BELOW) = 1.0 - cfrac(ABOVE)

        ! Determine nutrients available to plants for growth.
        soilm = 0.0
        ! available nutrients for Grass/crop are in the top claypg layers, cak 01/29/03
        soilm(1:nelem) = sum(minerl(1:claypg, 1:nelem), DIM=1, &
                     MASK= minerl(1:claypg, 1:nelem) .gt. toler)

        ! save the soil total minerals; savannas change the availability KLK - 03/27/2009
        availm = soilm
!==================================================================
        ! Calculate savanna available fractions
        if (cursys .eq. SAVSYS) then
          if(otfrac .ge. 0) then
!           availm(iel)= (availm(iel)-fertnet(iel)) * (1.-otfrac) +  &
!                        fertnet(iel) * (1. - otffrc)

            ! add in designated fraction of the fertilizer in soil
            availm(1:nelem) = availm(1:nelem) * (1.-otfrac) + &
                              fertnet(1:nelem) * (otfrac - otffrc)

          else
            ! calculate the soil N availability based on tree basal area
            ! Century : tm = availm(N);  if(dofert) tm = tm - feramt(N);  tm = min(tm, 1.5)
            tm = MIN(availm(N), 1.5)
            gnfrac = MAX(0.0, MIN(1.0, & ! Limit GNFRAC between 0 and 1
               exp(-1.664*exp(-.00102*tm*sitpot)*basfc2*trbasl)))

            ! reduce the available nutrient by the crop fraction
            availm(1:nelem) = availm(1:nelem) * gnfrac
          endif
        endif
!==================================================================
! ..... Determine actual production values, restricting the C/E ratios
        ! Added crop favail as a parameter to the restrp call   KLK 30Dec13
        ! Added eup as a parameter to the restrp call   KLK May 2015
        call restrp(elimit, nelem, availm, favail(1:3,CRPSYS), cercrp, &
                    CPARTS-1, cfrac, pcropc, rimpct, crpstg, &
                    snfxmx(CRPSYS), cprodcdy, eprodcdy, uptake, eup, &
                    cropNfix, relyld)

! ..... If growth occurs...
        if (cprodcdy .gt. 0.) then

! ....... If the carbohydrate storage pool falls below a critical value
! ....... add a minimal amount of carbon from the csrsnk to allow plant
! ....... growth.  This should only occur when the plants are small.
          if (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) .lt. 15.0) then
!            write(*,*) 'Warning, carbostg pool below minimal in growth'
!            write(*,*) 'time = , carbostg = ', time,
!     &                  carbostg(CRPSYS,UNLABL)+carbostg(CRPSYS,LABELD)
            carbostg(CRPSYS,UNLABL) = carbostg(CRPSYS,UNLABL) + (15.0 * (1.0 - fraclblstg))
            carbostg(CRPSYS,LABELD) = carbostg(CRPSYS,LABELD) + (15.0 * fraclblstg)
            csrsnk(UNLABL) = csrsnk(UNLABL) - (15.0 * (1.0 - fraclblstg))
            csrsnk(LABELD) = csrsnk(LABELD) - (15.0 * fraclblstg)
          endif

! ....... Increment the counter that is tracking the number of days to
! ....... growth for the current growing season, cak - 03/11/2010
          cgrwdys = cgrwdys + 1
          ! On Nov 30, 2014, at 1:40 PM, Melannie.Hartman@colostate.edu wrote:
          !"1) Should clsgres be checked instead of flsgres (flsgres is the late
          !    season growth restriction on forests)?
          ! 2) Should the if statement check for perennials (frtcindex = 1 or 3)?
          ! All other equations that use clsgres seem to be restricted to perennials."
          ! Standardize the test for late season growth   KLK Jan 2015
          ! latgrow = (cgrwdys .gt. curgdys) .and. ((frtcindx .le. 1) .or. (frtcindx .eq. 3))
          latgrow = (clsgres .gt. 0.0) .and. (cgrwdys .gt. curgdys) .and.  &
              ((frtcindx .le. 1) .or. (frtcindx .eq. 3))

! ....... Calculations for symbiotic N fixation accumulators moved
! ....... from nutrlm subroutine, cak - 10/17/02
! ....... Compute N fixation which actually occurs and add to the
! ....... N fixation accumulator.
          nfix = nfix + cropNfix
          snfxac(CRPSYS) = snfxac(CRPSYS) + cropNfix
! ....... Add computation for nfixac -mdh 1/16/02
          nfixac = nfixac + cropNfix

! ....... C/N ratio for production
          tcnpro = cprodcdy/eprodcdy(N)

          ! Calculate production for each grass/crop part
          ! instead of repeating the fractions reorder and algebraically simplify
          agfrac = agp/tgprod
          mcprd(ABOVE)  = cprodcdy * agfrac
          mcprd(BELOWM) = cprodcdy * mrtfrac * (1.0 - agfrac)
          mcprd(BELOWJ) = cprodcdy - (mcprd(ABOVE) + mcprd(BELOWM))

          ! If we restrict production late in the growing season update mcprd
          ! Moved from the end of the routine so corrected mcprd is available
          ! for calculations rather than duplicating the if  KLK Jan 2014
          if (latgrow) mcprd(ABOVE:BELOWM) = mcprd(ABOVE:BELOWM) * (1.0 - clsgres)

! ....... Use NPP as an estimate of GPP to produce carbon that is stored
! ....... in the carbohydrate storage pool
          ! totprod = sum(mcprd(1:CPARTS))
          ! but mcprd is part allocation of cprodcdy so the sum must be cprodcdy
          if (cprodcdy .gt. 0.0) then
            totprod = cprodcdy * npp2cs(CRPSYS)
            call csched(totprod,cisofr,1.0, &
                        csrsnk(UNLABL),carbostg(CRPSYS,UNLABL), &
                        csrsnk(LABELD),carbostg(CRPSYS,LABELD), &
                        1.0,accum)
          endif

! ....... Crop/grass growth
! ....... All growth comes from the carbohydrate pool, cak - 08/12/2009
! ....... For a perennial, if we have reached the late growing season
! ....... use the late season growth restriction parameter value to
! ....... determine how much carbohydrate to flow out of the crop/grass
! .....,, carbohydrate storage pool, cak - 03/11/2010

! ....... Growth of shoots
          ! mcprd has been corrected for late growth
          call csched(mcprd(ABOVE),fraclblstg,1.0, &
                      carbostg(CRPSYS,UNLABL),aglcis(UNLABL), &
                      carbostg(CRPSYS,LABELD),aglcis(LABELD), &
                      1.0,agcisa)
! ....... Growth of juvenile roots
          call csched(mcprd(BELOWJ),fraclblstg,1.0, &
                      carbostg(CRPSYS,UNLABL),bglcisj(UNLABL), &
                      carbostg(CRPSYS,LABELD),bglcisj(LABELD), &
                      1.0,bgcisja)
! ....... Growth of mature roots
          call csched(mcprd(BELOWM),fraclblstg,1.0, &
                      carbostg(CRPSYS,UNLABL),bglcism(UNLABL), &
                      carbostg(CRPSYS,LABELD),bglcism(LABELD), &
                      1.0,bgcisma)

          ! Growth respiration;  mcprd has been corrected for late growth
          cgrspdyflux(ABOVE:BELOWM) = mcprd(ABOVE:BELOWM) * cgresp(ABOVE:BELOWM)

          ! Growth respiration is subtracted from the carbohydrate storage pool.
          grspdyflux(CRPSYS) = cgrspdyflux(ABOVE) + cgrspdyflux(BELOWJ) + cgrspdyflux(BELOWM)
          call csched(grspdyflux(CRPSYS),fraclblstg,1.0, &
                      carbostg(CRPSYS,UNLABL),csrsnk(UNLABL), &
                      carbostg(CRPSYS,LABELD),csrsnk(LABELD), &
                      1.0,cautoresp)

! ....... Actual uptake
          do iel = 1, nelem
            ! NOTE: sum(euf) = 1 sum(eup(:,iel)) = eprodcdy(iel) and
            !       sum(uptake(:,iel)) =  eprodcdy(iel)
            euf(ABOVE) = eup(ABOVE,iel) / eprodcdy(iel)
            euf(BELOWJ) = (eup(BELOW,iel) * (1.0 - mrtfrac)) / eprodcdy(iel)
            euf(BELOWM) = (eup(BELOW,iel) * mrtfrac) / eprodcdy(iel)
!      write(*,*) 'eup',eup(:,iel),sum(eup(:,iel))/eprodcdy(iel)
!      write(*,*) '  uptake', uptake(:,iel),sum(uptake(:,iel))/ eprodcdy(iel)
!      write(*,*) '  euf', euf,sum(euf)
            afrtup(iel) = uptake(EFERT,iel)
            hrvafert(iel) = hrvafert(iel) +  uptake(EFERT,iel) *afue
            ! remove soil uptake (NOT EFERT, that is autofert) from net fertilizer budget.
            fertnet(iel)  = max(fertnet(iel) - uptake(ESOIL,iel), 0.0)

! ......... Take up nutrients from internal storage pool
            ! Don't allow uptake from storage if crpstg is negative, cak 07/21/03
            ! If we have reached the late growing season use the late
            ! season growth restriction parameter value to determine how
            ! much nutrients to flow out of the crop/grass nutrient
            ! storage pool, cak - 03/11/2010
            if (crpstg(iel) .gt. 0.0) then
              if (latgrow) then
                amt = uptake(ESTOR,iel) * euf(ABOVE) * (1.0 - clsgres)
                call flow(crpstg(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt

                amt = uptake(ESTOR,iel) * euf(BELOWJ) * (1.0 - clsgres)
                call flow(crpstg(iel),bglivej(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt

                amt = uptake(ESTOR,iel) * euf(BELOWM) * (1.0 - clsgres)
                call flow(crpstg(iel),bglivem(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = uptake(ESTOR,iel) * euf(ABOVE)
                call flow(crpstg(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt

                amt = uptake(ESTOR,iel) * euf(BELOWJ)
                call flow(crpstg(iel),bglivej(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt

                amt = uptake(ESTOR,iel) * euf(BELOWM)
                call flow(crpstg(iel),bglivem(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
            endif

! ......... Take up nutrients from soil
            ! Nutrients for uptake are available in the top claypg layers,
            ! cak 01/29/03
            ! If we have reached the late growing season flow nutrients to
            ! the crop/grass storage pool rather than to the component
            ! nutrient pools based on the crop/grass late season growth
            ! restriction parameter value, cak - 03/11/2010
            do lyr = 1, claypg
              if (minerl(lyr,iel) .gt. toler) then
                fsol = 1.0
                if (iel .eq. P) then
                  fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
                endif
                call cmpnfrac(lyr,ammonium,nitrate,minerl, frac_nh4,frac_no3)
                ! Changed availm to soilm for the correct normalization, KLK - 03/27/2009
                calcup = uptake(ESOIL,iel) * minerl(lyr,iel) * fsol / soilm(iel)  !  / availm(iel)
                ! Aboveground live
                namt = 0.0
                if (latgrow) then
                  amt = calcup * euf(ABOVE) * clsgres
                  namt = namt + amt
                  call flow(minerl(lyr,iel),crpstg(iel),time,amt)
                  eupaga(iel) = eupaga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                  amt = calcup * euf(ABOVE) * (1.0 - clsgres)
                  namt = namt + amt
                  call flow(minerl(lyr,iel),aglive(iel),time,amt)
                  eupaga(iel) = eupaga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                else
                  amt = calcup * euf(ABOVE)
                  namt = namt + amt
                  call flow(minerl(lyr,iel),aglive(iel),time,amt)
                  eupaga(iel) = eupaga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                endif
                if (iel .eq. N) then
                  namt = -1.0*namt
                  call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                    ammonium, nitrate, subname)
                endif
                ! Juvenile fine roots
                namt = 0.0
                if (latgrow) then
                  amt = calcup * euf(BELOWJ) * clsgres
                  namt = namt + amt
                  call flow(minerl(lyr,iel),crpstg(iel),time,amt)
                  eupbga(iel) = eupbga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt

                  amt = calcup * euf(BELOWJ) * (1.0 - clsgres)
                  namt = namt + amt
                  call flow(minerl(lyr,iel),bglivej(iel),time,amt)
                  eupbga(iel) = eupbga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                else
                  amt = calcup * euf(BELOWJ)
                  namt = namt + amt
                  call flow(minerl(lyr,iel),bglivej(iel),time,amt)
                  eupbga(iel) = eupbga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                endif
                if (iel .eq. N) then
                  namt = -1.0*namt
                  call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                    ammonium, nitrate, subname)
                endif
                ! Mature fine roots
                namt = 0.0
                if (latgrow) then
                  amt = calcup * euf(BELOWM) * clsgres
                  namt = namt + amt
                  call flow(minerl(lyr,iel),crpstg(iel),time,amt)
                  eupbga(iel) = eupbga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                  amt = calcup * euf(BELOWM) * (1.0 - clsgres)
                  namt = namt + amt
                  call flow(minerl(lyr,iel),bglivem(iel),time,amt)
                  eupbga(iel) = eupbga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                else
                  amt = calcup * euf(BELOWM)
                  namt = namt + amt
                  call flow(minerl(lyr,iel),bglivem(iel),time,amt)
                  eupbga(iel) = eupbga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                endif
                if (iel .eq. N) then
                  namt = -1.0*namt
                  call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                    ammonium, nitrate, subname)
                endif
              endif
            end do

            ! Take up nutrients from nitrogen fixation
            if (iel .eq. N .and. cropNfix .gt. 0) then
              ! If we have reached the late growing season flow nutrients to
              ! the crop/grass storage pool rather than to the component
              ! nutrient pools based on the crop/grass late season growth
              ! restriction parameter value, cak - 03/11/2010
              if (latgrow) then
                amt = uptake(ENFIX,iel) * euf(ABOVE) * clsgres
                call flow(esrsnk(iel),crpstg(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(ENFIX,iel) * euf(BELOWJ) * clsgres
                call flow(esrsnk(iel),crpstg(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(ENFIX,iel) * euf(BELOWM) * clsgres
                call flow(esrsnk(iel),crpstg(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt

                amt = uptake(ENFIX,iel) * euf(ABOVE) * (1.0 - clsgres)
                call flow(esrsnk(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(ENFIX,iel) * euf(BELOWJ) * (1.0 - clsgres)
                call flow(esrsnk(iel),bglivej(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(ENFIX,iel) * euf(BELOWM) * (1.0 - clsgres)
                call flow(esrsnk(iel),bglivem(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = uptake(ENFIX,iel) * euf(ABOVE)
                call flow(esrsnk(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(ENFIX,iel) * euf(BELOWJ)
                call flow(esrsnk(iel),bglivej(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(ENFIX,iel) * euf(BELOWM)
                call flow(esrsnk(iel),bglivem(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
            endif

! ......... Take up nutrients from automatic fertilizer
            if (((aufert .ne. 0) .or. eviFlag) .and. &
                (uptake(EFERT,iel) .gt. 0.0)) then
              ! Automatic fertilizer added to plant pools
              ! If we have reached the late growing season flow nutrients to
              ! the crop/grass storage pool rather than to the component
              ! nutrient pools based on the crop/grass late season growth
              ! restriction parameter value, cak - 03/11/2010
              if (latgrow) then
                amt = uptake(EFERT,iel) * euf(ABOVE) * clsgres
                call flow(esrsnk(iel),crpstg(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(BELOWJ) * clsgres
                call flow(esrsnk(iel),crpstg(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(BELOWM) * clsgres
                call flow(esrsnk(iel),crpstg(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt

                amt = uptake(EFERT,iel) * euf(ABOVE) * (1.0 - clsgres)
                call flow(esrsnk(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(BELOWJ) * (1.0 - clsgres)
                call flow(esrsnk(iel),bglivej(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(BELOWM) * (1.0 - clsgres)
                call flow(esrsnk(iel),bglivem(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(ABOVE)
                call flow(esrsnk(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(BELOWJ)
                call flow(esrsnk(iel),bglivej(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(BELOWM)
                call flow(esrsnk(iel),bglivem(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
! ........... Automatic fertilizer added to mineral pool moved to simsom
            endif
          end do
        else
          ! Not enough nutrients for production this time step
          cprodcdy = 0.0
          eprodcdy(1:nelem) = 0.0
        endif

      else
        ! Else no production this time step
        cprodcdy = 0.0
        eprodcdy(1:nelem) = 0.0
      endif

! ... Accumulate monthly output
      mrspflux(CRPSYS) = mrspflux(CRPSYS) + mrspdyflux(CRPSYS)
      grspflux(CRPSYS) = grspflux(CRPSYS) + grspdyflux(CRPSYS)
      cmrspflux        = cmrspflux  + cmrspdyflux
      mrspmth(CRPSYS)  = mrspmth(CRPSYS) + sum(cmrspdyflux)
      cgrspflux        = cgrspflux  + cgrspdyflux
      grspmth(CRPSYS)  = grspmth(CRPSYS) + sum(cgrspdyflux)
      srspmth(CRPSYS)  = srspmth(CRPSYS) + &
                         cmrspdyflux(BELOWJ) + cmrspdyflux(BELOWM) + &
                         cgrspdyflux(BELOWJ) + cgrspdyflux(BELOWM)

! ... Accumulate annual output
      mrspann(CRPSYS) = mrspann(CRPSYS) + sum(cmrspdyflux)
      grspann(CRPSYS) = grspann(CRPSYS) + sum(cgrspdyflux)
      srspann(CRPSYS) = srspann(CRPSYS) + &
                        cmrspdyflux(BELOWJ) + cmrspdyflux(BELOWM) + &
                        cgrspdyflux(BELOWJ) + cgrspdyflux(BELOWM)


      ! Ignore the update if cprodcdy is Zero
      if (cprodcdy .eq. 0.) return

      ! If we restrict late season production season update production output
      if (latgrow) then
        cprodcdy = cprodcdy * (1.0 - clsgres)
        eprodcdy(1:nelem) = eprodcdy(1:nelem) * (1.0 - clsgres)
      endif

      ! Sum the daily production variables for output to the monthly *.bin file
      cprodc = cprodc + cprodcdy
      eprodc(1:nelem) = eprodc(1:nelem) + eprodcdy(1:nelem)

      return
      contains
        include 'line.f'
      end subroutine growth
