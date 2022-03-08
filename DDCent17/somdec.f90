
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... SOMDEC.F

      subroutine somdec(amovdly, dtm, newminrl, anerdcmp, minrlavail, &
                        couts2, eins2, couts3, eins3)
      use calflow;

      implicit none
      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'zztim.inc'

! ... Argument declarations
      real   :: amovdly(CMXLYR)
      real   :: dtm
      real   :: anerdcmp(4,2)
      real, dimension(MAXIEL)  :: minrlavail
      real                      :: couts3
      real, dimension(MAXIEL)   :: eins3
      real, dimension(2)        :: couts2
      real, dimension(2,MAXIEL) :: eins2

      integer, parameter :: met = 1, strc = 2, som23 = 3, wood = 4
      double precision newminrl

! ... Soil Organic Matter Decomposition                written by vek, 04/91
! ... Decompose SOM1 (surface and soil), SOM2, and SOM3.
! ... defac = decomposition factor based on water and
! ...         temperature computed in prelim and in cycle
! ... dtm   = time step (dt/ntspm)

! ... Modified to turn on flows for surface SOM2, cak - 6/14/05
      ! Minor changes KLK Nov 2012
      ! Use agdrat for som2 surface decomposition. -MDH 9/24/2012
      ! terminated do loops with enddo instead of continue
      ! comments now start with '!' instead of '^c'

! ... Function declarations
      real     :: agdrat, bgdrat
      logical  :: candec
      external :: agdrat, bgdrat, candec

! ... Local variables
      integer :: iel
      real    :: accum(ISOS), cfs1s2, cfs1s3, cfs2s1, cfs2s3
      real    :: cfs3s1, cfsfs2, cleach, co2los, linten, mnrflo, orgflow
      real    :: pheff, rceof1, rceto1(MAXIEL), rceto2(MAXIEL)
      real    :: rceto3(MAXIEL), tcflow, ntcflw
      real    :: tval, einflo
      real    :: mix

! ... *******************************************************************
! ... Initialize ACCUM even though it is not really used
      accum = 0.0
      ! somc = 0.0

      ! C/E ratios for flows to SOM3; do now rather than repeat as needed
      do iel=1,nelem
        rceto3(iel) = bgdrat(minrlavail,varat3,iel)
      end do

! ... Surface SOM1 decomposes to SOM2 with CO2 loss
! ... Now goes to SURFACE SOM2 not SOIL SOM2, cak - 06/14/05
      if (som1c(SRFC) .gt. 1.e-07) then

        ! Determine C/E ratios for flows to surface som2
        do iel=1,nelem
          rceto2(iel) = agdrat(minrlavail,varat21,iel)
        enddo

! ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

! ..... Compute total C flow out of surface microbes.
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = som1c(SRFC) * defac * dec3(SRFC) * dtm
        tcflow = som1c(SRFC) * agdefac * dec3(SRFC) * anerdcmp(met,SRFC) * dtm * pheff
! ..... where
! .....   som1c(SRFC) =  unlabeled and labeled C in surface microbes
! .....                  (som1ci(SRFC,1)+som1ci(SRFC,2))
! .....   dec3(SRFC)  =  intrinsic decomposition rate of surface microbes

! ..... If decomposition can occur, schedule flows associated with respiration
! ..... and decomposition
        if (candec(nelem,minrlavail,som1c(SRFC),som1e,2,SRFC,rceto2)) then

! ....... CO2 loss - Compute and schedule respiration flows.
          co2los = tcflow * p1co2(SRFC)
! ....... where
! .......   p1co2(SRFC)  = set to p1co2a(SRFC) in prelim.
! .......   p1co2a(SRFC) = a fixed parameter;
! .......                  intercept parameter which controls flow from soil
! .......                  organic matter with fast turnover to CO2 (fraction
! .......                  of carbon lost to CO2 when there is no sand in the
! .......                  soil)
! ....... Changed csrsnk to s11c2 (10/92)

          call respir(co2los,2,SRFC,som1c,som1ci,s11c2,resp, &
                      som1e,minerl,gromin,s1mnr,newminrl)

! ....... Decompose Surface SOM1 to SOM2

! ....... cfsfs2 is C Flow from SurFace som1 to Som2
          cfsfs2 = tcflow - co2los
          som11tosom21 = som11tosom21 + cfsfs2

! ....... Partition and schedule C flows by isotope
          call csched(cfsfs2,som1ci(SRFC,LABELD),som1c(SRFC), &
                      som1ci(SRFC,UNLABL),som2ci(SRFC,UNLABL), &
                      som1ci(SRFC,LABELD),som2ci(SRFC,LABELD), &
                      1.0,accum)

! ....... Compute and schedule N, P, and S flows.

! ....... Update mineralization accumulators.
          do iel=1,nelem
            call esched(cfsfs2,som1c(SRFC),rceto2(iel), &
                        som1e(SRFC,iel),som2e(SRFC,iel), &
                        minerl(SRFC,iel),mnrflo,einflo)
            call mnracc(mnrflo,gromin(iel),s1mnr(SRFC,iel))
            eins2(1,iel) = eins2(1,iel) + einflo
! ......... newminrl should be updated only for nitrogen
! ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
          enddo
        endif
      endif

! ... End of SOM1 (surface layer) Decomposition

! ... *******************************************************************

! ... Soil SOM1 decomposes to soil SOM2 and SOM3 with CO2 loss and
! ... possible leaching of organics.
      if (som1c(SOIL) .gt. 1.e-07) then

! ..... Determine C/E ratios for flows to soil som2
        do iel=1,nelem
          rceto2(iel) = bgdrat(minrlavail,varat22,iel)
        enddo

! ..... Compute total C flow out of soil microbes.

! ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.8, 0.5, 1.14, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

! ..... Added impact of soil anaerobic conditions -rm 12/91
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = som1c(SOIL) * defac * dec3(SOIL) * cltfac(1) * eftext *
!     &           anerb * dtm
        ntcflw = som1c(SOIL) * dec3(SOIL) * eftext * bgdefac * &
                 anerdcmp(met,SOIL) * dtm * pheff
        tcflow = ntcflw * cltfac(1)

! ..... where
! .....   som1c(SOIL) = unlabeled and labeled C in soil microbes
! .....                 (som1ci(SOIL,1)+som1ci(SOIL,2))
! .....   dec3(SOIL)  = intrinsic decomposition rate of
! .....                 soil microbes
! .....   cltfac(1)   = cultivation factor for som1
! .....                 (set in cycle)
! .....   eftext      = effect of soil texture on the soil microbe
! .....                 decomposition rate (computed in prelim)

! ..... If soil som1 can decompose to som2, it will also go to som3.
! ..... If it can't go to som2, it can't decompose at all.

! ..... If decomposition can occur,
        if (candec(nelem,minrlavail,som1c(SOIL),som1e,2,SOIL,rceto2)) then

! ....... CO2 Loss - Compute and schedule respiration flows
          co2los = tcflow * p1co2(SOIL)
          sdco2sum = sdco2sum + co2los
          ntdco2sm = ntdco2sm + ntcflw * p1co2(SOIL)
! ....... where
! .......   p1co2(SOIL) is computed in prelim as a function of fixed parameters
! .......   p1co2a(SOIL) and p1co2b(SOIL) and soil texture.
! ....... Changed csrsnk to s21c2 (10/92)
! ....... Changed s21c2 to s12c2, cak - 08/30/2007
          call respir(co2los,2,SOIL,som1c,som1ci,s12c2,resp, &
                      som1e,minerl,gromin,s1mnr,newminrl)

! ....... Decompose Soil SOM1 to SOM3
! ....... The fraction of tcflow that goes to SOM3 is a function of
! ....... clay content (fps1s3 is computed in prelim).
          cfs1s3 = tcflow * fps1s3 * (1.0 + animpt * (1.0 -anerdcmp(met,SOIL)))
          som12tosom3 = som12tosom3 + cfs1s3

! ....... Partition and schedule C flows by isotope
          call csched(cfs1s3,som1ci(SOIL,LABELD),som1c(SOIL), &
                      som1ci(SOIL,UNLABL),som3ci(UNLABL), &
                      som1ci(SOIL,LABELD),som3ci(LABELD), &
                      1.0,accum)

! ....... Compute and schedule N, P, and S flows and update mineralization accumulators.
          do iel=1,nelem
            call esched(cfs1s3,som1c(SOIL),rceto3(iel), &
                        som1e(SOIL,iel),som3e(iel), &
                        minerl(SRFC,iel),mnrflo,einflo)
            eins3(iel) = eins3(iel) + einflo
            call mnracc(mnrflo,gromin(iel),s1mnr(SOIL,iel))
! ......... newminrl should be updated only for nitrogen  akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
          enddo

! ....... Leaching of Organics
! ....... This only occurs when the water flow out of water layer 2
! ....... exceeds a critical value.  Use the same C/N, C/P, and C/S
! ....... ratios as for the flow to SOM3.

! ....... Removed organic leaching sink and replaced it with the
! ....... stream flows. -rm 2/92
          if(amovdly(2) .gt. 0.0) then
            linten = min(1.0-(omlech(3)-amovdly(2))/omlech(3), 1.0)
            cleach = tcflow * orglch * linten
! ......... Partition and schedule C flows by isotope
            call csched(cleach,som1ci(SOIL,LABELD),som1c(SOIL), &
                        som1ci(SOIL,UNLABL),strm5u, &
                        som1ci(SOIL,LABELD),strm5l, &
                        1.0,accum)

! ......... Compute and schedule N, P, and S flows and update mineralization
! ......... accumulators.
            do iel=1,nelem

! ........... Need to use the ratio for som1 for organic leaching     -rm 3/92
!              rceof1 = som1c(SOIL) / som1e(SOIL,iel)
! ........... Dissolved organic matter leaching rates for N and P contents
! ........... were too high, add adjustment to make these rates be
! ........... consistent with observed data, cak - 04/07/03
              if (iel .eq. P) then
                rceof1 = (som1c(SOIL) / som1e(SOIL,iel)) * 35.0
              else
                rceof1 = (som1c(SOIL) / som1e(SOIL,iel)) * 2.0
              endif
              orgflow = cleach / rceof1
              call flow(som1e(SOIL,iel),stream(iel+5),time,orgflow)
            enddo

! ....... No leaching this time step
          else
            cleach = 0.
          endif

! ....... Decompose Soil SOM1 to SOM2.
! ....... SOM2 gets what's left of tcflow.
          cfs1s2 = tcflow - co2los - cfs1s3 - cleach
          som12tosom22 = som12tosom22 + cfs1s2

! ....... Partition and schedule C flows by isotope
          call csched(cfs1s2,som1ci(SOIL,LABELD),som1c(SOIL), &
                      som1ci(SOIL,UNLABL),som2ci(SOIL,UNLABL), &
                      som1ci(SOIL,LABELD),som2ci(SOIL,LABELD), &
                      1.0,accum)

! ....... Compute and schedule N, P, and S flows and update mineralization
! ....... accumulators.
          do iel=1,nelem
            call esched(cfs1s2,som1c(SOIL),rceto2(iel), &
                        som1e(SOIL,iel),som2e(SOIL,iel), &
                        minerl(SRFC,iel),mnrflo,einflo)
            call mnracc(mnrflo,gromin(iel),s1mnr(SOIL,iel))
            eins2(2,iel) = eins2(2,iel) + einflo
! ......... newminrl should be updated only for nitrogen
! ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
          enddo
        endif
      endif

! ... End of Soil SOM1 decomposition

! ... *****************************************************************

! ... Soil SOM2 decomposes to soil SOM1 and SOM3 with CO2 loss
      if (som2c(SOIL) .gt. 1.e-07) then

! ..... Determine C/E ratios for flows to soil SOM1
        do iel=1,nelem
          rceto1(iel) = bgdrat(minrlavail,varat12,iel)
        enddo

! ..... Compute total C flow out of SOM2C

! ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

! ..... Added impact of soil anaerobic conditions -rm 12/91
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = som2c(SOIL) * defac * dec5(SOIL) * cltfac(2) * anerb *
!     &           dtm
!        tcflow = som2c(SOIL) * dec5(SOIL) * cltfac(2) *
!     &           bgdefac * anerb * dtm * pheff
        ntcflw= som2c(SOIL) * dec5(SOIL) * bgdefac * anerdcmp(som23,SOIL) * dtm *pheff
        tcflow = ntcflw * cltfac(2)

! ..... where
! .....   dec5(SOIL) = intrinsic decomposition rate of soil som2
! .....   cltfac(2)  = cultivation factor for som2 (set in cycle)

! ..... If som2 can decompose to som1, it will also go to som3.
! ..... If it can't go to som1, it can't decompose at all.

! ..... If decomposition can occur,
        if (candec(nelem,minrlavail,som2c(SOIL),som2e,2,SOIL,rceto1)) then
          couts2(2) = couts2(2) + tcflow

! ....... CO2 loss - Compute and schedule respiration flows
          co2los = tcflow * p2co2(SOIL)
          sdco2sum = sdco2sum + co2los
          ntdco2sm = ntdco2sm + ntcflw * p2co2(SOIL)

! ....... Changed csrsnk to s2c2 (10/92)
! ....... Changed s2c2 to s22c2, cak - 08/30/2007
          call respir(co2los,2,SOIL,som2c,som2ci,s22c2,resp, &
                      som2e,minerl,gromin,s2mnr,newminrl)

! ....... Rearranged the order of calculation.  There is a calculated
! ....... impact of clay on the decomposition of som2 to som3.  -rm 12/91

! ....... Decompose SOM2 to SOM3, SOM3 gets what's left of tcflow.
! ....... Added impact of soil anaerobic conditions -rm 12/91
          cfs2s3 = tcflow * fps2s3 * (1.0 + animpt * (1.0 - anerdcmp(som23,SOIL)))
          som22tosom3 = som22tosom3 + cfs2s3

! ....... Partition and schedule C flows by isotope
          call csched(cfs2s3,som2ci(SOIL,LABELD),som2c(SOIL), &
                      som2ci(SOIL,UNLABL),som3ci(UNLABL), &
                      som2ci(SOIL,LABELD),som3ci(LABELD), &
                      1.0,accum)

! ....... Compute and schedule N, P, and S flows and update mineralization accumulators.
          do iel=1,nelem
            call esched(cfs2s3,som2c(SOIL),rceto3(iel), &
                        som2e(SOIL,iel),som3e(iel), &
                        minerl(SRFC,iel),mnrflo,einflo)
            call mnracc(mnrflo,gromin(iel),s2mnr(SOIL,iel))
            eins3(iel)   = eins3(iel)   + einflo
! ......... newminrl should be updated only for nitrogen
! ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
          enddo

! ....... Decompose SOM2 to SOM1

! ....... Added impact of soil anaerobic conditions -rm 12/91
          cfs2s1 = tcflow  - co2los - cfs2s3
          som22tosom12 = som22tosom12 + cfs2s1

! ....... Partition and schedule C flows by isotope
          call csched(cfs2s1,som2ci(SOIL,LABELD),som2c(SOIL), &
                      som2ci(SOIL,UNLABL),som1ci(SOIL,UNLABL), &
                      som2ci(SOIL,LABELD),som1ci(SOIL,LABELD), &
                      1.0,accum)

! ....... Compute and schedule N, P, and S flows and update mineralization
! ....... accumulators.
          do iel=1,nelem
            call esched(cfs2s1,som2c(SOIL),rceto1(iel), &
                        som2e(SOIL,iel),som1e(SOIL,iel), &
                        minerl(SRFC,iel),mnrflo,einflo)
            call mnracc(mnrflo,gromin(iel),s2mnr(SOIL,iel))
! ......... newminrl should be updated only for nitrogen
! ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
          enddo
        endif
      endif

! ... End of soil SOM2 Decompositon

! ... *****************************************************************

! ... Surface SOM2 decomposes to surface SOM1 with CO2 loss
! ... This section added by pulliam

      if (som2c(SRFC) .gt. 1.e-07) then

! ..... Determine C/E ratios for flows to surface SOM1
        do iel=1,nelem
          ! Use agdrat for above ground decomposition. -MDH 9/24/2012
          rceto1(iel) = agdrat(minrlavail,varat11,iel)
        enddo

! ..... Compute total C flow out of surface SOM2C

! ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

! ..... No impact of soil anaerobic conditions or cultivation pulliam 09/95
!        tcflow = som2c(SRFC) * defac * dec5(SRFC) * dtm
        tcflow = som2c(SRFC) * agdefac * dec5(SRFC) * anerdcmp(som23,SRFC) * dtm * pheff
! ..... where
! .....   dec5(SRFC) = intrinsic decomposition rate of surface som2

! ..... No flow to som3 from surface.

! ..... If decomposition can occur,
        if (candec(nelem,minrlavail,som2c(SRFC),som2e,2,SRFC, rceto1)) then
           couts2(1) = couts2(1) + tcflow

! ....... CO2 loss - Compute and schedule respiration flows
          co2los = tcflow * p2co2(SRFC)

! ....... Changed csrsnk to s2c2 (10/92)
! ....... Changed s2c2 to s21c2, cak - 08/30/2007
          call respir(co2los,2,SRFC,som2c,som2ci,s21c2,resp, &
                      som2e,minerl,gromin,s2mnr,newminrl)

! ....... Decompose surface SOM2 to surface SOM1

          cfs2s1 = tcflow  - co2los
          som21tosom11 = som21tosom11 + cfs2s1

! ....... Partition and schedule C flows by isotope
          call csched(cfs2s1,som2ci(SRFC,LABELD),som2c(SRFC), &
                      som2ci(SRFC,UNLABL),som1ci(SRFC,UNLABL), &
                      som2ci(SRFC,LABELD),som1ci(SRFC,LABELD), &
                      1.0,accum)

! ....... Compute and schedule N, P, and S flows and update mineralization
! ....... accumulators.
          do iel=1,nelem
            call esched(cfs2s1,som2c(SRFC),rceto1(iel), &
                        som2e(SRFC,iel),som1e(SRFC,iel), &
                        minerl(SRFC,iel),mnrflo,einflo)
            call mnracc(mnrflo,gromin(iel),s2mnr(SRFC,iel))
! ......... newminrl should be updated only for nitrogen
! ......... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
          enddo
        endif
      endif

! ... End of surface SOM2 Decompositon

! ... **********************************************************************

! ... SOM3 decomposes to soil SOM1 with CO2 loss.
      if (som3c .gt. 1.e-07) then

! ..... Determine C/E ratios for flows to SOM1.
        do iel=1,nelem
          rceto1(iel) = bgdrat(minrlavail,varat12,iel)
        enddo

! ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 3.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

! ..... Compute total C flow out of SOM3C
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = som3c * defac * dec4 * cltfac(3) * anerb * dtm
!        tcflow= som3c * dec4 * cltfac(3) * bgdefac * anerb * dtm *pheff

        ntcflw = som3c * dec4 * bgdefac * anerdcmp(som23,SOIL) * dtm *pheff
        tcflow = ntcflw * cltfac(3)

! ..... where
! .....   dec4      = intrinsic decomposition rate of som3
! .....   cltfac(3) = cultivation factor for som3 (set in cycle)

! ..... If decomposition can occur,
        if (candec(nelem,minrlavail,som3c,som3e,1,1,rceto1)) then
          couts3  = couts3 + tcflow

! ....... CO2 loss - Compute and schedule respiration flows.
          co2los = tcflow * p3co2 * anerdcmp(som23,SOIL)
          sdco2sum = sdco2sum + co2los
          ntdco2sm = ntdco2sm + ntcflw * p3co2 * anerdcmp(som23,SOIL)

! ....... Changed csrsnk to s3c2 (10/92)
          !somc(1) = som3c
          tval = sum(s3c2)
          call respir(co2los,1,SRFC,som3c,som3ci,s3c2,resp, &
                      som3e,minerl,gromin,s3mnr,newminrl)

! ....... Decompose SOM3 to soil SOM1
          cfs3s1 = tcflow - co2los
          som3tosom12 = som3tosom12 + cfs3s1

! ....... Partition and schedule C flows by isotope
          call csched(cfs3s1,som3ci(LABELD),som3c, &
                      som3ci(UNLABL),som1ci(SOIL,UNLABL), &
                      som3ci(LABELD),som1ci(SOIL,LABELD), &
                      1.0,accum)

! ....... Compute and schedule N, P, and S flows and update mineralization accumulators.
          do iel=1,nelem
            call esched(cfs3s1,som3c,rceto1(iel), &
                        som3e(iel),som1e(SOIL,iel), &
                        minerl(SRFC,iel),mnrflo,einflo)
            call mnracc(mnrflo,gromin(iel),s3mnr(iel))
            ! newminrl updated only for nitrogen .... akm via cak 07/31/01
            if (iel .eq. N) then
              newminrl = newminrl + mnrflo
            endif
          end do
        endif
      endif

! ... End of SOM3 Decomposition

! ... **********************************************************************

! ... Transfer of surface SOM2 to soil SOM2, mixing (added by pulliam)

      if (som2c(SRFC) .gt. 1.0E-7) then

        if (CURSYS .eq. FORSYS) then
          mix = tmix
        else
          mix = cmix
        endif

! ..... No pH effect on this mixing flow
        tcflow = som2c(SRFC) * mix * agdefac * anerdcmp(som23,SRFC) * dtm
        som21tosom22 = som21tosom22 + tcflow
        couts2(1) = couts2(1) + tcflow

! ..... Partition and schedule C flows by isotope
        call csched(tcflow,som2ci(SRFC,LABELD),som2c(SRFC), &
                    som2ci(SRFC,UNLABL),som2ci(SOIL,UNLABL), &
                    som2ci(SRFC,LABELD),som2ci(SOIL,LABELD), &
                    1.0,accum)

! ..... Compute and schedule N, P, and S flows
        do iel=1,nelem
          if(som2e(SRFC,iel) .gt. 0) then
            call esched(tcflow,som2c(SRFC),(som2c(SRFC)/som2e(SRFC,iel)), &
                        som2e(SRFC,iel),som2e(SOIL,iel), &
                        minerl(SRFC,iel),mnrflo,einflo)
            eins2(2,iel)  = eins2(2,iel)  + einflo
          endif
        end do
      endif

! ... **********************************************************************

      return
      contains
       include 'catanf.f'
      end
