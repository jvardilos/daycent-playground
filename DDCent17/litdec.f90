
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... LITDEC.F

      subroutine litdec(dtm, newminrl, anerdcmp, minrlavail, eins2)

      implicit none
      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'

! ... Argument declarations
      real      dtm
      real, dimension(MAXIEL)   :: minrlavail
      real, dimension(2,MAXIEL) :: eins2
      real      anerdcmp(4,2)
      integer, parameter :: met = 1, strc = 2, som23 = 3, wood = 4
      double precision newminrl

! ... Litter Decomposition
! ... Decompose structural and metabolic material for surface and soil.
! ... written by vek 04/91
      ! separate surface and soil C/E ratios MDH Nov 2012

! ... Function declarations
      logical  candec
      real     agdrat, bgdrat
      external agdrat, bgdrat, candec

! ... Local variables
      integer  iel, lyr
      real     accum(ISOS), biocnv, cfmes1, co2los, mnrflo, pheff
      real     rceto1(3), tcflow, einflo

! ... Factor to convert C to biomass is 2.5 for everything but wood.
      parameter (biocnv = 2.5)

      accum = 0.0

! ... Compute pH effect on decomposition
      pheff = max(min(carctanf(ph, 4.0, 0.5, 1.1, 0.7), 1.0), 0.0)

! ... Surface STRUCTURAL Material
      if(strucc(SRFC) .gt. 1.e-07) then

! ..... Compute total C flow out of structural in layer SRFC
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = min(strucc(SRFC),strmax(SRFC)) * defac * dec1(SRFC) *
!     &               exp(-pligst(SRFC)*strlig(SRFC)) * dtm
        tcflow = min(strucc(SRFC),strmax(SRFC)) * agdefac * dec1(SRFC) * &
                 exp(-pligst(SRFC)*strlig(SRFC)) * anerdcmp(strc,SRFC) * dtm * pheff
! ..... where
! .....   tcflow       = grams C
! .....   strucc(SRFC) = the current value for total structural C
! .....                  (grams C) in layer SRFC
! .....   strmax(SRFC) = the maximum amount of structural C in
! .....                  layer SRFC that will decompose (grams C).
! .....   agdefac      = the aboveground decomposition factor based on
! .....                  water and temperature computed in prelim and
! .....                  in cycle
! .....   dec1(SRFC)   = the intrinsic decomposition rate of
! .....                  structural C (a fixed parameter)
! .....   pligst       = a fixed parameter that represents the effect
! .....                  of lignin-to-structural-ratio on structural
! .....                  decomposition
! .....   strlig(SRFC) = the lignin content of surface structural
! .....                  residue (grams lignin/grams biomass)
! .....   dtm          = the time step in years
! .....   pheff        = pH effect on decomposition

! ..... Decompose structural into som1 and som2 with CO2 loss.
! ..... Changed csrsnk to st1c2 (10/92)
        ! use surface decomposition ratios MDH Nov 2012
        call declig(minrlavail,strlig(SRFC),SRFC,nelem,2,ps1co2,ratnew1, &
                    rsplig,tcflow,strucc,st1c2,strcis,struce, &
                    gromin,minerl,strmnr,resp,som1ci,som1e,som2ci, &
                    som2e,struc1tosom11,struc1tosom21,newminrl,eins2(SRFC,:))

      endif

! ... Soil STRUCTURAL Material
      if (strucc(SOIL) .gt. 1.e-07) then

! ..... Compute total C flow out of structural in layer SOIL
! ..... Added impact of soil anerobic conditions -rm 12/91
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = min(strucc(SOIL),strmax(SOIL)) * defac * dec1(SOIL) *
!     &           exp(-pligst(SOIL)*strlig(SOIL)) * cltfac(4) *anerb*dtm
        tcflow = min(strucc(SOIL),strmax(SOIL)) * dec1(SOIL) * &
                 exp(-pligst(SOIL)*strlig(SOIL)) * cltfac(4) * &
                     bgdefac * anerdcmp(strc,SOIL) * dtm * pheff

! ..... where
! .....   tcflow       = grams C
! .....   strucc(SOIL) = the current value for total structural C
! .....                  (grams C) in layer SOIL
! .....   strmax(SOIL) = the maximum amount of structural C in
! .....                  layer SOIL that will decompose (grams C).
! .....   bgdefac      = the belowground decomposition factor based on
! .....                  water and temperature computed in prelim and
! .....                  in cycle
! .....   dec1(SOIL)   = the intrinsic decomposition rate of
! .....                  structural C (a fixed parameter)
! .....   pligst       = a fixed parameter that represents the effect
! .....                  of lignin-to-structural-ratio on structural
! .....                  decomposition
! .....   strlig(SOIL) = the lignin content of soil structural
! .....                  residue (grams lignin/grams biomass)
! .....   cltfac(4)    = the cultivation factor for soil
! .....                  structural material (set in cycle)
! .....   anerdcmp(strc,SOIL) = impact of soil anaerobic conditions on decomposition
! .....   dtm          = the time step in years (= dt/ntspm)
! .....   pheff        = pH effect on decomposition

! ..... Decompose structural into som1 and som2 with CO2 loss.
! ..... Changed csrsnk to st2c2 (10/92)
        ! use soil decomposition ratios MDH Nov 2012
        call declig(minrlavail,strlig(SOIL),SOIL,nelem,2,ps1co2,ratnew2, &
                    rsplig,tcflow,strucc,st2c2,strcis,struce, &
                    gromin,minerl,strmnr,resp,som1ci,som1e,som2ci, &
                    som2e,struc2tosom12,struc2tosom22,newminrl,eins2(SOIL,:))

      endif
! ... End of Structural Decomposition

! ... METABOLIC Material

! ... Process each layer
      do 30 lyr = SRFC, SOIL

        if (metabc(lyr) .gt. 1.e-07) then

! ....... Determine C/E ratios for flows to SOM1
          do 10 iel=1,nelem

! ......... Compute ratios for surface metabolic residue
            if (lyr .eq. SRFC) then
              rceto1(iel) = agdrat(minrlavail,varat11,iel)

! ......... Compute ratios for soil metabolic residue
            else
              rceto1(iel) = bgdrat(minrlavail,varat12,iel)
            endif
10        continue

! ....... Compute pH effect on decomposition
          pheff = max(min(carctanf(ph, 4.8, 0.5, 1.14, 0.7), 1.0), 0.0)

! ....... Compute total C flow out of metabolic in layer lyr
! ....... Add pH effect on decomposition to calculation, cak - 08/02/02
!          tcflow = metabc(lyr) * defac * dec2(lyr) * dtm
          if (lyr .eq. SRFC) then
            tcflow = metabc(lyr) * agdefac
          else
            tcflow = metabc(lyr) * bgdefac
          endif
! ....... Added impact of soil anerobic conditions -rm 12/91
          tcflow = tcflow * dec2(lyr) * dtm * pheff * anerdcmp(met,lyr)

! ....... where:
! .......   tcflow      = grams C
! .......   metabc(lyr) = the current value for total metabolic C
! .......                 (grams C) in layer lyr
! .......   agdefac     = the aboveground decomposition factor
! .......   bgdefac     = the belowground decomposition factor
! .......   dec2(lyr)   = the intrinsic decomposition rate of
! .......                 metabolic C
! .......   dtm         = the time step in years
! .......   pheff       = pH effect on decomposition
! .......   anerb       = impact of soil anaerobic conditions on decomposition

! ....... Make sure metab does not go negative.
          if (tcflow .gt. metabc(lyr)) then
            tcflow = metabc(lyr)
          endif

! ....... If decomposition can occur,
          if (candec(nelem,minrlavail,metabc(lyr),metabe,2,lyr,rceto1)) then

! ......... CO2 loss
            co2los = tcflow * pmco2(lyr)
            if (lyr .ne. SRFC) then
              ! Add these values in for the correct ratios
              sdco2sum = sdco2sum + co2los
              ntdco2sm = ntdco2sm + co2los ! no decomposition factor
            endif

! ......... Changed csrsnk to mt1c2, mt2c2 (10/92)
            if (lyr .eq. SRFC) then
              call respir(co2los,2,lyr,metabc,metcis,mt1c2, &
                          resp,metabe,minerl,gromin,metmnr,newminrl)
            else
              call respir(co2los,2,lyr,metabc,metcis,mt2c2, &
                          resp,metabe,minerl,gromin,metmnr,newminrl)
            endif

! ......... Decompose metabolic into som1
            cfmes1 = tcflow - co2los
            if (lyr .eq. 1) then
              metc1tosom11 = cfmes1
            else
              metc2tosom12 = cfmes1
            endif

! ......... Partition and schedule C flows by isotope
            call csched (cfmes1,metcis(lyr,LABELD),metabc(lyr), &
                         metcis(lyr,UNLABL),som1ci(lyr,UNLABL), &
                         metcis(lyr,LABELD),som1ci(lyr,LABELD), &
                         1.0,accum)

! ......... Compute and schedule N, P, and S flows and update
! ......... mineralization accumulators.
            do 20 iel = 1, nelem

              call esched(cfmes1,metabc(lyr),rceto1(iel), &
                          metabe(lyr,iel),som1e(lyr,iel), &
                          minerl(SRFC,iel),mnrflo,einflo)
              call mnracc(mnrflo,gromin(iel),metmnr(lyr,iel))
! ........... newminrl should be updated only for nitrogen
! ........... akm via cak 07/31/01
              if (iel .eq. N) then
                newminrl = newminrl + mnrflo
              endif
20          continue
          endif
        endif

! ..... End of Metabolic Decomposition

! ... Next layer
30    continue

      return
      contains
       include 'catanf.f'
      end
