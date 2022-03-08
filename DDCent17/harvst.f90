
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine harvst(month, pltlig, curday)
      use calflow;

      implicit none
      include 'const.inc'
      include 'fertil.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

! ... Argument declarations
      integer   month, curday
      real      pltlig(3)

! ... Harvest the crop
      ! modified Nov 12 KLK
      ! moved storage removed during root harvest to harvest sections so it can be
      ! added to the harvested E. Transfer storage to metabolic litter during other
      ! root death. Helps fix a bug where crpstg is non-zero with no live crop/grass.


! ... Local variables
      integer   iel, mm, hmonth(2)
      real      accum(ISOS), addsdc, addsde(MAXIEL), bgd
      real      cisbgdj(ISOS), cisbgdm(ISOS), cstraw, ctubesj, ctubesm
      real      etubes, etubesj(MAXIEL), etubesm(MAXIEL), fr14
      real      recres(MAXIEL), recresj(MAXIEL), recresm(MAXIEL)
      real      resid, reside(MAXIEL), sumpttr, sumtran, harv_volpl
      real      stdstraw, srfclittrj, srfclittrm, soillittrj
      real      soillittrm, esrfclittrj(MAXIEL), esoillittrj(MAXIEL)
      real      esrfclittrm(MAXIEL), esoillittrm(MAXIEL)
      real      estraw(MAXIEL), estdstraw(MAXIEL)
      real      rfracm, storhrv, stordth(MAXIEL), storsrf

      ! functions
      character*10 :: crpchar

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

! ... Initialization
      cgrain = 0.0
      crmvst = 0.0
      cstraw = 0.0
      stdstraw = 0.0
      addsdc = 0.0
      resid = 0.0
      ctubesj = 0.0
      ctubesm = 0.0
      etubes = 0.0
      srfclittrj = 0.0
      srfclittrm = 0.0
      soillittrj = 0.0
      soillittrm = 0.0
        egrain = 0.0
        ermvst = 0.0
        estraw = 0.0
        estdstraw = 0.0
        addsde = 0.0
        reside = 0.0
        etubesj = 0.0
        etubesm = 0.0
        esrfclittrj = 0.0
        esrfclittrm = 0.0
        esoillittrj = 0.0
        esoillittrm = 0.0

! ... Check that there is material to harvest
      if (aglivc .le. 0.001 .and. (flghrv .eq. 1  .or.  &
          (flghrv .eq. 0  .and.  stdedc .le. 0.001))) return

      ! mature root fraction; used to allocate stored E   KLK 2Nov12
      rfracm = 0.
      if(bglivcm .gt. 0.) rfracm = bglivcm/(bglivcm + bglivcj)

! ... Carbon

! ... Grain
! ... (Alister's new way of calculating:)
      if (flghrv .eq. 1) then
        if (frtcindx .lt. 5) then
          sumtran = 0
          sumpttr = 0
          hmonth(1) = month - himon(1)
          hmonth(2) = month - himon(2)
          if (hmonth(1) .lt. 1) then
            hmonth(1) = hmonth(1) + MONTHS
          endif
          if (hmonth(2) .lt. 1)  then
            hmonth(2) = hmonth(2) + MONTHS
          endif
          if (hmonth(2) .ge. hmonth(1)) then
            do 10 mm = hmonth(1), hmonth(2)
              sumtran = sumtran + htran(mm)
              sumpttr = sumpttr + hpttr(mm)
10          continue
          else
            do 15 mm = hmonth(1), MONTHS
              sumtran = sumtran + htran(mm)
              sumpttr = sumpttr + hpttr(mm)
15          continue
            do 16 mm = 1, hmonth(2)
              sumtran = sumtran + htran(mm)
              sumpttr = sumpttr + hpttr(mm)
16          continue
          endif
          if (sumpttr .eq. 0.0) then
            sumpttr = 0.0001
          endif
          hi = himax * (1.0 - hiwsf * (1.0 - (sumtran / sumpttr)))
        else if (frtcindx .ge. 5) then
          if (grnfldys .gt. 0) then
            hi = himax * (1.0 - hiwsf * (1.0 - (gwstress/grnfldys)))
          else
            hi = 0.0
          endif
        endif
        cgrain = hi * aglivc * (1.0 - aglrem)
        call csched(cgrain,aglcis(LABELD),aglivc, &
                    aglcis(UNLABL),csrsnk(UNLABL), &
                    aglcis(LABELD),csrsnk(LABELD), &
                    1.0,cisgra)
! ..... Straw
        cstraw = aglivc * (1.0 - aglrem) - cgrain
! ..... Straw removal
        crmvst = rmvstr * cstraw
        accrst = accrst + crmvst
!      write(*,*) "hi",time,hi, aglivc,cgrain,cstraw, crmvst, accrst,    &
!     &            sumtran,sumpttr,sumtran/sumpttr
        call csched(crmvst,aglcis(LABELD),aglivc, &
                    aglcis(UNLABL),csrsnk(UNLABL), &
                    aglcis(LABELD),csrsnk(LABELD), &
                    1.0,accum)
! ..... Some straw will remain as standing dead
        addsdc = remwsd * (cstraw-crmvst)
        call csched(addsdc,aglcis(LABELD),aglivc, &
                    aglcis(UNLABL),stdcis(UNLABL), &
                    aglcis(LABELD),stdcis(LABELD), &
                    1.0,accum)

! ... Non-grain harvest
      else
        cgrain = 0.0
! ..... Straw can come from aboveground live carbon and/or standing
! ..... dead carbon
        if (aglivc .gt. 0.001) then
          cstraw = aglivc * (1.0 - aglrem)
          call csched(cstraw,aglcis(LABELD),aglivc, &
                      aglcis(UNLABL),csrsnk(UNLABL), &
                      aglcis(LABELD),csrsnk(LABELD), &
                      1.0,accum)
        else
          cstraw = 0.0
        endif
        if (stdedc .gt. 0.001) then
          stdstraw = stdedc * (1.0 - aglrem)
          call csched(stdstraw,stdcis(LABELD),stdedc, &
                      stdcis(UNLABL),csrsnk(UNLABL), &
                      stdcis(LABELD),csrsnk(LABELD), &
                      1.0,accum)
        else
          stdstraw = 0.0
        endif
        crmvst = cstraw + stdstraw
        accrst = accrst + crmvst
      endif

! ... Other elements
      do 20 iel = 1, nelem

! ..... Grain
        if (flghrv .eq. 1) then
          egrain(iel) = efrgrn(iel) * aglive(iel) * (1.0 - aglrem) * sqrt(hi/himax)
          call flow(aglive(iel),esrsnk(iel),time,egrain(iel))
! ....... Volatilization of N from plants
! ....... Use the local variable harv_volpl so that volatilization that
! ....... occurs at harvest and senescence can both be tracked, see dshoot
! ....... cak - 01/02
          if (iel .eq. N) then
            harv_volpl = vlossp * aglive(iel)
            call flow(aglive(iel),esrsnk(iel),time,harv_volpl)
            volpl = volpl + harv_volpl
            volpla = volpla + harv_volpl
            volpac = volpac + harv_volpl
! ......... N/C ratio in straw
            recres(iel) = ((aglive(iel) - harv_volpl) * (1.0 - aglrem) - &
                            egrain(iel)) / cstraw
          else
! ......... P/C, or S/C ratio in straw
            recres(iel) = (aglive(iel) * (1.0 - aglrem) - egrain(iel)) / cstraw
          endif
! ....... Straw removal
          ermvst(iel) = crmvst * recres(iel)
          call flow(aglive(iel),esrsnk(iel),time,ermvst(iel))
! ....... Some straw remains as standing dead
          addsde(iel) = addsdc * recres(iel)
          call flow(aglive(iel),stdede(iel),time,addsde(iel))
        else
          egrain(iel) = 0.0
! ....... E/C ratio in aboveground live removed as straw/hay
          if (aglivc .gt. 0.001) then
            recres(iel) = aglive(iel) / aglivc
! ......... Straw/hay removal from aboveground live
            estraw(iel) = cstraw * recres(iel)
            call flow(aglive(iel),esrsnk(iel),time,estraw(iel))
          endif
! ....... E/C ratio in standing dead removed as straw/hay
          if (stdedc .gt. 0.001) then
            recres(iel) = stdede(iel) / stdedc
! ......... Straw/hay removal from standing dead
            estdstraw(iel) = stdstraw * recres(iel)
            call flow(stdede(iel),esrsnk(iel),time,estdstraw(iel))
          endif
          ermvst(iel) = estraw(iel) + estdstraw(iel)
          accrste(iel) = accrste(iel) + ermvst(iel)
        endif
20    continue

      if (flghrv .eq. 1) then
! ..... Partition c, n, p, and s in remaining straw into top layer
! ..... of structural and metabolic
        resid = cstraw - crmvst - addsdc
        if (aglivc .gt. 0.0) then
          fr14 = aglcis(LABELD)/aglivc
        else
          fr14 = 0.0
        endif
        call partit(resid,recres,1,aglcis,aglive,pltlig(ABOVE),fr14)
        do 25 iel = 1, nelem
          reside(iel) = resid * recres(iel)
25      continue
      endif

! ... Below ground removal (root harvest) -lh 8/91
! ... Harvest both juvenile and mature roots, cak - 05/18/2007

      stordth(1:nelem)  = crpstg(1:nelem) * (1.0 - bglrem) ! crop storage lost during root death

! ... Juvenile root harvest
      if (bglivcj .gt. 0.001) then
        ctubesj = hibg * bglivcj * (1.0 - bglrem)
        cgrain = cgrain + ctubesj
        call csched(ctubesj,bglcisj(LABELD),bglivcj, &
                    bglcisj(UNLABL), csrsnk(UNLABL), &
                    bglcisj(LABELD), csrsnk(LABELD), &
                    1.0,cisgra)
        ! removed if (bglivcj .gt. 0.0001); previous if requires
        ! bglivcj > 0.001  to get here KLK Nov 2012
        cisbgdj(LABELD) = (bglivcj * (1.0 - bglrem) - ctubesj) * &
                          (bglcisj(LABELD) / bglivcj)
        cisbgdj(UNLABL) = (bglivcj * (1.0 - bglrem) - ctubesj) - cisbgdj(LABELD)
        do iel = 1, nelem
          etubesj(iel) = hibg * (1.0 - bglrem) * bglivej(iel)
          call flow(bglivej(iel), esrsnk(iel), time, etubesj(iel))
          ! add stored E to the harvest  KLK 2 Nov 12
          storhrv    = hibg * stordth(iel) * (1.-rfracm)
          call flow(crpstg(iel), esrsnk(iel), time, storhrv)
          egrain(iel) = egrain(iel) + etubesj(iel) + storhrv
          recresj(iel) = bglivej(iel) / bglivcj
        end do
      else
        cisbgdj = 0.0
        recresj = 0.0
      endif

! ... Mature root harvest
      if (bglivcm .gt. 0.001) then
        ctubesm = hibg * bglivcm * (1.0 - bglrem)
        cgrain = cgrain + ctubesm
        call csched(ctubesm,bglcism(LABELD),bglivcm, &
                    bglcism(UNLABL), csrsnk(UNLABL), &
                    bglcism(LABELD), csrsnk(LABELD), &
                    1.0,cisgra)
        ! removed if(bglivcm .gt. 0.0001); again bglivcj > 0.001 to get here
        ! KLK Nov 2012
        cisbgdm(LABELD) = (bglivcm * (1.0 - bglrem) - ctubesm) * &
                          (bglcism(LABELD) / bglivcm)
        cisbgdm(UNLABL) = (bglivcm * (1.0 - bglrem) - ctubesm) - cisbgdm(LABELD)
        do 35 iel = 1, nelem
          etubesm(iel) = hibg * (1.0 - bglrem) * bglivem(iel)
          call flow(bglivem(iel), esrsnk(iel), time, etubesm(iel))
          ! add stored E to the harvest  KLK 2 Nov 12
          storhrv      = hibg * stordth(iel) * rfracm
          call flow(crpstg(iel), esrsnk(iel), time, storhrv)
          egrain(iel) = egrain(iel) + etubesm(iel) + storhrv
          recresm(iel) = bglivem(iel) / bglivcm
35      continue
      else
        cisbgdm = 0.0
        recresm = 0.0
      endif

! ... Calculation of accumulator for grain production
      cgracc = cgracc + cgrain
      egracc(1:nelem) = egracc(1:nelem)+egrain(1:nelem)

! ... Partition c, n, p, and s in remaining roots into bottom layer of
! ... structural and metabolic
! ... Juvenile roots
      bgd = cisbgdj(LABELD) + cisbgdj(UNLABL)
      if(bgd > 0) then
        if (bglivcj .gt. 0.0) then
          fr14 = bglcisj(LABELD)/bglivcj
        else
          fr14 = 0.0
        endif
        ! A fraction of the dead roots are transferred to the surface
        ! litter layer, the remainder goes to the soil litter layer cak - 05/14/2007
        srfclittrj = bgd * rdsrfc
        soillittrj = bgd - srfclittrj
        call partit(srfclittrj, recresj, SRFC, bglcisj, bglivej, &
                   pltlig(BELOWJ), fr14)
        call partit(soillittrj, recresj, SOIL, bglcisj, bglivej, &
                 pltlig(BELOWJ), fr14)
        esrfclittrj(1:nelem) = srfclittrj * recresj(1:nelem)
        esoillittrj(1:nelem) = soillittrj * recresj(1:nelem)
      endif

! ... Mature roots
      bgd = cisbgdm(LABELD) + cisbgdm(UNLABL)
      if(bgd > 0) then
        if (bglivcm .gt. 0.0) then
          fr14 = bglcism(LABELD)/bglivcm
        else
          fr14 = 0.0
        endif
       ! A fraction of the dead roots are transferred to the surface
       ! litter layer, the remainder goes to the soil litter layer  cak - 05/14/2007
        srfclittrm = bgd * rdsrfc
        soillittrm = bgd - srfclittrm
        call partit(srfclittrm, recresm, SRFC, bglcism, bglivem, &
                   pltlig(BELOWM), fr14)
        call partit(soillittrm, recresm, SOIL, bglcism, bglivem, &
                   pltlig(BELOWM), fr14)
        esrfclittrm(1:nelem) = srfclittrm * recresm(1:nelem)
        esoillittrm(1:nelem) = soillittrm * recresm(1:nelem)
      endif

      ! storage E in remaining dead roots flows to metabolic
      ! determine fraction of roots going to surface
      ! correct the denominator to sum soil mature and juvenile litter  KLK 23 June 2015
      storsrf = min((srfclittrj + srfclittrm)/ &
                    ((srfclittrj + srfclittrm) + soillittrj + soillittrm), 1.0)
      ! loop for flow commands
      do iel = 1, nelem
        stordth(iel) = stordth(iel) * (1.0-hibg)  ! storage in non harvested kill
        call flow(crpstg(iel), metabe(SRFC,iel), time, stordth(iel) * storsrf) ! surface litter
        call flow(crpstg(iel), metabe(SOIL,iel), time, stordth(iel) * (1.0-storsrf)) ! soil litter
      end do

! ... Write output to the harvest.csv file
      if (flghrv .eq. 1) then
        cstraw = crmvst
        estraw(N) = ermvst(N)
        estraw(P) = ermvst(P)
        estraw(S) = ermvst(S)
      endif
!        afertapp(1:nelem) = gafertot(1:nelem) - safertot(1:nelem)
!        fertapp(1:nelem) = gfertot(1:nelem) - sfertot(1:nelem)
!        omaeapp(1:nelem) = gomaetot(1:nelem) - somaetot(1:nelem)
!        safertot(1:nelem) = gafertot(1:nelem)
!        sfertot(1:nelem) = gfertot(1:nelem)
!        somaetot(1:nelem) = gomaetot(1:nelem)
!      irrapp = girrtot - sirrtot
!      omadapp = gomadtot - somadtot
!      sirrtot = girrtot
!      somadtot = gomadtot
      call wrtharvest(time, curday, trim(crpchar(crpval))//char(0), &
                      agcacc, bgcjacc, bgcmacc, &
                      cgrain, egrain(N), egrain(P), egrain(S), &
                      crmvst, ermvst(N), ermvst(P), ermvst(S), &
                      cstraw, estraw(N), estraw(P), estraw(S), &
                      stdstraw, estdstraw(N), estdstraw(P), estdstraw(S), &
                      addsdc, addsde(N), addsde(P), addsde(S), &
                      resid, reside(N), reside(P), reside(S), &
                      hrvirr, hrvfert(N), hrvfert(P), hrvfert(S), &       !  irrapp, fertapp(N), fertapp(P), fertapp(S),
                              hrvafert(N), hrvafert(P), hrvafert(S), &    !  afertapp(N), afertapp(P), afertapp(S),
                      hrvomadc, hrvomade(N), hrvomade(P), hrvomade(S), &  !  omadapp, omaeapp(N), omaeapp(P), omaeapp(S),
                      strmac(1), strmac(2), strmac(3), strmac(4), &
                      strmac(5), strmac(6), strmac(7), strmac(8), &
                      cgracc, egracc(N), egracc(P), egracc(S), &
                      accrst, accrste(N), accrste(P), accrste(S), &
                      ctubesj, etubesj(N), etubesj(P), etubesj(S), &
                      ctubesm, etubesm(N), etubesm(P), etubesm(S), &
                      srfclittrj, esrfclittrj(N), esrfclittrj(P), esrfclittrj(S), &
                      soillittrj, esoillittrj(N), esoillittrj(P), esoillittrj(S), &
                      srfclittrm, esrfclittrm(N), esrfclittrm(P), esrfclittrm(S), &
                      soillittrm, esoillittrm(N), esoillittrm(P), esoillittrm(S))
      ! clear the harvest specific accumulators
      ! NOTE: should't these be in the plot
        hrvirr   = 0.0 ! clear harvest irrigation
        hrvfert  = 0.0 ! clear harvest regular fert
        hrvafert = 0.0 ! clear harvest auto fert
        hrvomadc = 0.0 ! clear harvest omad c
        hrvomade = 0.0 ! clear harvest omad minerals

! ... Update state variables
      call flowup(time)

! ... Check status of live and standing dead to eliminate any roundoff residuals.
      ! Simplified clears using single if for each pool and FORTRAN arrays.
      ! Added the cleared residuals to sink.
      ! Changed the threshold to a fraction of pre flow totals; this works
      ! because flow update the components but the sums updated are updates
      ! in the later call to sumcar.     KLK Nov 2012

      if (sum(aglcis) .lt. aglivc*1.e-05) then
        csrsnk = csrsnk + sum(aglcis); ! add roundoff C to sink
        aglcis = 0.0                   ! clear array
        esrsnk = esrsnk + aglive;      ! add roundoff E to sink
        aglive = 0.0                   ! clear array
      endif

      ! combined root check since they shouldn't have been removed separately
      if (sum(bglcisj)+sum(bglcism) .lt. (bglivcj+bglivcm)*1.e-05) then
        csrsnk = csrsnk + (sum(bglcisj)+sum(bglcism)); ! add roundoff C to sink
        bglcisj = 0.0                                  ! clear array
        bglcism = 0.0                                  ! clear array
        esrsnk = esrsnk + (bglivej + bglivem + crpstg);! add roundoff E to sink
        bglivej = 0.0                                  ! clear array
        bglivem = 0.0                                  ! clear array
        crpstg = 0.0                                   ! clear array
      endif

      if (stdcis(UNLABL)+stdcis(LABELD) .lt. 1.e-05)then
        csrsnk = csrsnk + sum(stdcis); ! add roundoff C to sink
        stdcis = 0.0                   ! clear array
        esrsnk = esrsnk + stdede;      ! add roundoff E to sink
        stdede = 0.0                   ! clear array
      endif

      call sumcar ! moved sumcar until after we are through mucking with the components

      return
      end
