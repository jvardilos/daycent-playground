
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... WDEATH

      subroutine wdeath (tavewk, bgwfunc, tfrac, avgstemp)
      use calflow;

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'forrem.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'zztim.inc'

! ... Argument declarations
      real tavewk, bgwfunc, tfrac, avgstemp

! ... Death of leaves, fine branches, large wood, fine roots, and coarse roots.
! ... Modifications:
! ...
! ... Corrected a bug in the death of lableled fine branches, coarse wood,
! ... and coarse roots. Use the actual labled fraction instead of the growth
! ... ratio, cisotf. Prevents the removal of to much labeled material when
! ... non equilibrium, under ratio, wood dies.  7/2007  K. Killian

! ... Function declarations
      real      gpdf, maxswpot
      external  gpdf, maxswpot

! ... Local variables
      integer iel
      integer, save :: drpdys =0
      logical, save :: drpdlv = .FALSE.
      real    accum(ISOS), ctodie, etodie, fr14, recres(MAXIEL)
      real    tostore, srfclittr, soillittr
      real    rtdh, agetempeff, tempeff, watreff, cturn, eturn, temp
      real, parameter :: Clowlim = 0.00001 ! fparts 0 C level to prevent underflow

      recres = 0.0
      accum  = 0.0

! ... Death of leaves
! ... NOTE:  WOODDR(1)   - the death rate in fall for deciduous forests
! ...        LEAFDR(MTH) - the monthly death rate for leaves in every
! ...                      case except for fall in deciduous forests.
      if (rleavc .ne. 0.0) then
        if (decid .ge. 1) then

! ....... Deciduous forest
! ....... If the daylight hours are increasing - it must be spring
          if (hrsinc .and. (tavewk .gt. tmplff)) drpdlv = .FALSE.

! ....... If daylight hours are decreasing and the temperature is low
! ....... enough drop leaves for fall season.
! ....... Add check for number of daylight hours to conditional for
! ....... determining if leaf drop should occur, cak - 06/30/03
! ....... If leaf drop has not occurred by the time the winter solstice
! ....... is reached force leaf drop to occur, cak - 10/28/04
          if (decid .eq. 1) then
            if (.not. drpdlv .and.                                       &
     &          (((tavewk .lt. tmplff) .and. (.not. hrsinc) .and.        &
     &            (dayhrs .lt. 12.0)) .or.                               &
     &          ((sitlat .ge. 0) .and. (month .eq. 12)) .or.             &
     &          ((sitlat .lt. 0) .and. (month .eq. 6)))) then
! ........... Allow leaf drop to occur over a 30 day period, dropping
! ........... all of the remaining leaves on the last day.
              if (drpdys .lt. 30) then
                ctodie = rleavc * wooddr(LEAF)*tfrac
                drpdys = drpdys + 1
                decidgrow = .FALSE.
               else
                ctodie = rleavc * wooddr(LEAF)
                drpdlv = .TRUE.
                drpdys = 0
              endif
            else
              ctodie = rleavc * leafdr(month)*tfrac
            endif
          elseif (decid .eq. 2) then
            ! Compute death for drought deciduous forests
            ctodie = rleavc * (1. - bgwfunc) * wooddr(LEAF)*tfrac
          endif
        else
! ....... Continuous forest; Use leaf death rate multiplier from EACHYR
          ctodie = rleavc * leafdr(month)*tfrac * ldrmlt
        endif

        if (ctodie .ne. 0. .and. rleavc .lt. Clowlim) ctodie = rleavc ! zero a negligible pool

! ..... Compute E/C ratios
        do iel = 1, nelem
          recres(iel) = rleave(iel) / rleavc

! ....... Compute flow to retranslocation storage
          tostore = recres(iel) * ctodie * forrtf(iel)
          call flow(rleave(iel), forstg(iel), time, tostore)

! ....... Decrease E/C by the amount that is retranslocated
          recres(iel) = recres(iel) * (1 - forrtf(iel))
        end do

! ..... If evntyp is greater than 1 the leaves go to the source/sink
! ..... rather than to the litter, cak - 02/07/2006
        fr14 = rlvcis(LABELD) / rleavc
        if (evntyp .lt. 2) then
          call partit(ctodie, recres, 1, rlvcis, rleave, wdlig(LEAF), fr14)
        else
          call csched(ctodie, fr14, 1.0,               &
                      rlvcis(UNLABL), csrsnk(UNLABL),  &
                      rlvcis(LABELD), csrsnk(LABELD),  &
                      1.0, accum)
          do iel = 1, nelem
            etodie = ctodie * (rleave(iel) / rleavc)
            call flow(rleave(iel), esrsnk(iel), time, etodie - tostore)
          enddo
        endif
      endif

! ... Add code to age fine roots, juvenile fine roots age to the mature
! ... fine root pool.  Modify this subroutine so that the death rate of
! ... roots is a function of soil water potential and soil temperature,
! ... cak - 06/28/2007
      ! See: 'A Model of Production and Turnover of Roots in Shortgrass Prairie
      !       Parton, Singh, and Coleman/, 1978, Journal of Applied Ecology
      temp = avgstemp
      ! Cap the temperature effect on fine roots at -2 and +28 degrees C
      if (temp .gt. 28.0) then
        temp = 28.0
      else if (temp .lt. -2.0) then
        temp = -2.0
      endif

! ... Soil temperature effect on root death rate
      tempeff = min((temp - 10.0)**2 / 4.0 * 0.00175 + 0.1, 0.5)
! ... Soil water potential effect on root death rate
      watreff = maxswpot(tlaypg)
      watreff = carctanf(watreff, 35.0, 0.5, 1.0, 0.05)
      rtdh = max(tempeff, watreff)

! ... removal of juvenile fine roots
      ! combine death and maturation into a single loop
      if (frootcj .ne. 0.0) then
        agetempeff = gpdf(temp, 37.0, 0.0, 3.0, 3.0) ! Soil temperature effect on juvenile root aging

        recres(1:nelem) = frootej(1:nelem) / frootcj ! local copy of frootej
        fr14 = frtcisj(LABELD) / frootcj             ! E/C ratio

        ! to prevent double booking root C, give maturing roots a 1 step reprieve from death
        ctodie = frootcj * wooddr(FROOTJ) * rtdh * tfrac
        if (ctodie .ne. 0. .and.  frootcj .lt. Clowlim) then
          ctodie = frootcj   ! zero a negligible pool
          cturn  = 0.        ! dead roots don't mature
          soillittr = ctodie ! send the residual to soil to eliminate round-off error
          srfclittr = 0.
        else
          cturn = (frootcj - ctodie) * tmxturn * agetempeff * tfrac
          ! Schedule juvenile root death
          ! A fraction of the dead roots are transferred to surface litter
          ! the remainder goes to soil litter cak - 05/14/2007
          srfclittr = ctodie * wrdsrfc
          soillittr = ctodie - srfclittr
        endif

        if(srfclittr > 0.) call partit(srfclittr, recres, SRFC, frtcisj,  &
                           frootej, wdlig(FROOTJ), fr14)
        call partit(soillittr, recres, SOIL, frtcisj, frootej, wdlig(FROOTJ), fr14)

        ! Schedule the root maturation
        call csched(cturn, fr14, 1.0,                  &
                    frtcisj(UNLABL), frtcism(UNLABL),  &
                    frtcisj(LABELD), frtcism(LABELD),  &
                    1.0, accum)
        do iel = 1, nelem
          eturn = cturn * (frootej(iel) / frootcj)
          call flow(frootej(iel), frootem(iel), time, eturn)
        end do
      endif

! ... Death of mature fine roots
      ! don't need to worry about maturing roots; they show up in the next time step
      if (frootcm .ne. 0.0) then
        ctodie = frootcm * wooddr(FROOTM) * rtdh * tfrac
        if (ctodie .ne. 0. .and. frootcm .lt. Clowlim) ctodie = frootcm ! zero a negligible pool
        recres(1:nelem) = frootem(1:nelem) / frootcm
        fr14 = frtcism(LABELD) / frootcm
        ! A fraction of the dead roots are transferred to the surface litter
        ! layer, the remainder goes to the soil litter layer cak - 05/14/2007
        srfclittr = ctodie * wrdsrfc
        soillittr = ctodie - srfclittr
        call partit(srfclittr, recres, SRFC, frtcism, frootem, wdlig(FROOTM), fr14)
        call partit(soillittr, recres, SOIL, frtcism, frootem, wdlig(FROOTM), fr14)
      endif

! ... Fine Branches, Large Wood, and Coarse Roots go to the dead wood
! ... compartments: WOOD1, WOOD2, WOOD3

! ... Death of fine branches
      if (fbrchc .gt. 0.0) then
        ctodie = fbrchc * wooddr(FBRCH)*tfrac
        if (ctodie .ne. 0. .and. fbrchc .lt. Clowlim) ctodie = fbrchc ! zero a negligible pool
        ! remove labled fraction instead of growth ratio, cisotf.  KLK 7/2007
        call csched(ctodie, fbrcis(LABELD), fbrchc,  &
                    fbrcis(UNLABL), wd1cis(UNLABL),  &
                    fbrcis(LABELD), wd1cis(LABELD),  &
                    1.0, accum)

        do iel = 1, nelem
          etodie = ctodie * (fbrche(iel) / fbrchc)
          call flow(fbrche(iel), wood1e(iel), time, etodie)
        end do
      endif

! ... Death of large wood
      if (rlwodc .gt. 0.0) then
        ctodie = rlwodc * wooddr(LWOOD)*tfrac
        if (ctodie .ne. 0. .and. rlwodc .lt. Clowlim) ctodie = rlwodc ! zero a negligible pool
        ! remove labled fraction instead of growth ratio, cisotf. KLK 7/2007
        call csched(ctodie, rlwcis(LABELD), rlwodc,  &
                    rlwcis(UNLABL), wd2cis(UNLABL),  &
                    rlwcis(LABELD), wd2cis(LABELD),  &
                    1.0, accum)

        do iel = 1, nelem
          etodie = ctodie * (rlwode(iel) / rlwodc)
          call flow(rlwode(iel), wood2e(iel), time, etodie)
        end do
      endif

! ... Death of coarse roots
      if (crootc .gt. 0.0) then
        ctodie = crootc * wooddr(CROOT)*tfrac
        if (ctodie .ne. 0. .and. crootc .lt. Clowlim) ctodie = crootc ! zero a negligible pool
        ! remove labled fraction instead of growth ratio, cisotf. KLK 7/2007
        call csched(ctodie, crtcis(LABELD), crootc,  &
                    crtcis(UNLABL), wd3cis(UNLABL),  &
                    crtcis(LABELD), wd3cis(LABELD),  &
                    1.0, accum)

        do iel = 1, nelem
          etodie = ctodie * (croote(iel) / crootc)
          call flow(croote(iel), wood3e(iel), time, etodie)
        end do
      endif


!...Death of fruits/nuts
      ! fraction of C to be removed
      ! the die (drop) fraction depends on whether we are growing or dropping
      if(dofngw) then  ! if(frnutc.gt.0.  .and.  fndday.lt.fngddl(1)) then
        ! this is the normal fall rate during growing season
        ctodie = frnutc * wooddr(FRNUT)

      elseif(fnftrm .gt. 0) then
       ! fall rate after growth; instead of constant fraction, the rate is
       ! inversely proportional to the fruit/nut fall time remaining.
       ! take it all if all the fruit can all die in this time step
        ctodie = min(frnutc/fnftrm, frnutc)
       ! Decrement the fall time remaining
        fnftrm = max(fnftrm -1, 0.)
      else
        ctodie = 0.0
      endif

      ! based on the leaf death code without the minerals movement to storage.
      if (ctodie .gt. 0.0) then

        ! Compute E/C ratios
        recres(1:nelem) = frnute(1:nelem) / frnutc

!       move the fruit into into structural and metabolic
        fr14 = frncis(LABELD) /frnutc
        call partit(ctodie, recres, 1,frncis,frnute,wdlig(FRNUT), fr14)
      endif

      return
      contains
       include 'catanf.f'
      end subroutine wdeath
