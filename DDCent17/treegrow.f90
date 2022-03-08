
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine treegrow (tfrac, tavedly)
      use calflow;

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'dynam.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'timvar.inc'
      include 'zztim.inc'

! ... Argument declarations
      real             tfrac, tavedly

! ... Simulate forest production associated with tree growth.
! ... This function was a part of TREES.F.
! ...   tfrac   - time step;  fraction of month (0-1)
! ...   tavedly - mean air temperature over production period (deg C)

      ! Remove avgstemp from the growth parameter list.
      ! plant root respiration rates are controlled by the aboveground photosynthesis
      ! so use air temperature to effect on root maintenance respiration CAK 31Aug2011

      ! Cindy Keough May 8, 2014,
      ! When the carbohydrate storage pool goes to zero set the labeled fraction for
      ! new material entering this pool for the current tree.

! ... Function declarations
      real      leafa, maxswpot, rtimp, lacalc
      external  leafa, maxswpot, rtimp, lacalc

! ... Local variables
      integer :: iel, ipart, iptr, lyr
      integer,save :: grwthyr = -16777215
      real    :: accum(ISOS), availm(MAXIEL), amt
      real    :: calcup, cprodfLeft
      real    :: euf(FPARTS)
      real    :: namt
      real    :: remCfrac, rimpct
      real    :: sum_cfrac, totCup
      real    :: uptake(4,MAXIEL)
      real    :: treeNfix
      real    :: cprodfdy, eprodfdy(MAXIEL)
      real    :: rleavc_opt, mrspReduce, fraclblstg, cstor
      real    :: totprod
      real    :: eup(FPARTS, MAXIEL)
      character (len=32) :: bufr
      character (len=10) :: subname = 'treegrow  '
      double precision frac_nh4, frac_no3
      real, parameter :: toler = 1.0E-30

      if (tfrac .lt. 0.0 .or. tfrac .gt. 1.0) then
        write(bufr,*) tfrac
        call abortrun('tfrac = '//trim(bufr)//' out of bounds in treegrow')
      endif
      ! if (tavedly .le. -999.0) then
      !   write(bufr,*) tavedly
      !   call abortrun('Error in treegrow, tavedly = '//trim(bufr)
      ! endif

      ! initialize arrays
      accum = 0.0
      uptake = 0.0
      cprodfdy = 0.0
      eprodfdy = 0.0

!     reset fruit growing degree day accumulator if growth cycle is complete
      if(dofngw) then
         fnftrm = 0 ! growing clears fall time remaining
      else  !     .not. dofngw
         fndday = 0  ! not growing; reset fruit GDD
!     start the post growing season fall after manually stopping fruit growth
!     Conditions: no fruit growth, hanging fruit, fruit fall rate >0, and fall
!     time remaining of zero. Since fruit fall SHOULD zero fruit C, this should
!     prevent extraneous fall conditions except possibly at run start
         if(frnutc .gt. 0. .and. fnftim .gt. 0. .and. fnftrm .eq. 0.) fnftrm = fnftim
      endif

      ! labeled fraction in the forest carbohydrate storage pool
      cstor = sum(carbostg(FORSYS,:)) ! intermediate value; total C in forest storage
      if (cstor .gt. 0.0) then
        fraclblstg = carbostg(FORSYS,LABELD) / cstor
      else
        fraclblstg = cisotf ! set to new material fraction
      endif

! ... Determine nutrients available to plants for growth.
      ! Nutrients available to trees are in the top tlaypg layers, cak 01/29/03
      ! Previous flowup, in trees, should have updated mineral pools.
      ! F95 note: The nested loops, looped vector sum, and full sum are equivalent
      availm(1:nelem) = sum(minerl(1:tlaypg,1:nelem), DIM= 1, &
                   MASK= minerl(1:tlaypg,1:nelem) .GT. toler)

! ... Determine old or new forest
       ! iptr is the forest carbon allocation index; new (iptr = 1) or mature (iptr = 2)
       ! Switch from new forest to old forest allocation fractions when frstage
       ! > swold where frstage is time - tree start  KLK 3/2004 - 3/2017
       ! this fixes a bug that makes the fractions dependent on absolute date
       ! NOTE: a float frstage is "poor practice" but integers <16,777,216 are OK
      if (frstage .le. swold) then  ! time since plant instead of absolute date
        iptr = 1  ! juvenile forest C allocation fractions
      else
        iptr = 2  ! mature forest C allocation fractions
      endif

! ... Temperature effect on maintenance respiration for aboveground components
      ! 0<= maintenance respiration temperature effect <= 1.0, cak - 09/16/02
      mrspTempEffect(FORSYS,SRFC) = max(0.0, min(1.0, 0.1 * exp(0.07 * tavedly)))

! ... Temperature effect on maintenance respiration for belowground components
      mrspTempEffect(FORSYS,SOIL) = max(0.0, min(1.0, 0.1 * exp(0.07 * tavedly)))

! ... Add a soil water term to fine root maintenance respiration, cak - 06/27/2007
      ! Calculate the soil water potential for the wettest soil layer in tree rooting zone
      mrspWaterEffect(FORSYS) = maxswpot(tlaypg)
      if (mrspWaterEffect(FORSYS) .le. 76.0) then
        mrspWaterEffect(FORSYS) = 1.0 - mrspWaterEffect(FORSYS) / 80.0
      else
        mrspWaterEffect(FORSYS) = 0.05
      endif

      ! Linearly decrease maintenance respiration as carbohydrate in the storage pool,
      ! based on the C required for optimal leaf area, gets smaller  cak - 10/20/2006
      ! Use two line functions for this calculation depending on the optimal leaf
      ! carbon value, cak - 08/13/2009
      ! Calculate C in the leaves as if the tree is at optimal LAI, cak - 10/20/2006
      rleavc_opt = lacalc(fbrchc, rlwodc, maxlai, klai) / (2.5 * btolai)
      if (rleavc_opt > 0.000000001) then
        ! limit mrspReduce below the minimum value
        if (cstor .le. fmrsplai(1) * rleavc_opt) then         ! minimum value
          mrspReduce = fmrsplai(2)

        elseif (cstor .ge. fmrsplai(5) * rleavc_opt) then     ! maximum value
          mrspReduce = fmrsplai(6)

        elseif (cstor .lt. fmrsplai(3) * rleavc_opt) then     ! first segment

          mrspReduce = line(cstor, fmrsplai(1) * rleavc_opt, fmrsplai(2), &
                                   fmrsplai(3) * rleavc_opt, fmrsplai(4))

        else                                                  ! second segment
          mrspReduce = line(cstor, fmrsplai(3) * rleavc_opt, fmrsplai(4), &
                                   fmrsplai(5) * rleavc_opt, fmrsplai(6))
        endif

! ...   Maintenance respiration flux calculation added, mdh - 9/4/01
        fmrspdyflux(LEAF) = fkmrspmx(LEAF) * &
                            mrspTempEffect(FORSYS,SRFC) * rleavc * &
                            tfrac * mrspReduce
        fmrspdyflux(FROOTJ) = fkmrspmx(FROOTJ) * &
                              mrspTempEffect(FORSYS,SOIL) * &
                              mrspWaterEffect(FORSYS) * frootcj * tfrac * &
                              mrspReduce
        fmrspdyflux(FROOTM) = fkmrspmx(FROOTM) * &
                              mrspTempEffect(FORSYS,SOIL) * &
                              mrspWaterEffect(FORSYS) * frootcm * tfrac * &
                              mrspReduce
        fmrspdyflux(FBRCH) = fkmrspmx(FBRCH) * &
                             mrspTempEffect(FORSYS,SRFC) * fbrchc * &
                             tfrac * mrspReduce
        fmrspdyflux(LWOOD) = fkmrspmx(LWOOD) * &
                             mrspTempEffect(FORSYS,SRFC) * rlwodc * &
                             tfrac * mrspReduce
        fmrspdyflux(CROOT) = fkmrspmx(CROOT) * &
                             mrspTempEffect(FORSYS,SOIL) * crootc * &
                             tfrac * mrspReduce
        ! added maintenance respiration for fruit/nuts similar to leaf
        fmrspdyflux(FRNUT) = fkmrspmx(FRNUT) * &
                             mrspTempEffect(FORSYS,SRFC) * frnutc * &
                             tfrac * mrspReduce

! ...   Maintenance respiration comes from carbohydrate storage, cak - 08/12/2009
        mrspdyflux(FORSYS) = min(sum(fmrspdyflux), cstor)
        call csched(mrspdyflux(FORSYS),fraclblstg,1.0, &
                    carbostg(FORSYS,UNLABL),csrsnk(UNLABL), &
                    carbostg(FORSYS,LABELD),csrsnk(LABELD), &
                    1.0,fautoresp)
      else
        fmrspdyflux = 0.0
      endif

! ... If growth can occur
      if (forgrw .eq. 1 .and. pforc .gt. 0.0) then
        if(grwthyr .ne. cyear) frstage = frstage +1; ! increment forest age with minimum round-off
        grwthyr = cyear

! ..... Calculate actual production values
! ..... Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, frootcj+frootcm)
! ..... Determine actual production values, restricting the C/E ratios
! ..... When calling restrp we are only looking at allocation to fine roots
! ..... and leaves, cak - 07/02/02

! ..... Determine actual production values, restricting the C/E ratios
        ! Added tree favail as a parameter to the restrp call   KLK 30Dec13
        ! Added eup as a parameter to the restrp call   KLK May 2015
        call restrp(elimit, nelem, availm, favail(1:3,FORSYS), ccefor, &
                    2, tree_cfrac, pforc, rimpct, forstg, snfxmx(FORSYS), &
                    cprodfdy, eprodfdy, uptake, eup, treeNfix, relyld)
      else
        cprodfdy = 0.0
        eprodfdy = 0.0
      endif

! ... If growth occurs...
      if (cprodfdy .gt. 0.) then

        ! If the carbohydrate storage pool falls below a critical value add carbon from
        !  csrsnk to allow plant growth.  This should only occur when the tree is small.
        if (carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD) .lt. 15.0) then
!          write(*,*) 'Warning, adding C to carbostg in treegrow time= ',doy,'/',cyear, &
!          ' carbostg = ', cstor
          carbostg(FORSYS,UNLABL) = carbostg(FORSYS,UNLABL) + (15.0 * (1.0 - fraclblstg))
          carbostg(FORSYS,LABELD) = carbostg(FORSYS,LABELD) + (15.0 * fraclblstg)
          csrsnk(UNLABL) = csrsnk(UNLABL) - (15.0 * (1.0 - fraclblstg))
          csrsnk(LABELD) = csrsnk(LABELD) - (15.0 * fraclblstg)
        endif

        ! Increment counter tracking the number of days to allow woody growth, cak - 10/20/2006
        fgrwdys = fgrwdys + 1

! ..... Compute carbon allocation fractions for each tree part
        ! C use by roots
        cprodfLeft = cprodfdy * (1.0 - tree_cfrac(FROOTJ))

        ! C allocation to leaves, up to a optimal LAI
        tree_cfrac(LEAF) = leafa(rleavc_opt, rleavc, cprodfLeft, cprodfdy)
        remCfrac = 1.0 - (tree_cfrac(FROOTJ) + tree_cfrac(LEAF))

        ! allocate C to fruit/nut growth if it is warm, we have fruit growth and the leaves aren't blooming
        if(tavedly .gt. 0.0  .and.  dofngw .and. remCfrac .gt. 0.0) then
          tree_cfrac(FRNUT) = remCfrac* fcfrac(FRNUT, iptr) ! set a non-zero fruit growth potential
          ! remCfrac = 1.0 - (tree_cfrac(FROOTJ) + tree_cfrac(LEAF) + tree_cfrac(FRNUT))
          remCfrac = remCfrac - tree_cfrac(FRNUT)
          ! add daily increment to fruit growing season
          fndday   = fndday + max(min(tavedly,fngddl(3))-fngddl(2), 0.0)

          ! growth exceeds limit Stop fruit growth; force harvest if it is flagged
          if (fndday .ge. fngddl(1) .and. dofngw) then
            dofngw = .false.   ! no more fruit growth after today
            fndday = fngddl(1) ! fruit Growing degree days equals limit
            fnftrm = fnftim    ! initialize the time for fruit fall
          endif
        endif

! ..... If we have leftover carbon allocate it to the woody plant parts
! ..... using a weighted average
        if (remCfrac .lt. 1.0E-05) then
          tree_cfrac(FBRCH:CROOT) = 0.0 ! zero allocation to other parts FBRCH, LWOOD, and CROOT ...
        else
! ....... for FBRCH, LWOOD, and CROOT ...
          totCup = sum(fcfrac(FBRCH:CROOT, iptr))
          if (totCup .gt. 0.0) then
! ....... for FBRCH, LWOOD, and CROOT ...
            tree_cfrac(FBRCH:CROOT) = fcfrac(FBRCH:CROOT, iptr) * (remCfrac / totCup)
          else
            write(bufr,*) totCup
            call abortrun('wood growth fraction fcfrac(FBRCH)+fcfrac(LWOOD)+fcfrac(CROOT)'// &
                          trim(bufr)//' <= 0 in treegrow')
          endif
        endif

! ... Error checking
        sum_cfrac = 0.0
        do ipart = 1, FPARTS-1
          if (tree_cfrac(ipart) .lt. 0.0) then
            write(bufr,'(i2)') ipart;
            call abortrun('tree_cfrac('//bufr(:2)//') < 0 in treegrow')
          else if (tree_cfrac(ipart) .gt. 1.0) then
            write(bufr,'(i2)') ipart;
            call abortrun('tree_cfrac('//bufr(:2)//') > 1 in treegrow')
          else
            sum_cfrac = sum_cfrac + tree_cfrac(ipart)
          endif
        end do
        if (abs(1.0 - sum_cfrac) .gt. 0.001) then
          write(bufr,*) sum_cfrac
          call message('Error in treegrow, unnormalized carbon allocation fractions '//trim(bufr))
        endif

! ..... Recalculate actual production values with updated C-allocation
! ..... fractions, restricting the C/E ratios -mdh 5/11/01
! ..... Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, frootcj+frootcm)

! ..... Determine actual production values, restricting the C/E ratios
        ! Added tree favail as a parameter to the restrp call   KLK 30Dec13
        call restrp(elimit, nelem, availm, favail(1:3,FORSYS), ccefor, &
                    FPARTS-1,tree_cfrac, pforc, rimpct, forstg, &
                    snfxmx(FORSYS), cprodfdy, eprodfdy, uptake, eup, &
                    treeNfix, relyld)

! ..... Calculations for symbiotic N fixation accumulators moved from
! ..... nutrlm subroutine, cak - 10/17/02
! ..... Compute N fixation which actually occurs and add to the
! ..... N fixation accumulator.
        nfix = nfix + treeNfix
        snfxac(FORSYS) = snfxac(FORSYS) + treeNfix
! ..... Add computation for nfixac -mdh 1/16/02
        nfixac = nfixac + treeNfix

! ..... C/N ratio for production
        if (eprodfdy(N) .eq. 0.0) call abortrun('eprodfdy(N) = 0.0 in treegrow')
        tcnpro = cprodfdy/eprodfdy(N)

! ..... Calculate production for each tree part
! ..... New variable MFPRD added for gridded output - 6/96 rm
        do ipart = 1, FPARTS
          if (ipart .eq. FROOTJ) then
            mfprd(ipart) = tree_cfrac(FROOTJ) * cprodfdy * (1.0 - wmrtfrac)
          else if (ipart .eq. FROOTM) then
            mfprd(ipart) = tree_cfrac(FROOTJ) * cprodfdy * wmrtfrac
          else
            mfprd(ipart) = tree_cfrac(ipart) * cprodfdy
          endif
        end do

! ..... Use NPP as an estimate of GPP to produce carbon that is stored
! ..... in the carbohydrate storage pool
        totprod = sum(mfprd(1:FPARTS))
        if (totprod .gt. 0.0) then
          totprod = totprod * npp2cs(FORSYS)
          call csched(totprod,cisotf,1.0, &
                      csrsnk(UNLABL),carbostg(FORSYS,UNLABL), &
                      csrsnk(LABELD),carbostg(FORSYS,LABELD), &
                      1.0,accum)
        endif

! ..... Forest Growth
        ! All growth comes from the carbohydrate pool, cak - 08/12/2009
        ! Growth of leaves split into labeled & unlabeled parts
        call csched(mfprd(LEAF),fraclblstg,1.0, &
                    carbostg(FORSYS,UNLABL),rlvcis(UNLABL), &
                    carbostg(FORSYS,LABELD),rlvcis(LABELD), &
                    1.0,alvcis)
        ! Growth juvenile fine roots; split into labeled & unlabeled parts
        call csched(mfprd(FROOTJ),fraclblstg,1.0, &
                    carbostg(FORSYS,UNLABL),frtcisj(UNLABL), &
                    carbostg(FORSYS,LABELD),frtcisj(LABELD), &
                    1.0,afrcisj)
        ! Growth mature fine roots; split into labeled & unlabeled parts
        call csched(mfprd(FROOTM),fraclblstg,1.0, &
                    carbostg(FORSYS,UNLABL),frtcism(UNLABL), &
                    carbostg(FORSYS,LABELD),frtcism(LABELD), &
                    1.0,afrcism)
        ! Growth fruit/nut; split into labeled & unlabeled parts
        call csched(mfprd(FRNUT),fraclblstg,1.0, &
                    carbostg(FORSYS,UNLABL),frncis(UNLABL), &
                    carbostg(FORSYS,LABELD),frncis(LABELD), &
                    1.0,afncis)
        ! Growth respiration
        fgrspdyflux(LEAF) = mfprd(LEAF) * fgresp(LEAF)
        fgrspdyflux(FROOTJ) = mfprd(FROOTJ) * fgresp(FROOTJ)
        fgrspdyflux(FROOTM) = mfprd(FROOTM) * fgresp(FROOTM)
        fgrspdyflux(FRNUT) = mfprd(FRNUT) * fgresp(FRNUT)
! ..... If the tree is decidious and the time alloted to grow the
! ..... woody components of the trees has passed use the late season
! ..... growth restriction parameter value to determine how much
! ..... carbohydrate to flow out of the forest carbohydrate storage
! ..... pool for the woody components, cak - 03/11/2010
        if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys)) then
          ! Growth of fine branches; split into labeled & unlabeled parts
          call csched(mfprd(FBRCH)*(1.0-flsgres),fraclblstg,1.0, &
                      carbostg(FORSYS,UNLABL),fbrcis(UNLABL), &
                      carbostg(FORSYS,LABELD),fbrcis(LABELD), 1.0,afbcis)
          ! Growth of large wood; split into labeled & unlabeled parts
          call csched(mfprd(LWOOD)*(1.0-flsgres),fraclblstg,1.0, &
                      carbostg(FORSYS,UNLABL),rlwcis(UNLABL), &
                      carbostg(FORSYS,LABELD),rlwcis(LABELD), 1.0,alwcis)
          ! Growth of coarse roots; split into labeled & unlabeled parts
          call csched(mfprd(CROOT)*(1.0-flsgres),fraclblstg,1.0, &
                      carbostg(FORSYS,UNLABL),crtcis(UNLABL), &
                      carbostg(FORSYS,LABELD),crtcis(LABELD), 1.0,acrcis)
          ! Growth respiration
          fgrspdyflux(LWOOD) = mfprd(LWOOD) * (1.0 - flsgres) * fgresp(LWOOD)
          fgrspdyflux(FBRCH) = mfprd(FBRCH) * (1.0 - flsgres) * fgresp(FBRCH)
          fgrspdyflux(CROOT) = mfprd(CROOT) * (1.0 - flsgres) * fgresp(CROOT)
        else
          ! Growth of fine branches; split into labeled & unlabeled parts
          call csched(mfprd(FBRCH),fraclblstg,1.0, &
                      carbostg(FORSYS,UNLABL),fbrcis(UNLABL), &
                      carbostg(FORSYS,LABELD),fbrcis(LABELD), 1.0,afbcis)
          ! Growth of large wood; split into labeled & unlabeled parts
          call csched(mfprd(LWOOD),fraclblstg,1.0, &
                      carbostg(FORSYS,UNLABL),rlwcis(UNLABL), &
                      carbostg(FORSYS,LABELD),rlwcis(LABELD), 1.0,alwcis)
          ! Growth of coarse roots; split into labeled & unlabeled parts
          call csched(mfprd(CROOT),fraclblstg,1.0, &
                      carbostg(FORSYS,UNLABL),crtcis(UNLABL), &
                      carbostg(FORSYS,LABELD),crtcis(LABELD), 1.0,acrcis)
          ! Growth respiration
          fgrspdyflux(FBRCH) = mfprd(FBRCH) * fgresp(FBRCH)
          fgrspdyflux(LWOOD) = mfprd(LWOOD) * fgresp(LWOOD)
          fgrspdyflux(CROOT) = mfprd(CROOT) * fgresp(CROOT)
        endif

! ..... Growth respiration is subtracted from the carbohydrate storage pool.
        grspdyflux(FORSYS) = sum(fgrspdyflux(1:FPARTS))
        call csched(grspdyflux(FORSYS),fraclblstg,1.0, &
                    carbostg(FORSYS,UNLABL),csrsnk(UNLABL), &
                    carbostg(FORSYS,LABELD),csrsnk(LABELD), &
                    1.0,fautoresp)

! ..... Actual Uptake
        do iel = 1, nelem
          ! remove any soil uptake (NOT EFERT, that is autofert) from net fertilizer budget.
          fertnet(iel)  = max(fertnet(iel) - uptake(ESOIL,iel), 0.0)

          if (eprodfdy(iel) .eq. 0.0) call abortrun('Divide by zero in treegrow, eprodfdy(iel) = 0')
          ! split eup(froot)  into mature and juvenile fractions
          amt             = eup(FROOT,iel)
          eup(FROOTM,iel) = amt * wmrtfrac;
          eup(FROOTJ,iel) = amt - eup(FROOTM,iel)
          euf = eup(:,iel) / eprodfdy(iel)

! ....... Take up nutrients from internal storage pool
! ....... Don't allow uptake from storage if forstg is negative -mdh 8/8/00
! ....... If the tree is decidious and the time alloted to grow the
! ....... woody components of the trees has passed use the late season
! ....... growth restriction parameter value to determine how much
! ....... nutrients to flow out of the forest nutrient storage pool for
! ....... the woody components, cak - 03/11/2010
          if (forstg(iel) .gt. 0.0) then
            amt = uptake(ESTOR,iel) * euf(LEAF)
            call flow(forstg(iel),rleave(iel),time,amt)
            eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
            eupacc(iel) = eupacc(iel) + amt

            amt = uptake(ESTOR,iel) * euf(FROOTJ)
            call flow(forstg(iel),frootej(iel),time,amt)
            eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt

            amt = uptake(ESTOR,iel) * euf(FROOTM)
            call flow(forstg(iel),frootem(iel),time,amt)
            eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt

            amt = uptake(ESTOR,iel) * euf(FRNUT)
            call flow(forstg(iel),frnute(iel),time,amt)
            eupprt(FRNUT,iel) = eupprt(FRNUT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
            if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys)) then
              amt = uptake(ESTOR,iel) * euf(FBRCH) * (1.0 - flsgres)
              call flow(forstg(iel),fbrche(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt

              amt = uptake(ESTOR,iel) * euf(LWOOD) * (1.0 - flsgres)
              call flow(forstg(iel),rlwode(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt

              amt = uptake(ESTOR,iel) * euf(CROOT) * (1.0 - flsgres)
              call flow(forstg(iel),croote(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
            else
              amt = uptake(ESTOR,iel) * euf(FBRCH)
              call flow(forstg(iel),fbrche(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              amt = uptake(ESTOR,iel) * euf(LWOOD)
              call flow(forstg(iel),rlwode(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              amt = uptake(ESTOR,iel) * euf(CROOT)
              call flow(forstg(iel),croote(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
            endif
          endif

! ....... Take up nutrients from soil
! ....... Nutrients for uptake are available in the top tlaypg layers,
! ....... cak 01/29/03
! ....... If the tree is decidious and the time alloted to grow the
! ....... woody components of the trees has passed flow nutrients to
! ....... the forest storage pool rather than to the component nutrient
! ....... pools based on the forest late season growth restriction
! ....... parameter value, cak - 03/11/2010
          do lyr = 1, tlaypg
            if (minerl(lyr,iel) .gt. toler) then
              call cmpnfrac(lyr,ammonium,nitrate,minerl,frac_nh4,frac_no3)
              ! The P fsol isn't needed here for the weighted average, cak - 04/05/02
!              fsol = 1.0; if (iel .eq. P) fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
              calcup = uptake(ESOIL,iel)*minerl(lyr,iel)/availm(iel) ! * fsol

              ! Leaves
              amt = calcup * euf(LEAF)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                  ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),rleave(iel),time,amt)
              eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              ! Juvenile fine roots
              amt = calcup * euf(FROOTJ)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                  ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),frootej(iel),time,amt)
              eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              ! Mature fine roots
              amt = calcup * euf(FROOTM)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                  ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),frootem(iel),time,amt)
              eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              ! Fruit Nut
              amt = calcup * euf(FRNUT)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                  ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),frnute(iel),time,amt)
              eupprt(FRNUT,iel) = eupprt(FRNUT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              ! Fine branch
              namt = 0.0
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = calcup * euf(FBRCH) * flsgres
                namt = namt + amt
                call flow(minerl(lyr,iel),forstg(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = calcup * euf(FBRCH) * (1.0 - flsgres)
                namt = namt + amt
                call flow(minerl(lyr,iel),fbrche(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = calcup * euf(FBRCH)
                namt = namt + amt
                call flow(minerl(lyr,iel),fbrche(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
              if (iel .eq. N) then
                namt = -1.0*namt
                call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                  ammonium, nitrate, subname)
              endif
              ! Large wood
              namt = 0.0
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = calcup * euf(LWOOD) * flsgres
                namt = namt + amt
                call flow(minerl(lyr,iel),forstg(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = calcup * euf(LWOOD) * (1.0 - flsgres)
                namt = namt + amt
                call flow(minerl(lyr,iel),rlwode(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = calcup * euf(LWOOD)
                namt = namt + amt
                call flow(minerl(lyr,iel),rlwode(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
              if (iel .eq. N) then
                namt = -1.0*namt
                call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                  ammonium, nitrate, subname)
              endif
              ! Coarse roots
              namt = 0.0
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = calcup * euf(CROOT) * flsgres
                namt = namt + amt
                call flow(minerl(lyr,iel),forstg(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = calcup * euf(CROOT) * (1.0 - flsgres)
                namt = namt + amt
                call flow(minerl(lyr,iel),croote(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = calcup * euf(CROOT)
                namt = namt + amt
                call flow(minerl(lyr,iel),croote(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
              if (iel .eq. N) then
                namt = -1.0*namt
                call update_npool(lyr, namt, frac_nh4, frac_no3, &
                                  ammonium, nitrate, subname)
              endif
            endif
          end do

! ....... Take up nutrients from nitrogen fixation
! ....... If the tree is decidious and the time alloted to grow the
! ....... woody components of the trees has passed flow nutrients to
! ....... the forest storage pool rather than to the component nutrient
! ....... pools based on the forest late season growth restriction
! ....... parameter value, cak - 03/11/2010
          if (iel .eq. N .and. treeNfix .gt. 0) then
! ......... Leaves
            amt = uptake(ENFIX,iel) * euf(LEAF)
            call flow(esrsnk(iel),rleave(iel),time,amt)
            eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
! ......... Juvenile fine roots
            amt = uptake(ENFIX,iel) * euf(FROOTJ)
            call flow(esrsnk(iel),frootej(iel),time,amt)
            eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
! ......... Mature fine roots
            amt = uptake(ENFIX,iel) * euf(FROOTM)
            call flow(esrsnk(iel),frootem(iel),time,amt)
            eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
! ......... fruit Nut
            amt = uptake(ENFIX,iel) * euf(FRNUT)
            call flow(esrsnk(iel),frnute(iel),time,amt)
            eupprt(FRNUT,iel) = eupprt(FRNUT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
! ......... Fine branch
            if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
              amt = uptake(ENFIX,iel) * euf(FBRCH) * flsgres
              call flow(esrsnk(iel),forstg(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              amt = uptake(ENFIX,iel) * euf(FBRCH) * (1.0 - flsgres)
              call flow(esrsnk(iel),fbrche(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(FBRCH)
              call flow(esrsnk(iel),fbrche(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
            endif
! ......... Large wood
            if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
              amt = uptake(ENFIX,iel) * euf(LWOOD) * flsgres
              call flow(esrsnk(iel),forstg(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              amt = uptake(ENFIX,iel) * euf(LWOOD) * (1.0 - flsgres)
              call flow(esrsnk(iel),rlwode(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(LWOOD)
              call flow(esrsnk(iel),rlwode(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
            endif
! ......... Coarse roots
            if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
              amt = uptake(ENFIX,iel) * euf(CROOT) * flsgres
              call flow(esrsnk(iel),forstg(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              amt = uptake(ENFIX,iel) * euf(CROOT) * (1.0 - flsgres)
              call flow(esrsnk(iel),croote(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(CROOT)
              call flow(esrsnk(iel),croote(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
            endif
          endif

! ....... Take up nutrients from automatic fertilizer
! ....... If the tree is decidious and the time alloted to grow the
! ....... woody components of the trees has passed flow nutrients to
! ....... the forest storage pool rather than to the component nutrient
! ....... pools based on the forest late season growth restriction
! ....... parameter value, cak - 03/11/2010
          if (aufert .ne. 0) then
            if (uptake(EFERT,iel) .gt. 0.) then
              ! Automatic fertilizer added to plant pools

              ! Moved the automatic fertilizer addition to mineral pool, update_npool
              ! and accumulators to simsom so one component of an orchard, fertilized
              ! savanna, can't strip the other's minerl addition.   KLK 3 Mar 2014
              afrtup(iel) = afrtup(iel) + uptake(EFERT,iel)

              ! Leaves
              amt = uptake(EFERT,iel) * euf(LEAF)
              call flow(esrsnk(iel),rleave(iel),time,amt)
              eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              ! Juvenile fine roots
              amt = uptake(EFERT,iel) * euf(FROOTJ)
              call flow(esrsnk(iel),frootej(iel),time,amt)
              eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              ! Mature fine roots
              amt = uptake(EFERT,iel) * euf(FROOTM)
              call flow(esrsnk(iel),frootem(iel),time,amt)
              eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt

              ! Fruit/Nut
              amt = uptake(EFERT,iel) * euf(FRNUT)
              call flow(esrsnk(iel),frnute(iel),time,amt)
              eupprt(FRNUT,iel) = eupprt(FRNUT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt

              ! Fine branch
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = uptake(EFERT,iel) * euf(FBRCH) * flsgres
                call flow(esrsnk(iel),forstg(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(FBRCH) * (1.0 - flsgres)
                call flow(esrsnk(iel),fbrche(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(FBRCH)
                call flow(esrsnk(iel),fbrche(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
              ! Large wood
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = uptake(EFERT,iel) * euf(LWOOD) * flsgres
                call flow(esrsnk(iel),forstg(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(LWOOD) * (1.0 - flsgres)
                call flow(esrsnk(iel),rlwode(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(LWOOD)
                call flow(esrsnk(iel),rlwode(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
              ! Coarse roots
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = uptake(EFERT,iel) * euf(CROOT) * flsgres
                call flow(esrsnk(iel),forstg(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                amt = uptake(EFERT,iel) * euf(CROOT) * (1.0 - flsgres)
                call flow(esrsnk(iel),croote(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(CROOT)
                call flow(esrsnk(iel),croote(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
              endif
              ! Automatic fertilizer to be added to mineral pool
              afrtup(iel) = afrtup(iel) + uptake(EFERT,iel)
            endif
          endif
        end do

! ... Else there is no production this time step due to nutrient limitation
      else
        cprodfdy = 0.0 ! vector set
        eprodfdy = 0.0 ! vector set
      endif

! ... Accumulate monthly output
      mrspflux(FORSYS) = mrspflux(FORSYS) + mrspdyflux(FORSYS)
      grspflux(FORSYS) = grspflux(FORSYS) + grspdyflux(FORSYS)

      amt       = sum(fmrspdyflux(1:FPARTS))
      sumrsp    = sumrsp + amt
      mrspmth(FORSYS) = mrspmth(FORSYS) + amt
      mrspann(FORSYS) = mrspann(FORSYS) + amt ! annual output
      fmrspflux       = fmrspflux + fmrspdyflux

      amt       = sum(fgrspdyflux(1:FPARTS))
      grspmth(FORSYS) = grspmth(FORSYS) + amt
      grspann(FORSYS) = grspann(FORSYS) + amt ! annual output
      fgrspflux = fgrspflux + fgrspdyflux

! ... Accumulate monthly and annual output for soil respiration
      srspmth(FORSYS) = srspmth(FORSYS) + &
                        fmrspdyflux(FROOTJ) + fmrspdyflux(FROOTM) + fmrspdyflux(CROOT) + &
                        fgrspdyflux(FROOTJ) + fgrspdyflux(FROOTM) + fgrspdyflux(CROOT)
      srspann(FORSYS) = srspann(FORSYS) + &
                        fmrspdyflux(FROOTJ) + fmrspdyflux(FROOTM) + fmrspdyflux(CROOT) + &
                        fgrspdyflux(FROOTJ) + fgrspdyflux(FROOTM) + fgrspdyflux(CROOT)

! ... If the tree is decidious and the time alloted to grow the
! ... woody components of the trees has passed and the late season
! ... growth restriction parameter value is used to determine how much
! ... carbohydrate to flow out of the forest carbohydrate storage
! ... pool for the woody components reset the wood production
! ... components so the output accumulators are tracking correctly
      if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys)) then
        mfprd(FBRCH) = mfprd(FBRCH) * (1.0 - flsgres)
        mfprd(LWOOD) = mfprd(LWOOD) * (1.0 - flsgres)
        mfprd(CROOT) = mfprd(CROOT) * (1.0 - flsgres)
        cprodfdy = cprodfdy * (1.0 - flsgres)
        eprodfdy(1:nelem) = eprodfdy(1:nelem) * (1.0 - flsgres)
      endif

! ... Sum the daily production variables for output to the monthly *.bin file
      cprodf          = cprodf + cprodfdy
      eprodf(1:nelem) = eprodf(1:nelem) + eprodfdy(1:nelem)

      return
      contains
        include 'line.f'
      end subroutine treegrow
