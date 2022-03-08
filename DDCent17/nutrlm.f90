
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... NUTRLM.F

      subroutine nutrlm(elimit, nelem, nparts, cfrac, eavail, &
                        maxeci, mineci, snfxmx, eup, cprodl, eprodl, &
                        plantNfix, afert, aufert, eviFlag)

      implicit none
      include 'const.inc'
      include 'timvar.inc'

! ... Argument declarations
      logical eviFlag
      integer nelem, nparts
      real    cfrac(nparts), cprodl, elimit
      real    eprodl(MAXIEL), eavail(MAXIEL)
      real    maxeci(FPARTS-1,MAXIEL), mineci(FPARTS-1,MAXIEL)
      real    plantNfix, snfxmx, afert(MAXIEL), aufert
      real :: eup(nparts, MAXIEL)

      ! Nutrient limitation for plants is based on demand

      ! Changes:
      ! make autofert local, limit auto fertilization to 1.0 and don't set the
      !  persistent aufert with eviflag
      !  changed fixation check to better describe reduced fixation totale > demand
      !  restrict autofert to nutrient limited growth; prevents problems if aufert > 1
      !  fixed a bug that added N fixation to other minerals     KLK 23 Feb 2016
      ! Not sure if fixation really limits that fast but this may be just a matter of degree
      ! removed unused common block parfs    KLK May 2015
        ! added eup to the parameter list
      ! Correct an auto fertilization/N fixation bug    KLK 24 March 2015
        ! removed unneeded local variable maxNfix
        ! ensure fixation is not reduced if auto fertilized N is > 0. Corrects a bug that
        ! could convert fixation to fertilizer but seemed to double count fertilizer
        ! as fixation
      ! changed to f90 format KLK 22 Aug 2014
        ! converted several smaller loops to f90 array calculations
        ! move the demand calculation to nutrlm
        ! simplified the limiting element calculation and N fixation calculations
        ! removed a warning loop that changed eprodl but not the actual uptake eup
        ! changed fixation reported so it remains dependent on the potential growth
        !  (available C source) and not nutrient limited growth. It caused N balance
        !  errors and, if allowed to iterate, would suppress N fixation entirely

      ! Local Variables
      ! NOTE:  Local variables cannot have adjustable array size.  ECFOR
      !        is set to the largest array size which may occur.
      integer   iel
      real      demand, ecfor(FPARTS-1,MAXIEL)
      real      totale, ratio, minratio, aufertl
      character bfr*80

! ... Definitions of Local variables
! ...   demand  - E demand based on maximum E/C ratio
! ...   ecfor   - Actual E/C ratio by part limited by element
! ...   ratio   - ratio of available to potential demand with nutrient limitation
! ...   totale  - E available

      ! remove excess spaces from error messages KLK 2014
      if ((nelem .le. 0) .or. (nelem .gt. 3)) then
        write(bfr,*) nelem
        call abortrun('nelem = '//trim(bfr)//' out of bounds in nutrlm')
      endif
      if (nparts .gt. FPARTS-1) then
        write(bfr,*) nparts
        call abortrun('nparts = '//trim(bfr)//' out of bounds in nutrlm')
      endif

! ..... Automatic fertilizer explanation
        !  aufert 0.0 to 1.0
        !    automatic fertilizer may be applied to remove some nutrient
        !    stress so that relyld is at least the value of aufert
        !   EVI - there will be no nutrient limitation when running EVI
        !   BUG: if (eviFlag) aufert = 1.0
        !        eviFlag shouldn't set the persistant aufert since EVI doesn't clear it
        !        This will hold until a fert event regardless of whether evi is avalable.
        !        the alternative if (eviFlag) aufertl = 1.0 would switch
        !        with missing data which also may be problematic
      if (eviFlag) aufert = 1.0 ! (eviFlag  .and.  aufertl .eq. 0) aufertl = 1.0
      aufertl  = min(aufert, 1.0)
      elimit   = 0.0      ! limiting element
      minratio = 1.0      ! max ratio
      plantNfix  = max(snfxmx * cprodl, 0.)  ! N FIXATION
      eup      = 0.

! ... Compute production limitation
      do iel = 1, nelem

        ! DEMAND based on the maximum E/C ratio.
        ! load eup with the maximum E/C
        eup(1:nparts,iel)   = cfrac(1:nparts) * maxeci(1:nparts,iel)
        demand = sum(eup(1:nparts,iel)) * cprodl !  maxdemand(iel)
        if (demand .eq. 0.0) call abortrun('demand = 0.0 in nutrlm')

        totale = eavail(iel)
        if (iel .eq. N) totale = totale + plantNfix

! ..... New E/C ratios by part based on E available.
        if (totale .ge. demand) then
          ! Not nutrient limited; Use the max E/C ratio.
          ecfor(1:nparts,iel) = maxeci(1:nparts,iel)
          eprodl(iel) = demand       ! set nutrient production limit to maximum demand
          ratio = 1                  ! unlimited
        else
          ! New E/C ratios by part based on limited E available.
          ecfor(1:nparts,iel) = mineci(1:nparts,iel) + &
            (maxeci(1:nparts,iel) - mineci(1:nparts,iel)) * (totale / demand)

          ! Adjust ratios so N/P ratio of the leaves does not exceed the
          ! observed critical value of 13.5, cak - 04/05/02
          ! if (nelem .eq. P  .and.  ecfor(LEAF,N)/ecfor(LEAF,P) .gt. maxnp) &
          !          ecfor(LEAF,N) = ecfor(LEAF,P) * maxnp

          ! replace eup vector with the limited E/C
          eup(1:nparts,iel)   = cfrac(1:nparts) * ecfor(1:nparts,iel)

          ! Calculate the potential demand based on revised E/C ratios
          eprodl(iel) = sum(eup(1:nparts,iel)) * cprodl
          ratio = totale / eprodl(iel) ! limited growth ratio with reduced E/C

          ! Record the limiting element; this is still the limit even with aufert>0
          if (ratio .lt. minratio) then
            minratio = ratio
            elimit = real(iel)
          endif

          if ((aufertl .gt. 0.0) .and. (ratio .lt. aufertl)) then
            ! Fraction of E demand that needs to be added = aufertl - ratio(iel)
            afert(iel) = eprodl(iel) * (aufertl - ratio)
          endif
        endif
      end do

! ... Recompute EPRODL and fixation based on new production values
      if(minratio .lt. aufertl) minratio = aufertl
      cprodl                = cprodl * minratio
      eprodl(1:nelem)       = eprodl(1:nelem) * minratio
      eup(1:nparts,1:nelem) = eup(1:nparts,1:nelem) * cprodl

      ! limit fixation if there is more N than required for full growth
      ! QUESTION should the plant really strip ALL available N
      if(plantNfix .gt. 0  .and.  afert(N) .eq. 0.  .and.  eprodl(N) .lt. eavail(N)+plantNfix) then
         plantNfix = max(eprodl(N) -eavail(N), 0.0)
      end if

      return
      end
