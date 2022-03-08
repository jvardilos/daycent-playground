
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine treeDynC(pforc, tree_a2drat, tavemth, tree_cfrac)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot3.inc'

! ... Argument declarations
      real, intent(INOUT) :: pforc, tree_a2drat(MAXIEL)
      real, intent(in)    :: tavemth
      real, intent(out)   :: tree_cfrac(FPARTS-1)

      ! Set initial carbon allocation fractions for tree LEAF and FROOT
      ! The fine root fraction is calculated using available N, demand and leaf
      ! greenup. The remainder is given to the leaves.
      ! RBRCH, LWOOD, and CROOT, FRNUT are set to zero at this point.
      ! this FULLY initializes tree_cfrac

      ! For seasonal deciduous/conifer forest system, reapportion growth
      ! This applies only to non-drought systems.
      ! From Bill Parton (April 2000):
      ! "For the drought decidious plants the main impact is to have leaves
      ! drop in response to drought stress. You don't want leaf growth
      ! initiated the same way as traditional decidious plants. Drought
      ! decidious plants have their leaf growth start when it get wet enough
      ! and is not controlled by temperature directly like the traditional
      ! decidious plants."

      ! Changes
      ! Jan 2017 KLK
      !   remove potent.inc; pass pforc, tree_a2drat as parameteras.
      ! May 2015 KLK
      !   removed potforc, it is just an alias for pforc in potent.
      ! July 2015 KLK
      !   changed to use more F90 syntax
      !   initial Fruit/nut changes
      ! March 2015 KLK
      !   Updated routine description
      !   Incorporate more F90 matrix syntax
      !   Simplify allocation code so the IF block sets the root growth then allocates
      !     the rest to the leaves. Once minimum root growth was set (cak - 03/27/2007)
      !     this was what the routine did but it was harder to read.

! ... Function declarations
      integer   grochk
      real      froota, rtimp
      external  froota, grochk, rtimp

! ... Local variables
      integer, save :: greenUpCnt = 0
      logical, save :: greenUpErr = .false.
      real          :: availm(MAXIEL), eavail(MAXIEL)
      real          :: rimpct, fracrc
      real, parameter :: toler = 1.0E-30


! ... Estimate the fraction of carbon going to the roots
      fracrc = sum(tfrtcw + tfrtcn)/4.0

! ... Estimate fine root and leaf production from total potential production
      ! rootprod = pforc * fracrc
      ! leafprod = pforc - rootprod  ! = pforc * (1.0 - fracrc)

      ! Parallel sum of nutrients available for plant growth.
      ! Nutrients available to trees are in the top tlaypg layers, cak 01/29/03
      ! RESULT = SUM(ARRAY, DIM[, MASK])
      availm(1:nelem) = sum(minerl(1:tlaypg, 1:nelem), 1, (minerl(1:tlaypg, 1:nelem).gt.0))

      ! Calculate impact of root biomass on available nutrients
      rimpct = rtimp(riint, rictrl, frootcj+frootcm)

      ! Calculate soil available nutrients, based on a maximum fraction
      ! (favail) and the impact of root biomass (rimpct), adding storage.
      ! converted to use tree favail   KLK 30Jan13
      eavail(1:nelem) = (availm(1:nelem) * favail(1:nelem,FORSYS) * rimpct) + forstg(1:nelem)
      eavail(N) = eavail(N) + snfxmx(FORSYS) * pforc  ! increase available N from fixation

      ! Estimate the demand based on the maximum E/C ratio.
      !   demand = leafprod / cerfor(IMIN,LEAF,iel) + rootprod / cerfor(IMIN,FROOTJ,iel)
      ! New calculation -mdh 5/10/01
      !   tree_a2drat(iel) = min(1.0, eavail(iel)/ demand)
      !   tree_a2drat(iel) = max(0.0, tree_a2drat(iel))
      !   tree_a2drat(iel) = max(0.0, min(1.0, eavail(iel)/ &
      !                      (leafprod / cerfor(IMIN,LEAF,iel) + rootprod / cerfor(IMIN,FROOTJ,iel))))
      tree_a2drat(1:nelem) = max(0.0, min(1.0, eavail(1:nelem)/ &
                         (pforc * ((1.0 - fracrc)/cerfor(IMIN,LEAF,1:nelem) + fracrc/cerfor(IMIN,FROOTJ,1:nelem)))))

      ! assume minimal fine root growth during leaf out, cak - 03/27/2007
      tree_cfrac(FROOT) = 0.1

      if (decid .eq. 1) then      ! Decidious forest
        ! Within the greenup period, deciduous trees allocate only a small root growth
        if (grochk(tavemth) .eq. 1) then
          greenUpCnt = greenUpCnt + 1
          ! All allocation can go to leaves during greenup
          ! Allow some fine root growth during leaf out, cak - 03/27/2007
          tree_cfrac(FROOT) = 0.1
          if (greenUpCnt .gt. 30  .and. .not. greenUpErr) then
            call message('Warning: in treeDynC, tree green up, greenUpCnt, > 30 days')
            greenUpErr = .true.
          !  call abortrun('treeDynC, greenUpCnt > 30')
          endif
        elseif (decidgrow) then
          ! If we are in the period between leaf out and leaf drop, allocate
          ! to fine roots first and then to leaves
          greenUpCnt = 0
          tree_cfrac(FROOT) = froota(tree_a2drat,h2ogef(2),FORSYS)
        else
          ! No growth occurs this time period
          pforc = 0.0
          tree_cfrac(FROOT) = 0
          goto 999
        endif

      else if (decid .eq. 2) then      ! Drought decidious forest
! ..... Determine if we are within the greenup period for drought deciduous trees
        if (hrsinc .and. h2ogef(2) .gt. 0.5) then
! ....... All allocation can be to leaves during greenup
! ....... Allow some fine root growth during leaf out, cak - 03/27/2007
          tree_cfrac(FROOT) = 0.1
        else
! ....... Allocate to fine roots first and then to leaves
          tree_cfrac(FROOT) = froota(tree_a2drat,h2ogef(2),FORSYS)
        endif

      else if (decid .eq. 0) then      ! Evergreen forest
        ! Use froota to allocate to fine roots first
        tree_cfrac(FROOT) = froota(tree_a2drat,h2ogef(2),FORSYS)
      else
        call abortrun("Invalid forest type in treeDynC!")
      endif

      tree_cfrac(LEAF)  = 1.0 - tree_cfrac(FROOT)
      tree_cfrac(FBRCH) = 0.0
      tree_cfrac(LWOOD) = 0.0
      tree_cfrac(CROOT) = 0.0
      tree_cfrac(FRNUT) = 0.0

999   continue

      return
      end
