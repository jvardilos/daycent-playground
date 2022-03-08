
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine calciv(swdbgf)
      use calflow;

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'fertil.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'zztim.inc'

! ... Calculate initial values for temperature, water and live root
! ... carbon variables.
! ... Called from detiv.
! ... Note that variables which are functions of changeable parameters
! ... (i.e. read from 'site'.par) should be computed in prelim instead
! ... of calciv.

      ! Cindy Keough May 8, 2014,
      ! Correct the double accounting of labeled C while initializing C pools
      ! 31Aug11  CAK
      ! Modifications to Burke equation SOM allocation to the som2 and som3
      ! pools based on the clay content for the grassland and forest systems

      integer :: swdbgf(2) ! site verbose/debug flags

! ... Local variables
      integer   iel
      real      avtemp, arain, dumye(MAXIEL), dumyc(ISOS), fraclabl
      real      k2, reclt1(MAXIEL), reclt2(MAXIEL), tcg
!      real      storFrac
      real      som21frac, som22frac, som3frac
      real      tcl
      real      ramp
      character string*80, char1*1
      character subname*10

! ... Initialize soil C pools using Burke's equations.
! ...   ivauto = 0  the user has supplied the initial values
! ...   ivauto = 1  initialize using the grassland soil parameters
! ...   ivauto = 2  initialize using the crop soil parameters
! ...   ivauto = 3  initialize soil C pools for a forest soil

      subname = 'calciv    '

! ... Initialize dumyc and dumye variables.
      dumyc(LABELD) = 1000.
      dumyc(UNLABL) = 1000.
      dumye(N) = 100.
      dumye(P) = 100.
      dumye(S) = 100.

! ... Initialize the variable for computing the fraction of labeled
! ... material based on the type of labeling being done
      if (labtyp .eq. 1) then
! ..... C14 labeling
!        fraclabl = 0.0011
! ..... Calculate a fraction that is equivalent to a zero delta 14C
! ..... value, cak - 05/03/2007
        k2 = FRAC_C14 * (1.0 + (1.0 / 1000.0))
        fraclabl = k2 / (1.0 + k2)
      elseif (labtyp .eq. 2) then
! ..... C13 labeling
        fraclabl = 0.011
      else
! ..... No labeling
        fraclabl = 0.0
      endif

! ... Initialize irrtot, accumulator in irrigt.  -mdh 12/9/96
      irrtot = 0.0

! ... Compute mean annual temperature (avtemp) and mean annual
! ... precipitation (arain)
      arain = sum(precip)
      avtemp = sum(tmn2m + tmx2m)/24. ! sum((tmn2m + tmx2m)/2)/12.
      if (avtemp .gt. 23.) then
        avtemp = 23.
      endif
      if (arain .gt. 120.) then
        arain = 120.
      endif

! ... Initialize soil C pools for cultivated soils
      if (ivauto .eq. 2) then

! ..... tcg = total soil carbon in grams (som1c + som2c + som3c)
        ! limit initial soil carbon to not less than 500 g/m^2, cak - 03/22/02
        tcg = max(((2.10E-02 * avtemp - 7.50E-01) * avtemp + &
              arain * (-4.58E-04 * arain + 5.81E-02 +        &
              4.94E-02 * silt + 5.82E-02 * clay) + 5.15) * 1000., 500.)

! ..... Assign a fixed value to surface som1.   vek  08-91
        som1ci(SRFC,UNLABL) = 10.
        som2ci(SRFC,UNLABL) = 0.0

! ..... Burke's equations only apply to soil compartments. vek  08-91
        som1ci(SOIL,UNLABL) = tcg * .02
        som2ci(SOIL,UNLABL) = tcg * .54
        som3ci(UNLABL) = tcg * .44
        stdcis(UNLABL) = 20.

        stdede(N) = .40
        stdede(P) = .075
        stdede(S) = .075
        clittr(SRFC,UNLABL) = 10.

      elseif(ivauto .eq. 1  .or.  ivauto .eq. 3) then

! ..... tcg = total soil carbon in grams (som1c + som2c + som3c)
        ! limit initial soil carbon to not less than 500 g/m^2, cak - 03/22/02
        tcg = max(((2.24E-02 * avtemp - 8.27E-01) * avtemp + &
              arain * (1.27E-01 - 9.38E-04 * arain +         &
              silt * 8.99E-02 + clay * 6.00E-02) + 4.09) *1000., 500.)

        ! use clay to determine the som2 and som3 C ratio
        som3frac = ramp(clay, 0.08, 0.47, 0.50, 0.57)

        ! Assign a fixed value to surface som1.   vek  08-91
        som1ci(SRFC,UNLABL) = 10.

        ! Assign initial values to the labeled pools as well as unlabeled
        ! pools using a percentage of the value from the unlabeled pool, cak - 03/22/02

        if(ivauto .eq. 1) then ! soil C pools for a grassland soil
          som22frac = (0.98 - som3frac)       ! * 1.0
          ! On Nov 2, 2015, at 1:02 PM, Parton,William wrote:
          ! For grasslands ... use a 10 % rule where surface som2 is equal to 10 % of mineral som2.
          som21frac = 0.1 * som22frac         ! (0.98 - som3frac) * 0.0
        else                    ! soil C pools for a forest soil
          som21frac = (0.98 - som3frac) * 0.2
          som22frac = (0.98 - som3frac) * 0.8
        endif

        ! som2 pool has been split into surface and soil pools,
        ! initialize both separately, cak - 06/14/05
        som2ci(SRFC,UNLABL) = tcg * som21frac

! ..... Burke's equations only apply to soil compartments.
        som1ci(SOIL,UNLABL) = tcg * .02
        som2ci(SOIL,UNLABL) = tcg * som22frac

        som3ci(UNLABL) = tcg * som3frac

        stdcis(UNLABL) = 80.
        stdede(N) = 1.6
        stdede(P) = .3
        stdede(S) = .3

        bglcisj(UNLABL) = 100.
        bglcisj(LABELD) = bglcisj(UNLABL) * fraclabl
        bglcisj(UNLABL) = bglcisj(UNLABL) - bglcisj(LABELD)
        bglivej(N) = 1.5
        bglivej(P) = .25
        bglivej(S) = .25
        bglcism(UNLABL) = 100.
        bglcism(LABELD) = bglcism(UNLABL) * fraclabl
        bglcism(UNLABL) = bglcism(UNLABL) - bglcism(LABELD)
        bglivem(N) = 1.5
        bglivem(P) = .25
        bglivem(S) = .25

        clittr(SRFC,UNLABL) = 100.
      endif
      ! the unlabeled pools were not being set Cindy Keough May 8, 2014, at 1:17 PM
      ! Done here since the code for all options was the same
      som1ci(:,LABELD) = som1ci(:,UNLABL) * fraclabl
      som1ci(:,UNLABL) = som1ci(:,UNLABL) - som1ci(:,LABELD)
      som2ci(:,LABELD) = som2ci(:,UNLABL) * fraclabl
      som2ci(:,UNLABL) = som2ci(:,UNLABL) - som2ci(:,LABELD)
      som3ci(LABELD) = som3ci(UNLABL) * fraclabl
      som3ci(UNLABL) = som3ci(UNLABL) - som3ci(LABELD)
      stdcis(LABELD) = stdcis(UNLABL) * fraclabl
      stdcis(UNLABL) = stdcis(UNLABL) - stdcis(LABELD)
      clittr(:,LABELD) = clittr(:,UNLABL) * fraclabl
      clittr(:,UNLABL) = clittr(:,UNLABL) - clittr(:,LABELD)

! ... End of soil C pool initialization

! ... Starting values for nitrogen, phosphorous, and sulfur depend on
! ... carbon values and the ratios of carbon to each other element.
! ... Initialize structural and metabolic pools C, N, P, and S.
! ... First set them to zero and calculate N/C, P/C, & S/C ratios.
      strcis = 0.
      metcis = 0.
      struce = 0.
      metabe = 0.

! ... Compute N/C, P/C, and S/C ratios from C/N, C/P, and C/S.
! ... This is for use in partit.
! ... Added the conditional set to zero if rcelit <= 0 -rm 7/98
      do iel = 1, MAXIEL
        if (rcelit(SRFC, iel) .gt. 0.) then
          reclt1(iel) = 1. / rcelit(SRFC, iel)
        else
          reclt1(iel) = 0.0
        endif

        if (rcelit(SOIL, iel) .gt. 0.) then
          reclt2(iel) = 1. / rcelit(SOIL, iel)
        else
          reclt2(iel) = 0.0
        endif
      end do

! ... Sum carbon isotopes for use in partit.
      call sumcar

! ... Split litter C content into structural/metabolic based upon
! ... litter C and litter lignin content and compute structural and
! ... metabolic N, P, & S based upon amount of C and the ratios
! ... computed above.
      if (curcrp .ne. ' ' .and. curtre .ne. ' ') then
        pltlig(ABOVE) = (wdlig(LEAF)+fligni(INTCPT,ABOVE) +     &
                        fligni(SLOPE,ABOVE) * arain) / 2.0
        pltlig(BELOWJ) = (wdlig(FROOTJ)+fligni(INTCPT,BELOWJ) + &
                         fligni(SLOPE,BELOWJ) * arain) / 2.0
        pltlig(BELOWM) = (wdlig(FROOTM)+fligni(INTCPT,BELOWM) + &
                         fligni(SLOPE,BELOWM) * arain) / 2.0
      else if (curcrp .ne. ' ') then
        pltlig(ABOVE) = fligni(INTCPT,ABOVE)+fligni(SLOPE,ABOVE) * arain
        pltlig(BELOWJ) = fligni(INTCPT,BELOWJ)+fligni(SLOPE,BELOWJ) * arain
        pltlig(BELOWM) = fligni(INTCPT,BELOWM)+fligni(SLOPE,BELOWM) * arain
      else if (curtre .ne. ' ') then
        pltlig(ABOVE) = wdlig(LEAF)
        pltlig(BELOWJ) = wdlig(FROOTJ)
        pltlig(BELOWM) = wdlig(FROOTM)
      endif

! ... Total C in litter
      tcl = clittr(SRFC,UNLABL)+clittr(SRFC,LABELD)
      if(tcl .gt. 0) call partit(tcl,reclt1,1,dumyc,dumye,pltlig(SRFC),clittr(SRFC,LABELD)/tcl)
      tcl = clittr(SOIL,UNLABL)+clittr(SOIL,LABELD)
      if(tcl .gt. 0) call partit(tcl,reclt2,2,dumyc,dumye,pltlig(SOIL),clittr(SOIL,LABELD)/tcl)

      call flowup(time)
      call sumcar

      if(swdbgf(1).gt.0) call showminrl(nlayer,minerl,ammonium,nitrate,subname)

! ... If the C/E ratio as read from the site file is <= 0.0 throw an
! ... error message, cak - 06/14/05
      do iel=1,MAXIEL
! ..... Compute N, P, and S for surface and soil som1, as well as for
! ..... som2 and som3.   vek  08-91
        if (rces1(SRFC,iel) .gt. 0.) then
          som1e(SRFC,iel)=som1c(SRFC)/rces1(SRFC,iel)
        else
          call message('Warning: bad input from <site>.100.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES1(1,' // char1 // ') is <= 0.0'
          call message(string)
        endif
        if (rces1(SOIL,iel) .gt. 0.) then
          som1e(SOIL,iel)=som1c(SOIL)/rces1(SOIL,iel)
        else
          call message('Warning: bad input from <site>.100.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES1(2,' // char1 // ') is <= 0.0'
          call message(string)
        endif
! ..... som2 pool has been split into surface and soil pools, compute N,
! ..... P, and S for both pools, cak - 06/14/05
        if (rces2(SRFC,iel) .gt. 0.) then
          som2e(SRFC,iel)=som2c(SRFC)/rces2(SRFC,iel)
        else
          call message('Warning: bad input from <site>.100.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES2(1,' // char1 // ') is <= 0.0'
          call message(string)
        endif
        if (rces2(SOIL,iel) .gt. 0.) then
          som2e(SOIL,iel)=som2c(SOIL)/rces2(SOIL,iel)
        else
          call message('Warning: bad input from <site>.100.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES2(2,' // char1 // ') is <= 0.0'
          call message(string)
        endif
        if (rces3(iel) .gt. 0.) then
          som3e(iel)=som3c/rces3(iel)
        else
          call message('Warning: bad input from <site>.100.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES3(' // char1 // ') is <= 0.0'
          call message(string)
        endif
      end do

      if (curtre .ne. ' ') then
        do iel = 1, MAXIEL
          if (cerfor(IVAL,FBRCH,iel) .gt. 0.) then
            wood1e(iel)=wood1c/cerfor(IVAL,FBRCH,iel)
          else
            call message('Warning: bad input from <site>.100.')
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   cerfor(1,3,' // char1 // ') is <= 0.0'
            call message(string)
          endif
          if (cerfor(IVAL,LWOOD,iel) .gt. 0.) then
            wood2e(iel)=wood2c/cerfor(IVAL,LWOOD,iel)
          else
            call message('Warning: bad input from <site>.100.')
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   cerfor(1,4,' // char1 // ') is <= 0.0'
            call message(string)
          endif
          if (cerfor(IVAL,CROOT,iel) .gt. 0.) then
            wood3e(iel)=wood3c/cerfor(IVAL,CROOT,iel)
          else
            call message('Warning: bad input from <site>.100.')
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   cerfor(1,5,' // char1 // ') is <= 0.0'
            call message(string)
          endif
        end do
      endif

! ... Surface temperature and soil temperature
      tave = (tmn2m(1) + tmx2m(1)) / 2.0
      stemp = tave

! ... Make sure there is N, P, and S for roots
      if (bglcisj(UNLABL)+bglcisj(LABELD) .gt. 0.0) then
        do iel = 1, nelem
          if (bglivej(iel) .le. 0.) then
            string = 'Initial bglivej value ('//char(ichar(char(iel)) + ichar('0'))// &
                     ') must be greater than 0.'
            call abortrun(string)
          endif
        end do
      endif
      if (bglcism(UNLABL)+bglcism(LABELD) .gt. 0.0) then
        do iel = 1, nelem
          if (bglivem(iel) .le. 0.) then
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = 'Initial bglivem value ('//char1//') must be greater than 0.'
            call abortrun(string)
          endif
        end do
      endif

! ... Initialize grain pools
      cgrain = 0.0
      egrain = 0.0

! ... Initialize crop and tree carbohydrate storage to 20-30% of live C
! ... (Bill Parton - 11/29/01)
!      storFrac = 0.25
!      carbostg(CRPSYS,UNLABL) = storFrac *
!     &                          (aglcis(UNLABL) + bglcis(UNLABL))
!      carbostg(CRPSYS,LABELD) = storFrac *
!     &                          (aglcis(LABELD) + bglcis(LABELD))
!      carbostg(FORSYS,UNLABL) = storFrac *
!     &                          (rlvcis(UNLABL) + frtcis(UNLABL) +
!     &                           crtcis(UNLABL) + rlwcis(UNLABL) +
!     &                           fbrcis(UNLABL))
!      carbostg(FORSYS,LABELD) = storFrac *
!     &                          (rlvcis(LABELD) + frtcis(LABELD) +
!     &                           crtcis(LABELD) + rlwcis(LABELD) +
!     &                           fbrcis(LABELD))

! ... To grow vegetation up from zero set an initial values for the
! ... crop and tree carbohydrate storage pools, cak - 02/23/2007
      if (cursys .eq. CRPSYS) then
        carbostg(CRPSYS,UNLABL) = 250.0
      else if (cursys .eq. FORSYS) then
        carbostg(FORSYS,UNLABL) = 250.0
      else if (cursys .eq. SAVSYS) then
        carbostg(CRPSYS,UNLABL) = 250.0
        carbostg(FORSYS,UNLABL) = 250.0
      endif
      carbostg(:,LABELD) = carbostg(:,UNLABL) * fraclabl
      carbostg(:,UNLABL) = carbostg(:,UNLABL) - carbostg(:,LABELD)

! ... a2drat - available E to plant demand for E.  -mdh 8/25/00
! ... Create separate a2drat arrays for crops and trees, mdh 5/11/01
      crop_a2drat = 1.0
      tree_a2drat = 1.0

      return
      end
