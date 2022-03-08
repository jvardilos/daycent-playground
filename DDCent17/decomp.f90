
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... DECOMP.F

      subroutine decomp(decodt,decsys,amovdly,newminrl, bgwfunc, &
                        agdefacsum, bgdefacsum, avgstemp, avgwfps, pptday, rprpet,petdly, &
                        time, cyear, month, doy, aceqcnt)

      USE ISO_C_Binding
      implicit none

      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'sitsoil.inc'

! ... Argument declarations
      real    :: decodt, time
      real    :: bgwfunc, agdefacsum, bgdefacsum, avgwfps, pptday, rprpet, petdly
      integer :: decsys, cyear, month, doy, aceqcnt
      real    :: amovdly(CMXLYR)
      real    :: anerdcmp(4,2)
      double precision newminrl

! ... Decomposition Submodel (rewritten by vek 04/91)
      ! include the full dailymoist.f90 decomposition loop into this routine

! ... Function declarations
      real     agdrat, bgdrat, fsfunc
      external agdrat, bgdrat, fsfunc

! ... Local variables
      integer       :: iel, kts

      real                          :: tfunc, avgstemp
      real, dimension(MAXIEL)       :: minrlavail
      real, dimension(2,MAXIEL)     :: eins2
      real, dimension(2)            :: cins2, couts2
      real, dimension(MAXIEL)       :: eins3
      real                          :: couts3

! ... Zero output variables tracking the daily carbon decomposition flows.
      metc1tosom11  = 0.0
      metc2tosom12  = 0.0
      struc1tosom11 = 0.0
      struc1tosom21 = 0.0
      struc2tosom12 = 0.0
      struc2tosom22 = 0.0
      som11tosom21  = 0.0
      som12tosom22  = 0.0
      som12tosom3   = 0.0
      som21tosom11  = 0.0
      som21tosom22  = 0.0
      som22tosom12  = 0.0
      som22tosom3   = 0.0
      som3tosom12   = 0.0
      wood1tosom11  = 0.0
      wood1tosom21  = 0.0
      wood2tosom11  = 0.0
      wood2tosom21  = 0.0
      wood3tosom12  = 0.0
      wood3tosom22  = 0.0

      ! equilibrium acceleration
      eins3  = 0.0;      couts3  = 0.0;
      eins2  = 0.0;      couts2  = 0.0;


! ... Combined effects of temperature and moisture on decomposition; returns anerdcmp
      anerdcmp = 1;
      call calcdefac(texture, tfunc, bgwfunc, agdefac, bgdefac, avgwfps, teff, &
                     rprpet, idef, pptday, snow, avgstemp)

! ... Calculate the effect of anerobic conditions on decomposition
      anerdcmp(:,2) = anerob(aneref,drain,rprpet,petdly,0) ! decomposition under anerobic conditions
      ! call thawdecp(doy, anerdcmp);  ! multiply anerdcmp by the thaw pulse modifications

      ! for output, average the soil anerobic conditions on decomposition
      anerb = anerb + anerdcmp(3,2) ! anerob(aneref,drain,rprpet,petdly,0)

! ... calculate defacm(month) in subroutine simsom. -mdh 10/94
      agdefacsum = agdefacsum + agdefac
      bgdefacsum = bgdefacsum + bgdefac

      ! run the cultivation effect model
      call tilleff(month, doy, bgwfunc)

      ! Scale pH values if necessary
      if (phsys .gt. 0) ph = phstart * pHscalar(month)


      do kts = 1, ntspm

        ! minrlavail is the mineral N, P, and S available for decomposition
        ! minerl contains the current value of mineral N, P, and S by layer.
        !minrlavail = minerl
        !if (nelem .ge. P) minrlavail(P) = minrlavail(P) * fsfunc(minerl(1,P), pslsrb, sorpmx)

        minrlavail = minerl(1,:)
        if (nelem .ge. P) minrlavail(P) = minrlavail(P) * fsfunc(minerl(1,P), pslsrb, sorpmx)

        aminrl(1:nelem) = aminrl(1:nelem) + minrlavail(1:nelem)

        ! -----old decomp subroutine -----------------------------------------
        ! ***** call original decomp(decodt,decsys,amovdly,newminrl) ****

        ! Create ratnew1 and ratnew2 for surface and soil. -MDH 9/24/2012
        ! ratnew1(iel,1) - the C/E ratio for new material created when a
        !                  lignin component decomposes to SRFC som1.
        ! ratnew1(iel,2) - the C/E ratio for new material created when a
        !                  lignin component decomposes to SRFC som2.
        ! ratnew2(iel,1) - the C/E ratio for new material created when a
        !                  lignin component decomposes to SOIL som1.
        ! ratnew2(iel,2) - the C/E ratio for new material created when a
        !                  lignin component decomposes to SOIL som2.

        ! Determine C/E ratios for flows from structural material to
        ! surface som1 and surface som2
        do iel=1,nelem
          ! ratnew1: SRFC som1 and som2
          ratnew1(iel,1) = agdrat(minrlavail,varat11,iel)
          ratnew1(iel,2) = agdrat(minrlavail,varat21,iel)
          ! ratnew2: SOIL som1 and som2
          ratnew2(iel,1) = bgdrat(minrlavail,varat12,iel)
          ratnew2(iel,2) = bgdrat(minrlavail,varat22,iel)
        end do

        !             ********** LITTER **********
        ! Decompose structural and metabolic components for surface and soil.
        call litdec(decodt, newminrl, anerdcmp, minrlavail, eins2)

        !             *********** WOOD ***********
        ! If the system is a forest or savanna...
        ! Decompose dead fine branches, large wood, and coarse roots.
        ! Dead fine roots are in the soil structural compartment.
        if (decsys .eq. FORSYS) then
          call woodec(decodt, newminrl, anerdcmp, minrlavail, eins2)
        endif

        !             ***** SOIL ORGANIC MATTER *****
        ! Decompose som1 and som2 (surface and soil) and som3.
        ! Added amovdly parameter for daily version. -mdh 10/10/94
        call somdec(amovdly, decodt, newminrl, anerdcmp, &
                    minrlavail, couts2, eins2, couts3, eins3)
        ! -----end old decomp subroutine -------------------------------------

        !             ***** update portion of dailymoist loop *****
        if (nelem .ge. P) call pschem(decodt)

        ! Update decomposition and nitrogen fixation flows.
        call flowup(time)
        call sumcar

      end do

      ! Accumulate output variables tracking decomposition carbon flows  cak - 11/09/2010
      ametc1tosom11  = ametc1tosom11  + metc1tosom11
      ametc2tosom12  = ametc2tosom12  + metc2tosom12
      astruc1tosom11 = astruc1tosom11 + struc1tosom11
      astruc1tosom21 = astruc1tosom21 + struc1tosom21
      astruc2tosom12 = astruc2tosom12 + struc2tosom12
      astruc2tosom22 = astruc2tosom22 + struc2tosom22
      asom11tosom21  = asom11tosom21  + som11tosom21
      asom12tosom22  = asom12tosom22  + som12tosom22
      asom12tosom3   = asom12tosom3   + som12tosom3
      asom21tosom11  = asom21tosom11  + som21tosom11
      asom21tosom22  = asom21tosom22  + som21tosom22
      asom22tosom12  = asom22tosom12  + som22tosom12
      asom22tosom3   = asom22tosom3   + som22tosom3
      asom3tosom12   = asom3tosom12   + som3tosom12
      awood1tosom11  = awood1tosom11  + wood1tosom11
      awood1tosom21  = awood1tosom21  + wood1tosom21
      awood2tosom11  = awood2tosom11  + wood2tosom11
      awood2tosom21  = awood2tosom21  + wood2tosom21
      awood3tosom12  = awood3tosom12  + wood3tosom12
      awood3tosom22  = awood3tosom22  + wood3tosom22

      ! Accumulate flows for accelerated equilibrium
      if(aceqcnt > 0) then
        cins2 = [struc1tosom21 + som11tosom21 + wood1tosom21 + wood2tosom21, &
                 struc2tosom22 + som12tosom22 + som21tosom22 + wood3tosom22]
        call EqAclratCE(som12tosom3 + som22tosom3, eins3, couts3, &
                        cins2, eins2, couts2)
      endif

      return
      contains

        subroutine tilleff(month, curday, wfunc)
          integer month, curday
          real wfunc

          include 'dovars.inc'    ! docult
          include 'parcp.inc'     ! clteff
          include 'parfx.inc'     ! maxcltef, xefcltef, cfita, cfitb
          include 'plot1.inc'     ! bgdefac
          include 'plot2.inc'     ! cltfac

        !...Local variables
          character mssg*132
          integer i

        !... tilleff is the routine that calculates the cultivation effect on
          ! decomposition. This code was moved from cycle so that the more
          ! complex model can be called in both the microcosm and regular modes.

          !  Effect of cultivation on decomposition (used in decomp routine)
          !  This code has been moved to the simsom subroutine, cak - 04/17/03
          !  Determine effect of cultivation on decomposition. vk 03-13-91
          !  cltfac is this month's effect of cultivation on decomposition
          !  of som1, som2, som3, and structural.  It is set to clteff
          !  in months when cultivation occurs; otherwise it equals 1.
          !  clteff is the effect of cultivation on decomposition read from the cult.100 file

          ! Modifications
          !...8/2003  K. Killian
          !   Two Modifications to the mixing to CLTEFF model. First, defined the
          !    intercept = 1.0 to ensure the function is continous at zero.
          !    Second, restricted the equation to a general linear or quadratic fit.
          !    Removed upper limit on wfunc that is extraneous with current equations
          !...9/2001  K Killian
          !    Added a method for adding cultivations in consecutive months for the
          !    deferred cultivation effect model based on the relationship between mixing
          !    and cultef developed by M Sperow. Cultivations are converted to mixing
          !    then the mixings are added by cultivating the appropriate fraction
          !    of the unmixed soil.
          !
          !  cltfit is the equation to convert cultef to mixing.
          !  cltfit   = (cfita * tmix + cfitb) *tmix + cfitc
          !
          !  with the coefficients defined as:
          !  cfita  = 0.63759*(maxcltef-1)/1.7594
          !  cfitb  = 1.0232 *(maxcltef-1)/1.7594
          !  cfitc  = maxcltef-(cltfita+cltfitb)-0.01350368 *(maxcltef-1)/1.7594
          !         = maxcltef-(0.63759+1.0232+0.01350368)*(maxcltef-1)/1.7594
          !         = maxcltef-1.6742936751*(maxcltef-1)/1.7594
          !
          !  mxfit is the equation to convert cultef to mixing
          !  mixfit   =(-cfitb + SQRT(cfitb**2-4*cfita*(cfitc-CULTEF)))/(2*cfita)
          !
          ! ------------- original Month by month code -------------------
          !        if (docult) then
          !          do 32 i = 1, 4
          !            cltfac(i) = abs(clteff(i))
          !32        continue
          !        else
          !          do 33 i = 1, 4
          !            cltfac(i) = 1.0
          !33        continue
          !        endif
          ! -------------------------------------------------------------

          ! mixing variables
          real cltef, cltefn, cfac, mx1, mx2, tmix, cltfc
          real :: coef
          coef = 31. * xefcltef * log(2.)
          ! mixfit   =-cfitb/(2*cfita) + SQRT((cfitb/(2*cfita))**2 - (1.0-CULTEF)/cfita))



          ! ------ DayCent Simulated Month code -------------------
          if(xefcltef .eq. 0) then
            if (docult .and. (cultday .eq. curday) .or. &
              (cultcnt .gt. 0. .and. cultcnt .lt. 31)) then
              cltfac = abs(clteff)
              if (cultday .eq. curday) cultcnt = 0
              cultcnt = cultcnt + 1
            else
              cltfac = 1.0
              cultcnt = 0
            endif
            return
          endif

    !...Effect of cultivation on decomposition (used in decomp routine)
    !...Determine effect of cultivation on decomposition. vk 03-13-91
    !...clteff is the effect of cultivation on decomposition read from
    !     the cult.100 file
    !
    !    if clteff < 0 print out the debug status prints.
    !
    ! New cultivation effect model   KLK
    !     cltfac is the effect of cultivation on decomposition
    !     of som1, som2, som3, and structural.
    !...xefcltef  is greater than zero:
    !     clteff is set to the maximum of cltfac or the existing effect.
    !     The enhancement decays by a fraction of calculated bgdefac
    !...xefcltef  equals zero
    !     It is a single month's effect on decomposition.
    !     The effect is clteff in cultivation months; otherwise it equals 1.
    !
    ! -----------------
    !  should change the function to an exponential sum
    !      (Exp[-f1*bgdefac*days] - Exp[-f2*bgdefac*days])*bgdefac
    !     where
    !           x  = decomposition days, bgdefac *days since since CULT event
    !           f1 = 0.173287 (assume 4 day half life  xefcltef=0.13147 )
    !           f2 = 9/8
    !
    !     The integral is 1.   The function peaks at x = 2.
    !     The function can be approximated as:
    !     cltef = (t <= 3.931)? 0.494009 : 0.999284 exp(-t*ln(2)/4);
    !     This both simplifies the function and preserves the integral.
    ! -----------------

            if (docult) then
              do i = 1, 4
                ! Take care of the new model where we add the cultivation effects
                if (maxcltef .gt. 0.) then
                  ! mxfit is the equation to convert cultef to mixing
                  cfac  = cltfac(i) ! save for debug output

                  if(cltfac(i) .eq. 1.) then
                    cltefn   = abs(clteff(i))
                  else
                    ! convert cltfac back to a value a clteff value so we can
                    ! use the fit to convert it to a mix fraction
                    cltef = 1+ (cltfac(i)-1.)/(coef*(xefcltef+1.))

                    if(cfita .eq. 0.) then
                      mx1= (cltef         -1.)/cfitb
                      mx2= (abs(clteff(i))-1.)/cfitb
                    else
                      ! use tmix as a temporary variable to hold the quadratic b/(2a)
                      tmix= cfitb/(2.*cfita)
                      mx1= SQRT(tmix*tmix +(cltef         -1.)/cfita) -tmix
                      mx2= SQRT(tmix*tmix +(abs(clteff(i))-1.)/cfita) -tmix
                    endif

                    ! now add the mixing fractions
                    if(mx1 .gt. 1.0  .and.  mx2 .gt. 1.0) then
                    ! Calculate and addition for overmixed values
                      tmix   = mx1 * mx2
                    ! add the fractions.
                    else
                      tmix   = mx1+(1.-mx1)*mx2
                    end if

                   ! Calculate the effective clteff from the mixing
                    cltefn   = (cfita * tmix + cfitb) *tmix + 1.0
                  endif

                 ! Convert the effective clteff to CLTFAC using xefcltef
                  cltfac(i)= 1. + (cltefn-1.)*(coef*(xefcltef +1.))/ &
                                  (1.+bgdefac*xefcltef)

                else
                  cltfac(i)= abs(clteff(i))
                endif

                ! verbose print out if clteff(i) < 0
                if (clteff(i) .lt. 0) then
                  write(mssg,'(a,i2,i4,2x,11f9.4)') 'CULT event ',month, &
                  curday,bgdefac,xefcltef,cfac,cltef,mx1, &
                        abs(clteff(2)),mx2,tmix,cltefn,cltfc,cltfac(i)
                  call message(mssg)
                endif
              end do
            else
              do i = 1, 4
                if (cltfac(i) .eq. 1.0) then
                elseif (cltfac(i) .le. 1.01) then
                  !if (clteff(i).lt.0 .and. cltfac(2).gt.1.0) then
                  if (clteff(i) .lt. 0.) then
                    write(mssg,'(i5,i4,3f10.3,a)') month,curday, &
                       bgdefac,xefcltef,cltfac(i),' ->    0.'
                    call message('clear cltfac  '//trim(mssg))
                  end if
                  cltfac(i) = 1.0
                else
                  cfac = cltfac(i)
                  cltfac(i) = 1. + (cltfac(i)-1.)/(1.+bgdefac*xefcltef)
                  if (clteff(i) .lt. 0.) then
                    write(mssg,'(a,i5,i4,3f10.3,a,f8.3)') 'change cltfac ', &
                       month,curday,bgdefac,xefcltef,cfac,' -> ',cltfac(i)
                    call message(mssg)
                  end if
                endif
              enddo

            endif

          return
        end subroutine tilleff

        real function anerob(aneref, drain, rprpet, pet, micosm)

         implicit none

   ! ... Argument declarations
         real     aneref(3), drain, rprpet, pet
         integer  micosm

   ! ... Calculates the impact of soil anaerobic conditions on decomposition.
   !     It returns a multiplier 'anerob' whose value is 0-1.

   ! ... If a microcosm is being simulated, return the value for aneref(3)

   !     Declaration explanations:
   !       aneref(1) - ratio RAIN/PET with maximum impact
   !       aneref(2) - ratio RAIN/PET with minimum impact
   !       aneref(3) - minimum impact
   !       drain     - percentage of excess water lost by drainage
   !       newrat    - local var calculated new (RAIN+IRRACT+AVH2O(3))/PET ratio
   !       pet       - potential evapotranspiration
   !       rprpet    - actual (RAIN+IRRACT+AVH2O(3))/PET ratio

        ! Local variables
         real      newrat, slope, xh2o

         anerob = 1.0

         ! Check if simulating a microcosm
         if (micosm .eq. 1) then
           anerob = aneref(3)

         ! Determine if RAIN/PET ratio is greater than the ratio with
         ! maximum impact.
         ! If it is winter time the value for rprpet has been calculated
         ! using snow that is melting into the soil profile, cak - 10/21/02
         else if (rprpet .gt. aneref(1)) then
           xh2o = (rprpet - aneref(1)) * pet * (1.0 - drain)
           if (xh2o .gt. 0) then
             newrat = aneref(1) + (xh2o / pet)
             slope = (1.0 - aneref(3)) / (aneref(1) - aneref(2))
             anerob = 1.0 + slope * (newrat - aneref(1))
           endif

           anerob = max(anerob, aneref(3))
         endif

         return
        end function anerob

      end subroutine decomp
