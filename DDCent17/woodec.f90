
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine woodec (dtm, newminrl, anerdcmp, minrlavail, eins2)

      implicit none
      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

! ... Argument declarations
      real             dtm
      real, dimension(MAXIEL)   :: minrlavail
      real, dimension(2,MAXIEL) :: eins2(2,MAXIEL)
      real             anerdcmp(4,2)
      integer, parameter :: met = 1, strc = 2, som23 = 3, wood = 4
      double precision newminrl

! ... Wood decomposition       written by vek 04/91

! ... defac  = decomposition factor based on water and temperature
! ...          (computed in prelim and in cycle)
! ... pligst = fixed parameter that represents the effect of
! ...          of lignin-to-structural-ratio on structural
! ...          decomposition
      ! separate surface and soil C/E ratios MDH Nov 2012

! ... Local variables
      real      tcflow, pheff

! ... FINE BRANCHES
! ...   wood1c       = C in dead fine branch component of forest system (g/m2)
! ...   decw1        = intrinsic rate of decomposition of dead fine branches
! ...   wdlig(FBRCH) = lignin fraction for fine branches

      if (wood1c .gt. 1.e-07) then

! ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

! ..... Compute total C flow out of fine branches
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = wood1c * defac * decw1 * exp(-pligst(SRFC) *
!     &           wdlig(FBRCH)) * dtm
        tcflow = wood1c * agdefac * decw1 * exp(-pligst(SRFC) * &
                 wdlig(FBRCH)) * anerdcmp(wood,SRFC) * dtm * pheff

! ..... Decompose fine branches into som1 and som2 with CO2 loss.
        ! use surface decomposition ratios MDH Nov 2012
        call declig(minrlavail,wdlig(FBRCH),SRFC,nelem,1,ps1co2,ratnew1, &
                    rsplig,tcflow,wood1c,wd1c2,wd1cis,wood1e, &
                    gromin,minerl,w1mnr,resp,som1ci,som1e,som2ci,som2e, &
                    wood1tosom11,wood1tosom21,newminrl,eins2(SRFC,:))
      endif

! ... LARGE WOOD
! ...   wood2c       = C in dead large wood component of forest system (g/m2)
! ...   decw2        = intrinsic rate of decomposition of dead large wood
! ...   wdlig(LWOOD) = lignin fraction for large wood

      if (wood2c .gt. 1.e-07) then

! ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

! ..... Compute total C flow out of large wood
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = wood2c * defac * decw2 * exp(-pligst(SRFC) *
!     &           wdlig(LWOOD)) * dtm
        tcflow = wood2c * agdefac * decw2 * exp(-pligst(SRFC) * &
                 wdlig(LWOOD)) * anerdcmp(wood,SRFC) * dtm * pheff

! ..... Decompose large wood into som1 and som2 with CO2 loss.
        ! use surface decomposition ratios MDH Nov 2012
        call declig(minrlavail,wdlig(LWOOD),SRFC,nelem,1,ps1co2,ratnew1, &
                    rsplig,tcflow,wood2c,wd2c2,wd2cis,wood2e, &
                    gromin,minerl,w2mnr,resp,som1ci,som1e,som2ci,som2e, &
                    wood2tosom11,wood2tosom21,newminrl,eins2(SRFC,:))
      endif

! ... COARSE ROOTS
! ...   wood3c       = C in dead coarse root component of forest system (g/m2)
! ...   decw3        = intrinsic rate of decomposition of dead coarse roots
! ...   wdlig(CROOT) = lignin fraction for coarse roots

      if (wood3c .gt. 1.e-07) then

! ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

! ..... Compute total C flow out of coarse roots.
! ..... Add pH effect on decomposition to calculation, cak - 08/02/02
!        tcflow = wood3c * defac * decw3 * exp(-pligst(SOIL) *
!     &           wdlig(CROOT)) *  anerb * dtm
        tcflow = wood3c * bgdefac * decw3 * exp(-pligst(SOIL) * &
                 wdlig(CROOT)) *  anerdcmp(wood,SOIL) * dtm * pheff

! ..... Decompose coarse roots into som1 and som2 with CO2 loss.
        ! use soil decomposition ratios MDH Nov 2012
        call declig(minrlavail,wdlig(CROOT),SOIL,nelem,1,ps1co2,ratnew2, &
                    rsplig,tcflow,wood3c,wd3c2,wd3cis,wood3e, &
                    gromin,minerl,w3mnr,resp,som1ci,som1e,som2ci,som2e, &
                    wood3tosom12,wood3tosom22,newminrl,eins2(SOIL,:))

      endif

      return
      contains
       include 'catanf.f'
      end subroutine woodec
