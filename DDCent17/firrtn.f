
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

! ... FIRRTN

      subroutine firrtn()
      use calflow;

      implicit none
      include 'const.inc'
      include 'fertil.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

! ... Elemental return from a fire event.

! ... Called from:  frem

! ... Local Variables
      integer   iel, clyr
      real      egain(MAXIEL)
      character subname*10
      double precision frac_nh4, frac_no3

! ... LITTER BURN

      subname = 'firrtn    '

! ... Only burn litter if a forest system.  Litter is burned in
! ... grem.f for the savanna system.

      egain = 0.0 ! vector implied do

! ... Litter is being burned for the forest systems in grem.f as well,
! ... cak - 08/23/02
!      if (dofire(FORSYS)) then
!        call litburn(egain)
!      endif

! ... Return from TREE compartments
! ... Carbon return is usually 0.  It is ignored since it
! ... would be returned as charcoal.  N, P, and S returns
! ... go to the top layer of minerl.  EGAIN will contain
! ... the total returns for N, P, and S across pools.
! ... No longer returning elements from the dead fine branch
! ... and dead large wood forest components since these
! ... components are no longer burned during a TREM event,
! ... cak - 01/02

      do 20 iel = 1, nelem
        egain(iel) = egain(iel) +
     &               remf(1) * retf(1,iel+1) * rleave(iel) +
     &               remf(2) * retf(2,iel+1) * fbrche(iel) +
     &               remf(3) * retf(3,iel+1) * rlwode(iel) +
     &               remf(3) * retf(3,iel+1) * forstg(iel) +
     &               remf(6) * retf(4,iel+1) * frnute(iel) ! remove fruit/nuts at fine branch rate
!     &               remf(4) * retf(2,iel+1) * wood1e(iel) +
!     &               remf(5) * retf(3,iel+1) * wood2e(iel) +

        frac_nh4 = 0.5
        frac_no3 = 0.5
        if (iel .eq. N) then
          clyr = 1
          call update_npool(clyr, egain(iel), frac_nh4, frac_no3,
     &                      ammonium, nitrate, subname)
        endif
        call flow(esrsnk(iel),minerl(1,iel),time,egain(iel))
20    continue

      return
      end
