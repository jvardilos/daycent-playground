
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... CUTRTN
      subroutine cutrtn(accum)
      use calflow;

! ... Elemental return from a cutting event.

! ... Called from:  frem

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

! ... Argument declarations
      real      accum(ISOS)

! ... Local Variables
      integer   iel
      real      cgain, egain(MAXIEL), frc14, recres(MAXIEL)

! ... LEAVES are returned to LITTER
      if (rleavc .gt. 0.001) then
        cgain = remf(1) * retf(1,1) * rleavc
        if (cgain .gt. 0.0) then
          do iel = 1, nelem
            egain(iel) = remf(1) * retf(1,iel+1) * rleave(iel)
            recres(iel) = egain(iel) / cgain
          enddo
          frc14 = rlvcis(LABELD) / rleavc
          call partit(cgain,recres,1,csrsnk,esrsnk,wdlig(LEAF),frc14)
        endif
      endif

! ... FINE BRANCHES go to DEAD FINE BRANCHES
      if (fbrchc .gt. 0.001) then
        cgain = remf(2) * retf(2,1) * fbrchc
        call csched(cgain,fbrcis(LABELD),fbrchc,
     &              csrsnk(UNLABL),wd1cis(UNLABL),
     &              csrsnk(LABELD),wd1cis(LABELD),
     &              1.0,accum)
        do iel = 1, nelem
          egain(iel) = remf(2) * retf(2,iel+1) * fbrche(iel)
          call flow(esrsnk(iel),wood1e(iel),time,egain(iel))
        end do
      endif

! ... LARGE WOOD goes to DEAD LARGE WOOD
      if (rlwodc .gt. 0.001) then
        cgain = remf(3) * retf(3,1) * rlwodc
        call csched(cgain,rlwcis(LABELD),rlwodc,
     &              csrsnk(UNLABL),wd2cis(UNLABL),
     &              csrsnk(LABELD),wd2cis(LABELD),
     &              1.0,accum)
        do iel = 1, nelem
          egain(iel) = remf(3) * retf(3,iel+1) * rlwode(iel)
          call flow(esrsnk(iel),wood2e(iel),time,egain(iel))
        end do
      endif

! ... Add STORAGE back
      do iel = 1, nelem
        egain(iel) = remf(3) * retf(3,iel+1) * forstg(iel)
        call flow(esrsnk(iel),metabe(SRFC,iel),time,egain(iel))
      end do

      ! Fruits/Nuts are returned to LITTER
      cgain = remf(6) * retf(4,1) * frnutc
      if (cgain .gt. 0.0) then
        do iel = 1, nelem
          egain(iel) = remf(6) * retf(4,iel+1) * frnute(iel)
!         recres(iel) = egain(iel) / cgain
          recres(iel) = (retf(1,iel+1) * frnute(iel))/
     &                  (retf(4,1) * frnutc)
        end do
        frc14 = frncis(LABELD) / frnutc
        call partit(cgain,recres,1,csrsnk,esrsnk,wdlig(LEAF),frc14)
      endif

      return
      end
