
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... LIVREM

      subroutine livrem(accum)
      use calflow;

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'timvar.inc'
      include 'zztim.inc'

! ... Argument declarations
      real      accum(ISOS)

! ... Removal of live biomass due to cutting or fire in a forest.

! ... Called from:  frem

! ... Local variables
      integer :: iel
      real    :: closs, eloss, afdeath, alfc
      real    :: liveCtotal, liveCremoved
      real    :: mRespStorage
      real,dimension(2)  :: Cstgloss
      character :: ebffr*40

      liveCremoved = 0.0
      liveCtotal = fbrchc + rlwodc
      afdeath = 0 ! total above ground forest removed by this event.
      alfc   = (rlwodc + fbrchc + rleavc) ! save current aboveground C

! ... Remove live LEAVES

      if (rleavc .gt. 0.0) then
        closs = remf(1) * rleavc
        afdeath = afdeath + closs
        call csched(closs,rlvcis(LABELD),rleavc, &
                    rlvcis(UNLABL),csrsnk(UNLABL), &
                    rlvcis(LABELD),csrsnk(LABELD), &
                    1.0,accum)

        do iel = 1, nelem
          eloss = closs * (rleave(iel) / rleavc)
          terem(iel) = terem(iel) + eloss
          call flow(rleave(iel),esrsnk(iel),time,eloss)
        end do
      endif

!...Remove live FRUIT/NUT

      if (frnutc .gt. 0) then
        closs = remf(6) * frnutc
        afdeath = afdeath + closs
        call csched(closs,frncis(LABELD),frnutc, & ! should be fruit labeled frncis not leaf rlvcis  KLK
                     frncis(UNLABL),csrsnk(UNLABL), &
                     frncis(LABELD),csrsnk(LABELD), &
                     1.0,accum)

        do iel = 1, nelem
          eloss = closs * (frnute(iel) / frnutc)
          terem(iel) = terem(iel) + eloss
          call flow(frnute(iel),esrsnk(iel),time,eloss)
       end do
      endif

! ... Remove live FINE BRANCHES

      if (fbrchc .gt. 0.0) then
        closs = remf(2) * fbrchc
        afdeath = afdeath + closs
! ..... Add calculation of liveCremoved for maintenance respiration, mdh - 05/01
        liveCremoved = liveCremoved + closs
        call csched(closs,fbrcis(LABELD),fbrchc, &
                    fbrcis(UNLABL),csrsnk(UNLABL), &
                    fbrcis(LABELD),csrsnk(LABELD), &
                    1.0,accum)

        do iel = 1, nelem
          eloss = closs * (fbrche(iel) / fbrchc)
          terem(iel) = terem(iel) + eloss
          call flow(fbrche(iel),esrsnk(iel),time,eloss)
        end do
      endif

! ... Remove live LARGE WOOD

      if (rlwodc .gt. 0.0) then
        closs = remf(3) * rlwodc
        afdeath = afdeath + closs
! ..... Add calculation of liveCremoved for maintenance respiration, mdh - 05/01
        liveCremoved = liveCremoved + closs
        call csched(closs,rlwcis(LABELD),rlwodc, &
                    rlwcis(UNLABL),csrsnk(UNLABL), &
                    rlwcis(LABELD),csrsnk(LABELD), &
                    1.0,accum)

        do iel = 1, nelem
          eloss = closs * (rlwode(iel) / rlwodc)
          terem(iel) = terem(iel) + eloss
          call flow(rlwode(iel),esrsnk(iel),time,eloss)
        end do

! ..... Remove C from storage pool based on fraction of live wood removed, mdh - 7/9/01
        if (carbostg(FORSYS,UNLABL) .lt. 0.0  .or.  carbostg(FORSYS,LABELD) .lt. 0.0) then
          write(ebffr,*) carbostg(FORSYS,:)
          call abortrun('carbostg(FORSYS,:) = '//trim(ebffr)//' < 0. in livrem')
        endif
        ! remove lower limit on storage so we don't leave residual stored C
        if (liveCremoved .gt. 0.0  .and.  liveCtotal .gt. 0.0) then
          ! underflow control; replace intermediate csched calculations with resulting flow calls.
          if(liveCremoved >= liveCtotal) then;
            Cstgloss = carbostg(forsys,:);
          else;
            Cstgloss = carbostg(forsys,:) * (liveCremoved / liveCtotal);
          endif;
          call flow(carbostg(FORSYS,UNLABL), csrsnk(UNLABL), time, Cstgloss(UNLABL))
          call flow(carbostg(FORSYS,LABELD), csrsnk(LABELD), time, Cstgloss(LABELD))
          accum = accum + Cstgloss

          ! mRespStorage = sum(carbostg(FORSYS,:))
          ! mrspstgLoss = mRespStorage * (liveCremoved / liveCtotal)
          ! call csched(mrspstgLoss, carbostg(FORSYS,LABELD),mRespStorage, &
          !             carbostg(FORSYS,UNLABL), csrsnk(UNLABL), &
          !             carbostg(FORSYS,LABELD), csrsnk(LABELD), &
          !             1.0, accum)
        endif

! ..... Remove from STORAGE pool based on fraction of large wood removed.

        do iel = 1, nelem
          eloss = MAX(remf(3) * forstg(iel), 0.0)
          call flow(forstg(iel),esrsnk(iel),time,eloss)
        end do
      endif

      tcrem = tcrem + afdeath ! save the C removed to the accumulator
      if(afdeath .gt. 0.8*alfc) frstage = 0 ! reset orchard age if we remove to much

      return
      end


      subroutine trehrvst(farvfrac,sumharvc, sumharve)
      use calflow;

      implicit none

! ... Fruit nut harvest routine. It returns a specified fraction of fruit/nut
!     component to the sinks. This is a specialized live tree removal
!     event/routine and is included in that source file.
!
!     This duplicates the live removal found in liverem. It may be possible
!     to replace that code with this routine since this inputs the fraction
!     and sum as formal parameters
!
!
! Modifications
!   06/2007  K. Killian
!     written  based on leaf removal code from in livrem

      include 'const.inc'
      include 'parfs.inc'
      include 'param.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

! ... Argument declarations
      real    farvfrac, sumharvc, sumharve(MAXIEL)

! ... Local variables
      integer iel
      real    closs, eloss, accum(ISOS)

!     clear the dummy accum; Prevent possible IEEE floating errors
      accum = 0.;

!...Remove live FRUIT/NUT
        closs = frnutc
        if(closs .gt. 0.0001) closs = farvfrac * frnutc

        if (closs .gt. 0) then
          sumharvc = sumharvc + closs
          call csched(closs,frncis(LABELD),frnutc, &
                      frncis(UNLABL),csrsnk(UNLABL), &
                      frncis(LABELD),csrsnk(LABELD), &
                      1.0,accum)

          do iel = 1, nelem
            eloss = closs * (frnute(iel) / frnutc)
            sumharve(iel) = sumharve(iel) + eloss
            call flow(frnute(iel),esrsnk(iel),time,eloss)
          end do
        endif

      end
