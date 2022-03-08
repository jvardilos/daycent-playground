
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine cultiv(doy, pltlig)
      use calflow;

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

      ! Argument declarations
      integer :: doy
      real      pltlig(3)

! ... Implement cultivation option
      ! modified Nov  13 KLK:     storage removal uses existing variables.
      !    Corrected a bug where uninitialized intermediate variables blew up
      !    the metabolic E
      ! modified Sept 13 KLK:     commented storage removal from aglivc
      ! modified Nov 12 KLK
      ! changed storage removal from aglivc loss to root loss to be consistant
      ! with convention in harvst. This returns stored E to both surface and soil
      ! metabolic. Helps fix a bug where crpstg is non-zero with no live crop/grass.

      ! Local variables
      integer   iel
      real   :: accum(ISOS), fr14, recres(MAXIEL), tagsfc, tagsoi
      real   :: trans, tsdsfc, tsdsoi, surfc
      real   :: gclivtc ! grass/crop live before cultivation
      real, parameter :: clowlim = 1.0e-5;

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

      ! Update flows and sum carbon isotopes
      call flowup(time)
      call sumcar

      gclivtc = 0 ! how much live to we have

      ! standing dead going into surface litter
      tsdsfc = stdedc * cultra(4)
      if (tsdsfc .gt. 0.0) then
        recres(1:nelem) = stdede(1:nelem)/stdedc
        fr14 = stdcis(LABELD)/stdedc
        call partit(tsdsfc,recres,1,stdcis,stdede,pltlig(ABOVE),fr14)
      endif

      ! surface litter going into the top soil layer.

      ! Structural
      trans = strucc(SRFC) * cultra(6)
      if (trans .gt. 0.) then
        call csched(trans,strcis(SRFC,LABELD),strucc(SRFC), &
                    strcis(SRFC,UNLABL),strcis(SOIL,UNLABL), &
                    strcis(SRFC,LABELD),strcis(SOIL,LABELD), &
                    1.0,accum)

! ..... Recompute lignin fraction in structural soil C
        call adjlig(strucc(SOIL),strlig(SRFC),trans,strlig(SOIL))

        do 20 iel = 1, nelem
          trans = struce(SRFC,iel) * cultra(6)
          call flow(struce(SRFC,iel),struce(SOIL,iel),time,trans)
20      continue
      endif

      ! Metabolic
      trans = metabc(SRFC) * cultra(6)
      if (trans .gt. 0.) then
        call csched(trans,metcis(SRFC,LABELD),metabc(SRFC), &
                    metcis(SRFC,UNLABL),metcis(SOIL,UNLABL), &
                    metcis(SRFC,LABELD),metcis(SOIL,LABELD), &
                    1.0,accum)
        do 30 iel = 1, nelem
          trans = metabe(SRFC,iel) * cultra(6)
          call flow(metabe(SRFC,iel),metabe(SOIL,iel),time,trans)
30      continue
      endif

      ! Surface SOM1
      trans = som1c(SRFC) * cultra(6)
      if (trans .gt. 0.) then
        call csched(trans,som1ci(SRFC,LABELD),som1c(SRFC), &
                    som1ci(SRFC,UNLABL),som1ci(SOIL,UNLABL), &
                    som1ci(SRFC,LABELD),som1ci(SOIL,LABELD), &
                    1.0,accum)
        do 40 iel = 1, nelem
          trans = som1e(SRFC,iel) * cultra(6)
          call flow(som1e(SRFC,iel),som1e(SOIL,iel),time,trans)
40      continue
      endif

      ! Surface SOM2
      trans = som2c(SRFC) * cultra(6)
      if (trans .gt. 0.) then
        call csched(trans,som2ci(SRFC,LABELD),som2c(SRFC), &
                    som2ci(SRFC,UNLABL),som2ci(SOIL,UNLABL), &
                    som2ci(SRFC,LABELD),som2ci(SOIL,LABELD), &
                    1.0,accum)
        do 45 iel = 1, nelem
          trans = som2e(SRFC,iel) * cultra(6)
          call flow(som2e(SRFC,iel),som2e(SOIL,iel),time,trans)
45      continue
      endif

      ! standing dead going to the top soil layer.
      tsdsoi = stdedc * cultra(5)
      if (tsdsoi .gt. 0.0) then
        call partit(tsdsoi,recres,2,stdcis,stdede,pltlig(ABOVE),fr14)
      endif

      ! above ground live going into surface litter
      if (aglivc .gt. 0.0) then
        tagsfc = aglivc * cultra(2)
        gclivtc = gclivtc + tagsfc
        recres(1:nelem) = aglive(1:nelem)/aglivc
        fr14 = aglcis(LABELD)/aglivc
        call partit(tagsfc,recres,1,aglcis,aglive,pltlig(ABOVE),fr14)
      endif

!      Disable the aglivc code since storage pool is now associated with roots  KLK
!      if (aglivc .gt. 0.0) then
!        do iel = 1, nelem
!c         Some storage pool going to metabolic surface pool
!          trans = crpstg(iel) * cultra(2)
!          call flow(crpstg(iel),metabe(SRFC,iel),time,trans)
!c         Some storage pool going to metabolic soil pool
!          trans = crpstg(iel) * cultra(3)
!          call flow(crpstg(iel),metabe(SOIL,iel),time,trans)
!        enddo
!      endif

      ! above ground live going to the top soil layer.
      if (aglivc .gt. 0.0) then
        tagsoi = aglivc * cultra(3)
        gclivtc = gclivtc + tagsoi
        call partit(tagsoi,recres,2,aglcis,aglive,pltlig(ABOVE),fr14)
      endif

      ! Live roots go to the top soil layer.
      ! Juvenile fine roots
      trans = bglivcj * cultra(7)
      if (trans .gt. 0.0) then
        gclivtc = gclivtc + trans
        do iel = 1, nelem
          recres(iel) = bglivej(iel)/bglivcj
        enddo
        fr14 = bglcisj(LABELD)/bglivcj
! ..... A fraction of the live roots are transferred to the surface litter,
        ! the remainder going to the soil litter layer, cak - 05/14/2007
        surfc = trans * rdsrfc
        call partit(surfc, recres,SRFC,bglcisj, &
                    bglivej,pltlig(BELOWJ),fr14)
        call partit(trans-surfc , recres,SOIL,bglcisj, &
                    bglivej,pltlig(BELOWJ),fr14)
      endif
      ! Mature fine roots
      trans = bglivcm * cultra(7)
      if (trans .gt. 0.0) then
        gclivtc = gclivtc + trans
        recres(1:nelem) = bglivem(1:nelem)/bglivcm
        fr14 = bglcism(LABELD)/bglivcm
! ..... A fraction of the live roots are transferred to the surface litter,
        ! the remainder going to the soil litter    cak - 05/14/2007
        surfc = trans * rdsrfc
        call partit(surfc, recres,SRFC,bglcism, &
                    bglivem,pltlig(BELOWM),fr14)
        call partit(trans-surfc, recres,SOIL,bglcism, &
                    bglivem,pltlig(BELOWM),fr14)
      endif

      ! allocate storage in the dead roots to surface and soil metabolic
      if(bglivcj+bglivcm .gt. 0) then
        ! loop for flow commands
        do iel = 1, nelem
          trans = crpstg(iel) * cultra(7)
          surfc = trans * rdsrfc
          call flow(crpstg(iel), metabe(SRFC,iel), time, surfc)
          call flow(crpstg(iel), metabe(SOIL,iel), time, trans-surfc)
        end do
      endif

      ! above ground live going to standing dead
      trans = aglivc * cultra(1)
      if (trans .gt. 0.0) then
        gclivtc = gclivtc + trans
        call csched(trans,aglcis(LABELD),aglivc, &
                    aglcis(UNLABL),stdcis(UNLABL), &
                    aglcis(LABELD),stdcis(LABELD), &
                    1.0,accum)
        do 90 iel = 1, nelem
          trans = aglive(iel) * cultra(1)
          call flow(aglive(iel),stdede(iel),time,trans)
90      continue
      endif

      ! State variables and accumulators
      call flowup(time)

      ! Check status of carbon values to make sure that everything
      ! has been reset correctly; if carbon = 0, reset elements to 0 as well
      if (sum(aglcis) .lt. clowlim) then
        aglcis = 0.0;  ! aglivc = 0. ! C sum
        aglive = 0.0
      endif
      if (sum(bglcisj) .lt. clowlim) then
        bglcisj = 0.0; ! bglivcj = 0. ! C sum
        bglivej = 0.0
      endif
      if (sum(bglcism) .lt. clowlim) then
        bglcism = 0.0; ! bglivcm = 0. ! C sum
        bglivem = 0.0
      endif

      if (sum(stdcis) .lt. clowlim) then
        stdcis = 0.0
        stdede = 0.
      endif

      call sumcar ! do the summation AFTER we touched the components

      ! zombie control; stop growth if live is < 1.0e-5
      !  first ignore plant dqtes, we want the plant to grow and non-killing cultivations
      if(.not. (dofrst .or. doplnt .or.  gclivtc .lt. 1.0e-5)  &
         .and. &
          aglivc .le. 1.e-5 .and. (bglivcj + bglivcm) .le. 1.e-5) then ! no live remaining
        crpgrw = 0
        dolast = .true.
        lastday = doy
      endif

      return
      end
