
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      logical function candec(nelem,minrlavail,tca,elstva,nlr,lyr,rcenew)

      implicit none
      include 'const.inc'

! ... Argument declarations
      integer   nelem, nlr, lyr
      real      minrlavail(MAXIEL), tca, elstva(nlr,MAXIEL), rcenew(MAXIEL)

! ... Determine if decomposition can occur.

! ... Input:
! ...   nelem       = number of elements
! ...   minrlavail  = mineral available for decomposition before uptake by plants
! ...   tca         = total C in Box A
! ...   elstva      = N, P, and S in Box A by layer and element
! ...   nlr         = 1st dimension of elstva... elstva(nlr,3)
! ...   lyr         = layer of Box A (1=surface, 2=soil)
! ...   rcenew(iel) = C/N, C/P, and C/S ratios of new material
! ...                 being added to Box B (iel=1,3)

! ... Output:
! ...   candec      = true if Box A can decompose to Box B,
! ...                 otherwise false

! ... When candec is called from declig, the array passed into rcenew
! ... has 6 elements, but candec only needs the first 3:  (1,1), (2,1),
! ... and (3,1)).  When candec is called from somdec, a 3-element array
! ... is passed into rcenew.

! ... Local variables
      integer  iel
      logical  cando(MAXIEL)

! ... Initialize cando to true for each array location.
      cando = .true.

! ... For each element (N, P, and S)
      do 10 iel = 1, nelem

! ..... If there is no available mineral E
        if (minrlavail(iel) .lt. 1.e-07) then

! ....... Compare the C/E of new material to the C/E of Box A if C/E of
! ....... Box A > C/E of new material
          if (tca/elstva(lyr,iel) .gt. rcenew(iel)) then

! ......... Immobilization is necessary and the stuff in Box A can't
! ......... decompose to Box B.
            cando(iel) = .false.
          endif
        endif
10    continue

! ... If there is some available mineral E, decomposition can
! ... proceed even if mineral E (minerl) is driven negative in
! ... the next time step.

      candec = (cando(N) .and. cando(P) .and. cando(S))

      return
      end
