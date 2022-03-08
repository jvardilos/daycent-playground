
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine esched(cflow,tca,rcetob,anps,bnps,labile,mnrflo,etob)
      use calflow;

      implicit none
      include 'zztim.inc'

! ... Argument declarations
      real     cflow, tca, rcetob, anps, bnps, labile, mnrflo, etob

! ... Schedule N, P, or S flow and associated mineralization or immobilization
      ! for decomposition from Box A to Box B.
      ! written by vek 05/91
      ! added outofa to simplify equilibrium acceleration

! ... Input:
      ! cflow  = C flow from Box A to Box B
      ! tca    = total C (unlabeled + labeled) in Box A
      ! rcetob = C/N, C/P, or C/S ratio of new material being added
      !          to Box B
      ! time   = simulation time, passed in /zztim/

! ... Transput:
      ! anps   = N, P, or S state variable for Box A
      ! bnps   = N, P, or S state variable for Box B
      ! labile = minerl(1,iel) where iel indicates N, P, or S

! ... Output:
      ! mnrflo = amount of N, P, or S that is mineralized or immobilized.
      !          Positive value indicates mineralization negative for immobilization.
      ! outofa = N, P, or S removed from box A.
      !          flow into B is outofa -mnrflo

! ... Local variables
      real     outofa

! ... Compute and schedule N, P, and S flows

! ... N, P, or S flowing out of Box A is proportional to C flow.
      outofa = anps * (cflow/tca) ! tca > 0 because of restrictions in calling routines
      etob   = cflow/rcetob       !
      mnrflo = outofa - etob

      ! prevent a 0/0 error from which happens with Microcosm option  (mse 2/95)
      ! Change from .AND. to .OR. could still get an underflow. -rm 1/96
      if (cflow .le. 0.0  .or.  outofa .le. 0.0                          &
      ! Removed code that turns of decomposition E flow when labile would E fall below zero, cak - 01/26/2011
      !    .or. labile + mnrflo .lt. 0.0 ! prevent immobilization if there is not enough mineral E, cak - 09/10/02
     &   ) then
        mnrflo = 0.0
        outofa = 0.0
        etob   = 0.0


      ! else if () then
      !
      !   mnrflo = 0.0
      !   goto 999
      ! endif

      ! IMMOBILIZATION occurs If C/E of Box A > C/E of material entering Box B
      else if (mnrflo < 0.0) then ! was (cflow/outofa .gt. rcetob)
        ! the amount of E required < 0. The extra E will be immobilized from the mineral pool,

        call flow(anps,bnps,time,outofa)    ! flow from Box A to Box B (outofa)

        call flow(labile,bnps,time,-mnrflo) ! flow from mineral pool to Box B (-mnrflo)


      else
! ..... MINERALIZATION occurs
        call flow(anps,bnps,time,etob)      ! flow from Box A to Box B (outofa)

        call flow(anps,labile,time,mnrflo)  ! flow from Box A to mineral pool

      endif

      return
      end
