
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... Forest removal variables
!       evntyp = type of event
!                0 for a cutting event (this includes storms)
!                1 for a fire event
!
!       fd(ipool) = fraction of fine and coarse roots that die
!
!       remf(ipool) = removal fraction
!                     ipool indicates which state variable
!
!       retf(ipool,iel) = return fraction
!                         ipool indicates which state variable
!                         iel indicates C, N, P, or S
!                         updated to include the fruit return

      integer evntyp
      real    fd(2), remf(6), retf(4,4)

      common /forrem/ fd, remf, retf, evntyp
      save /forrem/
