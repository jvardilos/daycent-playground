
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

!     Removed eup, It didn't need to be global.  KLK May 2015
!       Growth and treegrow used it for  E allocation during growth but both
!        routines cleared it; communicating no values between them that way.

      common/potent/agp,tgprod,pcropc,pforc,crop_a2drat,tree_a2drat

      real agp,tgprod,pcropc,pforc
      real crop_a2drat(3),tree_a2drat(3)

      save /potent/
