
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      common/seq/cursys, decsys, otfrac

      integer cursys, decsys
      real    otfrac

! ... cursys = 1 grass/crop
! ...        = 2 forest
! ...        = 3 savanna

! ... decsys = 1 No wood is present to decompose
! ...        = 2 Wood is present to decompose

!...otfrac  sets a fixed savanna tree fraction for orchards  0 < 1.0

      save /seq/
