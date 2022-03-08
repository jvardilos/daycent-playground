
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


!   Variables used to track carbon flows from decomposition
!   identifying source and destination.
!
!   ametc1tosom11  - annual accumulator for carbon flow from surface
!                    metabolic pool to fast surface organic matter pool
!   ametc2tosom12  - annual accumulator for  carbon flow from soil
!                    metabolic pool to fast soil organic matter pool
!   asom11tosom21  - annual accumulator for carbon flow from fast
!                    surface organic matter pool to intermediate
!                    surface organic matter pool
!   asom12tosom22  - annual accumulator for carbon flow from fast soil
!                    organic matter pool to intermediate soil organic
!                    matter pool
!   asom12tosom3   - annual accumulator for carbon flow from fast soil
!                    organic matter pool to slow soil organic matter pool
!   asom21tosom11  - annual accumulator for carbon flow from
!                    intermediate surface organic matter pool to fast
!                    surface organic matter pool
!   asom21tosom22  - annual accumulator for carbon flow from
!                    intermediate surface organic matter pool to
!                    intermediate soil organic matter pool
!   asom22tosom12  - annual accumulator for carbon flow from
!                    intermediate soil organic matter pool to fast soil
!                    organic matter pool
!   asom22tosom3   - annual accumulator for carbon flow from
!                    intermediate soil organic matter pool to slow
!                    soil organic matter pool
!   asom3tosom12   - annual accumulator for carbon flow from slow soil
!                    organic matter pool to fast soil organic matter pool
!   astruc1tosom11 - annual accumulator for carbon flow from surface
!                    structural pool to fast surface organic matter pool
!   astruc1tosom21 - annual accumulator for carbon flow from surface
!                    structural pool to intermediate surface organic matter pool
!   astruc2tosom12 - annual accumulator for carbon flow from soil
!                    structural pool to fast soil organic matter pool
!   astruc2tosom22 - annual accumulator for carbon flow from soil
!                    structural pool to intermediate soil organic matter pool
!   awood1tosom11  - annual accumulator for carbon flow from dead
!                    fine branch pool to fast surface organic matter pool
!   awood1tosom21  - annual accumulator for carbon flow from dead fine branch
!                    pool pool to intermediate surface organic matter pool
!   awood2tosom11  - annual accumulator for carbon flow from dead large wood
!                    pool to fast surface organic matter pool
!   awood2tosom21  - annual accumulator for carbon flow from dead large wood
!                    pool to intermediate surface organic matter pool
!   awood3tosom12  - annual accumulator for carbon flow from dead coarse
!                    root pool to fast soil organic matter pool
!   awood3tosom22  - annual accumulator for carbon flow from dead coarse
!                    root pool to intermediate soil organic matter pool
!   metc1tosom11   - carbon flow from surface metabolic pool to fast
!                    surface organic matter pool
!   metc2tosom12   - carbon flow from soil metabolic pool to fast soil
!                    organic matter pool
!   som11tosom21   - carbon flow from fast surface organic matter pool
!                    to intermediate surface organic matter pool
!   som12tosom22   - carbon flow from fast soil organic matter pool
!                    to intermediate soil organic matter pool
!   som12tosom3    - carbon flow from fast soil organic matter pool to
!                    slow soil organic matter pool
!   som21tosom11   - carbon flow from intermediate surface organic
!                    matter pool to fast surface organic matter pool
!   som21tosom22   - carbon flow from intermediate surface organic
!                    matter pool to intermediate soil organic matter pool
!   som22tosom12   - carbon flow from intermediate soil organic matter
!                    pool to fast soil organic matter pool
!   som22tosom3    - carbon flow from intermediate soil organic matter
!                    pool to slow soil organic matter pool
!   som3tosom12    - carbon flow from slow soil organic matter pool to
!                    fast soil organic matter pool
!   struc1tosom11  - carbon flow from surface structural pool to fast
!                    surface organic matter pool
!   struc1tosom21  - carbon flow from surface structural pool to
!                    intermediate surface organic matter pool
!   struc2tosom12  - carbon flow from soil structural pool to fast
!                    soil organic matter pool
!   struc2tosom22  - carbon flow from soil structural pool to
!                    intermediate soil organic matter pool
!   wood1tosom11   - carbon flow from dead fine branch pool to fast
!                    surface organic matter pool
!   wood1tosom21   - carbon flow from dead fine branch pool to
!                    intermediate surface organic matter pool
!   wood2tosom11   - carbon flow from dead large wood pool to fast
!                    surface organic matter pool
!   wood2tosom21   - carbon flow from dead large wood pool pool to
!                    intermediate surface organic matter pool
!   wood3tosom12   - carbon flow from dead coarse root pool to fast
!                    soil organic matter pool
!   wood3tosom22   - carbon flow from dead coarse root pool to
!                    intermediate soil organic matter pool



      common/cflows/ametc1tosom11, ametc2tosom12,                        &
     &  asom11tosom21, asom12tosom22, asom12tosom3, asom21tosom11,       &
     &  asom21tosom22, asom22tosom12, asom22tosom3, asom3tosom12,        &
     &  astruc1tosom11, astruc1tosom21, astruc2tosom12, astruc2tosom22,  &
     &  awood1tosom11, awood1tosom21, awood2tosom11, awood2tosom21,      &
     &  awood3tosom12, awood3tosom22,                                    &
     &  metc1tosom11,  metc2tosom12,                                     &
     &  som11tosom21,  som12tosom22, som12tosom3, som21tosom11,          &
     &  som21tosom22,  som22tosom12, som22tosom3, som3tosom12,           &
     &  struc1tosom11, struc1tosom21, struc2tosom12, struc2tosom22,      &
     &  wood1tosom11,  wood1tosom21, wood2tosom11, wood2tosom21,         &
     &  wood3tosom12,  wood3tosom22

      real ametc1tosom11, ametc2tosom12
      real asom11tosom21, asom12tosom22, asom12tosom3,  asom21tosom11
      real asom21tosom22, asom22tosom12, asom22tosom3,  asom3tosom12
      real astruc1tosom11,astruc1tosom21,astruc2tosom12,astruc2tosom22
      real awood1tosom11, awood1tosom21, awood2tosom11, awood2tosom21
      real awood3tosom12, awood3tosom22
      real metc1tosom11,  metc2tosom12
      real som11tosom21,  som12tosom22,  som12tosom3,   som21tosom11
      real som21tosom22,  som22tosom12,  som22tosom3,   som3tosom12
      real struc1tosom11, struc1tosom21, struc2tosom12, struc2tosom22
      real wood1tosom11,  wood1tosom21,  wood2tosom11,  wood2tosom21
      real wood3tosom12,  wood3tosom22

      save /cflows/
