!               Copyright 1993 Colorado State University
!                         All Rights Reserved



      common /plot1/ vals1
      common /plot2/ vals2
      common /plot3/ vals3
      common /dummy/ vals4
      common /dummy/ vals5


      real :: vals1(467)
      real :: vals2(224)
      real :: vals3(159)
      real :: vals4(1)
      real :: vals5(1)
      save  /plot1/, /plot2/, /plot3/
