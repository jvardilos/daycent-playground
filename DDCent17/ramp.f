
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... RAMP.F

      real function ramp(x, x1, y1, x2, y2)

      implicit none

      real      x, x1, y1, x2, y2 ! arguments

      !  This function models a "ramp":
      !              /-----
      !             /
      !       -----/

      if (x .le. x1) then
        ramp = y1
      else if (x .ge. x2) then
        ramp = y2
      else
        ramp = (y2 - y1) / (x2 - x1) * (x - x2) + y2
      endif

      return
      end function ramp
