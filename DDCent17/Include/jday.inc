!************************************************************************
!  jday.inc
!  For use with the daily version of century (in getwth & simsom)
!
!  Melannie D. Hartman
!  made dysimo a real; its always used as a real klk 2018/2/11
!
!************************************************************************

      common /jdays/  frstdy, lstdy, dysimo

      integer :: frstdy(12), lstdy(12)
      real    :: dysimo(12)
!      real    :: sradadj(12)

      integer, parameter :: ilstdy(12)  = (/31,59,90,120,151,181,212,       &
     &                                      243,273,304,334,365/)
      integer, parameter :: ifrstdy(12) = (/1,32,60,91,121,152,182,213,     &
     &                                      244,274,305,335/)

      save /jdays/

!      integer function idaymnth(m)
!        integer :: m
!        if(m=2) then; idaymnth = 28; else;
!                      idaymnth = 30+int(modulo(m,2) + modulo(abs(m-2)/6,2));
!        endif
!      end function idaymnth
