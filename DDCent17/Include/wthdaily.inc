      common /wthdaily/ srad, avgtemp, tempmax, tempmin, ppt, solrad,    &
     &            rhumid, windsp, leapyr, yrmatch, wstrtyr !, usexdrvrs

! leapyr moved to weather data since it is determined by the weather NOT the clock
! yrmatch match weather file to clock
! wstrtyr clock start year
      integer :: wstrtyr !, usexdrvrs
      logical :: leapyr, yrmatch
      double precision srad(NDAY+1)

      real avgtemp(NDAY+1), tempmin(NDAY+1), tempmax(NDAY+1)
      real ppt(NDAY+1), solrad(NDAY+1), rhumid(NDAY+1), windsp(NDAY+1)

      save /wthdaily/
