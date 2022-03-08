
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

!       Modifications
!
!
!   moved forstrt = to frstage in plot 3    04/2017 K. Killian
!  12/2013 K. Killian
!   added forstrt = the time that forrest growth started for the Century orchard upgrade
!   converted tend, and blktnd to integer
!
!  6/2011  K. Killian
!   added an integer year counter to stabilize the time calculations
!   included t0par variables. the first 3 variables maintain the T0 order
!      common /t0par/ tend, dtpl, dt  !    t0par order
!   Left zztime alone since time is everywhere.


!****  integer time variables.
    ! cyear   = current year
    ! month   = current month
    ! strtyr  = starting year of the block
    ! tend    = simulation end time (always used as an integer)
    ! blktnd    = ending time of the block

! ... Time variables not declared by Time-0
    ! dt   = model time step (currently one month)
    ! dtpl = how often (years) the output should be written out
    ! strplt    = time (years) for first plot record to be written to plot file
    ! tplt      = next time (years) when variables should be written to
    !             the plot file
    ! decodt    = time step used for the decomposition model
    ! daylength = amount of daylight hours (1..24)




      common/timvar/tend, dtpl, dt, cyear, month, doy, strtyr, blktnd,        &
     &              strplt, tplt, decodt, daylength

      integer cyear, month, doy, strtyr, blktnd, tend
      real    dtpl, dt, strplt, tplt, decodt
      real    daylength(366)

      save /timvar/
