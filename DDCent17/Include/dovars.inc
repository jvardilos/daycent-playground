
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


!     Scheduling option variables
!     doaceq      time for accelerated equilibrium normalization
!     docult      time for a cultivation event
!     doerod      time for an erosion event
!     dofert      time for a fertilization event
!     dofire      time for a fire event
!                 (3) - grass/crop,forest,savannah
!     doflst      time for a forest last month growth event
!     dofone      time for a forest first month growth event
!     dofrst      time for a crop first month growth event
!     dograz      time for a grazing event
!     dohrvt      time for a harvest event
!     domirri     monthly irrigation event
!     dolast      time for a crop last month growth event
!     doomad      time for a organic matter addition event
!     doplnt      time for a planting event
!     dosene      do senescence event; set .true. when the actions are triggered
!     dotrem      time for a forest removal event
!     dofngw      forest fruit nut growth flag
!     dothrv      automatic tree fruit nut harvest
!
!     cultday     scheduled day of year for cultivation event
!     erodday     scheduled day of year for erosion event
!     fertday     scheduled day of year for fertilization event
!     fireday     scheduled day of year for fire event
!     flstday     scheduled day of year for forest last month growth event
!     foneday     scheduled day of year for forest first month growth event
!     frstday     scheduled day of year for crop first month growth event
!     grazday     scheduled day of year for grazing event
!     hrvtday     scheduled day of year for harvest event
!     irriday     scheduled day of year for irrigation event
!     lastday     scheduled day of year for crop last growth month event
!     omadday     scheduled day of year for organic matter addition event
!     plntday     scheduled day of year for planting event
!     seneday     scheduled day of year for senescence event
!     tremday     scheduled day of year for forest removal event
!
!     aceqcnt     number of accelerated equilibrium events possible (set by first call)
!     cultcnt     number of days that cultivation effect on decomposition
!                 has occurred
!     erodcnt     number of days that the erosion event has occurred
!     fertcnt     number of days since the fertilization event has occurred
!     grazcnt     number of days that the grazing event has occurred
!     irrcnt      number of days remaining for the irrigation event
!     irintvl     IRrigation interval in days
!                 > 0   timed irrigation in progress, 28 month default, 1 day default
!                   0   no irrigation
!                 < 0   growing season termination
!     plntcnt     number of days since the planting event has occurred
!     senecnt     number of days since the senescence event has occurred
!     frstschd    set to .true. when a FRST event is scheduled, .false. when
!                 the FRST event occurs
!     harvschd    set to .true. when a HARV event is scheduled, .false. when
!                 the harvest event occurs
!     plntschd    set to .true. when a PLNT event is scheduled, .false. when
!                 the planting event occurs
!     senmschd    set to .true. when a SENM event is scheduled, .false. when
!                 the senensence event occurs
!     savefrstday scheduled day of year for crop first month growth event as
!                 read from the schedule file
!     saveplntday scheduled day of year for planting event as read from
!                 schedule file
!     watertable  flag to determine if the watertable condition is off
!                 (watertable = 0) or on (watertable = 1)

      common/dovars/cultday, erodday, fertday, fireday, flstday,         &
     &              foneday, frstday, grazday, hrvtday, irriday,         &
     &              lastday, omadday, plntday, seneday, tremday,         &
     &              aceqcnt, cultcnt, erodcnt, fertcnt, grazcnt,         &
     &              irrcnt, irintvl, plntcnt, senecnt, savefrstday,      &
     &              saveplntday, watertable,                             &
     &              frstschd, harvschd, plntschd, senmschd,              & !logicals
     &              doaceq, docult, doerod, dofert, dofire, doflst,      &
     &              dofone, dofrst, dograz, dohrvt, domirri, dolast,     &
     &              doomad, doplnt, dosene, dotrem, dofngw, dothrv

      integer :: cultday, erodday, fertday, fireday, flstday
      integer :: foneday, frstday, grazday, hrvtday, irriday
      integer :: lastday, omadday, plntday, seneday, tremday

      integer :: aceqcnt ! ! don't run Eq sums unless this is >0
      integer :: cultcnt, erodcnt, fertcnt, grazcnt, irrcnt, irintvl
      integer :: plntcnt, senecnt, savefrstday, saveplntday
      integer :: watertable

      logical :: frstschd, harvschd, plntschd, senmschd
      logical :: doaceq, docult, doerod, dofert, dofire(3), doflst
      logical :: dofone, dofrst, dograz, dohrvt, domirri, dolast
      logical :: doomad, doplnt, dosene, dotrem, dofngw, dothrv

      save /dovars/
