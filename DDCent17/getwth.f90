!               Copyright 1993 Colorado State University
!                       All Rights Reserved


!***********************************************************************
!**
!**  FILE:     getwth.f
!**
!**  PURPOSE:  Retrieve a day's worth of weather for the weather file.
!**            Compute weekly average temperature, weekly min and max
!**            temperature, weekly pet, weekly pecip, and weekly soil
!**            surface temperature using circular arrays over a 7 day
!**            period.  Compute monthly average temperature using a
!**            circular array over a 30 day period.
!**
!**  This routine was developed for the RAMS / Daily Century linkage
!**     weather year = wthrstrt + MOD((yr - strtyr),(weather end - weather start)+1)
!**
!**  History
!**  written Melannie D. Hartman 12/5/96
!**  Add more robust checking for valid weather data values. CAK - 04/09/01
!**  most average array calculations use Fortran array intrinsics KLK Dec2011
!**  changed error messages to use message system KLK Dec2011
!**  moved EOF processing after return so you don't have to jump over it KLK Dec2011
!**  Added independent EVI file  KLK Dec2011
!**  Corrected EVI read bug;  EVI file supersedes EVi data in weather file  KLK Dec2012
!**  add more date information to error messages to better identify line KLK 13 Mar 2014
!**  Prevent an infinite loop trying to read an empty weather file  KLK 29 Sept 2015
!**  allow synching of file to simulation and changed input to a single read and
!**    a value parsing.  KLK 11 Mar 2016
!**
!**  INPUTS:
!**     curday       - the day of the year to read weather from the file
!**     cyear        - current simulation year for year search
!**     fwloss       - as read from the fix.100 file
!**     month        - current month of the year (1..12)
!**     precscalar   - monthly precipitation scalar values
!**     sitlat       - site latitude in decimal degrees as read from
!**                    <site>.100 file
!**     snow         - current snowpack (equiv. cm H2O)
!**     tmaxscalar   - monthly maximum temperature scalar values
!**     tminscalar   - monthly minimum temperature scalar values
!**     tmn2m        - average minimum air temperature for the month
!**                    as read from <site>.100 file (deg C - 2m)
!**     tmx2m        - average maximum air temperature for the month
!**                    as read from <site>.100 file (deg C - 2m)
!**     wthinput     - flag to indicate weather scalars, if any, to be
!**                    applied to data read from weather data file:
!**                      0 - No scalars used
!**                      1 - Use scalars for minimum temperature only
!**                      2 - Use scalars for maximum temperature only
!**                      3 - Use scalars for both minimum and maximum
!**                          temperatures
!**                      4 - Use scalars for precipitation only
!**                      5 - Use scalars for minimum and maximum
!**                          temperatures and precipitation
!**
!**  OUTPUTS:
!**     (From weather file):
!**     tempmax - maximum air temperature for the day (deg C)
!**     tempmin - minimum air temperature for the day (deg C)
!**     avgtemp - average air temperature for the day (deg C)
!**     ppt     - precipitation for the day (cm)
!**     solrad  - total incoming shortwave radiation (langleys/day)
!**     srad    - total incoming shortwave radiation (W/m^2/day)
!**     rhumid  - average relative humidity for the day (% 1..100)
!**     windsp  - average daily windspeed at 2 meters (mph)
!**
!**     tmaxwk  - average of tempmax for the last 7 days (deg C)
!**     tminwk  - average of tempmin for the last 7 days (deg C)
!**     tavewk  - average of avgtemp for the last 7 days (deg C)
!**     pptwk   - total precip for the last 7 days (cm)
!**     petwk   - total PET for the last 7 days (cm H2O)
!**     tavemth - average of avgtemp for the last 30 days (deg C)
!**     petdly  - potential evapotranspiration rate for day (cm H2O)
!**
!**  Called by:  simsom.f
!**
!**  Calls:  getevi    a routine to read a stand alone EVI file
!**
!**
!** usexdrvrs   data set
!** 0    ndy,nmth,nyr,njday,tempmax,tempmin,ppt
!** 1    ndy,nmth,nyr,njday,tempmax,tempmin,ppt, solrad, rhumid, windsp
!** 2    ndy,nmth,nyr,njday,tempmax,tempmin,ppt, srad
!** 3    ndy,nmth,nyr,njday,tempmax,tempmin,ppt, solrad, rhumid, windsp, srad
!** 4    ndy,nmth,nyr,njday,tempmax,tempmin,ppt, srad,  evi

!***********************************************************************

    subroutine getwth(curday, month, cyear, tempmax, tempmin, avgtemp, ppt,   &
                      solrad, rhumid, windsp, tavewk, petdly, fwloss,  &
                      sitlat, snow, tmn2m, tmx2m, tavemth, wthinput,   &
                      precscalar, tmaxscalar, tminscalar, srad, & !, usexdrvrs
                      leapyr, yrmatch, wstrtyr)

    USE ISO_C_Binding

    implicit none
    include 'const.inc'
    include 'evivars.inc'
    include 'jday.inc'
    include 'sitsoil.inc'

    ! Formal parameters

    integer :: curday, month, cyear, wthinput, wstrtyr !, usexdrvrs
    real    :: petdly, tavemth, tavewk
    real    :: sitlat
    real    :: snow
    real, dimension(NDAY+1) :: tempmax,tempmin,avgtemp,ppt,solrad,rhumid,windsp
    real, dimension(4) :: fwloss
    real, dimension(NMONTH) :: tmn2m, tmx2m, precscalar, tmaxscalar, tminscalar
    double precision srad(NDAY+1)
    logical :: leapyr, yrmatch

    ! Local Variables

    integer          :: istat, lch
    integer, save    :: ndy, nyr, njday, nmth
    integer, save    :: arrayindx = 1, wkarrayindx = 1, rwindcnt = 0
    real             :: value, tmaxwk, tminwk, pptwk, petwk
    real, save       :: tmaxarray(7),  tminarray(7), tavearray(7)
    real, save       :: temparray(30), pptarray(7),  petarray(7)
    double precision :: dailyPrecip, maxTemp, minTemp, tdew
    logical, save    :: startofrun = .true.,  debug = .false.
    logical          :: rewound
    character (len=128) :: bffr
    character (len=12)  :: yearstr

    ! functions
    real getevi
    integer             :: ReadInt
    double precision    :: ReadReal

    ! If it is the start of the run initialize all of the values in the
    ! temparray and the arrays being used to calculate the weekly averages.
    if (startofrun) then
      temparray = (tmn2m(month) + tmx2m(month)) / 2.0
      tavearray = (tmn2m(month) + tmx2m(month)) / 2.0
      tmaxarray = tmx2m(month)
      tminarray = tmn2m(month)
      pptarray  = 0.0
      petarray  = 0.0
      startofrun = .false.
    endif

    eviFlag = .false. ! default assumption is there is no EVI data
    inquire(unit=9, OPENED = rewound)
    if(.not. rewound) call abortrun("rewinding weather file in getwthr")
    rewound = .false. ! this read did not trigger a file rewind

    ! 0  ndy,nmth,nyr,njday,tempmax,tempmin,ppt
    ! 1  ndy,nmth,nyr,njday,tempmax,tempmin,ppt, solrad, rhumid, windsp
    ! 2  ndy,nmth,nyr,njday,tempmax,tempmin,ppt, srad
    ! 3  ndy,nmth,nyr,njday,tempmax,tempmin,ppt, solrad, rhumid, windsp, srad
    ! 4  ndy,nmth,nyr,njday,tempmax,tempmin,ppt, srad,   evi

    ! Restart here if we found an end of the weather file

     if(yrmatch  .and.  cyear .eq. wstrtyr  .and.  curday .eq. 1) then
       call message('synchronizing weather to simulation year; rewind file')
       rewind(9)
       yrmatch = .false.
     endif

10  continue
    ! initialize unread extra drivers
    solrad(curday) = -999.0
    rhumid(curday) = -999.0
    windsp(curday) = -999.0
    srad(curday)   = -999.0
    eviday         = -999.0

    lch =0
    call READCLIN(9, bffr, lch, .TRUE.) !
    if(lch .lt. 0) go to 20

    ! Read the date fields
                      ndy   = ReadInt(0,bffr,lch,istat) ! store date data
     if(istat .gt. 0) nmth  = ReadInt(0,bffr,lch,istat) ! store date data
     if(istat .gt. 0) nyr   = ReadInt(0,bffr,lch,istat) ! store date data
     if(istat .gt. 0) njday = ReadInt(0,bffr,lch,istat) ! store date data
     if(istat .le. 0) call abortrun('parsing weather date fields:'//trim(bffr))

     if (njday .ne. curday) then
       write(bffr,'(i4,a,4i4)') curday,' got day ', ndy,nmth,nyr,njday
       call abortrun('Expected day '//trim(bffr)//' in weather file.')
     endif

     if (nmth .ne. month) then
       if(nmth .eq. 3  .and.  njday .eq. 60) then
         write(bffr,*) ndy,nmth,nyr,njday
         call abortrun('missing leap day for year '//trim(bffr)//' in weather file')
       else
         write(bffr,*) nyr,njday,' was in month ',month,' not ',nmth
         call abortrun('Day '//trim(bffr)//' in weather file')
       endif
     endif

     ! log a rewind
     if(rewound  .and.  rwindcnt .lt. 3) then
       write(yearstr,*) cyear; yearstr = adjustl(yearstr);
       call message('EOF on weather file reading '//trim(yearstr)//'; rewinding to '//bffr(:lch-1))
       rwindcnt = rwindcnt +1
       if(rwindcnt .eq. 3) call message('further rewinds will be ignored')
     endif

     ! Check for leap year; It's dependent on the weather data not the run date
     if (curday .eq. 1) then
       ! first reset calendar arrays if last year was a leap year
       if (leapyr) then
         dysimo(2)    = dysimo(2) -1
         frstdy(3:12) = frstdy(3:12) -1
         lstdy(2:12)  = lstdy(2:12) -1
       endif
       ! Determine if the weather date is a leap year
       if ((mod(nyr,400) .eq. 0)) then
         leapyr = .TRUE.
       else if ((mod(nyr,100) .eq. 0)) then
         leapyr = .FALSE.
       else if ((mod(nyr,4) .eq. 0)) then
         leapyr = .TRUE.
       else
         leapyr = .FALSE.
       endif
       if (leapyr) then
         dysimo(2)    = dysimo(2)+1
         frstdy(3:12) = frstdy(3:12)+1
         lstdy(2:12)  = lstdy(2:12)+1
       endif
     endif

    ! Read the required weather fields
                       tempmax(curday) = ReadReal(0,bffr,lch,istat)
      if(istat .gt. 0) tempmin(curday) = ReadReal(0,bffr,lch,istat)
      if(istat .gt. 0) ppt(curday)     = ReadReal(0,bffr,lch,istat)
      if(istat .le. 0) call abortrun('parsing weather temp maximum, minimum or precip: '//trim(bffr))

  ! Checks for valid weather data
    if (tempmax(curday) .le. -99.0) then
      if (debug) call message('Warning: missing maximum temperature, day '//trim(bffr))
      tempmax(curday) = tmx2m(month)
    endif
    if (tempmin(curday) .le. -99.0) then
      if (debug) call message('Warning: missing minimum temperature, day '//trim(bffr))
      tempmin(curday) = tmn2m(month)
    endif
    if (ppt(curday) .le. -99.0) then
      if (debug) call message('Warning:  missing precipitation data, day '//trim(bffr))
      ppt(curday) = 0
    endif
    if (tempmax(curday) .lt. tempmin(curday)) then
      call message('Warning: Correcting invalid weather, tmax < tmin, day '//trim(bffr))
      tempmax(curday) = tmx2m(month)
      tempmin(curday) = tmn2m(month)
    endif

    ! Use extra weather drivers for PET calculations
    if (usexdrvrs .eq. 1) then
                       solrad(curday)  = ReadReal(0,bffr,lch,istat)
      if(istat .gt. 0) rhumid(curday)  = ReadReal(0,bffr,lch,istat)
      if(istat .gt. 0) windsp(curday)  = ReadReal(0,bffr,lch,istat)
      if(istat .le. 0) call abortrun('parsing type 1 drivers solrad, rhumid or windsp: '//trim(bffr))
      if ((solrad(curday) .le. -99.0) .or. &
          (rhumid(curday) .le. -99.0) .or. &
          (windsp(curday) .le. -99.0)) then
        call abortrun('Invalid weather driver data: '//trim(bffr))
      endif

  ! Read srad
    elseif (usexdrvrs .eq. 2) then
      srad(curday)    = ReadReal(0,bffr,lch,istat)
      if(istat .le. 0) call abortrun('parsing type 2 drivers srad: '//trim(bffr))

  ! Use extra weather drivers for PET and read srad
    elseif (usexdrvrs .eq. 3) then
                       solrad(curday)  = ReadReal(0,bffr,lch,istat)
      if(istat .gt. 0) rhumid(curday)  = ReadReal(0,bffr,lch,istat)
      if(istat .gt. 0) windsp(curday)  = ReadReal(0,bffr,lch,istat)
      if(istat .gt. 0) srad(curday)    = ReadReal(0,bffr,lch,istat)
      if(istat .le. 0) call abortrun('parsing type 3 drivers solrad, rhumid, windsp, srad: '//trim(bffr))
      if ((solrad(curday) .le. -99.0) .or. &
          (rhumid(curday) .le. -99.0) .or. &
          (windsp(curday) .le. -99.0)) then
        call abortrun('Invalid weather driver data: '//trim(bffr))
      endif

  ! Read weather file with EVI and solar radiation
    else if (usexdrvrs .eq. 4) then
                       srad(curday)    = ReadReal(0,bffr,lch,istat)
      if(istat .gt. 0) eviday          = ReadReal(0,bffr,lch,istat)
      if(istat .le. 0) call abortrun('parsing included EVI srad, eviday: '//trim(bffr))
      srad(curday) = -999.0 ! For now ignore solar radiation values read from file

    else if (usexdrvrs .ne. 0) then
      ! No extra weather drivers used
      write(bffr,'(i3)') usexdrvrs
      call abortrun('unknown weather file drivers '//trim(bffr))
    endif

    ! is there an external EVI file?
    value =  getevi('')
    if (value .gt. 0) eviday = value
    eviFlag = (eviday .gt. 0.0)  ! set eviFlag if we read EVI


  ! Apply weather scalars to the climate values as indicated, cak - 10/18/05
    if (wthinput .eq. 4 .or. wthinput .eq. 5) then
      ppt(curday) = ppt(curday) * precscalar(month)
    endif
    if (wthinput .eq. 2 .or. wthinput .eq. 3 .or. wthinput .eq. 5) then
      tempmax(curday) = tempmax(curday) + tmaxscalar(month)
    endif
    if (wthinput .eq. 1 .or. wthinput .eq. 3 .or. wthinput .eq. 5) then
      tempmin(curday) = tempmin(curday) + tminscalar(month)
    endif
  ! Check that application of scalars to temperature values does
  ! not result in invalid weather values
    if (tempmin(curday) .gt. tempmax(curday)) then
      write(bffr,*) tempmin(curday), ' greater than scaled maximum: ',tempmax(curday), &
       ' on day, month, year: ',ndy,nmth,nyr
      call message('Warning: scaled minimum temperature: '//trim(bffr))
      call message('Setting maximum temperature to scaled minimum'// &
                   'and minimum temperature 1 degree lower.')

      tempmax(curday) = tempmin(curday)
      tempmin(curday) = tempmax(curday) - 1.0
    endif

    avgtemp(curday) = (tempmax(curday) + tempmin(curday)) / 2.0

  ! If the solar radiation value has not been read from the weather data file
  ! calculate the solar radiation for the day using Peter Thornton's code, cak - 06/18/2009
    if (srad(curday) .le. -99.0) then
      minTemp = tempmin(curday)
      maxTemp = tempmax(curday)
      dailyPrecip = ppt(curday)
      tdew = 0.0
      call calc_srad(maxTemp, minTemp, dailyPrecip, curday, tdew,srad(curday))
      ! Adjust the total daily radiation value returned from Peter Thornton's code
      ! to reflect the the transmission coefficient and cloud cover at the site.
      srad(curday) = srad(curday) * sradadj(month)
    endif

    call calcpet(curday, month, tempmin(curday), tempmax(curday), &
                 avgtemp(curday), solrad(curday), rhumid(curday), &
                 windsp(curday), snow, usexdrvrs, fwloss, sitlat, &
                 petdly)

  ! Compute weeky average values over the last 7 days, the last 7 days
  ! worth of values are stored in the circular arrays, cak - 11/07/2208
    tmaxarray(wkarrayindx) = tempmax(curday)
    tminarray(wkarrayindx) = tempmin(curday)
    tavearray(wkarrayindx) = avgtemp(curday)
    pptarray(wkarrayindx)  = ppt(curday)
    petarray(wkarrayindx)  = petdly
    wkarrayindx = wkarrayindx + 1
    if (wkarrayindx .gt. 7) wkarrayindx = 1

  ! Code added to compute average temperature over the last 30 days,
  ! store the last 30 days worth of temperatures in the circular
  ! array, cak - 06/10/02
    temparray(arrayindx) = avgtemp(curday)
    arrayindx = arrayindx + 1
    if (arrayindx .gt. 30) arrayindx = 1

  ! Compute weekly values over the last 7 days, cak - 11/07/2007
      tmaxwk = sum(tmaxarray) / 7.0
      tminwk = sum(tminarray) / 7.0
      tavewk = sum(tavearray) / 7.0
      pptwk  = sum(pptarray)
      petwk  = sum(petarray)

  ! Code added to compute average temperature over the last 30 days,
  ! this value is used in the grochk subroutine for determining when
  ! leaf out should start, cak - 06/10/02
    tavemth = sum(temparray) / 30.0

    return

    ! If necessary, start reading the weather file again from the beginning
20    rewind(9)
      ! abort if the this call has rewound the file already; prevents a rewind/reread infinite loop
      if(rewound) call abortrun('weather file appears empty;')
      rewound = .true.
      goto 10
    end


    real function getevi(evifil)
    !
    ! this routine is special code to read EVI only files
    !
    ! The routine looks for an EVI with the specified name.
    ! the null name, "", causes the routine to read the file
    !
    ! The EVI data is tied to specific dates. It does NOT rewind like
    ! the weather file. Missing dates will default to -99.9 the PRDX mode flag.
    !
    ! the format of these files is month, year and a month's worth of EVI
    ! As defined RAM for the inventory.
    ! Missing months or years will be filled in with missing data
    ! Short records will be filled in with missing data
    !
    ! Modified 19 Nov 13 KLK
    !    recoded to ensure that EVI doesn't start until proper date
    !    Ensure that missing data, prefile, missing records, or closed file returns
    !    a negative value indicating no data.

    include 'timvar.inc'

    character (len=*)   :: evifil

    logical, save       :: eviopn = .FALSE.      ! evi file open
    integer, save       :: eyear = 0, emonth = 0 ! evi dates in buffer
    integer, save       :: eday                  ! days read in year
    integer, save       :: lch, istat            ! file read status
    integer, parameter  :: elu=16
    character (len=324),save :: bffr

    ! Functions
    double precision ReadReal
    integer          ReadInt

    ! close the evi unit
    if(evifil .eq. 'close'  .and.  eviopn) then
      close(elu,IOSTAT=istat)
      call message('Closed EVI data file ') ! log close
      getevi = istat
      eviopn = .FALSE. ! EVI file closed
      return

    ! run a set up on the first call
    else if(evifil .ne. '') then
      open(unit=elu, file=evifil, status='OLD', IOSTAT=istat)
      getevi = istat
      if(istat .eq. 0) then
        call message('Opened EVI data file '//evifil) ! log open
        eviopn = .TRUE. ! we opened an EVI file
      endif
      eyear = cyear -1 ! set the evi year to force a read
      return
    endif

    getevi = -999. ! default to no data

    ! check to see if an EVI file is open
    if(.not.eviopn) then
      getevi = -1.0 ! return no file

    else
      ! Do we need to buffer a new record
      if(cyear.gt.eyear .or. (cyear.eq.eyear .and. month.gt.emonth)) then
        ! we have stepped beyond the data, read another reord
10      lch   = 0
        istat = 0
        eyear = ReadInt(elu, bffr, lch, istat)
        ! this line doesn't look like data
        if(istat .eq. 0) then
          ! skip a comment ('#!') or heading line [Yy]ear;
          if(index('#!Yy', bffr(lch:lch)) .gt. 0) goto 10 !cycle to read a new buffer

          ! something else failed
          close(elu)         ! close the evi file
          call abortrun('parsing EVI file '//bffr(lch:)) !log and abort

        elseif(istat .lt. 0) then                  ! end of file
          close(elu)                               ! close the evi file
          eviopn = .FALSE.                         ! revert back to PRDX mode
          write(bffr,'(I2,a,i5)') month,"/",cyear  ! log change
          call message('Warning: EOF reading EVI for '//bffr(:11))
          getevi = -1.0 ! return no file
          return
        endif

        emonth = ReadInt(0, bffr, lch, istat) ! parse the month field
! write(*,*) "getevi read buffer", cyear,':',month,"    eyear ",eyear,':',emonth, lch, istat
        ! write(*,*) trim(bffr)

        !cycle to read a new buffer if we have data from earlier than
        ! the current simulation month
        if(eyear .lt. cyear  .or.  (cyear.eq.eyear .and. emonth.lt.month)) goto 10

        if(emonth .eq. 1) eday = 0 ! reset the evi day count
      endif

      ! pull a data value if the year and month match
      if(cyear.eq.eyear  .and.  month.eq.emonth) then
        eday   = eday+1
        getevi = real(ReadReal(0, bffr, lch, istat))
        ! write(*,*) "getevi", eyear, emonth, eday, getevi, lch, istat
        ! record parse error fill out this month with missing data
        if(istat .le. 0) getevi = -999.
      ! else
      ! To get here, the data must be in the future (cyear<eyear or month<emonth)
      ! if(cyear.lt.eyear .or. (cyear.eq.eyear .and. month.lt.emonth))
      endif
    endif

    return
    end
