!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine wthini(cyear,doy)

      USE ISO_C_Binding

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'sitsoil.inc'   ! location of usexdrvrs instead of wthdaily
      include 'wthdaily.inc'

      integer             :: cyear,doy

      ! local variables
      character (len=256) :: bffr
      character (len=8)   :: cdate
      integer             :: lch, istat, i, ndci
      integer             :: ndom, nmth, nyr, ndoy
      logical             :: rewnd
      integer, parameter  :: xdrv(5) = (/0, 2, 4, 1, 3/)
      double precision    :: tmp

      ! functions
      integer             :: ReadInt
      double precision    :: ReadReal
      real                :: getevi


      if     (wthr .eq. 'C') return ! no nothing

! ... Open daily weather data file
      close(unit=9)                                      !original

      if(wthnam .eq. ' ') call abortrun('No weather file specified')

      open(unit=9,file=wthnam,status='OLD',IOSTAT=istat)     !original
      if(istat .ne. 0) then
        write(bffr, '(i4)') istat
        call abortrun('opening weather file error '//trim(bffr)//' file: '//wthnam)
      endif
      if(wthr .eq. 'Y') then
        write(cdate,'(i7)') cyear
        cdate = adjustl(cdate)
        call message('Searching for: '//trim(cdate)//' weather in '//trim(wthnam))
      else if(wthr .eq. 'F') then
        call message('Opened weather file: '//trim(wthnam))
      endif
      nyr = cyear-1
      ndoy    = -1
      rewnd   = .false.
      yrmatch = .FALSE.

 ! added a year search, 'Y' feature on opening the file     KLK 10 Mar 2016
 !   If the search fails, a warning is printed and the file will be rewound.
 !
 ! Originally, the remainder of this routine was just a rewind
 !           rewind 9
 ! The routine has been modified to skip comment headers then characterize the
 ! data based on the structure of the first line of weather data.
 ! This overrides weather driver input in sitepar.in    KLK  8Dec2011
 !
 ! call getevi with a file name similar to the weather file
    !   1 change the extension to .EVI
    !   2 append .evi to the existing name.
    ! It would be possible to modify this with a special name. It was not done
    ! now to simplify the DayCent interface.
 !
 !
 ! set usexdrvrs by counting the number of fields in the weather file
 ! usexdrvrs    format
 !  0           ndy,nmth,nyr,njday,tempmax,tempmin,ppt
 !  2           ndy,nmth,nyr,njday,tempmax,tempmin,ppt, srad
 !  4           ndy,nmth,nyr,njday,tempmax,tempmin,ppt, srad,   evi
 !  1           ndy,nmth,nyr,njday,tempmax,tempmin,ppt, solrad, rhumid, windsp
 !  3           ndy,nmth,nyr,njday,tempmax,tempmin,ppt, solrad, rhumid, windsp, srad


 !
 ! Here is code to skip leading comment lines. Commented until getwth does
 ! the same thing to prevent a rewind from blowing up.
100   continue
      lch   = 0
      istat = 0
      call READCLIN(9, bffr, lch, .TRUE.)
      if(lch .lt. 0) call abortrun('EOF on weather file :'//wthnam)
      lch = 1;

      do while (nyr .ne. cyear)
        ! read the date fields ndom, nmth, nyr, ndoy
        ndom = ReadInt(9,bffr,lch,istat) ! store date data in
        if(istat .le. 0) then
          if(ndoy .lt. -1  .or.  rewnd) then  ! first line. It's fatal if we can't parse this.
             call abortrun('parsing data in: '//trim(wthnam))
          elseif(ndoy >= 1) then ! probably an EOF; rewind and accept
             call message('I/O error (EOF?) reading weather file: rewinding')
             rewind(9)
             rewnd = .true.
             go to 100;
          endif
        endif

        nmth = ReadInt(9,bffr,lch,istat) ! store date data in
        if(istat .le. 0) call abortrun('parsing weather month:'//trim(wthnam))
        nyr  = ReadInt(9,bffr,lch,istat) ! store date data in
        if(ndoy .eq. -1) wstrtyr = nyr   ! save the file start year
        if(istat .le. 0) call abortrun('parsing weather year:'//trim(wthnam))
        ndoy = ReadInt(9,bffr,lch,istat) ! store date data in
        if(istat .le. 0) call abortrun('parsing weather DOY:'//trim(wthnam))

        ! look for the year if requested.
        if(wthr .eq. 'F' .or. rewnd) then
          nyr = cyear ! don't care; terminate the search loop
          if(wthr .eq. 'Y') call message('Warning weather for '//trim(cdate)// &
              ' not found; starting simulation from first year in: '//trim(wthnam))

        else if(wthr .eq. 'Y') then
          if    (ndoy .gt. 1) then ! Not the year start so grab another line
            lch = 0; cycle;        ! grab another year

          elseif(nyr < cyear) then ! file starts before block
            ! search forward
            read(9,'(363(/))', iostat = istat) ! skip 264 days of entries
            if(istat .ne. 0) call abortrun('searching weather dates:'//trim(wthnam))
            lch = 0; cycle;        ! grab another year

          elseif(nyr > cyear) then ! file starts before block
            yrmatch = .TRUE.       ! match year may be some time in the future
            exit;                  ! We have to start with the present data.
          endif
        endif
      enddo

      ! backup so the actual data read doesn't get confused
      backspace(9, iostat = istat)
      if(istat .ne. 0) call abortrun('reading weather:'//wthnam)


      ! determine extra drivers in the file
      ! step through common temp and precip tempmax, tempmin, ppt
      do i=1,3
        tmp = ReadReal(0,bffr,lch,istat)
        if(istat .le. 0) call abortrun('parsing weather data: '//trim(wthnam))
      end do

      ! count extra remaining data fields
      ! this is the critical part that determines usexdrvrs
      ndci = 1
      do while (lch .lt. len(bffr)  .and.  ndci .lt. 5)
        tmp = ReadReal(0,bffr,lch,istat)
        if(istat .le. 0) call abortrun('parsing extra weather drivers: '//trim(wthnam))
        ndci = ndci+1
      end do

      ! record the weather data description
      if(usexdrvrs .gt. 5  .or. usexdrvrs .lt. 0) then
        usexdrvrs = xdrv(ndci)    ! set extra drivers flag for fortran read
        write(bffr,'(i2)') usexdrvrs
        call message('setting weather drivers to type '//trim(bffr))
      else if(usexdrvrs .gt. 0) then
        i = 1
        do while (i .le. 5 .and. xdrv(i) .ne. usexdrvrs)
          i = i + 1
        end do
        if(i .gt. ndci) then
          write(bffr,'(i2)') usexdrvrs
          call abortrun('insufficient data for driver mode'// &
                        bffr(:2)//' in file'//trim(wthnam))
        end if
      end if

      if (usexdrvrs .eq. 1) then
        call message(' Using solrad, rhumid, and windsp for PET calc')
      elseif (usexdrvrs .eq. 2) then
        call message(' Using only air temp for PET calc; Reading srad')
      elseif (usexdrvrs .eq. 3) then
        call message(' Using solrad, rhumid, and windsp for PET calc; Reading srad')
      else
        call message(' Using only air temp for PET calc')
      endif

      ! open an EVI file related to the weather file
      i = index(wthnam,'.',.TRUE.) -1 ! look for an extension
      if(i.eq.0) i=len_trim(wthnam)   ! no extension then append
      tmp = dble(getevi(wthnam(:i)//'.evi')) !run the evi open code

      return
      end subroutine wthini
