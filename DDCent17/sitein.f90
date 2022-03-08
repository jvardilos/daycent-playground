
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... SITEIN.F

      subroutine sitein(extndopt,flgvrb, sitparread, soilread)
      USE ISO_C_Binding

! ... Read the parameter file in the form needed by time0.

      implicit none

      include 'const.inc'
      include 'doubles.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'sitsoil.inc'
      include 'wth.inc'

!     functions
      real       :: chkdata; external  chkdata
      logical    :: wthrstat

      ! Argument declarations
      integer    :: extndopt, flgvrb
      logical    :: sitparread, soilread

      ! Local variables
      integer    :: ii, jj
      real       :: stk(1), value
      ! pull parameters rather than load parcp.inc parfs.inc
      logical    :: sitpar=.TRUE.
      character (len=16)  :: varfnd
      character (len=132) :: buffr
      character (len=6), parameter :: callr = 'sitein'
      real :: rcesd(MAXIEL,4) =  RESHAPE( (/10.0,50.0,50.0, 10.0,50.0,50.0, 17.0,117.0,117.0,       &
                                   7.0,62.0,62.0/), (/MAXIEL,4/))

      stk(1) = 0

      ii = 0;
      read(7,*) buffr         ! skip the file header
      !read(7,'(a)') buffr    ! read climate section head
      call READCLIN(7,buffr,ii,.TRUE.) ! read with standard comment/blank line parsing

      ! check to see if we have read a binary file
      ! check for a '*** Climate' header or some variation there of
      if(buffr(:2) .ne. '**'  .and.  index(buffr,'Climate').eq.0) then
        ! missing: look for non printing characters
        jj = len_trim(buffr)
        do ii=1, jj
          if(LLT(buffr(ii:ii),' ') .or. LGT(buffr(ii:ii),'~')) call &
             FILRERR(7,0,'sitein reading non-text data?',buffr(:40))
        enddo
      endif

      ! chkdata to echo the echo flag
      if(flgvrb .gt. 0) value = chkdata(7,callr,'echo',0,stk)
      ! UPGRADE: do statistics from the site file
      ii = index(buffr, 'tatistics')
      if(ii .ne. 0) then
        ii = ii+9
        jj = scan(buffr(ii:), ' 	') ! search for space/tab
        buffr = adjustl(buffr(ii+jj:))

        call fndfilnam(buffr,ii) ! do we have a weather file name

        if(ii .gt. 0) then
          if(wthrstat(buffr, precip, prcstd, prcskw, tmn2m, tmx2m)) then
            if(flgvrb .gt. 0) then
              write(*,*) 'precip:',precip
              write(*,*) 'prcstd:',prcstd
              write(*,*) 'prcskw:',prcskw
              write(*,*) 'tmn2m: ',tmn2m
              write(*,*) 'tmx2m: ',tmx2m
            endif
            varfnd = ' '
            value = chkdata(7,callr,varfnd,0,stk)
            do while (index(varfnd, 'precip') .ne. 0 .or. index(varfnd, 'prcstd') .ne. 0 .or. &
                     index(varfnd, 'prcskw') .ne. 0 .or. &
                     index(varfnd, 'tmn2m') .ne. 0 .or. index(varfnd, 'tmx2m') .ne. 0)
              varfnd = ' '
              value = chkdata(7,callr,varfnd,0,stk)
            enddo
            value = chkdata(7,callr,'reread',0,stk)
          else
            varfnd = ' '
            value = chkdata(7,callr,varfnd,0,stk)
            ! FATAL error if the climate data is missing
            if(index(varfnd, 'precip') .eq. 0) call abortrun("missing site file climate data")
            value = chkdata(7,callr,'reread',0,stk)
            ii = 0
          endif
        endif
      endif

      if(ii .le. 0) then ! read from the site file if no weather file statistics
        do ii = 1, MONTHS
          precip(ii) = chkdata(7,callr,'precip',0,stk)
        end do

        do ii = 1, MONTHS
          prcstd(ii) = chkdata(7,callr,'prcstd',0,stk)
        end do

        do ii = 1, MONTHS
          prcskw(ii) = chkdata(7,callr,'prcskw',0,stk)
        end do

        do ii = 1, MONTHS
          tmn2m(ii) = chkdata(7,callr,'tmn2m',0,stk)
        end do

        do ii = 1, MONTHS
          tmx2m(ii) = chkdata(7,callr,'tmx2m',0,stk)
        end do
      end if

!      varfnd = ' '
!      value = chkdata(7,callr,varfnd,0,stk)
!      if(index(varfnd,'cldcov').gt.0) then
!        cldcov(1) = value
!        do ii = 2, MONTHS
!          cldcov(ii) = chkdata(7,callr,'cldcov',0,stk)
!        end do
!        sitpar = .TRUE.
!      else
!        value = chkdata(7,callr,'reread',0,stk)
!      endif

      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(index(varfnd,'sradj').gt.0) then
        sradadj(1) = value
        do ii = 2, MONTHS
          sradadj(ii) =  chkdata(7,callr,'sradj',0,stk)
        end do
      else
        sitpar = .FALSE.
        value = chkdata(7,callr,'reread',0,stk)
      endif

      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(varfnd .eq. 'rainhr') then
        hours_rain = value
        read(7,*) ! read the header line if we found the target

      else
       ! don't recycle the failed line since it should be a header
       ! value = chkdata(7,callr,'reread',0,stk)
        sitpar = .FALSE.
      endif

! ... sitepar.in variables added to this file KLK 4Oct11

      !***Site and control parameter block
      ! This block has been completely restructured.
      ! There is an attempt to read legacy files
      ivauto = int(chkdata(7,callr,'ivauto',0,stk))
      nelem = int(chkdata(7,callr,'nelem',0,stk))

      sitlat = chkdata(7,callr,'sitlat',0,stk)
      sitlng = chkdata(7,callr,'sitlng',0,stk)

      ! *** NEW variables for site specific daylength
      ! elevation is the critical variable
      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(varfnd(:4) .eq. 'elev') then
        elevation = value
      else
        sitpar = .FALSE.
        value = chkdata(7,callr,'reread',0,stk)
      endif

      ! site slope, aspect, east and west horizons are probably unknown
      ! allow them to be omitted and default to zero
      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(varfnd .eq. 'sitslp' .or. varfnd .eq. 'sslope') then
        sitslp = value
        aspect = chkdata(7,callr,'aspect',0,stk)
        ehoriz = chkdata(7,callr,'ehoriz',0,stk)
        whoriz = chkdata(7,callr,'whoriz',0,stk)
      else
        sitpar = .FALSE.
        value = chkdata(7,callr,'reread',0,stk)
      endif

      ! LEGACY SAND, SILT, CLAY, BULKDEN, PH, AWILT and AFIEL inputs calculated
      ! from soil layer information. Read the group if they are there; skip them
      ! if not but don't specify some and not others. NO ERROR CHECKS
      call skipnames('sand silt clay rock bulkd', '',&
                           ' replaced by layer data ')

      ! soil water block, layers, drain, deep flow, precip controls have been conserved
      nlayer = 0
      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(varfnd(:6) .eq. 'nlayer') then
        nlayer = min(int(value), 9)
        if (value .gt. 9) call message('Warning: resetting nlayer input > 9, reset to 9')
      else
        value = chkdata(7,callr,'reread',0,stk)
        nlayer = 0
      endif

      nlaypgst = int(chkdata(7,callr,'nlaypg',0,stk))
      nlaypg   = nlaypgst !set dynamic layers to site limit
      drain    = chkdata(7,callr,'drain',0,stk)
      basef    = chkdata(7,callr,'basef',0,stk)
      stormf   = chkdata(7,callr,'stormf',0,stk)
! ... Add precipitation amount required for runoff and precipitation fraction
! ... above that level which is lost via runoff, cak - 05/27/03
      precro   = chkdata(7,callr,'precro',0,stk)
      fracro   = chkdata(7,callr,'fracro',0,stk)
      swflag   = int(chkdata(7,callr,'swflag',0,stk))
      ! not sure what to do with swflag

      ! again more soil data
      call skipnames('awilt afiel ph', '',' replaced by layer data ')

      pslsrb = chkdata(7,callr,'pslsrb',0,stk)
      sorpmx = chkdata(7,callr,'sorpmx',0,stk)

! ignore usexdrvrs; We are going to set it from the weather file
      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(varfnd(:6) .eq. 'sublim') then
        sublimscale = value
        reflec      = chkdata(7,callr,'reflec',0,stk)
        albedo      = chkdata(7,callr,'albedo',0,stk)
        dmpflux     = chkdata(7,callr,'dmpflx',0,stk)
      else
        sitpar      = .FALSE.
        value = chkdata(7,callr,'reread',0,stk)
      endif

      ! check for the drain lag
      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(varfnd(:6) .eq. 'drnlag') then
        drainlag    = value
      else
        drainlag    = 0.
        value = chkdata(7,callr,'reread',0,stk)
      endif

      ! ignore usexdrvrs; We are going to set it from the weather file
      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(varfnd(:5) .eq. 'hpotd') then
        hpotdeep    = value
        ksatdeep    = chkdata(7,callr,'ksatd',0,stk)
        call skipnames('tbmin tbmax', 'unused soil boundary temperature value ','')
        dmp         = chkdata(7,callr,'dmpfct',0,stk)
        call skipnames('timlag', 'unused soil boundary temperature value ','')

        ! changes for WFPS effect on N2O flux Cindy K 7 Aug 2013
        Ncoeff        = chkdata(7,callr,'ncoeff',0,stk)
        Ncoeff        = sign(min(abs(Ncoeff),1.),Ncoeff) ! restrict Ncoeff to 1
        jdayStart     = chkdata(7,callr,'dnstrt',0,stk)
        jdayEnd       = chkdata(7,callr,'dnend',0,stk)
        N2Oadjust_fc  = chkdata(7,callr,'nadjfc',0,stk) ! 'nadjfc|nadj'
        N2Oadjust_wp  = chkdata(7,callr,'nadjwp',0,stk)
        MaxNitAmt     = chkdata(7,callr,'maxnit',0,stk)
        netmn_to_no3  = chkdata(7,callr,'mino3',0,stk)
        wfpsdnitadj   = chkdata(7,callr,'wfpsnip',0,stk)
        N2N2Oadj      = chkdata(7,callr,'n2n2oa',0,stk)
      else
        value = chkdata(7,callr,'reread',0,stk)
        sitpar = .FALSE.
        varfnd = ' '
      endif

      ! skip variables for methanogenesis model
      call skipnames('fn2n2o co2ch4 acferm mxch4f frexud aeh deh behfl behdr methzr', &
                           'deprecated Methane variable ','')

      ! skip the old water table
      call skipnames('wtrtbl','water table replaced by event ','')

      epnfa(INTCPT) = chkdata(7,callr,'epnfa',0,stk)
      epnfa(SLOPE) = chkdata(7,callr,'epnfa',0,stk)
      epnfs(INTCPT) = chkdata(7,callr,'epnfs',0,stk)
      epnfs(SLOPE) = chkdata(7,callr,'epnfs',0,stk)
      satmos(INTCPT) = chkdata(7,callr,'satmos',0,stk)
      satmos(SLOPE) = chkdata(7,callr,'satmos',0,stk)
      sirri = chkdata(7,callr,'sirri',0,stk)

      varfnd = ''
      afue = min(1.0,chkdata(7,callr,varfnd,0,stk))
      if(index(varfnd,'afue') .eq. 0) then
        ! default: reading file and variable not found.
        ! Don't reuse this value! This is the last variable
        ! tafue = chkdata(unit,callr,'reread',0,stk) ! reuse value
        afue = 0.5
        ! stack(m-1, is) = tafue
      else
        read(7,*) ! gat the target now skip the header line
      endif

      ! On a binary extend, do not read the initial conditions
      ! Extends used to return from here
      ! if (extndopt .ne. 0) return
      !  The extended format needs to read soil data from later in the file
      !  This was a mistake but it is to hard to correct now    KLK 20 Sep 2012
      ! This commented code will skip these values but it doesn't appear necessary
      !  since a binary extend overwrites these values
      !if (extndopt .ne. 0) then
      !  ! Check for the soil layer parameters
      !  varfnd = ' '
      !  do while (index(varfnd,'soil param') .eq. 0)
      !    varfnd = ' '
      !    value = chkdata(7,callr,varfnd,0,stk)
      !    if(varfnd .eq. 'EOF') return
      !  end do
      !else

      som1ci(SRFC,UNLABL) = chkdata(7,callr,'som1ci',0,stk)
      som1ci(SRFC,LABELD) = chkdata(7,callr,'som1ci',0,stk)
      som1ci(SOIL,UNLABL) = chkdata(7,callr,'som1ci',0,stk)
      som1ci(SOIL,LABELD) = chkdata(7,callr,'som1ci',0,stk)
      som1c = sum(som1ci, dim=2)

      som2ci(SRFC,UNLABL) = chkdata(7,callr,'som2ci',0,stk)
      som2ci(SRFC,LABELD) = chkdata(7,callr,'som2ci',0,stk)
      som2ci(SOIL,UNLABL) = chkdata(7,callr,'som2ci',0,stk)
      som2ci(SOIL,LABELD) = chkdata(7,callr,'som2ci',0,stk)
      som2c = sum(som2ci, dim=2)

      som3ci(UNLABL) = chkdata(7,callr,'som3ci',0,stk)
      som3ci(LABELD) = chkdata(7,callr,'som3ci',0,stk)
      som3c = sum(som3ci)

      do ii = SRFC, SOIL
        do jj = 1, MAXIEL
          rces1(ii,jj) = chkdata(7,callr,'rces1',0,stk)
          if (rces1(ii,jj).le.0) rces1(ii,jj) = rcesd(jj,ii)
          som1e(ii,jj)=som1c(ii)/rces1(ii,jj)
        end do
      end do

      do ii = SRFC, SOIL
        do jj = 1, MAXIEL
          rces2(ii,jj) = chkdata(7,callr,'rces2',0,stk)
          if (rces2(ii,jj).le.0) rces2(ii,jj) = rcesd(jj,3)
          som2e(ii,jj)=som2c(ii)/rces2(ii,jj)
        end do
      end do

      do ii = 1, MAXIEL
        rces3(ii) = chkdata(7,callr,'rces3',0,stk)
!        if (rces3(ii).le.0) rces3(ii) = rcesd(ii,4)
        if (rces3(ii).gt.0) som3e(ii)=som3c/rces3(ii)
      end do

      clittr(SRFC,UNLABL) = chkdata(7,callr,'clittr',0,stk)
      clittr(SRFC,LABELD) = chkdata(7,callr,'clittr',0,stk)
      clittr(SOIL,UNLABL) = chkdata(7,callr,'clittr',0,stk)
      clittr(SOIL,LABELD) = chkdata(7,callr,'clittr',0,stk)
! ... Add check for initial litter values to prevent divide by zero error
! ... in calciv subroutine
      if (ivauto .eq. 0) then
        if(clittr(SRFC,UNLABL) + clittr(SRFC,LABELD) .le. 0.0) then
          call abortrun('Site initial surface litter, clittr(1,*),'//' <= 0.0')
        endif
        if(clittr(SOIL,UNLABL) + clittr(SOIL,LABELD) .le. 0.0) then
          call abortrun('Site initial soil litter, clittr(2,*),'//' <= 0.0')
        endif
      endif

      do ii = SRFC, SOIL
        do jj = 1, MAXIEL
          rcelit(ii,jj) = chkdata(7,callr,'rcelit',0,stk)
        end do
      end do

      aglcis(UNLABL) = chkdata(7,callr,'aglcis',0,stk)
      aglcis(LABELD) = chkdata(7,callr,'aglcis',0,stk)

      do ii = 1, MAXIEL
        aglive(ii) = chkdata(7,callr,'aglive',0,stk)
      end do

! ... Read the fine root values into the juvenile fine root pool, these will
      ! be split between fine root and mature root based on the fraction
      ! of fine root production that went to mature roots based on the
      ! parameter read for the initial crop. cak - 05/24/2007
      bglcisj(UNLABL) = chkdata(7,callr,'bglcis',0,stk)
      bglcisj(LABELD) = chkdata(7,callr,'bglcis',0,stk)

      do ii = 1, MAXIEL
        bglivej(ii) = chkdata(7,callr,'bglive',0,stk)
      end do

      stdcis(UNLABL) = chkdata(7,callr,'stdcis',0,stk)
      stdcis(LABELD) = chkdata(7,callr,'stdcis',0,stk)

      do ii = 1, MAXIEL
        stdede(ii) = chkdata(7,callr,'stdede',0,stk)
      end do

      read(7,*)
      rlvcis(UNLABL) = chkdata(7,callr,'rlvcis',0,stk)
      rlvcis(LABELD) = chkdata(7,callr,'rlvcis',0,stk)

      do ii = 1, MAXIEL
        rleave(ii) = chkdata(7,callr,'rleave',0,stk)
      end do

      fbrcis(UNLABL) = chkdata(7,callr,'fbrcis',0,stk)
      fbrcis(LABELD) = chkdata(7,callr,'fbrcis',0,stk)

      do ii = 1, MAXIEL
        fbrche(ii) = chkdata(7,callr,'fbrche',0,stk)
      end do

      rlwcis(UNLABL) = chkdata(7,callr,'rlwcis',0,stk)
      rlwcis(LABELD) = chkdata(7,callr,'rlwcis',0,stk)

      do ii = 1, MAXIEL
        rlwode(ii) = chkdata(7,callr,'rlwode',0,stk)
      end do

! ... Read the fine root values into the juvenile fine root pool, these will be
! ... split between the fine root and mature root pools based on the fraction of
! ... fine root production that goes to mature roots parameter for the initial
! ... tree. cak - 05/24/2007
      frtcisj(UNLABL) = chkdata(7,callr,'frtcis',0,stk)
      frtcisj(LABELD) = chkdata(7,callr,'frtcis',0,stk)

      do ii = 1, MAXIEL
        frootej(ii) = chkdata(7,callr,'froote',0,stk)
      end do

      crtcis(UNLABL) = chkdata(7,callr,'crtcis',0,stk)
      crtcis(LABELD) = chkdata(7,callr,'crtcis',0,stk)

      do ii = 1, MAXIEL
        croote(ii) = chkdata(7,callr,'croote',0,stk)
      end do

      wd1cis(UNLABL) = chkdata(7,callr,'wd1cis',0,stk)
      wd1cis(LABELD) = chkdata(7,callr,'wd1cis',0,stk)
      wood1c = sum(wd1cis)
      wd2cis(UNLABL) = chkdata(7,callr,'wd2cis',0,stk)
      wd2cis(LABELD) = chkdata(7,callr,'wd2cis',0,stk)
      wood2c = sum(wd2cis)
      wd3cis(UNLABL) = chkdata(7,callr,'wd3cis',0,stk)
      wd3cis(LABELD) = chkdata(7,callr,'wd3cis',0,stk)
      wood3c = sum(wd3cis)


      ! Check for the forest age parameter
      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(index(varfnd,'frstage').gt.0) then
        frstage = value
         read(7,*) ! skip the header
      else
        frstage = 0.
        if(sum(rlwcis).ne.0) call message(&
           'restarting existing forest with juvenile growth');
      ! elseif(index(varfnd,'mineral').gt.0) then  ! test header
      ! else; value = chkdata(7,callr,'reread',0,stk) ! skip this since the next line is a header
      endif

      do ii = 1, MAXIEL
        do jj = 1, CMXLYR
          minerl(jj,ii) = chkdata(7,callr,'minerl',0,stk)
        end do
      end do

      do ii = 1, MAXIEL
        parent(ii) = chkdata(7,callr,'parent',0,stk)
      end do

! ... The secndy and occlud input values can now be double precision, cak - 03/20/02
      do ii = 1, MAXIEL
        secndy(ii) = chkdata(7,callr,'secndy',0,stk)
      end do

      occlud = chkdata(7,callr,'occlud',0,stk)
      read(7,*)

      do ii = 1, CMXLYR
        rwcf(ii) = chkdata(7,callr,'rwcf',0,stk)
      end do

      snlq = chkdata(7,callr,'snlq',0,stk)
      snow = chkdata(7,callr,'snow',0,stk)

      varfnd = '';
      value = chkdata(7,callr,varfnd,0,stk)
      if(index(varfnd,'snwins').gt.0) then
        SnowFlag = value
        sitpar = .TRUE.
      else
        value = chkdata(7,callr,'reread',0,stk)
      endif

      ! Check for the soil layer parameters
      varfnd = ' '
      value = chkdata(7,callr,varfnd,0,stk)
      if(index(varfnd,'Soil') + index(varfnd,'soil') .NE. 0) then
        soilread = .TRUE.
        do ii = 1, MXSWLYR
           varfnd = ''
           dpthmx(ii) = chkdata(7,callr,varfnd,0,stk)
           if(index(varfnd,'Enhanced') .ne. 0  .or. index(varfnd,'enhanced') .ne. 0 &
              .or. varfnd .eq. ' ') then
             dpthmx(ii:MXSWLYR)   = 0.
             lyblkd(ii:MXSWLYR)   = 0.
             fieldc(ii:MXSWLYR)   = 0.
             wiltpt(ii:MXSWLYR)   = 0.
             ecoeff(ii:MXSWLYR)   = 0.
             tcoeff(ii:MXSWLYR)   = 0.
             sandfrac(ii:MXSWLYR) = 0.
             clayfrac(ii:MXSWLYR) = 0.
             orgfrac(ii:MXSWLYR)  = 0.
             swclimit(ii:MXSWLYR) = 0.
             satcond(ii:MXSWLYR)  = 0.
             lyrpH(ii:MXSWLYR)    = 0.
             exit;
           else if(index(varfnd, 'sldpmx') .ne. 0) then

!             dpthmx(ii) = chkdata(7,callr,'sldpmx',0,stk)
             if(dpthmx(II).gt.0) numlyrs = ii ! record the deepest layer
             lyblkd(ii)   = chkdata(7,callr,'slblkd',0,stk)
             fieldc(ii)   = chkdata(7,callr,'slfldc',0,stk)
             wiltpt(ii)   = chkdata(7,callr,'slwltp',0,stk)
             ecoeff(ii)   = chkdata(7,callr,'slecof',0,stk)
             tcoeff(ii)   = chkdata(7,callr,'sltcof',0,stk)
             sandfrac(ii) = chkdata(7,callr,'slsand',0,stk)
             clayfrac(ii) = chkdata(7,callr,'slclay',0,stk)
             orgfrac(ii)  = chkdata(7,callr,'slorgf',0,stk)
             swclimit(ii) = chkdata(7,callr,'slclim',0,stk)
             satcond(ii)  = chkdata(7,callr,'slsatc',0,stk)
             lyrpH(ii)    = chkdata(7,callr,'slph'  ,0,stk)

           else
             call abortrun('site file soil input error: found '//trim(varfnd)//' not sldpmx')
           endif
        end do
        dpthmn(1) = 0 ! Set minimum depth
        dpthmn(2:MXSWLYR) = dpthmx(1:MXSWLYR-1)
      endif

      ! check the slclim input.
      ! A nieve copy of soils.in data will get this value wrong
      ! deltamin model gives 0 at depth but not likely that the soil will dry to 0
      !  thus assume a zero minimum is an input error
      if(minval(swclimit(1:numlyrs)) .eq. 0) then
        swclimit(1:numlyrs) = wiltpt(1:numlyrs) - swclimit(1:numlyrs)
        call message('WARNING: soil layer input error, swclimit = 0, converting deltamin to swclimit')
      endif
      ! redo the layer input check of the limit after a possible subtraction
      if(minval(swclimit(1:numlyrs)) .lt. 0) call abortrun('soil layer input error, layer swclimit < 0')

      if (extndopt .ne. 0) return ! return on a real extend

      if(index(varfnd,'Enhanced ext') .ne. 0) then
      !  We found the enhancement block; read the variables
      !      write(*,*) "varfnd=",trim(varfnd)

         varfnd = ' '
         metcis(SRFC,UNLABL) = chkdata(7,callr,varfnd,0,stk)
         metcis(SRFC,LABELD) = chkdata(7,callr,varfnd,0,stk)
         metcis(SOIL,UNLABL) = chkdata(7,callr,varfnd,0,stk)
         metcis(SOIL,LABELD) = chkdata(7,callr,varfnd,0,stk)

         if(varfnd(:6) .eq. 'rmetcs') then
           metcis = metcis * clittr
         else if(varfnd(:6) .ne. 'metcis') then
           call abortrun('Found'//"'"//trim(varfnd)//"' instead of RMETCS/METCIS in <site>.100")
         endif

         metabc = sum(metcis, dim=2) ! do this here so we can calculate metabe from rcemet
         strcis = clittr -metcis     !Set structural C as clittr - metcis
         strucc = sum(strcis, dim=2)

         varfnd = ''
         do jj = SRFC, SOIL
           value = sum(clittr(jj,:))
           do ii = 1, MAXIEL
             metabe(jj,ii) = chkdata(7,callr,varfnd,0,stk)
             if(metabe(jj,ii) .gt.0  .and. varfnd(:6) .eq. 'rcemet') metabe(jj,ii) = &
                 metabc(jj) / metabe(jj,ii)
!            Calculate the strctural C as clittr - metcis
             if(ii .le. NELEM  .and.  rcelit(jj,ii) >0) then
                struce(jj,ii) = value/rcelit(jj,ii) -metabe(jj,ii)
             else
                struce(jj,ii) = 0
             endif
           end do
         end do

         varfnd = ''
         strlig(SRFC) = chkdata(7,callr,varfnd,0,stk)
         strlig(SOIL) = chkdata(7,callr,varfnd,0,stk)
         if(varfnd(:6) .eq. 'rstrlg') then
           strlig = strlig * sum(strcis, dim=2)
         else if(varfnd(:6) .ne. 'strlig') then
           call abortrun('Found'//"'"//trim(varfnd)//"' instead of RSTRLG/STRLIG in <site>.100")
         endif

         do ii = 1, MAXIEL
           wood1e(ii) = chkdata(7,callr,'rcewd1',0,stk)
           if(wood1e(ii) .gt. 0.) wood1e(ii) = wood1c/wood1e(ii)
         end do

         do ii = 1, MAXIEL
           wood2e(ii) = chkdata(7,callr,'rcewd2',0,stk)
           if(wood2e(ii) .gt. 0.) wood2e(ii) = wood2c/wood2e(ii)
         end do

         do ii = 1, MAXIEL
           wood3e(ii) = chkdata(7,callr,'rcewd3',0,stk)
           if(wood3e(ii) .gt. 0.) wood3e(ii) = wood3c/wood3e(ii)
         end do

         varfnd = ' '
         value = chkdata(7,callr,varfnd,0,stk)
         if(index(varfnd,'cstg').gt.0) then
           !indices were reversed; Match the input order; 2018/11/20
           carbostg(CRPSYS,UNLABL) = value
           carbostg(CRPSYS,LABELD) = chkdata(7,callr,'cstg',0,stk)
           carbostg(FORSYS,UNLABL) = chkdata(7,callr,'cstg',0,stk)
           carbostg(FORSYS,LABELD) = chkdata(7,callr,'cstg',0,stk)
         else
           value = chkdata(7,callr,'reread',0,stk)
         endif
         do ii = 1, MAXIEL
           crpstg(ii) = chkdata(7,callr,'crpstg',0,stk)
         end do
         do ii = 1, MAXIEL
           forstg(ii) = chkdata(7,callr,'forstg',0,stk)
         end do

         varfnd = ' '
         value = chkdata(7,callr,varfnd,0,stk)
         if(index(varfnd,'mrtfr').gt.0) then
           ! save the mature root fraction in the C arrays for DETIV
           bglcism = [-1.0,  value]
           frtcism = [-1.0, chkdata(7,callr,'wmrtfr',0,stk)]
         else
           bglcism = 0
           frtcism = 0
           value = chkdata(7,callr,'reread',0,stk)
         endif

         ammonium = chkdata(7,callr,'ammoni',0,stk)
         do ii = 1, MXSWLYR
           nitrate(ii) = chkdata(7,callr,'nitrat',0,stk)
         end do

         do ii = 1, MXSWLYR
           swc(ii) = dble(chkdata(7,callr,'swcini',0,stk))
         end do

!** Signal an extend if some of the critical enhanced values are non zero.
!   prevents extend errors when a normal site is written with extended values = zero
!   8/3/2001 KLK
       if(extndopt.le.1  .and. &
!          abs(som1c(2)).ne.0 .or. abs(som2c(2)).ne.0 .or. abs(som3c).ne.0) extndopt = 1
          som1c(2).gt.0 .and. som2c(2).gt.0 .and. som3c.gt.0) extndopt = 1

!      else
!
!         !... Zero the site extend variables
!     no need to re zero strcis, metcis, struce, metabe, somXe, woodXe, crpstg, or strlig
!     They have been cleared and on any extend they will be properly set.
      end if


!     Check initialization of litter pools since some may not be read directly
!     paranoia says we SHOULD we do this with ALL the C and mineral pools.
      if(min( &
        strcis(SRFC,UNLABL),clittr(SRFC,UNLABL),metcis(SRFC,UNLABL),    &
        strcis(SOIL,UNLABL),clittr(SOIL,UNLABL),metcis(SOIL,UNLABL),    &
        strcis(SRFC,LABELD),clittr(SRFC,LABELD),metcis(SRFC,LABELD),    &
        strcis(SOIL,LABELD),clittr(SOIL,LABELD),metcis(SOIL,LABELD))    &
        .lt. 0.) call abortrun("site initial strucc or metabc, < 0")

      close(unit=7)

      sitparread = sitparread .or. sitpar

      return

      contains
        subroutine fndfilnam(buffr,ii)
          character (len=*) :: buffr
          integer :: ii, jj

          jj = index(buffr,'#',back=.true.)         ! search for comment character
          if(jj .gt. 0) buffr = buffr(:jj-1)        ! skip any comments
          jj = scan(buffr, '"'//"'", back = .true.) ! do we have a quote
          if(jj .gt. 0) then                        ! found a close quote
            ii = index(buffr(:jj-1), buffr(jj:jj), back =.true.)  ! find the open quote
            if(ii .gt. 0) buffr = buffr(ii+1:jj-1)  ! take just the quoted string
          endif
          ii = len_trim(buffr)                      ! length of file name
          return
        end subroutine fndfilnam

        subroutine skipnames(list, mssg1, mssg2)
          character (len=*) :: list, mssg1, mssg2
          integer           :: lvf
          character (len=8) :: wrnd

          varfnd = ' '
          do while (varfnd .eq.' ')
            value = chkdata(7,callr,varfnd,0,stk)
            ! get out of here if there is an EOF or a section break
            if(varfnd.eq.'EOF')    return; ! call abortrun('EOF in site file')
            if(varfnd(:2).eq.'**') return

            ! find the end of the variable name
            lvf = index(varfnd,'(')
            if(lvf .gt. 0) then
              lvf = lvf-1
            else
              lvf = len_trim(varfnd);
            endif

            if(index(list,varfnd(:lvf)) .gt. 0) then
              if(wrnd .ne. varfnd(:lvf)) call message('SITE NOTE: '//mssg1// &
                                         '"'//varfnd(:lvf)//'"'//mssg2//' value ignored')
              wrnd = varfnd(:lvf)
              varfnd = ' '
            else
             value = chkdata(7,callr,'reread',0,stk)
             exit
            endif
          end do
        end subroutine skipnames


!       subroutine rdsoilin(soilname)
!
!         USE ISO_C_Binding
!
!         implicit none;
!
!         ! routine to read the soil.in format soil file
!         !
!         ! This read was originally part of initsw -> initlyrs. However, if sitein
!         ! can read these fortran structure, common, it is easier to force sitein
!         ! to read either format. If the soil layers are omitted from the site
!         ! file then immediately trigger a read of soil.in.
!         ! Avoid the cross language issues by reading with fortran.
!         !
!
!         ! formal parameters
!         character (len=*) :: soilname;
!
!         include "const.inc"
!         include "sitsoil.inc"
!         ! real :: numlyrs;
!
!         ! local variables
!         integer ierr, lun
!         integer :: izero = ichar("0")
!         real :: deltamin
!         character (len=20) :: line
!
!         numlyrs = 0;        !initialize layer counter
!
!         open(newunit=lun, file=soilname, status='old', iostat=ierr)
!         if (ierr.ne.0) call abortrun("opening file: "//trim(soilname))
!
!         ierr = 0
!         do while (ierr.eq.0)
!           numlyrs = numlyrs+1
!           if (numlyrs > MXSWLYR) then
!             write(line,'(i4)') MXSWLYR
!             call abortrun("Number of soil layers exceeds maximum MXSWLYR, "//line(:4));
!           endif
!           read(lun,*, IOSTAT = ierr) dpthmn(numlyrs), dpthmx(numlyrs), &
!                 bulkd(numlyrs), fieldc(numlyrs), wiltpt(numlyrs), ecoeff(numlyrs), &
!                 tcoeff(numlyrs), sandfrac(numlyrs), clayfrac(numlyrs), &
!                 orgfrac(numlyrs), deltamin, satcond(numlyrs), pH(numlyrs)
!                 swclimit(numlyrs) = wiltpt(numlyrs) - deltamin;
!
!           if(ierr.ne.0) exit
!
!          ! write(*,'("soil.in:",i3,2("	",f5.1),"	",f4.2,2("	",f6.4),'// &
!          !           '6("	",f4.2),"	",f7.5,"	",f4.2)')  numlyrs,dpthmn(numlyrs), &
!          !      dpthmx(numlyrs), bulkd(numlyrs), fieldc(numlyrs), wiltpt(numlyrs), &
!          !      ecoeff(numlyrs), tcoeff(numlyrs), sandfrac(numlyrs), clayfrac(numlyrs), &
!          !      orgfrac(numlyrs), swclimit(numlyrs), satcond(numlyrs), pH(numlyrs);
!
!
!           swclimit(numlyrs) = wiltpt(numlyrs) - deltamin;
!           if (swclimit(numlyrs) < 0.0) then
!             if (deltamin > 0 .and. swclimit(numlyrs) > -0.0001) then
!               swclimit(numlyrs) = 0
!               write(line,'(i3,f6.4," - ",f6.4)') numlyrs,wiltpt(numlyrs), deltamin
!               call message("Warning  soils water layer "//line(:3)//" swclimit = "// &
!                             line(4:)//" rounded up to zero");
!             else
!               write(line,'(i3,f5.3," - ",f5.3)') numlyrs,wiltpt(numlyrs), deltamin
!               call abortrun("site soil water layer "//line(:3)//" swclimit = "// &
!                             line(4:)//" < 0.0.");
!             endif
!           endif
!
!         enddo
!         close(lun);
!         numlyrs = numlyrs-1
!
!         return
!       end subroutine rdsoilin
      end subroutine sitein


      subroutine sitepostp()
        USE ISO_C_Binding
        implicit none

        include 'const.inc'
        include 'doubles.inc'
        include 'param.inc'     ! precip, tmn2m, tmx2m
        include 'plot1.inc'
        include 'sitsoil.inc'
        include 'wth.inc'

        integer :: ii
        character (len=80) :: buffr

        ! Check that tmax is greater than tmin read from <site>.100, cak - 09/17/03
        !  SANITY CHECK data precip > 0 and tmax>tmin
        do ii = 1, MONTHS
          if      (tmn2m(ii) .gt. tmx2m(ii)) then
            write(buffr,'(i3,f7.2," > ",f7.2)') ii, tmn2m(ii), tmx2m(ii)
            buffr = 'site temperature minimum exceeds maximum month'//trim(buffr)
            call abortrun(buffr)
          else if (tmn2m(ii) .gt. tmx2m(ii)) then
            tmn2m(ii)= tmn2m(ii) - 0.05
          else if (precip(ii) .lt. 0.) then
            write(buffr,'(i3," =",f7.2)') ii,precip(ii)
            buffr = 'site precipitation for month'//trim(buffr)//' < 0.0'
            call abortrun(buffr)
          end if
        end do

        if ((nelem .le. 0) .or. (nelem .gt. 3)) then
          write(buffr,'(i5)') nelem
          call abortrun('site NELEM out of bounds = '//buffr(:5))
        endif

        maxt = maxval(tmx2m)   ! maxt  for nitrify added 4/03/98 - mdh
        precipyr = sum(precip) ! initialize precipyr with the average data

        secndy_double = DBLE(secndy)
        occlud_double = DBLE(occlud)

        dpthmn(1) = 0 ! Set minimum depth
        dpthmn(2:MXSWLYR) = dpthmx(1:MXSWLYR-1)

      end subroutine sitepostp


      logical function wthrstat (daywth, prec, pstd, pskw, tmn2m, tmx2m)
        implicit none;
        character*(*) :: daywth; ! weather file name
        real, dimension(12) :: prec, pstd, pskw, tmn2m, tmx2m

        ! local variables
        character*256 :: iline

        integer       :: i, mm, lch;
        integer       :: istat, iost
        integer       :: dom, m, y, doy;
        integer       :: fstyr, lasty;
        real          :: maxt, mint, precip;
        integer, dimension(12)          :: mc, mpc, txc, tnc
        real, dimension(12)             :: mprec
        double precision, dimension(12) :: sump, ssump, csump
        double precision :: var
        ! real       :: cp(12), yk, tk

        character*64  :: mssg

        integer             :: ReadInt
        double precision    :: ReadReal

        ! read the weather file for climate.
        wthrstat = .TRUE.
        prec  = 0.0;
        pstd  = 0.0;
        pskw  = 0.0;
        fstyr = -1048576;

           mprec  = -2**17;
           tmx2m  = -2**17;
           tmn2m  = -2**17;
           mpc    = 0;
           txc    = 0;
           tnc    = 0;
           !cp    = 0;

          prec = 0;
          pstd = 0;
          pskw = 0;

          sump = 0;
          ssump = 0;
          csump = 0;
          mc   = 0;

        iost = 0;
        open(unit=9, file=daywth, status='OLD', IOSTAT=iost)
        if(iost.ne.0) then
          call message('Warning: unable to open "'//trim(daywth)//'"') ! non fatal warning
          wthrstat = .FALSE. ! no weather statistics
          return
        endif
        call message('Calculating mean weather from: "'//trim(daywth)//'"')

        lasty = 0;
        do while (iost .eq. 0)
          read(9, '(a)', iostat=iost, iomsg= mssg) iline
          ! remove everything after a comment character
          i=index(iline, '#');
          if(i .eq. 0) i=index(iline, '!');
          if(i .gt. 0) iline(i:) = ' '
          if (len_trim(iline) .eq. 0) cycle

          lch = 1;
                           dom     = ReadInt(0,iline,lch,istat)
          if(istat .gt. 0) m       = ReadInt(0,iline,lch,istat)
          if(istat .gt. 0) y       = ReadInt(0,iline,lch,istat)
          if(istat .gt. 0) doy     = ReadInt(0,iline,lch,istat)
          if(istat .le. 0) call abortrun('parsing weather date fields:'//trim(iline))
                           maxt   = ReadReal(0,iline,lch,istat)
          if(istat .gt. 0) mint   = ReadReal(0,iline,lch,istat)
          if(istat .gt. 0) precip = ReadReal(0,iline,lch,istat)
          if(istat .le. 0) call abortrun('parsing weather temp maximum, minimum or precip: '//trim(iline))

          if(fstyr .eq. -1048576) then
            fstyr = y
            lasty = y
          endif

          ! process this years precipitation into the the statistics arrays
          if(y .ne. lasty) call monthprec(lasty, mc, mprec,mpc, sump,ssump,csump)
          lasty = y

          if (precip >= -98.0) then
            mpc(m)   = mpc(m)+1;
            if(mprec(m) .eq. -2**17) then
              mprec(m) = precip;
     !         cp(m) = 0
            else
     !         yk     = precip - cp(m)         ! So far, so good: c is zero.
     !         tk     = mprec(m) + yk       ! Alas, sum is big, y small, so low-order digits of y are lost.
     !         cp(m) = (tk - mprec(m)) - yk ! (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
     !         mprec(m) = tk                ! Algebraically, c should always be zero. Beware overly-aggressive optimizing
              mprec(m) = mprec(m) + precip;
            endif
          endif

          if (maxt >= -98.0) then
            txc(m)   = txc(m)+1;
            if(tmx2m(m) .eq. -2**17) then
              tmx2m(m) = maxt;
            else
              tmx2m(m) = tmx2m(m)+maxt;
            endif
          endif

          if (mint >= -98.0) then
            tnc(m)   = tnc(m)+1;
            if(tmn2m(m) .eq. -2**17) then
              tmn2m(m) = mint;
            else
              tmn2m(m) = tmn2m(m)+mint;
            endif
          endif
        end do
        close(9);

        call monthprec(lasty, mc, mprec,mpc, sump,ssump,csump) ! last years precipitation

        ! check for empty temperatures basically this is an empty file check
        if(sum(txc) .eq. 0  .or.  sum(tnc) .eq. 0) then
          call message("Warning: Missing temperature data in: "//trim(daywth))
          wthrstat = .FALSE.
          return
        endif
        where (txc /= 0.)
           tmx2m = tmx2m / txc;
        end where
        where (tnc /= 0.)
           tmn2m = tmn2m / tnc;
        end where

        ! precipitation statistics
        do mm= 1,12
          if(mc(mm) > 0) then
            prec(mm) = 0;
            pstd(mm) = 0;
            pskw(mm) = 0;
            if(sump(mm) .gt. 0) then
              prec(mm) = sump(mm)/mc(mm);
              var = ssump(mm) / mc(mm) - prec(mm)*prec(mm);
              if(var > 0) then
                pstd(mm) = sqrt(var);
                pskw(mm) = (csump(mm)/mc(mm) - 3* prec(mm) * var - prec(mm)**3) / pstd(mm)**3;
              endif
            endif
          else
            prec(mm) = -99.99;
            pstd(mm) = 0;
            pskw(mm) = 0;
          endif
        end do
        return;
      contains

       subroutine monthprec(lasty, mc, mprec,mpc, sump,ssump,csump)
         integer lasty
         real, dimension(12)    :: mprec
         integer, dimension(12) :: mpc, mc(12)
         double precision  :: sump(12), ssump(12), csump(12)

         real    :: dc, pr
         integer :: mm
         real, dimension(12) :: days = [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.];

         do mm = 1,12
           if(mprec(mm) .ge. 0.) then
             mc(mm) = mc(mm) +1;
             dc = days(mm)/mpc(mm)
             if((mm .eq. 2  .and.  mpc(mm) .eq. 29)  .or.  days(mm) .eq. mpc(mm)) then
               pr = mprec(mm)
             else
               pr = mprec(mm) * dc;
             endif
             sump(mm) = sump(mm) + pr;
             ssump(mm) = ssump(mm) + pr*pr;
             csump(mm) = csump(mm) + pr*pr*pr;
           endif
         end do
         mprec = 0
         mpc   = 0
       end subroutine monthprec
      end function wthrstat
