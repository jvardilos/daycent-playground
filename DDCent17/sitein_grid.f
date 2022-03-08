
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... SITEIN_GRID.F

      subroutine sitein_grid(ivopt)

      implicit none
      include 'const.inc'
      include 'doubles.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'wth.inc'

! ... Argument declarations
      integer   ivopt

! NOTE THIS ROUTINE IS SERIOUSLY OUT OF DATE WITH CURRENT SITE FILE REQUIREMENTS!!!
!  included ckdata as a contains routine since it is not used anywhere else anymore
!    KLK MAR 2015
! ... Read the parameter file in the form needed by time0.
! ... This version of sitein reads an extended <site>.100 parameter file
! ... created from a site.nc file from a Gridded Century run. - cak, 10/05/01

! ... Local variables
      integer     ii, jj
      character*6 name
      character*11, parameter ::  callr = 'sitein_grid'
      real        temp

      read(7,*)
      read(7,*)

      do 10 ii = 1, MONTHS
        read(7,*) precip(ii), name
        call ckdata(callr,'precip',name)
10    continue

      do 20 ii = 1, MONTHS
        read(7,*) prcstd(ii), name
        call ckdata(callr,'prcstd',name)
20    continue

      do 30 ii = 1, MONTHS
        read(7,*) prcskw(ii), name
        call ckdata(callr,'prcskw',name)
30    continue

      do 40 ii = 1, MONTHS
        read(7,*) tmn2m(ii), name
        call ckdata(callr,'tmn2m',name)
40    continue

! ... maxt calculation for nitrify added 4/03/98 - mdh
      maxt = -99.0
      do 50 ii = 1, MONTHS
        read(7,*) tmx2m(ii), name
        if (tmx2m(ii) .gt. maxt) maxt = tmx2m(ii)
        call ckdata(callr,'tmx2m',name)
50    continue

! ... Check to be sure that the tmax value read from the <site>.100
! ... file is greater than the tmin value read from the <site>.100
! ... file, cak - 09/17/03
      do 55 ii = 1, MONTHS
        if (tmx2m(ii) .lt. tmn2m(ii)) then
          write(*,*) 'ERROR: Invalid weather data in site file, ',
     &               'tmx2m < tmn2m for month ', ii
          write(*,*) 'tmx2m(', ii, ') = ', tmx2m(ii)
          write(*,*) 'tmn2m(', ii, ') = ', tmn2m(ii)
          STOP
        endif
55    continue

      read(7,*)
      read(7,*) temp, name
      ivauto = int(temp)
      call ckdata(callr,'ivauto',name)
      read(7,*) temp, name
      nelem = int(temp)
      call ckdata(callr,'nelem',name)

      read(7,*) sitlat, name
      call ckdata(callr,'sitlat',name)
      read(7,*) sitlng, name
      call ckdata(callr,'sitlng',name)

      read(7,*) sand, name
      call ckdata(callr,'sand',name)
      read(7,*) silt, name
      call ckdata(callr,'silt',name)
      read(7,*) clay, name
      call ckdata(callr,'clay',name)
! ... Add rock fraction, cak - 05/27/03
      read(7,*) rock, name
      call ckdata(callr,'rock',name)
      if (rock .gt. 0.90) then
        write(*,*) 'Rock fraction too large, rock = ', rock
        STOP
      endif
      read(7,*) bulkd, name
      call ckdata(callr,'bulkd',name)

      read(7,*) temp, name
      nlayer = int(temp)
      if (nlayer .gt. 9) then
        nlayer = 9
        call message('   Warning: nlayer value too large, reset to 9')
      endif
      call ckdata(callr,'nlayer',name)
      read(7,*) temp, name
      nlaypg = int(temp)
      call ckdata(callr,'nlaypg',name)
      read(7,*) drain, name
      call ckdata(callr,'drain',name)
      read(7,*) basef, name
      call ckdata(callr,'basef', name)
      read(7,*) stormf, name
      call ckdata(callr,'stormf', name)
! ... Add precipitation amount required for runoff and fraction of
! ... precipitation above amount required for runoff which is lost
! ... via runoff, cak - 05/27/03
      read(7,*) precro, name
      call ckdata(callr,'precro', name)
      read(7,*) fracro, name
      call ckdata(callr,'fracro', name)
      read(7,*) temp, name
      swflag = int(temp)
      call ckdata(callr,'swflag', name)

      do 60 ii = 1, CMXLYR
        read(7,*) awilt(ii), name
        call ckdata(callr,'awilt', name)
60    continue

      do 70 ii = 1, CMXLYR
        read(7,*) afiel(ii), name
        call ckdata(callr,'afiel', name)
70    continue

      read(7,*) ph, name
      call ckdata(callr,'ph',name)
! ... New phstart variable added for pH shift, cak - 08/02/02
      phstart = ph
      read(7,*) pslsrb, name
      call ckdata(callr,'pslsrb',name)
      read(7,*) sorpmx, name
      call ckdata(callr,'sorpmx',name)
      read(7,*)
      read(7,*) epnfa(INTCPT), name
      call ckdata(callr,'epnfa',name)
      read(7,*) epnfa(SLOPE), name
      call ckdata(callr,'epnfa',name)
      read(7,*) epnfs(INTCPT), name
      call ckdata(callr,'epnfs',name)
      read(7,*) epnfs(SLOPE), name
      call ckdata(callr,'epnfs',name)
      read(7,*) satmos(INTCPT), name
      call ckdata(callr,'satmos',name)
      read(7,*) satmos(SLOPE), name
      call ckdata(callr,'satmos',name)
      read(7,*) sirri, name
      call ckdata(callr,'sirri',name)

! ... If extending, do not read in initial conditions
      if (ivopt.gt.0) goto 999

      read(7,*)
      read(7,*) som1ci(SRFC,UNLABL), name
      call ckdata(callr,'som1ci',name)
      read(7,*) som1ci(SRFC,LABELD), name
      call ckdata(callr,'som1ci',name)
      read(7,*) som1ci(SOIL,UNLABL), name
      call ckdata(callr,'som1ci',name)
      read(7,*) som1ci(SOIL,LABELD), name
      call ckdata(callr,'som1ci',name)

      read(7,*) som2ci(SRFC,UNLABL), name
      call ckdata(callr,'som2ci',name)
      read(7,*) som2ci(SRFC,LABELD), name
      call ckdata(callr,'som2ci',name)
      read(7,*) som2ci(SOIL,UNLABL), name
      call ckdata(callr,'som2ci',name)
      read(7,*) som2ci(SOIL,LABELD), name
      call ckdata(callr,'som2ci',name)

      read(7,*) som3ci(UNLABL), name
      call ckdata(callr,'som3ci',name)
      read(7,*) som3ci(LABELD), name
      call ckdata(callr,'som3ci',name)

      do 90 ii = SRFC, SOIL
        do 80 jj = 1, MAXIEL
          read(7,*) rces1(ii,jj), name
          call ckdata(callr,'rces1',name)
80      continue
90    continue

      do 105 ii = SRFC, SOIL
        do 100 jj = 1, MAXIEL
          read(7,*) rces2(ii,jj), name
          call ckdata(callr,'rces2',name)
100     continue
105   continue

      do 110 ii = 1, MAXIEL
        read(7,*) rces3(ii), name
        call ckdata(callr,'rces3',name)
110   continue

! ... Set of extended data values from Gridded Century run, cak - 10/05/01
      do 112 ii = SRFC, SOIL
        do 114 jj = 1, MAXIEL
          read(7,*) som1e(ii,jj), name
          call ckdata(callr,'som1e',name)
114     continue
112   continue
      do 116 ii = 1, MAXIEL
        read(7,*) som2e(2,ii), name
        call ckdata(callr,'som2e',name)
116   continue
      do 118 ii = 1, MAXIEL
        read(7,*) som3e(ii), name
        call ckdata(callr,'som3e',name)
118   continue

      read(7,*) clittr(SRFC,UNLABL), name
      call ckdata(callr,'clittr',name)
      read(7,*) clittr(SRFC,LABELD), name
      call ckdata(callr,'clittr',name)
      read(7,*) clittr(SOIL,UNLABL), name
      call ckdata(callr,'clittr',name)
      read(7,*) clittr(SOIL,LABELD), name
      call ckdata(callr,'clittr',name)
! ... Add check for initial litter values to prevent divide by zero error
! ... in calciv subroutine
      if ((ivauto .eq. 0) .and.
     &    (clittr(SRFC,UNLABL) + clittr(SRFC,LABELD) .eq. 0.0)) then
        call message('   Initial surface litter values, clittr(1,*),')
        call message('   must be greater than zero.')
        STOP
      endif
      if ((ivauto .eq. 0) .and.
     &    (clittr(SOIL,UNLABL) + clittr(SOIL,LABELD) .eq. 0.0)) then
        call message('   Initial soil litter values, clittr(2,*),')
        call message('   must be greater than zero.')
        STOP
      endif

      do 130 ii = SRFC, SOIL
        do 120 jj = 1, MAXIEL
          read(7,*) rcelit(ii,jj), name
          call ckdata(callr,'rcelit',name)
120     continue
130   continue

! ... Set of extended data values from Gridded Century run, cak - 10/05/01
      read(7,*) strcis(SRFC,UNLABL), name
      call ckdata(callr,'strcis',name)
      read(7,*) strcis(SRFC,LABELD), name
      call ckdata(callr,'strcis',name)
      read(7,*) strcis(SOIL,UNLABL), name
      call ckdata(callr,'strcis',name)
      read(7,*) strcis(SOIL,LABELD), name
      call ckdata(callr,'strcis',name)
      read(7,*) metcis(SRFC,UNLABL), name
      call ckdata(callr,'metcis',name)
      read(7,*) metcis(SRFC,LABELD), name
      call ckdata(callr,'metcis',name)
      read(7,*) metcis(SOIL,UNLABL), name
      call ckdata(callr,'metcis',name)
      read(7,*) metcis(SOIL,LABELD), name
      call ckdata(callr,'metcis',name)

      read(7,*) aglcis(UNLABL), name
      call ckdata(callr,'aglcis',name)
      read(7,*) aglcis(LABELD), name
      call ckdata(callr,'aglcis',name)

      do 140 ii = 1, MAXIEL
        read(7,*) aglive(ii), name
        call ckdata(callr,'aglive',name)
140   continue

! ... Read the fine root values into the juvenile fine root pool, these
! ... values will be split between the fine root and mature root pool
! ... based on the fraction of fine root production that goes to mature
! ... roots parameter value as read from the crop.100 file in the detiv
! ... subroutine, cak - 05/24/2007
      read(7,*) bglcisj(UNLABL), name
      call ckdata(callr,'bglcis',name)
      read(7,*) bglcisj(LABELD), name
      call ckdata(callr,'bglcis',name)

      do 150 ii = 1, MAXIEL
        read(7,*) bglivej(ii), name
        call ckdata(callr,'bglive',name)
150   continue

      read(7,*) stdcis(UNLABL), name
      call ckdata(callr,'stdcis',name)
      read(7,*) stdcis(LABELD), name
      call ckdata(callr,'stdcis',name)

      do 160 ii = 1, MAXIEL
        read(7,*) stdede(ii), name
        call ckdata(callr,'stdede',name)
160   continue

      read(7,*)
      read(7,*) rlvcis(UNLABL), name
      call ckdata(callr,'rlvcis',name)
      read(7,*) rlvcis(LABELD), name
      call ckdata(callr,'rlvcis',name)

      do 170 ii = 1, MAXIEL
        read(7,*) rleave(ii), name
        call ckdata(callr,'rleave',name)
170   continue

      read(7,*) fbrcis(UNLABL), name
      call ckdata(callr,'fbrcis',name)
      read(7,*) fbrcis(LABELD), name
      call ckdata(callr,'fbrcis',name)

      do 180 ii = 1, MAXIEL
        read(7,*) fbrche(ii), name
        call ckdata(callr,'fbrche',name)
180   continue

      read(7,*) rlwcis(UNLABL), name
      call ckdata(callr,'rlwcis',name)
      read(7,*) rlwcis(LABELD), name
      call ckdata(callr,'rlwcis',name)

      do 190 ii = 1, MAXIEL
        read(7,*) rlwode(ii), name
        call ckdata(callr,'rlwode',name)
190   continue

! ... Read the fine root values into the juvenile fine root pool, these
! ... values will be split between the fine root and mature root pool
! ... based on the fraction of fine root production that goes to mature
! ... roots parameter value as read from the tree.100 file in the detiv
! ... subroutine, cak - 05/24/2007
      read(7,*) frtcisj(UNLABL), name
      call ckdata(callr,'frtcis',name)
      read(7,*) frtcisj(LABELD), name
      call ckdata(callr,'frtcis',name)

      do 200 ii = 1, MAXIEL
        read(7,*) frootej(ii), name
        call ckdata(callr,'froote',name)
200   continue

      read(7,*) crtcis(UNLABL), name
      call ckdata(callr,'crtcis',name)
      read(7,*) crtcis(LABELD), name
      call ckdata(callr,'crtcis',name)

      do 210 ii = 1, MAXIEL
        read(7,*) croote(ii), name
        call ckdata(callr,'croote',name)
210   continue

      read(7,*) wd1cis(UNLABL), name
      call ckdata(callr,'wd1cis',name)
      read(7,*) wd1cis(LABELD), name
      call ckdata(callr,'wd1cis',name)
      read(7,*) wd2cis(UNLABL), name
      call ckdata(callr,'wd2cis',name)
      read(7,*) wd2cis(LABELD), name
      call ckdata(callr,'wd2cis',name)
      read(7,*) wd3cis(UNLABL), name
      call ckdata(callr,'wd3cis',name)
      read(7,*) wd3cis(LABELD), name
      call ckdata(callr,'wd3cis',name)

!      read(7,*) w1lig, name
!      call ckdata(callr,'w1lig',name)
!      read(7,*) w2lig, name
!      call ckdata(callr,'w2lig',name)
!      read(7,*) w3lig, name
!      call ckdata(callr,'w3lig',name)

! ... Set of extended data values from Gridded Century run, cak - 10/05/01
      do 212 ii = 1, MAXIEL
        read(7,*) wood1e(ii), name
        call ckdata(callr,'wood1e',name)
212   continue
      do 214 ii = 1, MAXIEL
        read(7,*) wood2e(ii), name
        call ckdata(callr,'wood2e',name)
214   continue
      do 216 ii = 1, MAXIEL
        read(7,*) wood3e(ii), name
        call ckdata(callr,'wood3e',name)
216   continue

      read(7,*)

      do 230 ii = 1, MAXIEL
        do 220 jj = 1, CMXLYR
          read(7,*) minerl(jj,ii), name
          call ckdata(callr,'minerl',name)
220     continue
230   continue

! ... Set of extended data values from Gridded Century run, cak - 10/05/01
      do 232 ii = 1, MAXIEL
        read(7,*) aminrl(ii), name
        call ckdata(callr,'aminrl',name)
232   continue
      do 234 ii = 1, MAXIEL
        read(7,*) crpstg(ii), name
        call ckdata(callr,'crpstg',name)
234   continue

      do 240 ii = 1, MAXIEL
        read(7,*) parent(ii), name
        call ckdata(callr,'parent',name)
240   continue

! ... The secndy and occlud input values can now be double precision,
! ... cak - 03/20/02
      do 250 ii = 1, MAXIEL
        read(7,*) secndy_double(ii), name
        call ckdata(callr,'secndy',name)
250   continue

      read(7,*) occlud_double, name
      call ckdata(callr,'occlud',name)
      read(7,*)

! ... Save the double precision secndy and occlud variables read into their
! ... single precision counterparts, cak - 03/20/02
      secndy(1) = real(secndy_double(1))
      secndy(2) = real(secndy_double(2))
      secndy(3) = real(secndy_double(3))
      occlud = real(occlud_double)

      do 260 ii = 1, CMXLYR
        read(7,*) rwcf(ii), name
        call ckdata(callr,'rwcf',name)
260   continue

! ... Set of extended data values from Gridded Century run, cak - 10/05/01
      do 262 ii = 1, CMXLYR
        read(7,*) asmos(ii), name
        call ckdata(callr,'asmos',name)
262   continue
      do 264 ii = 1, MAXIEL
        read(7,*) avh2o(ii), name
        call ckdata(callr,'avh2o',name)
264   continue

      read(7,*) snlq, name
      call ckdata(callr,'snlq',name)
      read(7,*) snow, name
      call ckdata(callr,'snow', name)

! ... Set of extended data values from Gridded Century run, cak - 10/05/01
      read(7,*) prcann, name
      call ckdata(callr,'prcann',name)
      read(7,*) petann, name
      call ckdata(callr,'petann',name)

      close(unit=7)

999   continue

      return
      contains


         subroutine ckdata(callr,expect,found)

         implicit none

         ! Argument declarations
         character*(*) callr, expect, found

         ! Check a parameter name read by fixin or sitein.
         !   callr is the name of the calling routine.
         !   expect is the name that should have been read.
         !   found  is the name that was read from the data file.

         ! Local variables
         integer ii, loc
         character par*6

         ! Extract junk from the parameter name

         loc = index(found,'(')
         if (loc .eq. 0) then
           loc = index(found,'|')
           if (loc .eq. 0) loc = index(found,',')
         endif
         if (loc .eq. 0) then
           par = found( :6)
         else
           par = found( :loc-1)
         endif

         ! Convert to lower case if needed (MACHINE DEPENDANT)
         do ii = 1, 6
           if (par(ii:ii) .ge. 'A' .and. par(ii:ii) .le. 'Z') then
             par(ii:ii) = char( ichar(par(ii:ii)) + 32 )
           endif
         end do

         ! Test for expected name
         if (expect .ne. par) then
           !call message('   There is an error in your <grid site>.100 file.')
           call abortrun('Data for '// par //' was read when the '//
     &     'data for '// expect //' was expected in <grid site>.100.')
         endif

         return
         end subroutine ckdata

      end subroutine sitein_grid
