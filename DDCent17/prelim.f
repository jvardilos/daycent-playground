
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine prelim

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Initialize variables and parameters

c ... Local variables
      integer   lyr, ostat
      real      dely, delx, xslope, textur, yint
      real   :: adepsum
      ! real   :: soildepth, rockf, ompc ! field capacity and wilting point variables

c ... Time initializations -  time step is one month
      dt = 1.0/12.0

c     Location initialization
      if (sitlat .ge. 0.0) then
        dayhrs = 0.0
        hrsinc = .TRUE.
        decidgrow = .FALSE.
      else
        dayhrs = 24.0
        hrsinc = .FALSE.
        decidgrow = .TRUE.
      endif

c       Save the current DRAIN value in the undo array
        drain = max(0.0,min(1.0,drain))
        undodran = drain   ! array

c ... Initialization for growing degree day implementation,
c ... cak - 04/17/03
      accumdd = .false.
      thermunits = 0.0
      cultcnt = 0
      fertcnt = 0
      erodcnt = 0
      grazcnt = 0
      plntcnt = 0
      senecnt = 0
      frstschd = .false.
      harvschd = .false.
      plntschd = .false.
      senmschd = .false.

c ... Allow for time step < 1 month for running decomp
c ... ntspm is the number of time steps per month for decomp
c ... (read from the fix.100 file)

c ... decodt is the time step used in subroutine decomp
c ... decodt is assigned in subroutine dailymoist in the daily
c ... water budget version of Century
c      decodt = dt/real(ntspm)

c ... Initializations
      crpgrw = 0
      seedl = 0
      forgrw = 0
      falprc = 0

c ... Initialize volitalization accumulators
      volgma = 0.0
      volexa = 0.0
      volpla = 0.0

c ... Initialize erosion variables
      scloss = 0.0
      sclosa = 0.0

c ... Initialize accumulators
      call annacc

c ... Open the c14 data file
      if (labtyp .eq. 1) then
        open(unit=10,file='c14data',status='OLD', iostat=ostat)
        if(ostat .gt. 0) call abortrun('opening file "c14data"')
      endif

c ... Open the N input scalar data file
      if (Ninput .gt. 0) then
        open(unit=20,file='nscale.dat',status='OLD', iostat=ostat)
        if(ostat .gt. 0) call abortrun('opening file "nscale.dat"')
      endif

c ... Open the OMAD input scalar data file
      if (OMADinput .gt. 0) then
        open(unit=30,file='omadscale.dat',status='OLD', iostat=ostat)
        if(ostat .gt. 0) call abortrun('opening file "omadscale.dat"')
      endif

c ... Open the pH scalar data file
      if (phsys .gt. 0) then
        open(unit=40,file='phscale.dat',status='OLD', iostat=ostat)
        if(ostat .gt. 0) call abortrun('opening file "phscale.dat"')
      endif

c ... Open the precipitation input scalar data file
      if (wthinput .eq. 4 .or. wthinput .eq. 5) then
        open(unit=50,file='precscale.dat',status='OLD', iostat=ostat)
        if(ostat .gt. 0) call abortrun('opening file "precscale.dat"')
      endif

c ... Open the maximum temperature input scalar data file
      if (wthinput .eq. 2 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        open(unit=55,file='tmaxscale.dat',status='OLD', iostat=ostat)
        if(ostat .gt. 0) call abortrun('opening file "tmaxscale.dat"')
      endif

c ... Open the minimum temperature input scalar data file
      if (wthinput .eq. 1 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        open(unit=60,file='tminscale.dat',status='OLD', iostat=ostat)
        if(ostat .gt. 0) call abortrun('opening file "tminscale.dat"')
      endif

c ... Calculate C,N,P,S in lower horizon soil pools for use as soil
c ... replacement with erosion events
        lhzci(1,:) = som1ci(SOIL,:)*lhzf(1)
        lhzci(2,:) = som2ci(SOIL,:)*lhzf(2)
        lhzci(3,:) = som3ci(:)*lhzf(3)

        lhze(1,1:nelem) = som1e(SOIL,1:nelem)*lhzf(1)
        lhze(2,1:nelem) = som2e(SOIL,1:nelem)*lhzf(2)
        lhze(3,1:nelem) = som3e(1:nelem)*lhzf(3)

c ... Field capacity and wilting point.
      ! Set somsc using initial values read from the <site>.100 file CAK - 11/21/00
      somsc = som1ci(2,1) + som1ci(2,2) + som2ci(2,1) + som2ci(2,2) +
     &        som3ci(1) + som3ci(2)
!      disable the field capacity and wilting point calculations for the old
!      CENTURY soil layers since they are being set by soil water profile
!      THIS PREVENTS ANY CONFUSION AS TO HOW THIS IS SET     KLK 7/2014
!      if (swflag .ne. 0) then
!        soildepth =sum(adep(1:nlayer))
!
!        rockf = rock
!        if (soildepth .gt. 150) rockf = 0
!
!        ompc = somsc*1.724/(10000*bulkd*edepth)
!        do lyr = 1, nlayer
!          call swmodel(afiel(lyr), awilt(lyr), swflag, ompc, sand,       &
!     &                   silt, clay, bulkd, rock)
!          ompc = ompc * 0.85
!        end do
!      endif

c ... Added calculation for water content which will be used to
c ... determine plant production in POTGRS. 10-90 -rm
      wc = afiel(1)-awilt(1)

c ... Calculate available water holding capacity in top 30 cm of soil
c ... as this is the zone where the majority of plant roots occur
      awhc  = 0.0
      adepsum = 0.0
      do lyr = 1, nlayer
        adepsum = adepsum + adep(lyr)
        if (adepsum .le. 30.0) then
          awhc = awhc + (afiel(lyr) - awilt(lyr)) * adep(lyr)
        endif
      enddo

c ... Compute ORGLCH for use in SOMDEC.
      orglch = omlech(1) + omlech(2) * sand

c ... Intercept for the texture equation of secondary P depends upon
c ... pH input.
      if (ph .le. phesp(1)) then
        texesp(2) = phesp(2)
      else if (ph .ge. phesp(3)) then
        texesp(2) = phesp(4)
      else
        dely = phesp(4) - phesp(2)
        delx = phesp(3) - phesp(1)
        xslope = dely / delx
        yint = phesp(2) - (xslope*phesp(1))
        texesp(2) = (xslope*ph) + yint
      endif

      if (micosm .eq. 0) then
c ..... Preset array which will contain monthly values of defac
          agdefacm = -1.      ! array
          bgdefacm = -1.      ! array
        aagdefac = 0.
        abgdefac = 0.
        agdefac = 0.
        bgdefac = 0.
      else
        call message(' ')
        call message('Microcosms are not implemented in this version')
        call message('of Daily Century')
        STOP
      endif

c ... Effect of soil texture on the microbe decomposition rate
      eftext = peftxa+peftxb*sand

c ... Compute parameters which control decomposition of som1
c ... p1co2 must be computed for surface and soil.   vek  08/91
c ... Note that p1co2b(1) must equal 0 because there is no
c ... soil texture effect on the surface.
      p1co2(SRFC) = p1co2a(SRFC)
      p1co2(SOIL) = p1co2a(SOIL)+p1co2b(SOIL)*sand

c ... Decomposition of som1 to som3 is a function of clay content
c ... vek june90
      fps1s3 = ps1s3(1) + ps1s3(2) * clay
      fps2s3 = ps2s3(1) + ps2s3(2) * clay

      if (texepp(1) .eq. 1.0) then
c ..... Calculate pparmn(2)
c ..... Include effect of texture; weathering factor should be per year
c ..... Note that this code changes the value of a 'fixed' parameter
c ..... (pparmn(2))
        textur = clay + silt
        pparmn(2) = 12.0 * carctanf(textur, texepp(2), texepp(3),
     &                            texepp(4), texepp(5))
      endif

      if (texesp(1) .eq. 1.0) then
c ..... Calculate psecmn(2)
c ..... Include effect of texture
c ..... Note that this code changes the value of a 'fixed' parameter
c ..... (psecmn(2))
        psecmn(2) = 12.0 * (texesp(2) + texesp(3) * sand)
      endif

c ... Compute VLOSSG as a function of soil texture based on clay content
c ... vlossg_m is the VLOSSG parameter value as read from the fix.100 file,
c ... cak - 11/21/01
      if (clay .lt. 0.10) then
c        vlossg = 0.015
        vlossg = 0.03
c      else if (clay .gt. 0.40) then
      else if (clay .gt. 0.30) then
c        vlossg = 0.003
        vlossg = 0.01
      else
        vlossg = line(clay, 0.10, 0.03, 0.30, 0.01) ! (clay, 0.10, 0.015, 0.40, 0.003)
      endif
      vlossg = vlossg * vlossg_m

c ... Save initial values for printing or plotting
      call savarp

c ... Clear the flow stack.
      call floclr

      return
      contains
        include 'catanf.f'

        include 'line.f'
      end subroutine prelim


      subroutine swmodel (afiel, awilt, swflag, ompc, sand, silt,        &
     &                      clay, bulkd, rockf)

c ... swflag lets the model user choose between using actual data
c ... for awilt and afiel or equations from Gupta and Larson (1979)
c ... or Rawls et al (1982).
c ...
c ... swflag=0 Use actual data
c ... swflag=1 Use G&L for both awilt (-15 bar) and afiel (-0.33 bar)
c ... swflag=2 Use G&L for both awilt (-15 bar) and afiel (-0.10 bar)
c ... swflag=3 Use Rawls for both awilt (-15 bar) and afiel (-0.33 bar)
c ... swflag=4 Use Rawls for both awilt (-15 bar) and afiel (-0.10 bar)
c ... swflag=5 Use Rawls for afiel (-0.33 bar) and actual data for awilt
c ... swflag=6 Use Rawls for afiel (-0.10 bar) and actual data for awilt
c ...
c ...     swflag   1          2          3        4       5       6
      real, parameter :: fcsa(6) =[0.3075,    0.5018,   -0.20,           &
     &                             -0.30,  -0.19,   0.31]
      real, parameter :: fcsi(6) =[0.5886,    0.8548,    0.0,            &
     &                            0.0,    0.0,    0.0]
      real, parameter :: fccl(6) =[0.8039,    0.8833,    0.36,           &
     &                            0.23,   0.0,    0.0]
      real, parameter :: fcom(6) =[2.208E-03, 4.966E-03, 0.0299,         &
     &                            0.0317, 0.0210, 0.026]
      real, parameter :: fcbd(6) =[-0.1434,   -0.2423,    0.0,           &
     &                            0.0,    0.0,    0.0]
      real, parameter :: fcwp(6) =[0.0,       0.0,       0.0,            &
     &                            0.0,    0.72,   0.41]
      real, parameter :: fcin(6) =[0.0,       0.0,       0.2576,         &
     &                            0.4118, 0.2391, 0.4103]
      real, parameter :: wpsa(6) =[-0.0059,   -0.0059,    0.0,           &
     &                            0.0,    0.0,    0.0]
      real, parameter :: wpsi(6) =[0.1142,    0.1142,    0.0,            &
     &                            0.0,    0.0,    0.0]
      real, parameter :: wpcl(6) =[0.5766,    0.5766,    0.50,           &
     &                            0.50,   0.0,    0.0]
      real, parameter :: wpom(6) =[2.228E-03, 2.228E-03, 0.0158,         &
     &                            0.0158, 0.0,    0.0]
      real, parameter :: wpbd(6) =[0.02671,   0.02671,   0.0,            &
     &                            0.0,    0.0,    0.0]
      real, parameter :: wpwp(6) =[0.0,       0.0,       0.0,            &
     &                            0.0,    1.0,    1.0]
      real, parameter :: wpin(6) =[0.0,       0.0,       0.0260,         &
     &                            0.0260, 0.0,    0.0]
      real      ompc, sand, silt, clay, bulkd, rockf
      integer   swflag

      real      afiel, awilt


c ... Field capacity and wilting point.  Computations based on
c ... Gupta and Larson 1979, 'Estimating soil and water retention
c ... characteristics from particle size distribution, organic
c ... matter percent and bulk density'. Water Resources Research 15:1633
c ... or Rawls et al (1982) 'Estimation of soil water properties'
c ... Trans. ASAE ???:1316
c ... Field capacity options of -0.1 or -0.33 bar.
c ... Wilting point assumed to be water content at -15 bars.
c ... Calculate organic matter from initial conditions, ivauto or
c ... value at the beginning of an extend
c ... Note that Gupta and Larson and Rawls use % for texture
c ... but values here are fractions.
          afiel = fcsa(swflag)*sand  + fcsi(swflag)*silt +
     &            fccl(swflag)*clay  + fcom(swflag)*ompc +
     &            fcbd(swflag)*bulkd + fcwp(swflag)*awilt +
     &            fcin(swflag)
          awilt = wpsa(swflag)*sand  + wpsi(swflag)*silt +
     &            wpcl(swflag)*clay  + wpom(swflag)*ompc +
     &            wpbd(swflag)*bulkd + wpwp(swflag)*awilt +
     &            wpin(swflag)
          ! modifiy afiel and awilt according to fractional volume of rock -mdh 5/27/99
          ! For really deep soils, don't use rock fraction -mdh 6/29/99
          if(rockf .gt. 0) then
            afiel = afiel * (1.0 - rockf)
            awilt = awilt * (1.0 - rockf)
          endif
        return
      end subroutine swmodel
