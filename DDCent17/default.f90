
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine default()
      use calflow;

      implicit none
      include 'const.inc'
      include 'chrvar.inc'
      include 'comput.inc'
      include 'doubles.inc'
      include 'dovars.inc'
      include 'dynam.inc'
      include 'fertil.inc'
      include 'forrem.inc'
      include 'isovar.inc'
      include 'jday.inc'
      include 'ligvar.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'outval.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'
      include 'evivars.inc'

! ... This subroutine does a "brute force" initialization of all common block
! ... variables.  Some of these variables were being used in the code without
! ... initialization.  This does not cause a problem if the compiler
! ... initializes common block variables but we should not be depending on the
! ... compiler to initialize these variables.
!         numeric variables will be given an default value of 0
!         logical variables will be give an default value of .false.,
!         character variables will be given an default value of ' '.
!     - cak - 06/04/02
!
! Modification
!   This is overkill!
!   Recoded to replace explicit loops with Fortran array initialization.
!   This forces array declarations to determine  sizes rather than
!   hard coded loops which can get out of synch with definitions.
!   all such changes are commented with '! array'
!

! ... Fortran to C prototype
      INTERFACE

        SUBROUTINE floclr()
          !MS$ATTRIBUTES ALIAS:'_floclr' :: floclr
        END SUBROUTINE floclr

        SUBROUTINE floclr_double()
          !MS$ATTRIBUTES ALIAS:'_floclr_double' :: floclr_double
        END SUBROUTINE floclr_double

        SUBROUTINE floclr_double_in()
          !MS$ATTRIBUTES ALIAS:'_floclr_double_in' :: floclr_double_in
        END SUBROUTINE floclr_double_in

        SUBROUTINE floclr_double_out()
          !MS$ATTRIBUTES ALIAS:'_floclr_double_out' :: floclr_double_out
        END SUBROUTINE floclr_double_out

      END INTERFACE

! ... Local variables  (SHOULD BE NONE)

! ... chrvar common block
       cmdary = ' ' ! array
       typary = ' ' ! array

      curcrp = ' '
      curtre = ' '
      wthnam = ' '
      wthr = ' '

! ... comput common block
        agdefacm = 0.0 ! array
        bgdefacm = 0.0 ! array
      baseNdep = 0.0
            cercrp = 0.0 ! array
      eftext = 0.0
      fps1s3 = 0.0
      fps2s3 = 0.0
          lhzci = 0.0 ! array
          ! initialize both surface and soil ratios   MDH 9/24/2012
          ratnew1 = 0.0 ! array
          ratnew2 = 0.0 ! array
          lhze = 0.0 ! array
      orglch = 0.0
        p1co2 = 0.0 ! array
        h2ogef = 0.0 ! array
      wc = 0.0

! ... doubles common block
      occlud_double = 0.0
        secndy_double = 0.0 ! array

! ... dovars common block
      docult = .false.
      doerod = .false.
      dofert = .false.
        dofire = .false. ! array
      doflst  = .false.
      dofone  = .false.
      dofrst  = .false.
      dofngw  = .false.  ! reset in schedl
      dograz  = .false.
      dohrvt  = .false.
      domirri = .false.
      dolast  = .false.
      doomad  = .false.
      doplnt  = .false.
      dosene  = .false.
      dothrv  = .false.  ! reset in schedl
      dotrem  = .false.
      frstschd = .false.
      harvschd = .false.
      plntschd = .false.
      senmschd = .false.
      aceqcnt = -1      ! NO accelerated equilibrium
      cultday = 0
      erodday = 0
      fertday = 0
      fireday = 0
      flstday = 0
      foneday = 0
      frstday = 0
      grazday = 0
      hrvtday = 0
      irriday = 0
      lastday = 0
      omadday = 0
      plntday = 0
      seneday = 0
      tremday = 0
      cultcnt = 0
      fertcnt = 0
      erodcnt = 0
      grazcnt = 0
      irintvl = 0 ! no default irrigation
      irrcnt  = 0 ! no irrigation pending
      plntcnt = 0
      senecnt = 0
      savefrstday = 0
      saveplntday = 0

! ... dynam common block
        tree_cfrac = 0.0 ! array

! ... evivars common block
      eviday = 0.0
      eviFlag = .false.

! ... fertil common block
      aufert = 0.0
        feramt    = 0.0   ! array
        hrvafert = 0.0   ! array harvest autofert accumulator
        hrvfert  = 0.0   ! array harvest regular fert accumulator
        hrvomade = 0.0   ! array harvest omad minerals accumulator
      hrvirr = 0.0       !       harvest irrigation accumulator
      hrvomadc = 0.0     !       harvest omad C accumulator
      ninhib = 0.0
      ninhtm = 0
      nreduce = 1.0
      otffrc= 0.         ! orchard fertilizer split
        Nscalar = 1.0    ! array
        OMADscalar = 1.0 ! array
        fertnet   = 0.   ! array

! ... Initialize the flowstack and the number of flows variables
      call floclr()

! ... forrem common block
      evntyp = 0
        fd = 0.0 ! array
        remf = 0.0 ! array
        retf = 0.0 ! array

! ... isovar common block
      cisofr = 0.0
      cisotf = 0.0

! ... jdays common block
      dysimo = (/31,28,31,30,31,30,31,31,30,31,30,31/);        ! array
      lstdy  =  ilstdy        ! array
      frstdy = ifrstdy        ! array

! ... ligvar common block
        pltlig = 0.0          ! array

! ... monprd common block
        grspdyflux = 0.0      ! array
        mrspdyflux = 0.0      ! array
        mrspTempEffect = 0.0  ! array
        mrspWaterEffect = 0.0 ! array
        cgrspdyflux = 0.0     ! array
        cmrspdyflux = 0.0     ! array
        mcprd = 0.0           ! array
        fgrspdyflux = 0.0     ! array
        fmrspdyflux = 0.0     ! array
        mfprd = 0.0           ! array
      N2O_year = 0.0
      NO_year = 0.0
      N2_year = 0.0
      nit_amt_year = 0.0
      stempmth = 0.0
      annppt = 0.0
      N2O_month = 0.0
      NO_month = 0.0
      N2_month = 0.0
      nit_amt_month = 0.0
      pptmonth = 0.0
      cwstress = 0.0
      gwstress = 0.0
      twstress = 0.0

! ... npool common block
        nitrate = 0.0 ! array
      ammonium = 0.0
      frac_nh4_fert = 0.5 ! used in detiv subroutine, cak - 04/16/2007
      frac_no3_fert = 0.5 ! used in detiv subroutine, cak - 04/16/2007
!      texture = 0

! ... param common block
      afue = 0.5
        afiel = 0.0 ! array
        amov = 0.0 ! array
        awilt = 0.0 ! array
      basef = 0.0
      bulkd = 0.0
      cmix = 0.0
        autoresp1 = 0.0   ! array
        autoresp2 = 0.0   ! array
        co2ipr = 0.0      ! array
        co2irs = 0.0      ! array
        co2itr = 0.0      ! array
        co2tm = 0.0       ! array
        epnfa = 0.0       ! array
        epnfs = 0.0       ! array
        newautoresp = 0.0 ! array
        no3pref = 0.0     ! array
        npp2cs = 0.0      ! array
        prdx = 0.0        ! array
        satmos = 0.0      ! array
        snfxmx = 0.0      ! array
        co2ice = 0.0      ! array
      co2sys = 0.0
      drain = 0.0
      falprc = 0
      fracro = 0.0
        hpttr = 0.0    ! array
        htran = 0.0    ! array
        maxtmp = 0.0   ! array
        mintmp = 0.0   ! array
        pHscalar = 1.0 ! array
        prcskw = 0.0   ! array
        prcstd = 0.0   ! array
        precip = 0.0   ! array
      ivauto = 0
      labtyp = 0
      labyr = 0
        cmrspnpp = 0.0 ! array
        fgresp = 0.0   ! array
        fkmrspmx = 0.0 ! array
        fmrsplai = 0.0 ! array
      mctemp = 0.0
      micosm = 0
      nelem = 0
      Ninput = 0
      nlayer = 0
      nlaypg = 0
      Nstart = 0
      OMADinput = 0
      OMADstart = 0
      ph = 0.0
      phstart = 0.0
      phsys = 0
      phtm = 0
        ppdf = 0.0 ! array
      precro = 0.0
      psloss = 0.0
      pslsrb = 0.0
        rcelit = 0.0 ! array
        rces1 = 0.0 ! array
        rces2 = 0.0 ! array
        cgresp = 0.0 ! array
        ckmrspmx = 0.0 ! array
        rces3 = 0.0 ! array
      remwsd = 0.0
      rock = 0.0
      satmt = 0.0
      sirri = 0.0
      sorpmx = 0.0
      stamt = 0.0
      stormf = 0.0
      strm5l = 0.0
      strm5u = 0.0
      ststart = 0.0
      stsys = 0.0
      swflag = 0
      trbasl = 0.0
      claypg = 0
      claypg_const = 0
      tlaypg = 0
      tmix = 0.0

! ... parcp common block
      aglivb = 0.0
        astrec = 0.0 ! array
        crprtf = 0.0 ! array
        efrgrn = 0.0 ! array
        fecf = 0.0 ! array
        gret = 0.0 ! array
      aglrem = 0.0
      astgc = 0.0
      astlbl = 0.0
      astlig = 0.0
      auirri = 0
      awhc = 0.0
      basfc2 = 0.0
      bglrem = 0.0
      bioflg = 0
      biok5 = 0.0
      biomax = 0.0
        cfrtcn = 0.0 ! array
        cfrtcw = 0.0 ! array
        fnue = 0.0   ! array
        himon = 0    ! array
        clteff = 0.0 ! array
        fdfrem = 0.0 ! array
        fsdeth = 0.0 ! array
      cmxturn = 0.0
      crpgrw = 0
        cultra = 0.0 ! array
      eMax = 0.0
      fallrt = 0.0
      fawhc = 0.0
      fdgrem = 0.0
      feclig = 0.0
      flfrem = 0.0
      flghrv = 0
      flgrem = 0.0
        fligni = 0.0 ! array
        fret = 0.0   ! array
        frtc = 0.0   ! array
      frtcindx = 0
      frtsh = 0.0
      fulcan = 0.0
      gfcret = 0.0
      grwprc = 0.0
      grzeff = 0
      hibg = 0.0
      himax = 0.0
      hiwsf = 0.0
      irramt = 0.0
      mrtfrac = 0.0
      pltmrf = 0.0
        pramn = 0.0 ! array
        pramx = 0.0 ! array
        prbmn = 0.0 ! array
        prbmx = 0.0 ! array
      rdrj = 0.0
      rdrm = 0.0
      rdsrfc = 0.0
      rmvstr = 0.0
      rtdtmp = 0.0
      sdethc = 0.0
      seedl = 0
      sfclit = 0.0
      stdead = 0.0
      vlossp = 0.0

! ... parfs common block
      basfct = 0.0
      btolai = 0.0
        ccefor = 0.0  ! array
        cerfor = 0.0  ! array
      decid = 0
      decw1 = 0.0
      decw2 = 0.0
      decw3 = 0.0
        fcfrac = 0.0  ! array
      forgrw = 0
        forrtf = 0.0  ! array
      klai = 0.0
      laitop = 0.0
      ldrmlt = 0.0
        leafdr = 0.0  ! array
      maxlai = 0.0
      maxldr = 0.0
      maxnp = 0.0
      sapk = 0.0
      swold = 0.0
      tmxturn = 0.0
        tfrtcn = 0.0  ! array
        tfrtcw = 0.0  ! array
        wdlig = 0.0   ! array
        wooddr = 0.0  ! array
      wmrtfrac = 0.0
      woodb = 0.0
      wrdsrfc = 0.0
        cwscoef = (/ 0.378, 9.0 /) ! array

! ... parfx common block
      sdco2sum = 0.0
      ntdco2sm = 0.0
        adep = 0.0    ! array
        awtl = 0.0    ! array
      agppa = 0.0
      agppb = 0.0
        aneref = 0.0  ! array
        damrmn = 0.0  ! array
        lhzf   = 0.0  ! array
        omlech = 0.0  ! array
        pmnsec = 0.0  ! array
        pparmn = 0.0  ! array
        pprpts = 0.0  ! array
        psecmn = 0.0  ! array
        rcestr = 0.0  ! array
        texesp = 0.0  ! array
      animpt = 0.0
      bgppa = 0.0
      bgppb = 0.0
        co2ppm = 0.0  ! array
        dec1 = 0.0    ! array
        dec2 = 0.0    ! array
        dec3 = 0.0    ! array
        dec5 = 0.0    ! array
        p1co2a = 0.0  ! array
        p1co2b = 0.0  ! array
        p2co2 = 0.0   ! array
        pligst = 0.0  ! array
        pmco2 = 0.0   ! array
        ps1co2 = 0.0  ! array
        ps1s3 = 0.0   ! array
        ps2s3 = 0.0   ! array
        spl = 0.0     ! array
        strmax = 0.0  ! array
        tmelt = 0.0   ! array
      co2rmp = 0.0
        damr = 0.0    ! array
      dec4 = 0.0
      deck5 = 0.0
      dligdf = 0.0
      dresp = 0.0
      edepth = 0.0
      elitst = 0.0
      enrich = 0.0
        favail = 0.0  ! array
        fleach = 0.0  ! array
        texepp = 0.0  ! array
        fwloss = 0.0  ! array
        phesp = 0.0   ! array
        teff = 0.0    ! array
      fxmca = 0.0
      fxmcb = 0.0
      fxmxs = 0.0
      fxnpb = 0.0
      gremb = 0.0
      idef = 0
      minlch = 0.0
      nsnfix = 0
      ntspm = 0
      p3co2 = 0.0
      pabres = 0.0
        rad1p = 0.0     ! array
        varat11 = 0.0   ! array
        varat12 = 0.0   ! array
        varat21 = 0.0   ! array
        varat22 = 0.0   ! array
        varat3 = 0.0    ! array
      peftxa = 0.0
      peftxb = 0.0
      pmntmp = 0.0
      pmxbio = 0.0
      pmxtmp = 0.0
      psecoc1 = 0.0
      psecoc2 = 0.0
      rictrl = 0.0
      riint = 0.0
      rsplig = 0.0
      seed = 0
      vlosse = 0.0
      vlossg = 0.0
      vlossg_m = 0.0
        twscoef = (/ 0.378, 9.0 /) ! array

! ... pheno common block
      accumdd = .false.
      basetemp = 0.0    ! array
      cgrwdys = 0
      clsgres = 0.0
      curgdys = 0
      dayhrs = 0.0
      ddbase = 0.0
      decidgrow = .false.
      fgrwdys = 0
      flsgres = 0.0
      furgdys = 0
      grnfill = .false.
      grnfldys = 0
      grnhrvt = .false.
      hrsinc = .false.
      mnddhrv = 0.0
      mxddhrv = 0.0
      plntkill = .false.
      soiltavewk = 0.0
      thermunits = 0.0
      tmpgerm = 0.0
      tmpkill = 0.0
      tmplff = 0.0
      tmplfs = 0.0

! ... plot common blocks
      vals1 = 0.0  ! plot1
      vals2 = 0.0  ! plot2
      vals3 = 0.0  ! plot2

! ... potent common block
      agp = 0.0
      tgprod = 0.0
      pcropc = 0.0
      pforc = 0.0
        crop_a2drat = 0.0 ! array
        tree_a2drat = 0.0 ! array

! ... schvar common block
      evtptr = 0
      rptyrs = 0
!      fltary = 0.0   !  removed from program
      timary = 0      ! array
      ttlind = 0

! ... seq common block
      cursys = 0
      decsys = 0
      otfrac = -1

! ... site common block
      sitlat = 0.0
      sitlng = 0.0
      sitpot = 0.0
      sand = 0.0
      silt = 0.0
      clay = 0.0
      sitpot_m = 0.0


! ... timvar common block
      tend = 0
      dtpl = 0.0
      dt = 0.0
      strtyr = 0
      blktnd = 0
      strplt = 0.0
      tplt = 0.0
      decodt = 0.0
      daylength = 0.0    ! array

! ... wth common block
        prcurr = 0.0     ! array
        prcnxt = 0.0     ! array
        precscalar = 1.0 ! array
        tmaxscalar = 0.0 ! array
        tminscalar = 0.0 ! array
        tmn2m = 0.0      ! array
        tmx2m = 0.0      ! array
      maxt = 0.0
      wthinput = 0
      wthstart = 0

! ... wthdaily common block
        avgtemp = 0.0    ! array
        tempmax = 0.0    ! array
        tempmin = 0.0    ! array
        ppt = 0.0        ! array
        solrad = 0.0     ! array
        srad = 0.0       ! array
        rhumid = 0.0     ! array
        windsp = 0.0     ! array
        leapyr  = .false.
        yrmatch = .false.

! ... zztim common block
      time = 0.0

      return
      end
