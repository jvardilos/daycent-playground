
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... FIXIN.F

      include 'n2oparam.inc'

      subroutine fixin

      USE n2oparam

      implicit none

      include 'const.inc'
      include 'parfx.inc'
      include 'sitsoil.inc'

! ... Set values of fixed parameters and initial values.

! ... Local variables
      integer     ii, jj
      real        stk(1), tmp
      character (len=16)  :: varfnd
      character (len=16)  :: rawval
      character (len=32)  :: mssg

      integer, parameter ::   unit = 8
      character (len=6), parameter :: callr = 'fixin'

      ! unused parameters read as local variables
      real :: wfdryt

! ... functions
      real        chkdata

      call oldopen(unit,'fix.100') ! open the fix file  with library search

      read(unit,*)

      do ii = 1, CMXLYR
        adep(ii) = int(chkdata(unit,callr,'adep',0,stk))
      enddo

      agppa = chkdata(unit,callr,'agppa',0,stk)
      agppb = chkdata(unit,callr,'agppb',0,stk)

      ! look for aneref or the replacement satdec parameters.
      aneref = [1.5, 3., 1.] ! these values from the 2011+ inventory disable the model.
      !satdec = [1.0, 0.45, 0.25, 0.45]  ! parton's recommended saturation decomposition (unused)
      varfnd = 'agppb';
      ii = 0;
      do while(varfnd .ne. 'EOF'  .and.  varfnd .ne. '')
        varfnd = ' '
        tmp = chkdata(unit,callr,varfnd,0,stk)

        if(index(varfnd,'aneref').gt.0) then
          ii = 1; ! we read ZNAREF values
          aneref(1) = tmp
          aneref(2) = chkdata(unit,callr,'aneref',0,stk)
          aneref(3) = chkdata(unit,callr,'aneref',0,stk)

        else if(index(varfnd,'satdec').gt.0) then
          call message('NOTE: saturated decomposition model not implemented: values ignored')
          ! satdec(1) = tmp                        ! (1)
          tmp = chkdata(unit,callr,'satdec',0,stk) ! (2)
          tmp = chkdata(unit,callr,'satdec',0,stk) ! (3)

          ! Bill didn't require a wood value but added just-in-case
          varfnd = ' '
          tmp = chkdata(unit,callr,varfnd,0,stk)   ! (4)
          if(index(varfnd,'satdec').eq.0) then
            tmp = chkdata(7,callr,'reread',0,stk)
          endif

        else
          tmp = chkdata(7,callr,'reread',0,stk)
          exit;   ! terminate this search loop
        endif
      end do
      ! Warn user if the ANAREF parameters were missed
      if(ii .eq. 0) then
        write(mssg,"('(',2(f6.3,','),f6.3,')')") aneref
        call message('Warning: fix.100, using default anaerobic decomposition  ANEREF = '//trim(mssg))
      endif

      animpt = chkdata(unit,callr,'animpt',0,stk)

      do ii = 1, CMXLYR
        awtl(ii) = chkdata(unit,callr,'awtl',0,stk)
      enddo

      bgppa = chkdata(unit,callr,'bgppa',0,stk)
      bgppb = chkdata(unit,callr,'bgppb',0,stk)

      co2ppm(1) = chkdata(unit,callr,'co2ppm',0,stk)
      co2ppm(2) = chkdata(unit,callr,'co2ppm',0,stk)
      co2rmp = chkdata(unit,callr,'co2rmp',0,stk)

      do ii = SRFC, SOIl
        do jj = 1, MAXIEL
          damr(ii,jj) = chkdata(unit,callr,'damr',0,stk)
        enddo
      enddo

      do ii = 1, MAXIEL
        damrmn(ii) = chkdata(unit,callr,'damrmn',0,stk)
      enddo

      do ii = SRFC, SOIl
        dec1(ii) = chkdata(unit,callr,'dec1',0,stk)
      enddo
      do ii = SRFC, SOIl
        dec2(ii) = chkdata(unit,callr,'dec2',0,stk)
      enddo
      do ii = SRFC, SOIl
        dec3(ii) = chkdata(unit,callr,'dec3',0,stk)
      enddo
      dec4 = chkdata(unit,callr,'dec4',0,stk)
      do ii = SRFC, SOIl
        dec5(ii) = chkdata(unit,callr,'dec5',0,stk)
      enddo
      deck5 = chkdata(unit,callr,'deck5',0,stk)
      dligdf = chkdata(unit,callr,'dligdf',0,stk)
      dresp = chkdata(unit,callr,'dresp',0,stk)
      edepth = chkdata(unit,callr,'edepth',0,stk)
      elitst = chkdata(unit,callr,'elitst',0,stk)
      enrich = chkdata(unit,callr,'enrich',0,stk)

      fixfavail(1) = chkdata(unit,callr,'favail',0,stk)
      fixfavail(2) = 0. ! make sure the unread P value is initialized
      fixfavail(3) = chkdata(unit,callr,'favail',0,stk)
      fixfavail(4) = chkdata(unit,callr,'favail',0,stk)
      fixfavail(5) = chkdata(unit,callr,'favail',0,stk)
      fixfavail(6) = chkdata(unit,callr,'favail',0,stk)
      favail(:,1) = fixfavail  ! save the default favail values from fix.100
      favail(:,2) = fixfavail  ! save the default favail values from fix.100

      do ii = 1, 5
        fleach(ii) = chkdata(unit,callr,'fleach',0,stk)
      enddo

      do ii = 1, 4
        fwloss(ii) = chkdata(unit,callr,'fwloss',0,stk)
      enddo

      fxmca = chkdata(unit,callr,'fxmca',0,stk)
      fxmcb = chkdata(unit,callr,'fxmcb',0,stk)
      fxmxs = chkdata(unit,callr,'fxmxs',0,stk)
      fxnpb = chkdata(unit,callr,'fxnpb',0,stk)
      gremb = chkdata(unit,callr,'gremb',0,stk)
      idef = int(chkdata(unit,callr,'idef',0,stk))

      do ii = 1, 3
        lhzf(ii) = chkdata(unit,callr,'lhzf',0,stk)
      enddo

      minlch = chkdata(unit,callr,'minlch',0,stk)
      if (minlch .gt. 5.0) then
        write(mssg,*) minlch
        call abortrun('minlch (in fix.100) is too high  minlch = '//trim(mssg))
      endif
      nsnfix = int(chkdata(unit,callr,'nsnfix',0,stk))
      ntspm = int(chkdata(unit,callr,'ntspm',0,stk))

      do ii = 1, 3
        omlech(ii) = chkdata(unit,callr,'omlech',0,stk)
      enddo
      if (omlech(3) .gt. 2.0) then
        write(mssg,*) omlech(3)
        call abortrun('omlech(3) (in fix.100) is too high  omlech(3) = '//trim(mssg))
      endif

      do ii = SRFC, SOIL
        p1co2a(ii) = chkdata(unit,callr,'p1co2a',0,stk)
      enddo

      do ii = SRFC, SOIL
        p1co2b(ii) = chkdata(unit,callr,'p1co2b',0,stk)
      enddo

      do ii = SRFC, SOIL
        p2co2(ii) = chkdata(unit,callr,'p2co2',0,stk)
      enddo

      p3co2 = chkdata(unit,callr,'p3co2',0,stk)
      pabres = chkdata(unit,callr,'pabres',0,stk)

      peftxa = chkdata(unit,callr,'peftxa',0,stk)
      peftxb = chkdata(unit,callr,'peftxb',0,stk)

      do ii = 1, 4
        phesp(ii) = chkdata(unit,callr,'phesp',0,stk)
      enddo

      do ii = SRFC, SOIL
        pligst(ii) = chkdata(unit,callr,'pligst',0,stk)
      enddo

      do ii = SRFC, SOIL
        pmco2(ii) = chkdata(unit,callr,'pmco2',0,stk)
      enddo

      do ii = 1, MAXIEL
        pmnsec(ii) = chkdata(unit,callr,'pmnsec',0,stk)
      enddo

      pmntmp = chkdata(unit,callr,'pmntmp',0,stk)
      pmxbio = chkdata(unit,callr,'pmxbio',0,stk)
      pmxtmp = chkdata(unit,callr,'pmxtmp',0,stk)

      do ii = 1, MAXIEL
        pparmn(ii) = chkdata(unit,callr,'pparmn',0,stk)
      enddo

      do ii = 1, 3
        pprpts(ii) = chkdata(unit,callr,'pprpts',0,stk)
      enddo

      do ii = SRFC, SOIL
        ps1co2(ii) = chkdata(unit,callr,'ps1co2',0,stk)
      enddo

      ps1s3(1) = chkdata(unit,callr,'ps1s3',0,stk)
      ps1s3(2) = chkdata(unit,callr,'ps1s3',0,stk)
      ps2s3(1) = chkdata(unit,callr,'ps2s3',0,stk)
      ps2s3(2) = chkdata(unit,callr,'ps2s3',0,stk)

      do ii = 1, MAXIEL
        psecmn(ii) = chkdata(unit,callr,'psecmn',0,stk)
      enddo

! ... Add new variable to fix.100 file for computing backflow from occluded P
! ... to secondary P, cak - 03/20/02
      psecoc1 = chkdata(unit,callr,'psecoc1',0,stk)
      psecoc2 = chkdata(unit,callr,'psecoc2',0,stk)

      do ii = 1, MAXIEL
        do jj = 1, 3
          rad1p(jj,ii) = chkdata(unit,callr,'rad1p',0,stk)
        enddo
      enddo

      do ii = 1, MAXIEL
        rcestr(ii) = chkdata(unit,callr,'rcestr',0,stk)
      enddo

      rictrl = chkdata(unit,callr,'rictrl',0,stk)
      riint = chkdata(unit,callr,'riint',0,stk)
      rsplig = chkdata(unit,callr,'rsplig',0,stk)

      seed = int(chkdata(unit,callr,'seed',0,stk))
      if (seed .gt. 0) seed = -seed ! Make sure seed is negative

      spl(INTCPT) = chkdata(unit,callr,'spl',0,stk)
      spl(SLOPE) = chkdata(unit,callr,'spl',0,stk)

      do ii = SRFC, SOIl
        strmax(ii) = chkdata(unit,callr,'strmax',0,stk)
      enddo

      do ii = 1, 5
        texepp(ii) = chkdata(unit,callr,'texepp',0,stk)
      enddo

      texesp(1) = chkdata(unit,callr,'texesp',0,stk)
      texesp(3) = chkdata(unit,callr,'texesp',0,stk)

! ... Added extra teff(*) parameter from new maintenance respiration code. -mdh
      teff(1) = chkdata(unit,callr,'teff',0,stk)
      teff(2) = chkdata(unit,callr,'teff',0,stk)
      teff(3) = chkdata(unit,callr,'teff',0,stk)
      teff(4) = chkdata(unit,callr,'teff',0,stk)

      tmelt(1) = chkdata(unit,callr,'tmelt',0,stk)
      tmelt(2) = chkdata(unit,callr,'tmelt',0,stk)
      if (tmelt(2) .gt. 0.5) then
        write(mssg,*) tmelt(2)
        call abortrun('tmelt(2) (in fix.100) is too high  tmelt(2) = '//trim(mssg))
      endif

      do ii = 1, MAXIEL
        do jj = 1,3
          varat11(jj,ii) = chkdata(unit,callr,'varat11',0,stk)
        enddo
      enddo
      do ii = 1, MAXIEL
        do jj = 1,3
          varat12(jj,ii) = chkdata(unit,callr,'varat12',0,stk)
        enddo
      enddo

      do ii = 1, MAXIEL
        do jj = 1,3
          varat21(jj,ii) = chkdata(unit,callr,'varat21',0,stk)
        enddo
      enddo
      do ii = 1, MAXIEL
        do jj = 1,3
          varat22(jj,ii) = chkdata(unit,callr,'varat22',0,stk)
        enddo
      enddo

      do ii = 1, MAXIEL
        do jj = 1,3
          varat3(jj,ii) = chkdata(unit,callr,'varat3',0,stk)
        enddo
      enddo

      vlosse = chkdata(unit,callr,'vlosse',0,stk)
! ... vlossg is now computed as a function of clay content, the vlossg
! ... parameter value read from the fix.100 file is used as a multiplier,
! ... see prelim subroutine, cak - 11/21/01
      vlossg_m = chkdata(unit,callr,'vlossg',0,stk)

      !#############################################################
      ! **TEST** READ NITRIFICATION EQUATION PARAMETERS FOR RAM
          cDno3   = (/9.23,  1.556, 76.91, 0.00222/); ! - nitrate parameters from Del Grosso
          cNsoilt = (/35.0, -5.0,    4.5,  7.0/);     ! - temperature effect parameters
          cNph    = (/5.0, 0.56, 1.0, 0.45/);         ! - ph parameters
          cNwfps  = (/30.0, -9.0/);                   ! - Water Filled Pore Space parameters exponential
         ! set the default values for the cultivation effect parameters
                                ! monthly daycent defaults
          maxcltef = 0.         ! maxcltef = 5.74
          xefcltef = 0.         ! xefcltef = 0.2  ! 15.0 default decomposition half life is 5 days
          cfita   = 0.          ! cfita  = 0.;
          cfitb   = 0.          ! cfitb  = 5.

          floodN2delay     =   7     ! TEST INPUT: days to transition to flooded N2/N2O ratio
          flood_N2toN2O    = 100.0   ! N2/N2O ratio for flooded state (100.0) (-1 disable)
          CO2_to_CH4       =   0.5;  ! fraction of heterotrophic soil respiration to CH4
          frCH4emit        =   0.55; ! carbohydrate anaerobic fermentation reaction with methanogenesis
          Aeh              =   0.23; ! methane differential coefficient (Aeh)
          Deh              =   0.16; ! methane differential coefficient (Deh)
          Beh_flood        = -250.0; ! methane Eh lower-limit during flooding
          Beh_drain        =  300.0; ! methane Eh upper-limit during drainage
          wfdryt           =     14; ! wfunc drying period (days)
          zero_root_frac   =    0.7; ! methane fraction emitted via bubbles at 0 root mass
          frac_to_exudates =   0.45; ! fraction of root production to root exudates (crop)
          ch4rootlim       =    1.0; ! mrtblm root biomass that starts to reduce methane bubble formation (crop)

          ! resetable in crop. save fixin inputs as defaults
          frexud = frac_to_exudates;
          methzr = zero_root_frac;
          mrtblm = ch4rootlim;

      ! Read the cNph array
      varfnd = ' '
      tmp = chkdata(unit,callr,varfnd,0,stk)
      do while(varfnd .ne. 'EOF'  .and. varfnd .ne. '')
        call message('FIX OPTION: '//trim(varfnd)//' = '//trim(rawval()))
        if(index(varfnd,'cnph') .ne. 0) then
          cNph(1) = tmp
          ! use the same name for the rest of the array
          jj = index(varfnd,'('); if(jj .ne. 0) varfnd(jj:) = " "
          do ii = 2,4
            cNph(ii) = chkdata(unit,callr,varfnd,0,stk)
          enddo

        ! Read the cNsoilt array
        else if(index(varfnd,'cnsoilt') .ne. 0) then
          cNsoilt(1) = tmp
          jj = index(varfnd,'(')
          if(jj .ne. 0) varfnd(jj:) = " "
          do ii = 2,4
            cNsoilt(ii) = chkdata(unit,callr,varfnd,0,stk)
          enddo

        ! Read the cDno3 array
        else if(index(varfnd,'cdno3') .ne. 0) then
          cDno3(1) = tmp
          ! use the same name for the rest of the array
          jj = index(varfnd,'(')
          if(jj .ne. 0) varfnd(jj:) = " "
          do ii = 2,4
            cDno3(ii) = chkdata(unit,callr,varfnd,0,stk)
          enddo

        ! Read the cNwfps array
        else if(index(varfnd,'cnwfps') .ne. 0) then
          cNwfps(1) = tmp
          ! use the same name for the rest of the array
          jj = index(varfnd,'(')
          if(jj .ne. 0) varfnd(jj:) = " "
          cNwfps(2) = chkdata(unit,callr,varfnd,0,stk)

        ! first read the efolding length xefcltef:
        ! xefcltef = ln(2)/(half life (month) * average decomp)
        else if (index(varfnd,'xefcltef') .ne. 0) then
          ! If this is 'xefcltef', save the value and read the next line
          if(tmp.gt.0) then
            xefcltef = 2**(1./(tmp * 30.4375))-1
            call message('temp adjusted cultivation effects')
          endif

        ! Do we have the 100% mixing clteff for the addition routines
        else if (index(varfnd,'maxcltef') .ne. 0) then
          ! process the data if we haven't disabled the exponential model
          maxcltef = tmp
          cfitb    = (maxcltef-1)/0.95

        ! variables for methanogenesis model
        else if (index(varfnd,'co2ch4').gt.0) then
          CO2_to_CH4   = tmp
        else if (index(varfnd,'mxch4f').gt.0) then
          frCH4emit = tmp
        else if (index(varfnd,'fldn2d').gt.0) then
          floodN2delay = nint(tmp)
          if(floodN2delay <0) then
            floodN2delay = 0
          else if(floodN2delay >21) then
            floodN2delay = 21
          endif
        else if (index(varfnd,'fln2or').gt.0) then
          flood_N2toN2O = tmp
          if(flood_N2toN2O < 0) then
            write(mssg,'(f10.4)') flood_N2toN2O
            call message('Warning: using surface layers to determine flooded N2/N20 ratio:'//trim(mssg))
          endif
        else if (index(varfnd,'frexud').gt.0) then
          frexud = tmp
          if(frexud .lt. 0.0  .or.  frexud .gt. 1.0) then
             write(mssg,'(f10.4)') frexud
             call message('frexud = '//trim(mssg)//' out of bounds: ')
             frac_to_exudates = max(min(frexud, 1.0), 0.0)
          endif
          frac_to_exudates = frexud ! raw value read from fix file
        else if (index(varfnd,'aeh').gt.0) then
          Aeh              = tmp
        else if (index(varfnd,'deh').gt.0) then
          Deh              = tmp
        else if (index(varfnd,'behfl').gt.0) then
          Beh_flood        = tmp
        else if (index(varfnd,'behdr').gt.0) then
          Beh_drain        = tmp
        else if (index(varfnd,'methzr').gt.0) then
          methzr = tmp
          if(methzr .lt. 0.0  .or.  methzr .gt. 1.0) then
             write(mssg,'(f10.4)') methzr
             call message('methzr = '//trim(mssg)//' out of bounds: ')
             zero_root_frac = max(min(methzr, 1.0), 0.0)
          endif
          zero_root_frac = methzr
        else if (index(varfnd,'mrtblm').gt.0) then
          mrtblm      = tmp;
          if(mrtblm .lt. 1.0) then
             write(mssg,'(f10.4)')  mrtblm
             call message("Warning: Low input for  mrtblm = "//trim(mssg)//' reset to 1.0')
             mrtblm = 1.0
          endif
          ch4rootlim  = mrtblm;
        ! End variables for methanogenesis model

        ! add precipitation soil dry time to control the wfunc_pulse dry time
        ! <=1 means every rain event > 0.5 cm will trigger unless a pulse is in progress
        else if (index(varfnd,'wfdryt').gt.0) then
          wfdryt   = max(nint(tmp+1),0)

        ! N20 pulse parameters denitrification respiration restraint
        else if (inpulse(tmp,varfnd)) then  ! turn off respiration restraint for layers

        else
           call message('FIX error: Unknown optional parameter: '//trim(varfnd))
        end if

        varfnd = ' '
        tmp = chkdata(unit,callr,varfnd,0,stk)
      end do

      if(xefcltef .ne. 0  .and.  xefcltef .ne. 0.) then
        call message('using interacting cultivations')
      else
        maxcltef = 0.0 ! make sure maxcltef and cfitb are not set without xefcltef
        cfitb  = 0.
      endif
!#########################################################

      close(unit=unit)

      return
      contains
        logical function inpulse(tmp,varfnd)
          real       tmp
          character (len=*)  :: varfnd

           real(c_float) :: tmp1, tmp2, plmfdd, plxfdd;
           integer :: tmi;                   ! tnstd Thaw N Snow Total Depth
           ! real :: decrate, decpool(8)
           ! integer :: lch, istat, ReadInt;

           ! ... Fortran to C prototype

           INTERFACE
             subroutine ssnowtrig(tnsad, tnstd) BIND(C) ! set original freeze/thaw pulse parameters
               USE ISO_C_BINDING
               real(c_float),   VALUE, INTENT(IN) :: tnsad, tnstd
             END subroutine ssnowtrig

             subroutine sfddtrig(plmfdd, plxfdd, ptrlyr) BIND(C) ! set pulse length frozen degree day limits
               USE ISO_C_BINDING
               real(c_float),   VALUE, INTENT(IN) :: plmfdd, plxfdd; ! pulse length frozen degree limits
               integer,         VALUE, INTENT(IN) :: ptrlyr         ! soil temperature trigger layer
             END subroutine sfddtrig

             subroutine stdpsf(tnpsif) BIND(C) ! set denitrification N2/N2O pore space inflection
               USE ISO_C_BINDING
               real(c_double),  VALUE, INTENT(IN) :: tnpsif
             END subroutine stdpsf

             subroutine stdpsfdd(dpsfdm, dpsfdx) BIND(C) ! set frozen degree day shift freezing degree day limits
               USE ISO_C_BINDING
               real(c_float),   VALUE, INTENT(IN) :: dpsfdm, dpsfdx; ! denitrification pore space frozen degree limits
             END subroutine stdpsfdd

             subroutine setdrrl(tndrrl) BIND(C) ! setting denitrification respiration restraint layer
               USE ISO_C_BINDING
               integer,       INTENT(IN) :: tndrrl
             END subroutine setdrrl

             subroutine sthwplsd(tplen) BIND(C) ! setting denitrification respiration restraint layer
               USE ISO_C_BINDING
               integer,       INTENT(IN) :: tplen
             END subroutine sthwplsd

             subroutine setsatfr(tnsatf) BIND(C) ! set N2N2O ratio saturated fraction
               USE ISO_C_BINDING
               real(c_float), INTENT(IN) :: tnsatf
             END subroutine setsatfr

             ! subroutine sethawdec(decrate) BIND(C)
             !   USE ISO_C_BINDING
             !   real(c_float), INTENT(IN) :: decrate(8)
             ! END subroutine sethawdec

           END INTERFACE

           inpulse = .true. ! assume we are going to handle this read
           ! N20 pulse parameters denitrification respiration restraint
           tmp1 = tmp
           if (index(varfnd,'tnpsif').gt.0) then  ! turn off respiration restraint for layers
             call stdpsf(dble(tmp)); ! tnpsf = readReal(0,inlin(:clen),lch,istat);

           ! thaw pulse maximum length
           else if (index(varfnd,'tplen').gt.0) then
             call sthwplsd(nint(tmp));    ! maximum number of days

           ! thaw pulse snow depth trigger
           else if (index(varfnd,'tnsad').gt.0) then
             tmp2  = chkdata(unit,callr,'tnstd',0,stk);  ! readReal(0,inlin(:clen),lch,istat)
             call ssnowtrig(tmp1, tmp2);  ! trigger snow level

           ! Pulse freezing degree day (temperature) trigger
           else if (index(varfnd,'plmfdd').gt.0) then
             plmfdd = tmp
             plxfdd = chkdata(unit,callr,'plxfdd',0,stk);  !  readReal(0,inlin(:clen),lch,istat)
             varfnd = ' '
             tmi   = int(chkdata(unit,callr,varfnd,0,stk))
             if (index(varfnd,'ptrlyr').gt.0) then
               call sfddtrig(plmfdd,plxfdd, tmi);     ! set frozen degree day pulse limiting (tnsad, tnstd)
             else
               call sfddtrig(plmfdd,plxfdd, 0);     ! set frozen degree day pulse limiting (tnsad, tnstd)
             endif

           ! pulse thermal (fdd) limiting
           else if (index(varfnd,'dpsfdm').gt.0) then  ! turn off respiration restraint for layers
             tmp2 = chkdata(unit,callr,'dpsfdx',0,stk);  ! arm snow level     readReal(0,inlin(:clen),lch,istat)
             call stdpsfdd(tmp1, tmp2);

           ! turn off respiration restraint for layers
           else if (index(varfnd,'tndrrl').gt.0) then
             call setdrrl(nint(tmp)) ! set the pulse respiration restraint layer

           ! mix fraction saturated and unsaturated N2/N20 ratio
           else if (index(varfnd,'tnsatf').gt.0) then
             call setsatfr(max(min(tmp,1.),0.)) ! set the pulse N2/N2o saturated mixing rate

           ! ! mix fraction saturated and unsaturated N2/N20 ratio
           ! else if (index(varfnd,'decrat').gt.0) then
           !   decrate = tmp;
           !   lch = index(varfnd,'(')+1
           !   if(lch.gt.1) then
           !   decpool = 1.0
           !     istat = 1;
           !     do while (lch > 1  .and.  lch <= len_trim(varfnd)  .and. istat>0)
           !       istat = 0;
           !       tmi = ReadInt(0,varfnd,lch,istat)
           !       if(istat<0) call abortrun("bad freeze/thaw decomposition pool"//trim(varfnd))
           !       if (istat == 0) exit
           !       decpool(tmi) = decrate
           !     end do
           !     call sethawdec(decpool) ! set the
           !   endif

           else ! oops don't know what this is
             inpulse = .false. ! signal to keep looking
           endif

           return

        end function inpulse
      end subroutine fixin

!  ^([ \t]+)read\(8,\*\) (\S+), *name\r\1call ckdata\('fixin', *'(\S+)', *name\)\r
!  \1\2 = chkdata(unit,callr,'\3',0,stk)

!  ^([ \t]+)read\(8,\*\) temp, *name\r\1(\S+) *= *int\(temp\)\r\1call ckdata\('fixin', *'(\S+)', *name\)\r
!  \1\2 = int(chkdata(unit,callr,'\3',0,stk))

! ^([ \t]+)do (\d++) (.*)\r([\s\S]*?)\2\s+continue
! \1do \3\r\4\1enddo
