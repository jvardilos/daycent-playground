
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

! modifications
!   5 Sept 2013  K. Killian
!      Added input checking to catch and convert the water stress parameters
!      used by Bill and Cindy to the independent form.
!   1Feb13    K. Killian
!      added optional lines for water stress coefficients wscoef(2)
!      added optional lines for nutrient availability favail (5)
!   4Aug2011 K. Killian
!     introduced crop stack using the new stack management routine stakfind
!     Code the full crop name in crpval instead of just the first 5 characters

      subroutine cropin(tomatch)
      USE ISO_C_Binding
      implicit none

! ... Argument declarations
      character*28 tomatch

!...Read in the new crop type
!
      include 'chrvar.inc'
      include 'const.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'seq.inc'
      include 'sitsoil.inc'

! ... Local variables
      integer   :: ii, jj
      integer   :: m, is, unit
      real      :: del13c, tmpval
      character (len=132) :: string
      character (len=16)  :: varfnd

      ! functions
      integer (kind=4) :: crophash
!       INTERFACE
!         real function chkdata(unit,routin,expect,m,stk, elin)
!           integer   :: unit, m
!           character :: routin*(*), expect*(*)
!           real      :: stk(*)
!           logical,  optional :: elin
!         END  FUNCTION chkdata
!       END INTERFACE

      ! Number of lines to read for each crop type; (CROPLNS = 110)
      ! optional lines favail(6);                   (OPTLN = 5)
! ... stack variables
      ! Number of lines to read for each crop type; (CROPLNS = 109)
      ! optional lines favail(6), cwscoef(2) tmxbio;  (OPTLN = 9)
      integer, parameter       :: REQLNS = 109, OPTLN = 12, SDEPTH = 8
      real, save               :: stack(REQLNS+OPTLN, SDEPTH)
      integer, save            :: Tstack(SDEPTH) = 0
      character (len=28), save :: stackname(SDEPTH) = '######'
      integer                  :: stakfind
      real                     :: chkdata ! , value
      character (len=6), parameter :: callr = 'cropin'

      unit = 0

!...Check the stack for the data
      is = stakfind(tomatch,SDEPTH,Tstack,stackname)

      if(is.lt.0) then
        is = abs(is)
        stackname(is) = tomatch        ! record the name in the name array

        unit = 11
        call oldopen(unit,'crop.100')  ! open input file
        call findopt(unit,'cropin.100',tomatch,REQLNS)  ! find the option in the crop file
      endif

!         write(*,'(a,i5,a,7x,8(2x,a6))') 'entry ',is,' is '//tomatch  &
!     &           ,(stackname(max(1,Tstack(ii))), ii=1,SDEPTH)

        m = 1

        prdx(1) = chkdata(unit,callr,'prdx',m,stack(1,is))


        do ii = 1, 4
          ppdf(ii,1)= chkdata(unit,callr,'ppdf',m,stack(1,is))
        end do

        bioflg = int(chkdata(unit,callr,'bioflg',m,stack(1,is)))
        biok5  = chkdata(unit,callr,'biok5',m,stack(1,is))
        pltmrf = chkdata(unit,callr,'pltmrf',m,stack(1,is))
          call boundscheck(pltmrf, 'pltmrf', "planting reduction factor")
        fulcan = chkdata(unit,callr,'fulcan',m,stack(1,is))

! ..... Added frtcindx and frtc(4) for dynamic-C allocation -mdh 8/23/00
        frtcindx = nint(chkdata(unit,callr,'frtcindx|frtcin',m,stack(1,is)))
        ! Allow frtcindx  3, 4, 5, 6 for the growing degree day implementation
        if (frtcindx .lt. 0 .or. frtcindx .gt. 6) then
          write(string,'(i5)') frtcindx
          call abortrun('cropin option '//trim(tomatch)//' out of range frtcindx ='//string(:5))
        endif

! ..... Add frtc(5) for annual plants, cak - 09/12/03
        do ii = 1, 5
          frtc(ii)  = chkdata(unit,callr,'frtc',m,stack(1,is))
        end do
        if ((frtcindx .eq. 2) .or. (frtcindx .ge. 4)) then
          if ((frtc(1) .le. 0.0) .or. (frtc(1) .gt. 1.0)) then
            write(string,*) 'frtc(1) = ', frtc(1)
            call abortrun('crop '//trim(tomatch)//'input out of bounds  0.0 < frtc(1) <= 1.0'//trim(string))
          endif
          if ((frtc(2) .le. 0.0) .or. (frtc(2) .gt. 1.0)) then
            write(string,*) 'frtc(2) = ', frtc(2)
            call abortrun('crop '//trim(tomatch)//'input out of bounds  0.0 < frtc(2) <= 1.0'//trim(string))
          endif
          if (frtc(3) .le. 0.0) then
            write(string,*) 'frtc(3) = ', frtc(3)
            call abortrun('crop '//trim(tomatch)//'input out of bounds frtc(3) > 0.0 '//trim(string))
          endif
          if ((frtc(4) .le. 0.0) .or. (frtc(4) .gt. 1.0)) then
            write(string,*) 'frtc(4) = ', frtc(4)
            call abortrun('crop '//trim(tomatch)//' input out of bounds  0.0 < frtc(4) <= 1.0'//trim(string))
          endif
          if ((frtc(5) .le. 0.0) .or. (frtc(5) .gt. 1.0)) then
            write(string,*) 'frtc(5) = ', frtc(4)
            call abortrun('crop '//trim(tomatch)//'input out of bounds  0.0 < frtc(5) <= 1.0'//trim(string))
          endif
        endif

! ..... Added cfrtcn(*) and cfrtcw(*) for calculating water and nutrient
! ..... stress using unique parameter values, cak - 09/12/03
        cfrtcn(1) = chkdata(unit,callr,'cfrtcn',m,stack(1,is))
          call boundscheck(cfrtcn(1), 'cfrtcn(1)', "fraction of C allocated under nutrient stress")
        cfrtcn(2) = chkdata(unit,callr,'cfrtcn',m,stack(1,is))
          call boundscheck(cfrtcn(2), 'cfrtcn(2)', "fraction of C allocated with no nutrient stress")
        cfrtcw(1) = chkdata(unit,callr,'cfrtcw',m,stack(1,is))
          call boundscheck(cfrtcw(1), 'cfrtcw(1)', "fraction of C allocated under water stress")
        cfrtcw(2) = chkdata(unit,callr,'cfrtcw',m,stack(1,is))
          call boundscheck(cfrtcw(2), 'cfrtcw(2)', "fraction of C allocated with no water stress")

        if ((frtcindx .eq. 1) .or. (frtcindx .eq. 3)) then
          if ((cfrtcn(1) .le. 0.0) .or. (cfrtcn(2) .le. 0.0) .or.        &
     &        (cfrtcw(1) .le. 0.0) .or. (cfrtcw(2) .le. 0.0) .or.        &
     &        (cfrtcn(1) .gt. 1.0) .or. (cfrtcn(2) .gt. 1.0) .or.        &
     &        (cfrtcw(1) .gt. 1.0) .or. (cfrtcw(2) .gt. 1.0)) then
            write(string,'("cfrtcn =",2f6.2,"  cfrtcw =",2f6.2) ') cfrtcn, cfrtcw
            call abortrun('crop '//trim(tomatch)//' out of bounds '//     &
     &      'cfrtcn(*) & cfrtcw(*) must be > 0 and <=1 '//trim(string))
          endif
          if ((cfrtcn(2) .gt. cfrtcn(1))) then
            write(string,'("cfrtcn(1),",f4.2," must be >= cfrtcn(2),",f4.2) ') cfrtcn
            call abortrun('crop '//trim(tomatch)//' out of bounds '//trim(string))
          endif
          if ((cfrtcw(2) .gt. cfrtcw(1))) then
            write(string,'("cfrtcw(1),",f4.2," must be >= cfrtcw(2),",f4.2) ') cfrtcw
            call abortrun('crop '//trim(tomatch)//' out of bounds '//trim(string))
          endif
        endif

        biomax = chkdata(unit,callr,'biomax',m,stack(1,is))

        do ii = 1, 2
          do jj = 1, MAXIEL
            pramn(jj,ii) = chkdata(unit,callr,'pramn',m,stack(1,is))
          end do
        end do

        do ii = 1, 2
          do jj = 1, MAXIEL
            pramx(jj,ii) = chkdata(unit,callr,'pramx',m,stack(1,is))
          end do
        end do

        do ii = 1, 2
          do jj = 1, MAXIEL
            prbmn(jj,ii) = chkdata(unit,callr,'prbmn',m,stack(1,is))
          end do
        end do

        do ii = 1, 2
          do jj = 1, MAXIEL
            prbmx(jj,ii) = chkdata(unit,callr,'prbmx',m,stack(1,is))
          end do
        end do

        do ii = ABOVE, BELOWM
          fligni(INTCPT,ii) = chkdata(unit,callr,'fligni',m,stack(1,is))
          call boundscheck(fligni(INTCPT,ii), 'fligni(INTCPT,'//itoindx(ii)//')', &
                            "equation intercept - lignin content fraction")
          fligni(SLOPE,ii)  = chkdata(unit,callr,'fligni',m,stack(1,is))
          call boundscheck(fligni(INTCPT,ii), 'fligni(INTCPT,'//itoindx(ii)//')', &
                            "equation slope     - lignin content fraction")
        end do

        himax = chkdata(unit,callr,'himax',m,stack(1,is))
          call boundscheck(himax, 'himax',"harvest index maximum")
        hiwsf = chkdata(unit,callr,'hiwsf',m,stack(1,is))
        himon(1) = int(chkdata(unit,callr,'himon',m,stack(1,is)))
        himon(2) = int(chkdata(unit,callr,'himon',m,stack(1,is)))

        do ii = 1, MAXIEL
          efrgrn(ii) = chkdata(unit,callr,'efrgrn',m,stack(1,is))
          call boundscheck(efrgrn(ii), 'efrgrn('//itoindx(ii)//')', &
                           "aboveground E fraction going to grain")
        end do

        vlossp = chkdata(unit,callr,'vlossp',m,stack(1,is))

        do ii = 1, 4
          fsdeth(ii) = chkdata(unit,callr,'fsdeth',m,stack(1,is))
        end do

        fallrt = chkdata(unit,callr,'fallrt',m,stack(1,is))
        rdrj = chkdata(unit,callr,'rdrj',m,stack(1,is))
          call boundscheck(rdrj, 'rdrj', "maximum juvenile fine root death rate")
        rdrm = chkdata(unit,callr,'rdrm',m,stack(1,is))
          call boundscheck(rdrm, 'rdrm', "maximum mature fine root death rate")
        rdsrfc = chkdata(unit,callr,'rdsrfc',m,stack(1,is))
          call boundscheck(rdsrfc, 'rdsrfc', "dead fine roots transferred to surface litter")
        rtdtmp = chkdata(unit,callr,'rtdtmp',m,stack(1,is))

        do ii = 1, MAXIEL
          crprtf(ii) = chkdata(unit,callr,'crprtf',m,stack(1,is))
          call boundscheck(crprtf(ii), 'crprtf('//itoindx(ii)//')', "leaf E transferred to storage at death")
        end do

        mrtfrac = chkdata(unit,callr,'mrtfrac',m,stack(1,is))
          call boundscheck(mrtfrac, 'mrtfrac', "fine root production fraction")
        snfxmx(CRPSYS) = chkdata(unit,callr,'snfxmx',m,stack(1,is))
          call boundscheck(snfxmx(CRPSYS), 'snfxmx(CRPSYS)', "maximum symbiotic N fixation fraction")
        del13c = chkdata(unit,callr,'del13c',m,stack(1,is))
        co2ipr(CRPSYS) = chkdata(unit,callr,'co2ipr',m,stack(1,is))
        co2itr(CRPSYS) = chkdata(unit,callr,'co2itr',m,stack(1,is))

        do ii = IMIN, IMAX
          do jj = 1, MAXIEL
            co2ice(CRPSYS,ii,jj) = chkdata(unit,callr,'co2ice',m,stack(1,is))
          end do
        end do

        co2irs(CRPSYS) = chkdata(unit,callr,'co2irs',m,stack(1,is))

! ..... Added ckmrspmx parameters for maintenance respiration code,
! ..... mdh - 11/30/01
        ckmrspmx(ABOVE)  = chkdata(unit,callr,'ckmrspmx',m,stack(1,is))
        ckmrspmx(BELOWJ) = chkdata(unit,callr,'ckmrspmx',m,stack(1,is))
        ckmrspmx(BELOWM) = chkdata(unit,callr,'ckmrspmx',m,stack(1,is))

! ..... Added parameters for controlling the linear decrease in
! ....  maintenance respiration as the amount of carbohydrate stored in
! ..... the carbohydrate storage pool gets smaller, cak - 01/08/2010
        cmrspnpp(1) = chkdata(unit,callr,'cmrspnpp',m,stack(1,is))
        cmrspnpp(2) = chkdata(unit,callr,'cmrspnpp',m,stack(1,is))
        cmrspnpp(3) = chkdata(unit,callr,'cmrspnpp',m,stack(1,is))
        cmrspnpp(4) = chkdata(unit,callr,'cmrspnpp',m,stack(1,is))
        cmrspnpp(5) = chkdata(unit,callr,'cmrspnpp',m,stack(1,is))
        cmrspnpp(6) = chkdata(unit,callr,'cmrspnpp',m,stack(1,is))

! ..... Added cgresp parameters for the growth respiration code,
! ..... cak - 01/16/2007
        cgresp(ABOVE) = chkdata(unit,callr,'cgresp',m,stack(1,is))
        cgresp(BELOWJ) = chkdata(unit,callr,'cgresp',m,stack(1,is))
        cgresp(BELOWM) = chkdata(unit,callr,'cgresp',m,stack(1,is))

! ..... Added no3pref, mdh - 9/11/01
        no3pref(CRPSYS) = chkdata(unit,callr,'no3pref',m,stack(1,is))
          call boundscheck(no3pref(CRPSYS), 'no3pref(1)', "NO3 fraction of N update")
        if (no3pref(CRPSYS) .lt. 0. .or. no3pref(CRPSYS) .gt. 1.) then
          write(string,*) 'no3pref = ', no3pref
          call abortrun('crop.100 input out of bounds  0.0 < no3pref <= 1.0'//trim(string))
        endif

! ..... Added claypg, cak - 01/29/03
         claypg_const = chkdata(unit,callr,'claypg',m,stack(1,is))
         if(claypg_const .gt. nlayer .and. unit .ne.0) then
            write(varfnd,'(2I3)') claypg_const, nlayer
            call message('NOTICE: claypg for '//trim(tomatch)//','//varfnd(:3)// &
     &      ' exceeds the site soil depth  nlayer '//varfnd(4:6))
            claypg_const = nlayer
            stack(m-1,is) = nlayer
         end if
! ..... For an annual plant initialize the rooting depth to 1
        if ((frtcindx .eq. 2) .or. (frtcindx .ge. 4)) then
          claypg = 1
        else
          claypg = claypg_const
        endif

! ..... Added cmix, cak - 06/14/05
        cmix = chkdata(unit,callr,'cmix',m,stack(1,is))

! ..... Added tmpgerm, ddbase and tmpkill, cak - 04/17/03
        tmpgerm = chkdata(unit,callr,'tmpgerm',m,stack(1,is))
        ddbase = chkdata(unit,callr,'ddbase',m,stack(1,is))
        tmpkill = chkdata(unit,callr,'tmpkill',m,stack(1,is))

! ..... Added basetemp, mnddhrv, and mxddhrv, cak - 06/01/05
! ..... Change the basetemp variable to a two member array, cak - 05/21/08
        basetemp(1) = chkdata(unit,callr,'basetemp',m,stack(1,is))
        basetemp(2) = chkdata(unit,callr,'basetemp',m,stack(1,is))
        mnddhrv = chkdata(unit,callr,'mnddhrv',m,stack(1,is))
        mxddhrv = chkdata(unit,callr,'mxddhrv',m,stack(1,is))

! ..... Add parameters used to restrict production late in the growing
! ..... season, cak - 03/11/2010
        curgdys = int(chkdata(unit,callr,'curgdys',m,stack(1,is)))
        clsgres = chkdata(unit,callr,'clsgres',m,stack(1,is))
          call boundscheck(clsgres, 'clsgres', "late season growth restriction factor")

! ..... Added cmxturn, cak - 06/28/2007
        cmxturn = chkdata(unit,callr,'cmxturn',m,stack(1,is))
          call boundscheck(cmxturn, 'cmxturn', "maximum juvenile to mature fine roots fraction")

        ! Added water stress equation coefficents, cak - 07/10/2012
        varfnd = ''
        cwscoef(1) = chkdata(unit,callr,varfnd,m,stack(1,is))
        ! set the rest of wscoef from the file or the crop stack
        if(unit .eq. 0  .or.  index(varfnd,'wscoef') .gt. 0) then
          cwscoef(2) = chkdata(unit,callr,'wscoef',m,stack(1,is))
          if(cwscoef(1) .gt. 1.0  .and.  cwscoef(2) .lt. 0) then
            write(string,'(2f5.1)') cwscoef
            cwscoef(2) = -cwscoef(2)
            cwscoef(1) = log(cwscoef(1))/cwscoef(2)
            stack(m-2:m-1, is) = cwscoef   ! put revised values on stack
            write(varfnd,'(2f6.3)') cwscoef
            call message("converting deprecated "//trim(tomatch)// &
                 " WSCOEF values from"//trim(string)//" to "//varfnd)
          endif
        else  ! reading file;  wscoef wasn't found. (default file condition)
          cwscoef(1) = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
          ! set default conditions
          cwscoef = (/ 0.378, 9. /) ! (/ 0.5, 9.0 /) 2011 inventory
          stack(m-1:m, is) = cwscoef   ! put default values on the stack
          m = m+1
        endif
!     write(string,'(2f7.3)') cwscoef
!     call message('cwscoef Check'//trim(string)//" in crop "// trim(tomatch))

! ..... Added npp2cs, cak - 03/16/2011
        npp2cs(CRPSYS) = chkdata(unit,callr,'npp2cs',m,stack(1,is))

! ..... Added tmxbio, cak - 02/15/2012
        varfnd = ''
        tmxbio = chkdata(unit,callr,varfnd,m,stack(1,is))
        if(unit .ne. 0  .and.  index(varfnd,'tmxbio') .eq. 0) then
          ! default: reading file and variable not found.
          tmxbio = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
          tmxbio = 1200.0 ! 2.0
!        else
!          tmxbio = 1200.0
          stack(m-1, is) = tmxbio
        endif

        ! Depricated Added cafue, cak - 02/21/2011
        varfnd = ''
        tmpval = chkdata(unit,callr,varfnd,m,stack(1,is))
        if(unit .ne. 0  .and.  index(varfnd,'cafue') .eq. 0) then
          tmpval = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
          call message("Warning cafue input in "//trim(tomatch)//" is depricated")
        else
          ! default: reading file and variable not found.
        endif

        ! favail let the crop modify the favail fix parameter for delGgrosso Oct 2012
        ! the indexing on this is beastly. favail(2) is **NOT** input.
        ! thus look for favail(1) or read that position from the stack
        !  If found then set 2 to zero and input 3-6 from stack or file
        ! NOT FOUND. Then use the fix.100 values.
        !            The ugly part is we only put 1 and 3-6 on stack.
        ! Cindy added the N availability as a standalone so break this into 2 tests 30Jan13
        varfnd = ''
        favail(1,1) = chkdata(unit,callr,varfnd,m,stack(1,is))
        if(unit .ne. 0  .and.  index(varfnd,'favail') .eq. 0) then
          ! default: reading file and variable not found.
          favail(1,1) = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
          favail(1,1) = fixfavail(1)
          stack(m-1, is)   = fixfavail(1)   ! put favail(1) on the stack
          stack(m:m+3, is) = fixfavail(3:6) ! put favail(3:6) on the stack
          m = m+4
        else
          ! now look for the rest of the favail array
          varfnd = ''
          favail(3,1) = chkdata(unit,callr,varfnd,m,stack(1,is))
          ! set the rest of favail from the file or the crop stack
          if(unit .eq. 0  .or.  index(varfnd,'favail') .gt. 0) then
            favail(2,1) = 0 ! this will be initialized later
            favail(4,1) = chkdata(unit,callr,'favail',m,stack(1,is))
            favail(5,1) = chkdata(unit,callr,'favail',m,stack(1,is))
            favail(6,1) = chkdata(unit,callr,'favail',m,stack(1,is))
          else  ! reading file but the favail wasn't found. (default file condition)
            favail(3,1)        = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
            favail(2:6,1)      = fixfavail(2:6) ! use values from fix.100
            stack(m-1:m+2, is) = fixfavail(3:6) ! put favail(3:6) on the stack
            m = m+3
          endif
        endif

! ..... Added eMax for EVI calculation rm - 12/2009
        eMax = chkdata(unit,callr,'emax',m,stack(1,is))

        ! plant specific updates to methane values
        varfnd = ''
        frac_to_exudates = chkdata(unit,callr,varfnd,m,stack(1,is))
        if(unit .ne. 0  .and.  index(varfnd,'frexud') .eq. 0) then
          ! default: reading file and variable not found.
          tmpval = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
          frac_to_exudates = frexud
          stack(m-1, is) = frexud
        endif

        varfnd = ''
        ch4rootlim = chkdata(unit,callr,varfnd,m,stack(1,is))
        if(unit .ne. 0  .and.  index(varfnd,'mrtblm') .eq. 0) then
          ! default: reading file and variable not found.
          tmpval = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
          ch4rootlim = mrtblm
          stack(m-1, is) = mrtblm
        endif

        varfnd = ''
        zero_root_frac = chkdata(unit,callr,varfnd,m,stack(1,is))
        if(unit .ne. 0  .and.  index(varfnd,'methzr') .eq. 0) then
          ! default: reading file and variable not found.
          tmpval = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
          zero_root_frac = methzr
          stack(m-1, is) = methzr
        endif

!...Close the file
      if (unit.ne.0) close(unit)

!...Hold on to the current crop just read in
      curcrp = tomatch

      ! determine the numerical hash for curcrp to use as an output variable
      ! and store the bit pattern in the float
      crpval = transfer(crophash(curcrp), 0.0)

! ... Recalculate lignin
      call cmplig(cursys,fligni,wdlig)

      ! Calculate cisofr as 13C if 13C labeling
      if (labtyp .eq. 2) then
        cisofr = del13c * PEEDEE * 1.0e-03 + PEEDEE
        cisofr = cisofr / (1.0 + cisofr) ! =  1 / (1/cisofr + 1)
      endif

      return


      contains
        character(len=1) function itoindx(ii)
          integer :: ii
          itoindx = achar(iachar("0")+ii)
        end function itoindx

        subroutine boundscheck(val, var, mssg)
          real             :: val
          character(len=*) :: var, mssg

          if(val .lt. 0.0  .or.  val .gt. 1.0) then
             write(string,'(a," = ",f10.4)') var, val
             call abortrun('"'//trim(tomatch)//'" '//mssg//' out of bounds: '//trim(string))
          end if
        end subroutine boundscheck
      end
