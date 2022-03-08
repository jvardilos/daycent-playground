
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... TREEIN.F

      subroutine treein(tomatch)

      implicit none

! ... Read in the new forest type

      ! Modifications
      !  28Feb2013  K. Killian
      !  Added parameters for fruit/nut partition
      !  removed tafue in favor of a site level parameter
      !   5 Sept 2013  K. Killian
      !  Added input checking to catch and convert the water stress parameters
      !  used by Bill and Cindy to the independent form.
      ! 1Feb13    K. Killian
      !  added optional lines for water stress coefficients twscoef(2)
      !  added optional lines for nutrient availability     favail (5)
      !  made TAFUE an optional input (1)
      ! 28Nov12    K. Killian
      !  Changed error messages to use the message system instead

      include 'chrvar.inc'
      include 'const.inc'
      include 'dynam.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot3.inc'
      include 'pheno.inc'
      include 'site.inc'

! ... Argument declarations
      character (len = 28) tomatch

      ! Local variables
      integer   ::  ii, jj, kk, ir, jr, kr, istat, lch
      integer   ::  ms, me
      real      ::  del13c, stk(1), chkfrc
      real      ::  value, tmp
      character(len=140) :: string
      character(len=16)  :: varfnd

      ! functions
      real                     :: chkdata
      integer                  :: ReadInt

      ! Number of lines to read for each tree type
      integer, parameter       :: REQLNS = 146, OPTLN = 8, SDEPTH = 0
      integer, parameter       ::   unit = 11
      character (len=6), parameter :: callr = 'treein'

        stk = 0

        me = len_trim(tomatch)
        ms = index(tomatch,'(')

        ! find the tree option
        call oldopen(unit,'tree.100') ! open the tree file  with library search
        call findopt(unit,'tree.100',tomatch(:me),REQLNS) ! find the option in the file

        decid = int(chkdata(unit,callr,'decid',0,stk))
        if(decid.lt.0 .or. decid.gt.2) call abortrun &
            ("out of range deciduous flag, decid, for tree "//trim(tomatch))

        prdx(2) = chkdata(unit,callr,'prdx',0,stk)

        do ii = 1, 4
          ppdf(ii,2) = chkdata(unit,callr,'ppdf',0,stk)
        end do

        ! Read the cerfor array
        ir = 0;        jr = 0;        kr = 0;
        do ii = IMIN, IVAL      ! 1, 3
          part: do jj = 1, FPARTS-1   ! 1, 7-1      !v12  FPARTS-1
            element: do kk = 1, MAXIEL   ! 1, 3
              ! find the value in the line
              cerfor(ii, jj, kk) = 0
              varfnd = ' '
              value = chkdata(unit,callr,varfnd,0,stk)

              ! parse the indices actually in the file
              if (index(varfnd,'cerfor') .ne. 0) then
                lch = index(varfnd,'(')
                if(lch .ne. 0) then
                  lch = lch +1
                  istat = 0
                  ir = ReadInt(0, varfnd, lch, istat)
                  jr = ReadInt(0, varfnd, lch, istat)
                  kr = ReadInt(0, varfnd, lch, istat)
                elseif(ir .eq. 0) then
                  call abortrun('Value out of order '//varfnd//' in '//trim(tomatch))
                endif
                if(jj.eq.FRNUT  .and.  jr.eq.LEAF) then
                  tmp = chkdata(unit,callr,'reread',0,stk)
                  cerfor(ii, FRNUT, 1:MAXIEL) = cerfor(ii, LEAF, 1:MAXIEL)
                  ! load the ccefor array that is used for climate change
                  ccefor(ii, FRNUT, 1:MAXIEL)  = cerfor(ii, LEAF, 1:MAXIEL)
                  exit element ! skip the inner element loop
                elseif(ii.ne.ir .or. jj.ne.jr .or. kk.ne.kr) then
                  call abortrun("cerfor values out of order for tree "//trim(tomatch))
                endif

                cerfor(ii, jj, kk) = value
                varfnd = ' '

!               now load the ccefor array that is used for climate change
               if(ii.le.2) ccefor(ii, jj, kk) = cerfor(ii, jj, kk)

!             Check to see if we need to cover the (3,6,x) option
              elseif (ii.eq.IVAL .and. jj.eq.FPARTS-1) then
                cerfor(IVAL, FRNUT, 1:MAXIEL) = cerfor(1, FRNUT, 1:MAXIEL)
                tmp = chkdata(unit,callr,'reread',0,stk)
              else
                call abortrun('Missing CERFOR value. Found '//varfnd//' instead in tree '//trim(tomatch))
                tmp = chkdata(unit,callr,'reread',0,stk)
              endif
            end do element
          end do part
        end do

        decw1 = chkdata(unit,callr,'decw1',0,stk)
        decw2 = chkdata(unit,callr,'decw2',0,stk)
        decw3 = chkdata(unit,callr,'decw3',0,stk)

        do jj = NEWFOR, OLDFOR
          chkfrc = 0
          do ii = 1, FPARTS-1   !v12    FPARTS-1
            varfnd = ' '
            value = chkdata(unit,callr,varfnd,0,stk)
            if (index(varfnd,'fcfrac') .ne. 0) then
              lch = index(varfnd,'(')
              if(lch .ne. 0) then
                lch = lch +1
                istat = 0
                ir = ReadInt(0, varfnd, lch, istat)
              endif

              ! special conditions on the new fruit/nut allocations
              if(ii .eq. FRNUT  .and.  ir .ne. FRNUT) then
                ! reject the current value from the next column and set default
                tmp = chkdata(unit,callr,'reread',0,stk)
                value = 0
              endif
              fcfrac(ii, jj) = value

 !          don't error on reading another value for the last element
            else if(ii .eq. FRNUT) then
              tmp = chkdata(unit,callr,'reread',0,stk)
              fcfrac(ii, jj) = 0  ! set the default 0
            else
 !            unknown variable in the array
              call abortrun('Missing FCFRAC value. Found '//varfnd//' in '//trim(tomatch))
            endif
          end do
        end do

! ..... Calculate carbon fraction in each part. -mdh 11/13/00
! ..... Assume juvenile (iptr=1) tree for initial carbon fractions
        tree_cfrac = fcfrac(1:FPARTS-1, 1)

! ..... frfrac(2) Added for dynamic allocation. -mdh 11/15/00
! ..... Removed frfrac(2) and added tfrtcn(*) and tfrtcw(*) for calculating
! ..... water and nutrient stress using unique parameter values, cak - 09/12/03
        tfrtcn(1) = chkdata(unit,callr,'tfrtcn',0,stk)
        tfrtcn(2) = chkdata(unit,callr,'tfrtcn',0,stk)
        tfrtcw(1) = chkdata(unit,callr,'tfrtcw',0,stk)
        tfrtcw(2) = chkdata(unit,callr,'tfrtcw',0,stk)

        ! if ((frtcindx .eq. 1) .or. (frtcindx .eq. 3)) then
        if(minval(min(tfrtcn,tfrtcw)) .le. 0.0  .or. & ! minval min element; min compare array elements
           maxval(max(tfrtcn,tfrtcw)) .gt. 1.0) then
          write(string,'("tfrtcn =",2f6.2,"  tfrtcw =",2f6.2) ') tfrtcn, tfrtcw
          call abortrun('tfrtcn or tfrtcw out of bounds in tree '//trim(tomatch)//trim(string)//' must be > 0 and <=1')
        endif
        if ((tfrtcn(2) .gt. tfrtcn(1))) then
          write(string,'("tfrtcn(1),",f4.2," must be >= tfrtcn(2),",f4.2) ') tfrtcn
          call abortrun('tfrtcn out of order in tree '//trim(tomatch)//trim(string))
        endif
        if ((tfrtcw(2) .gt. tfrtcw(1))) then
          write(string,'("tfrtcw(1),",f4.2," must be >= tfrtcw(2),",f4.2) ') tfrtcw
          call abortrun('tfrtcw out of order in tree '//trim(tomatch)//trim(string))
        endif

        ! look for the fruit/nut fall time
        varfnd = ' '
        fnftim = max(0, nint(30*chkdata(unit,callr,varfnd,0,stk)))
        if (index(varfnd,'fnftim') .eq. 0) then
          fnftim = 30 ! default fall time is 30 days, one month
          tmp = chkdata(unit,callr,'reread',0,stk)
        endif


        ! look for the fruit/nut growing degree day limit
        ! Default values, no GDD limit; base Temp =0.0 C; Maximum Temp 100
        fngddl = (/0, 0, 100/)
        do ii = 1, 3
          varfnd = ' '
          value = chkdata(unit,callr,varfnd,0,stk)

!         index decode; test with ir, defaults to 1 if no index
          lch = index(varfnd,'(') +1
          if(lch .gt. 1) then
!           read the actual index from the file
            istat = 0
            ir = ReadInt(0, varfnd, lch, istat)
          endif

! ***COMMENTED TREE STACK CODE: Added  stack structure for reference
!         use values from the stack;
!         if (unit .eq. 0) then
!          accept the stack value
!           fngddl(i) = value
!         else if (index(varfnd,'fngddl') .ne. 0  .and.  ir .eq. i) then

          if (index(varfnd,'fngddl') .ne. 0) then
!           accept the value
            fngddl(ir) = value

          else
!           reread .. this isn't the value wanted
            value = chkdata(unit,callr,'reread',0,stk)
!           stack(m-1,is) = fngddl(i)
!           as long as the stack is not available exit the loop
            exit
          endif
        end do

        do ii = 1, MONTHS
          leafdr(ii) = chkdata(unit,callr,'leafdr',0,stk)
            call boundscheck(leafdr(ii), 'leafdr('//itoindx(ii)//')',"monthly leaf death rate")
        end do

        btolai = chkdata(unit,callr,'btolai',0,stk)
        if(btolai.le.0) call abortrun('treein btolai <= 0 in tree '//trim(tomatch))

        klai = chkdata(unit,callr,'klai',0,stk)
        laitop = chkdata(unit,callr,'laitop',0,stk)
        maxlai = chkdata(unit,callr,'maxlai',0,stk)
        maxldr = chkdata(unit,callr,'maxldr',0,stk)
          call boundscheck(maxldr, 'maxldr',"effect of N availability on leaf death rates")

        do ii = 1, MAXIEL
          forrtf(ii) = chkdata(unit,callr,'forrtf',0,stk)
          call boundscheck(forrtf(ii), 'forrtf('//itoindx(ii)//')',"forest translocation fraction")
        end do

        sapk = chkdata(unit,callr,'sapk',0,stk)
        swold = chkdata(unit,callr,'swold',0,stk)

!        do ii = 1, FPARTS     !v12   FPARTS
!          wdlig(ii) = chkdata(unit,callr,'wdlig',0,stk)
!        end do
        call addfRpool('wdlig', wdlig, wdlig(LEAF))
        do ii=1,FPARTS
          call boundscheck(wdlig(ii), 'wdlig('//itoindx(ii)//')',"forest lignin fraction")
        end do

!        do ii = 1, FPARTS        !v12       FPARTS
!          wooddr(ii) = chkdata(unit,callr,'wooddr',0,stk)
!        end do
        call addfRpool('wooddr', wooddr, 0.)
        do ii=1,FPARTS
          call boundscheck(wooddr(ii), 'wooddr('//itoindx(ii)//')',"tree death rate fraction")
        end do

        wrdsrfc = chkdata(unit,callr,'wrdsrfc|wrdsrf',0,stk)
          call boundscheck(wrdsrfc, 'wrdsrfc',"fraction of tree fine roots transferred surface litter")
        wmrtfrac = chkdata(unit,callr,'wmrtfrac|wmrtfr',0,stk)
          call boundscheck(wmrtfrac, 'wmrtfr',"symbiotic N fixation")
        snfxmx(FORSYS) = chkdata(unit,callr,'snfxmx',0,stk)
          call boundscheck(snfxmx(FORSYS), 'snfxmx(FORSYS)', &
                      "fraction of tree fine roots transferred surface litter")
        del13c = chkdata(unit,callr,'del13c',0,stk)

        co2ipr(FORSYS) = chkdata(unit,callr,'co2ipr',0,stk)
        co2itr(FORSYS) = chkdata(unit,callr,'co2itr',0,stk)

        do ii = IMIN, IMAX
          do jj = 1, MAXIEL
            co2ice(FORSYS,ii,jj) = chkdata(unit,callr,'co2ice',0,stk)
          end do
        end do

        co2irs(FORSYS) = chkdata(unit,callr,'co2irs',0,stk)

        basfc2 = chkdata(unit,callr,'basfc2',0,stk)
        basfct = chkdata(unit,callr,'basfct',0,stk)
! ..... sitpot is now computed as a function long term annual precipitation
! ..... the sitpot parameter value read from the tree.100 file is used as a
! ..... multiplier, see prelim subroutine, cak - 11/21/01
        sitpot_m = chkdata(unit,callr,'sitpot',0,stk)

! ..... Added new parameter, maximum N/P ratio, for phosphorus code, cak - 07/23/02
        maxnp = chkdata(unit,callr,'maxnp',0,stk)
! ..... Add a check for leaf N/P ratio, give a warning message if this value
! ..... falls outside of the maximum N/P, cak - 05/18/2009
        if (nelem .gt. 1) then
         ! restate condition so it requires NO division rather than 3
         ! second step assumes cerfor(1,1,1) is positive (which SHOULD be true)
         ! (((1.0/cerfor(1,1,1))/(1.0/cerfor(2,1,2))) .gt. maxnp) =
         ! (cerfor(2,1,2)/cerfor(1,1,1) .gt. maxnp)               =
          if (cerfor(2,1,2) .gt. maxnp*cerfor(1,1,1)) then
            call message('WARNING: N/P for leaves falls outside of '//   &
                         'maxnp for tree option: '//trim(tomatch))
          endif
        endif

! ..... Added fkmrspmx parameters for maintenance respiration code, mdh - 11/30/01
!        do ii = 1, FPARTS     !v12      FPARTS
!          fkmrspmx(ii) = chkdata(unit,callr,'fkmrsp',0,stk)
!        end do
        call addfRpool('fkmrspmx', fkmrspmx, fkmrspmx(LEAF))

! ..... Added parameters for controlling the linear decrease in
! ....  maintenance respiration as the amount of carbohydrate stored in
! ..... the carbohydrate storage pool gets smaller, cak - 08/13/2009
!        fmrsplai(1) = chkdata(unit,callr,'fmrspl',0,stk)
!        fmrsplai(2) = chkdata(unit,callr,'fmrspl',0,stk)
!        fmrsplai(3) = chkdata(unit,callr,'fmrspl',0,stk)
!        fmrsplai(4) = chkdata(unit,callr,'fmrspl',0,stk)
!        fmrsplai(5) = chkdata(unit,callr,'fmrspl',0,stk)
!        fmrsplai(6) = chkdata(unit,callr,'fmrspl',0,stk)    !v12 ends at 6
        call addfRpool('fmrsplai|fmrspl', fmrsplai, fmrsplai(LEAF))
        if (fmrsplai(FROOTJ) .lt. 0.0) then
          write(string,*) fmrsplai(FROOTJ)
          call abortrun('fmrsplai(2) = '//trim(string)//' < 0.0 in tree'//trim(tomatch))
        endif
        if (fmrsplai(LWOOD) .lt. fmrsplai(FROOTJ)) then
          write(string,*) 'large wood fmrsplai(4)',fmrsplai(LWOOD), &
          'must be >= juvenile fine root fmrsplai(2)', fmrsplai(FROOTJ)
          call abortrun(trim(string)//' in tree '//trim(tomatch))
        endif
        if (fmrsplai(FROOTM) .lt. fmrsplai(LWOOD)) then
          write(string,*) 'mature fine root fmrsplai(7)',fmrsplai(FROOTM), &
          'must be >= large wood fmrsplai(4)', fmrsplai(LWOOD)
          call abortrun(trim(string)//' in tree '//trim(tomatch))
        endif

! ..... Added fgresp parameters for the growth respiration code, cak - 01/16/2007
!        do ii = 1, FPARTS     !v12 FPARTS
!          fgresp(ii) = chkdata(unit,callr,'fgresp',0,stk)
!        end do
        call addfRpool('fgresp', fgresp, fgresp(LEAF))

! ..... Added no3pref, mdh - 9/11/01
        no3pref(FORSYS) = chkdata(unit,callr,'no3pref|no3pre',0,stk)
          call boundscheck(no3pref(FORSYS), 'no3pref(FORSYS)', "fraction of N update that is NO3")

! ..... Added tlaypg, cak - 01/29/03
        tlaypg = int(int(chkdata(unit,callr,'tlaypg',0,stk)))
        if (tlaypg .gt. nlayer) then
          write(string,*) ' tlaypg',tlaypg,'must be <= nlayer', nlayer
          call message('Warning: for '//trim(tomatch)//trim(string)//   &
              ' Resetting tlaypg to nlayer.')
          tlaypg = nlayer
        endif

! ..... Added tmix, cak - 06/14/05
        tmix = chkdata(unit,callr,'tmix',0,stk)

! ..... Added tmplff and tmplfs, cak - 02/21/03
        tmplff = chkdata(unit,callr,'tmplff',0,stk)
        tmplfs = chkdata(unit,callr,'tmplfs',0,stk)

! ..... Add parameters used to restrict production late in the growing
! ..... season, cak - 03/11/2010
        furgdys = int(int(chkdata(unit,callr,'furgdys|furgdy',0,stk)))
        flsgres = chkdata(unit,callr,'flsgres|flsgre',0,stk)

! ..... Added tmxturn, cak - 06/28/2007
        tmxturn = chkdata(unit,callr,'tmxturn|tmxtur',0,stk)

        ! Added water stress equation coefficents, cak - 07/10/2012
        varfnd = ''
        twscoef(1) = chkdata(unit,callr,varfnd,0,stk)
        if(unit .eq. 0  .or.  index(varfnd,'wscoef') .gt. 0) then
          twscoef(2) = chkdata(unit,callr,'wscoef',0,stk)
          if(twscoef(1) .gt. 1.0  .and.  twscoef(2) .lt. 0) then
            write(string,'(2f5.1)') twscoef
            twscoef(2) = -twscoef(2)
            twscoef(1) = log(twscoef(1))/twscoef(2)
            !stack(m-2:m-1, is) = twscoef   ! put revised values on stack
            write(varfnd,'(2f6.3)') twscoef
            call message("converting deprecated "//trim(tomatch)// &
                 " WSCOEF values from"//trim(string)//" to "//varfnd)
          endif
        else  ! reading file;  wscoef wasn't found. (default file condition)
          twscoef(1) = chkdata(unit,callr,'reread',0,stk) ! reuse value
          ! set default conditions
          twscoef = (/ 0.378, 9. /) ! (/ 0.5, 9.0 /) 2011 inventory
          !stack(m-1:m, is) = twscoef   ! put default values on the stack
          !m = m+1
        endif

! ..... Added npp2cs, cak - 03/16/2011
        npp2cs(FORSYS) = chkdata(unit,callr,'npp2cs',0,stk)

        ! favail let the crop modify the favail fix parameter for DelGrosso Oct 2012
        ! the indexing on this is beastly. favail(2) is **NOT** input.
        ! thus look for favail(1)
        !  If found then set 2 to zero and input 3-6 from file
        ! NOT FOUND. Then use the fix.100 values.
        ! Cindy added the N availability as a standalone so break this into 2 tests 30Jan13
        varfnd = ''
        favail(1,FORSYS) = chkdata(unit,callr,varfnd,0,stk)
        if(unit .eq. 0  .or.  index(varfnd,'favail') .gt. 0) then
          ! now look for the rest of the favail array
          varfnd = ''
          favail(3,FORSYS) = chkdata(unit,callr,varfnd,0,stk)
          if(unit .eq. 0  .or.  index(varfnd,'favail') .gt. 0) then
            ! entire array available
            favail(4,FORSYS) = chkdata(unit,callr,'favail',0,stk)
            favail(5,FORSYS) = chkdata(unit,callr,'favail',0,stk)
            favail(6,FORSYS) = chkdata(unit,callr,'favail',0,stk)
          else  ! read favail(1) but not the rest of the element data
            ! favail(3,FORSYS)   = chkdata(unit,callr,'reread',0,stk) ! reuse value
            favail(3:6,FORSYS) = fixfavail(3:6) ! use values from fix.100
            ! stack(m-1:m+2, is) = fixfavail(3:6) ! put favail(3:6) on the stack
            ! m = m+3
          endif
        else  ! reading file but variable not found.
          ! favail(1,FORSYS) = chkdata(unit,callr,'reread',0,stk) ! reuse value
          favail(:,FORSYS) = fixfavail(:) ! use values from fix.100
          ! stack(m-1, is)   = fixfavail(1)   ! put favail(1) on the stack
          ! stack(m:m+3, is) = fixfavail(3:6) ! put favail(3:6) on the stack
          ! m = m+4
        endif
        favail(2,FORSYS) = 0 ! this will be initialized later

! ..... Added tafue, cak - 02/21/2011 ! removed tafue 28Feb 14
!        ! Depricated Added tafue, cak - 02/21/2011
!        varfnd = ''
!        value = chkdata(unit,callr,varfnd,m,stack(1,is))
!        if(unit .ne. 0  .and.  index(varfnd,'tafue') .eq. 0) then
!          value = chkdata(unit,callr,'reread',m,stack(1,is)) ! reuse value
!          call message("Warning tafue input in "//trim(tomatch)//" is depricated")
!        else
!          ! default: reading file and variable not found.
!        endif

! ..... Close the file
        close(unit)

! ..... Hold on to the current tree just read in
        curtre = tomatch

! ... Calculate cisotf as 13C if 13C labeling
      if (labtyp .eq. 2) then
        cisotf = del13c * PEEDEE * 1.0e-03 + PEEDEE
        cisotf = cisotf / (1.0 + cisotf)  ! 1.0 / (1.0/cisotf + 1.0)
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
             call abortrun(mssg//'input out of bounds: '//trim(string)//' in tree '//trim(tomatch))
          end if
        end subroutine boundscheck


        subroutine addfRpool(var, farray, defalt)
          character (len=*)   :: var
          integer             :: ii
          real                :: farray(FPARTS)
          real                :: defalt, tmp

          do ii = 1, FPARTS-1
            farray(ii) = chkdata(unit,callr,var,0,stk)
          end do

          varfnd = ' '
          farray(FPARTS) = chkdata(unit,callr,varfnd,0,stk)
          !if(varfnd(1:1) .eq. '"' .or. varfnd(1:1) .eq. "'") then ! shouldn't Happen
          !  ii = scan(varfnd, varfnd(1:1), .true.) ! scan back to the ending quote
          !  if(ii.eq.0) call abortrun('bad tree input variable'//trim(varfnd))
          !  varfnd = varfnd(2:ii-1)
          !endif
          ii = SCAN(varfnd, '()') -1
          if(ii < 0) ii = len_trim(varfnd)
          if (index(var, varfnd(:ii)) .eq. 0) then
            farray(FROOTM) = farray(FPARTS-1)
            farray(FRNUT) = defalt ! default to the leaf ratio
            tmp = chkdata(unit,callr,'reread',0,stk)
          endif
        end subroutine addfRpool

      end subroutine treein



      subroutine trplntin(tomatch)

      implicit none

      character(len=28) tomatch

!...Read in and execute tree plant/initialization data.
!    The routine does 3 tasks.
!      1 It reads the parameters from the event description.
!      2 Clears any live tree, removes the aboveground tree and kills the root
!      3 Plants the requested tree biomass.
!
!    NOTE: This routine is called from schedule. Probably shouldn't be but was
!    easy and it causes the plant to occur before the time step's growth
!
!    There are two modes for inputing the planted tree. An immediate mode is
!    available that specifies the above ground size, uses the initial C/E ratios
!    for the tree and assumes root fractions. To provide the user more extensive
!    control of the tree plant values a new 'trplant.100' file was created.
!    In BOTH input modes, the fruit values are assumed to be zero.
!
!
!    Modifications:
!     6/2008   KLK    Initial coding
!      the immediate form
!        Immediate mode inputs the:
!        above ground Wood C     (W)
!        fine branch fraction    (B)
!        leaf fraction           (L)
!        tree age                (Y)
!        Large root is assumed to be half of the large wood,
!        fruit and fine roots, labeled C fractions are zero.
!        C/E ratios are taken from cerfor(3,x,x) values for the correct tree
!
!      definition of the external .100 file
!        tplttc        total C planted in tree plant event
!        flvcis(1)     C fraction planted in unlabeled leaf component
!        flvcis(1)     C fraction planted in labeled leaf
!        celeaf(1)     C/N ratio in planted leaf
!        celeaf(2)     C/P ratio in planted leaf
!        celeaf(3)     C/S ratio in planted leaf
!        ffbcis(1)     C fraction planted in unlabeled fine branch
!        ffbcis(2)     C fraction planted in  labeled fine branch
!        cefbr(1)      C/N ratio in planted fine branch
!        cefbr(2)      C/P ratio in planted fine branch
!        cefbr(3)      C/S ratio in planted fine branch
!        flwcis(1)     C fraction planted in unlabeled large wood
!        flwcis(2)     C fraction planted in   labeled large wood
!        celwod(1)     C/N ratio in planted large wood
!        celwod(2)     C/P ratio in planted large wood
!        celwod(3)     C/S ratio in planted large wood
!        ffrtcis(1)    C fraction planted in unlabeled fine root
!        ffrtcis(2)    C fraction planted in   labeled fine root
!        cefroot(1)    C/N ratio in planted fine root
!        cefroot(2)    C/P ratio in planted fine root
!        cefroot(3)    C/S ratio in planted fine root
!        fcrtcis(1)    C fraction planted in unlabeled coarse root
!        fcrtcis(2)    C fraction planted in   labeled coarse root
!        cecroot(1)    C/N ratio in planted coarse root
!        cecroot(2)    C/P ratio in planted coarse root
!        cecroot(3)    C/S ratio in planted coarse root
!        frstage        the age of the planted trees

      include 'chrvar.inc'
      include 'const.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'timvar.inc'
      include 'zztim.inc'

!     local variables
      integer :: sordr, m, me, ms, nv, istat, unit
      real    :: value, wt, mfprd(2,FPARTS), tplttc, stk, cisf(2)
      real, dimension(FPARTS-1,MAXIEL) :: CEf
!      character varfnd*16
      character (len=13),parameter :: keylist = "CSBY?csby";
      character (len=4), parameter :: evnt = "TPLT";

!     Declare functions
      real      chkdata
      double precision :: ReadReal, parsImedFld

!     Declare and define parameters
!...Number of lines to read for each tree type
      integer   TPLNTLNS
      parameter (TPLNTLNS = 30)

      unit = 11
      mfprd = 0.
      !cisf = [0.5-cisotf+0.5, cisotf]
      cisf = [1, 0]
      CEf  = cerfor(3,:,:)

        ms = 1
        me = len_trim(tomatch)

! write(*,*) " tomatch:",trim(tomatch)," at ",cyear, month
        if(tomatch(:1) .eq. '(') then
          if(tomatch(me:me) .ne. ')') call  &  ! OOOOPS  the end of the data was lost
                     abortrun('incomplete tree plant: '//tomatch)

          ! Make sure a tree is defined before this mode is used
          if(curtre .eq. " " .or. cerfor(3,1,1) .eq. 0) call abortrun( &
             'tree type must be defined before planting '//trim(tomatch))

          ! read in the total C, and fine branch and leaf fraction
          ! the default course root is the 1/2 the large wood.
          m  = 2
          me = len_trim(tomatch(:me-1))
          if(me<=m) call abortrun("empty immediate tree Plant: "//trim(tomatch))


          sordr=0 ! set the sort order to 0
          do while (m.le.me) ! the last character is the terminator

            !treat tab as a soft delimiter
            if(tomatch(m:m) .eq. '	') then; m=m+1; cycle; endif

            ms = m ! save the field start for error reporting
            value = parsImedFld(nv, tomatch, sordr,istat,m, ms,me,keylist,evnt)
            write(*,*) "Immediate:",evnt,value,'nv=',nv,istat,m,"  '",tomatch(m:me),"'"

            select case (nv)
              case (1)
                ! store input, ordered above ground C, fine branch frac, leaf frac
                tplttc = value
              case (2)
                mfprd(:,LWOOD)  = value *cisf                ! large wood
              case (3)
                mfprd(:,FBRCH)  = value *cisf                ! fine branch fraction
              case (4) ! Y  tree age in years
                frstage = nint(value)
              case DEFAULT
                call abortrun('Unknown '//evnt//' key '//tomatch(ms:m)//' in event '//trim(tomatch))
            end select
          end do

          mfprd(:,LEAF)  = (1. - sum(mfprd(:,LWOOD) + mfprd(:,FBRCH))) *cisf  ! leaf fraction
          if(sum(mfprd(:,LWOOD)) < 0.0) call abortrun('tree plant - leaf fraction < 0')
          mfprd(:,FROOT)  = ((1-sum(mfprd(:,LWOOD))) * 0.5) * cisf            ! fine root
          mfprd(:,CROOT)  = (mfprd(:,LWOOD) * 0.5)                            ! coarse roots
          mfprd(:,FRNUT)  = 0                                                 ! fruit
        else
!         Read all the parameters from the tree file
          call oldopen(unit,'tplt.100')       !  open the input file

          wt = 1
          if(index(tomatch,'*') .gt. 1) then
            ms = 1;  istat =0
            wt= ReadReal(0,tomatch,ms,istat)
            if(istat .le. 0) call abortrun('bad immediate tree plant weight: '//trim(tomatch))

            ms= ms+1
!            write(*,*) 'wt=',wt, ' ms=',ms,'    ',tomatch(m:),' ',istat
          endif

!         find the tree option
          call findopt(unit,'tplt.100',tomatch(ms:me), TPLNTLNS)

!         total C in plant
          tplttc        = wt * chkdata(unit,'tplt.100','tplttc',0,stk)
          if(tplttc .lt. 0.) call lvtrem()
          tplttc = abs(tplttc)

          frstage = chkdata(unit,'tplt.100','frstage',0,stk) ! Age of planted tree

          call readFpart('flvcis','celeaf', LEAF)  ! leaf
          call readFpart('ffbcis','cefbr', FBRCH)   ! fine branch
          call readFpart('flwcis','celwod', LWOOD)  ! large wood
          call readFpart('ffrtcis','cefroot', FROOT)
          call readFpart('fcrtcis','cecroot', CROOT) ! coarse root

          ! value   = sum(rlvcis) + sum(fbrcis) + sum(rlwcis)
          ! frstage = frstage * (value /(value + sum(crtcis)+sum(frtcism)))
        endif

        if(tplttc .lt. 0.) call lvtrem()
        tplttc = abs(tplttc)

        call addFpart(rlvcis, rleave, mfprd(:,LEAF), CEf(LEAF,:))  ! leaf
        call addFpart(fbrcis, fbrche, mfprd(:,FBRCH), CEf(FBRCH,:))   ! fine branch
        call addFpart(rlwcis, rlwode, mfprd(:,LWOOD), CEf(LWOOD,:))  ! large wood
        call addFpart(crtcis, croote, mfprd(:,CROOT), CEf(CROOT,:)) ! coarse root
        call addFpart(frtcisj, frootej, mfprd(:,FROOT) * (0.5 - wmrtfrac + 0.5), CEf(FROOT,:))
        call addFpart(frtcism, frootem, mfprd(:,FROOT) * wmrtfrac,               CEf(FROOT,:))

!        Update state variables and accumulators
         ! call flowup(time)
         call sumcar
        return

      contains

        subroutine readFpart(clbl,elbl, fp)
          integer, intent(IN)       :: fp
          character*(*), intent(IN) :: clbl,elbl
          real tmp
          character(len=16)  :: varfnd
!logical, optional, intent (in) :: header

          mfprd(1,fp) = chkdata(unit,'tplt.100',clbl,0,stk)
          varfnd = ' ';  tmp = chkdata(unit,'tplt.100',varfnd,0,stk)    ! labeled fraction
          if (index(varfnd,clbl) .ne. 0) then
            ! found the labeled C
            mfprd(2,fp) = tmp
            !
            varfnd = ' '
            tmp = chkdata(unit,'tplt.100',varfnd,0,stk)
          else
            ! mfprd(:,fp) = mfprd(1,fp) * cisf ! convert component total to label/unlabel
            mfprd(2,fp) = 0 ! set labeled fraction to 0
          endif

          ! look for the C/E ratios
          if (index(varfnd,elbl) .ne. 0) then
            CEf(fp,1) = tmp
            CEf(fp,2) = chkdata(unit,'tplt.100',elbl,0,stk)
            CEf(fp,3) = chkdata(unit,'tplt.100',elbl,0,stk)
          else
            tmp = chkdata(unit,'tplt.100','reread',0,stk)
          endif

       end subroutine readFpart


        subroutine addFpart(fpc, fpe, mfprd, CEfp)
          real, intent(INOUT)       :: fpc(2), fpe(MAXIEL)
          real, intent(IN)          :: mfprd(2), CEfp(MAXIEL)

          real          :: value(MAXIEL), tmpc

          value(1:2)  = mfprd * tplttc ! convert fraction to C
          tmpc        = sum(value(1:2))
          fpc         = fpc    + value(1:2)
          csrsnk      = csrsnk - value(1:2)

          value = 0
          where (CEfp(1:nelem) > 0) value(1:nelem)  = tmpc/CEfp(1:nelem)
          fpe(1:nelem)    = fpe(1:nelem) + value(1:nelem)
          esrsnk(1:nelem) = esrsnk(1:nelem) - value(1:nelem)
        end subroutine addFpart

      end subroutine trplntin


      subroutine lvtrem()

        include 'const.inc'
        include 'plot3.inc'

!       lvtrem Live tree removal
!              this encapsulates the code needed to remove live trees.
!              It does not effect dead wood

!       Local variables
        character(len=8) :: mssg
        real             :: accum(ISOS), value
!       set up a live wood removal fraction KLUDGE: it works for root death
        real, dimension(FPARTS) :: remfl = (/1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 /)

        value = rleavc + frootcj + fbrchc + rlwodc + crootc + frootcm + frnutc
        if(value .gt. 0.) then
          value = rleavc + fbrchc + rlwodc + frnutc
          write(mssg,'(f8.2,a)') value
          call message("Warning: TPLT is removing "//mssg// &
                       " gC/m2 aboveground live tree")
          value = frootcj + frootcj + crootc
          write(mssg,'(f8.2,a)') value
          call message("              and killing "//mssg// &
                       " gC/m2 live tree root")

!         Live Removal
          call livrem(accum, remfl)
          call killrt(accum, remfl)

          frstage =0.
        endif

      return
      end subroutine lvtrem

! LEAF   = 1   LEAF    forest leaf
! FROOTJ = 2   FROOTJ  juvenile (or total) forest fine root
! FBRCH  = 3   FBRCH   forest fine branch
! LWOOD  = 4   LWOOD   forest large wood
! CROOT  = 5   CROOT   forest coarse root
! FRNUT  = 6   FRNUT   forest tree fruit/nut
! FROOTM = 7   FROOTM  forest mature fine root

