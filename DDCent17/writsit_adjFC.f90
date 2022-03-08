
!              Copyright 1993 Colorado State University
!                       All Rights Reserved

!....   WRTLIN        ! Writes the individual site.100 lines to an open ASCII site file
MODULE wrtsitfil
        USE ISO_C_Binding

        ! 7 in word, plus paranthesis, comma and 2 indexes of 2 characters
        integer, parameter :: namlen=7, namstart=16;
!        character (len=44) :: sitlin!, sitword*44
        character (len=max(44,namlen+7+namstart)) :: sitlin!, sitword*44

       INTERFACE wrtlin
          MODULE PROCEDURE wrtRlin, wrtIntlin, wrtDPlin
       END INTERFACE
CONTAINS
        subroutine wrtRlin(osf, k, value)

!...Writes the individual site.100 lines to an open ASCII site file
!
         integer,       intent(IN)    :: osf
         integer,       intent(INOUT) :: k
         real,          intent(IN)    :: value

!
!**** create the key word for the appropriate index
10          sitlin = sitword(k)
            if(sitlin(1:1) .ne. '*')  then
              sitlin(namstart:) = sitlin;
              sitlin(:namstart-1) = ' ';
              call formout(sitlin,value)
            endif

            write(osf,'(a)') trim(sitlin)
            k = k+1
            if(sitlin(1:1) .eq. '*')  goto 10
          return
          end

        subroutine wrtIntlin(osf, k, ival) ! convert integer to fortran real (same thing)
          integer, intent(IN)    :: osf
          integer, intent(INOUT) :: k
          integer, intent(IN)    :: ival

          call wrtRlin(osf, k, real(ival))

        return
        end

        subroutine wrtDPlin(osf, k, DP) ! convert c_double to fortran real (same thing)
          integer,        intent(IN)    :: osf
          integer,        intent(INOUT) :: k
          real(c_double), intent(IN)    :: DP

          call wrtRlin(osf, k, real(DP))

        return
        end


        character*44 function sitword(lin)
        integer lin

        ! character addindx*16
        integer, save ::   l=0, k=0, i=1
        integer, parameter :: maxvar = 119

!       Modifications
!  4/2017  KLK
!    added FRSTAGE to forest section
!          KLK
!    removed recalculated soil profile parameters
!  7/2011  K. Killian
!    completely recoded to use an expanded variable list and a loop
!     instead of the brute force if-then statements
!     added Daily DayCent variables
!  5/2006  K. Killian
!    Removed AMINRL from the enhanced site parameters

        type sitvar
          character (len=namlen) :: varlbl
          integer                :: du
          integer                :: dl
        end type sitvar

        ! NOTE: sitvar MUST contain abbreviations of the section labels
        !       or the count gets out of kilter
        type (sitvar), dimension(maxvar) :: sitwrd = (/ &
           sitvar('x site',1,1),   sitvar('*Climat',1,1),  & !labels 1, 2
           sitvar('PRECIP(',1,12), sitvar('PRCSTD(',1,12), sitvar('PRCSKW(',1,12), &
           sitvar('TMN2M(', 1,12), sitvar('TMX2M(', 1,12), &
           sitvar('SRADJ(', 1,12), sitvar('RAINHR', 1,1),&

           sitvar('*contro', 1,1), & !label 3
           sitvar('IVAUTO', 1,1),  sitvar('NELEM',  1,1),  sitvar('SITLAT', 1,1), &
           sitvar('SITLNG', 1,1),  sitvar('ELEV',   1,1),  sitvar('SITSLP', 1,1), &
           sitvar('ASPECT', 1,1),  sitvar('EHORIZ', 1,1),  sitvar('WHORIZ', 1,1), &

!          sitvar('SAND', 1,1),    sitvar('SILT', 1,1),    sitvar('CLAY', 1,1),   &
!          sitvar('ROCK', 1,1),    sitvar('BULKD', 1,1),   sitvar('NLAYER', 1,1)

           sitvar('NLAYPG', 1,1),  sitvar('DRAIN', 1,1),   sitvar('BASEF', 1,1),  &
           sitvar('STORMF', 1,1),  sitvar('PRECRO', 1,1),  sitvar('FRACRO', 1,1), &
           sitvar('SWFLAG', 1,1),  &
         ! sitvar('AWILT(', 1,10), sitvar('AFIEL(', 1,10), sitvar('PH', 1,1), &

           sitvar('PSLSRB', 1,1),  sitvar('SORPMX', 1,1),  &
           sitvar('SUBLIM', 1,1),  sitvar('REFLEC', 1,1),  sitvar('ALBEDO', 1,1), &

           sitvar('DMPFLX', 1,1),  sitvar('DRNLAG', 1,1),  sitvar('HPOTD', 1,1),   sitvar('KSATD', 1,1),   &

           sitvar('DMPFCT',1,1), &
           sitvar('NCOEFF',1,1),   sitvar('DNSTRT',1,1),   sitvar('DNEND',1,1), &
           sitvar('NADJFC',1,1),   sitvar('NADJWP',1,1),   sitvar('MAXNIT',1,1), &
           sitvar('MINO3', 1,1),   sitvar('WFPSNIP',1,1),  sitvar('N2N2OA',1,1), &
!           sitvar('FN2N2O',1,1),  sitvar('CO2CH4',1,1),   sitvar('ACFERM',1,1), &
!           sitvar('FREXUD',1,1),   sitvar('AEH',   1,1),   sitvar('DEH', 1,1), &
!           sitvar('BEHFL', 1,1),   sitvar('BEHDR', 1,1),   sitvar('METHZR',1,1), &

           sitvar('*nutrie', 1,1), & !label 4
           sitvar('EPNFA(',1,2),   sitvar('EPNFS(', 1,2),  sitvar('SATMOS(',1,2), &
           sitvar('SIRRI',1,1),    sitvar('AFUE', 1,1),    &

           sitvar('*Organi', 1,1), & !label
           sitvar('SOM1CI(',2,2),  sitvar('SOM2CI(',2,2),  sitvar('SOM3CI(',1,2), &
           sitvar('RCES1(', 2,3),  sitvar('RCES2(', 2,3),  sitvar('RCES3(', 1,3), &
           sitvar('CLITTR(',2,2),  sitvar('RCELIT(',2,3),  sitvar('AGLCIS(',1,2), &
           sitvar('AGLIVE(',1,3),  sitvar('BGLCIS(',1,2),  sitvar('BGLIVE(',1,3), &
           sitvar('STDCIS(',1,2),  sitvar('STDEDE(',1,3),  &

           sitvar('*Forest', 1,1), & !label
           sitvar('RLVCIS(', 1,2), sitvar('RLEAVE(',1,3),  sitvar('FBRCIS(',1,2), &
           sitvar('FBRCHE(', 1,3), sitvar('RLWCIS(',1,2),  sitvar('RLWODE(',1,3), &
           sitvar('FRTCIS(', 1,2), sitvar('FROOTE(',1,3),  &
           sitvar('CRTCIS(',1,2),  sitvar('CROOTE(',1,3),  &
           sitvar('WD1CIS(',1,2),  sitvar('WD2CIS(',1,2),  sitvar('WD3CIS(',1,2), &
           sitvar('FRSTAGE',1,1), &

           sitvar('*Minera', 1,1), & !label
           sitvar('MINERL(',10,3), sitvar('PARENT(',1,3),  sitvar('SECNDY(',1,3), &
           sitvar('OCCLUD',1,1), &

           sitvar('*Water', 1,1), & !label
           sitvar('RWCF(',1,10),   sitvar('SNLQ', 1,1),    sitvar('SNOW',1,1), &
           sitvar('SNWINS', 1,1), &

           sitvar('*Soilin',1,1), & !label
           sitvar('SLDPMX(',1,21), sitvar('SLBLKD(',1,21), sitvar('SLFLDC(',1,21), &
           sitvar('SLWLTP(',1,21), sitvar('SLECOF(',1,21), sitvar('SLTCOF(',1,21), &
           sitvar('SLSAND(',1,21), sitvar('SLCLAY(',1,21), sitvar('SLORGF(',1,21), &
           sitvar('SLCLIM(',1,21), sitvar('SLSATC(',1,21), sitvar('SLPH(',  1,21),   &

           sitvar('*Extend',1,1), & !label
           sitvar('RMETCS(',2,2),  sitvar('RCEMET(',2,3),  sitvar('RSTRLG(',1,2), &
           sitvar('RCEWD1(',1,3),  sitvar('RCEWD2(',1,3),  sitvar('RCEWD3(',1,3), &
           sitvar('CSTG(',  2,2),  sitvar('CRPSTG(',1,3),  sitvar('FORSTG(',1,3), &
           sitvar('MRTFR',  1,1),  sitvar('WMRTFR', 1,1),  &
           sitvar('AMMONI', 1,1),  sitvar('NITRAT(',1,21), sitvar('SWCINI(',1,21) &
           /)

          ! reset line and variable counters if someone backs up
          if (lin < l) then
            k = 0
            i = 1
          endif

          ! search until the variable can handle the current line
          do while (lin - (k+sitwrd(i)%du*sitwrd(i)%dl) > 0 .and. i<=maxvar)
            k = k + sitwrd(i)%du*sitwrd(i)%dl
            i = i +1
          end do

        ! return 'stop' if you ask for to many lines
        if(i>maxvar) then
             sitword = 'stop'
             return
        endif

        l = lin

        ! construct the variable label
        sitword = sitwrd(i)%varlbl

        ! Expand the long for section labels
        SELECT CASE (sitword)
          CASE ('x site')
             sitword = 'X    Site file data from '
          CASE ('*Climat')
             sitword = '*** Climate parameters'
          CASE ('*contro')
             sitword = '*** Site and control parameters'
          CASE ('*nutrie')
             sitword = '*** External nutrient input parameters'
          CASE ('*Organi')
             sitword = '*** Organic matter initial values'
          CASE ('*Forest')
             sitword = '*** Forest organic matter initial parameters'
          CASE ('*Minera')
             sitword = '*** Mineral initial parameters'
          CASE ('*Water')
             sitword = '*** Water initial parameters'
          CASE ('*Soilin')
             sitword = '*** Soil layer parameters'
          CASE ('*Extend')
             sitword = '*** Enhanced extend parameters'
        END SELECT

        ! leave now if we found a section label
        if(sitword .ne. sitwrd(i)%varlbl) return


          if(sitwrd(i)%du*sitwrd(i)%dl > 1) then
            if(sitword(:6) .eq. 'MINERL') then
              sitword= addindx(sitword,1+mod(lin-k-1,sitwrd(i)%du))
              sitword= addindx(sitword,1+(lin-1-k)/sitwrd(i)%du)
            else
              if (sitwrd(i)%du > 1) sitword= addindx(sitword,1+(lin-k-1)/sitwrd(i)%dl)
              sitword= addindx(sitword,1+mod(lin-k-1,sitwrd(i)%dl))
            endif
          end if

        return
        end


        character*16 function addindx(sitlin,indx)
        character sitlin*(*)
        integer l,l2,indx

! find the word end, change existing ')' to ',' then add the closing paren
        addindx = sitlin
        l = index(addindx,' ')
        if(addindx(l-1:l-1).eq.')') addindx(l-1:l-1)=','
        addindx(l:l) = ')'

! insert the index digits; moving everything else to the right
        l2 = indx
10      if (l2 .eq. 0) return
          addindx(l:) = char(ichar('0')+mod(l2,10))//addindx(l:)
          l2 =l2/10
        go to 10
        end


!....   FORMOUT
        subroutine formout(sitlin,value)
        character*(*) sitlin
        real value
        integer j, je, jd
        character frmat*7

!       Modifications
!  1/2008  K. Killian
!    increased the output precision to 6 digits.
!     Recoded the zero removal loops

! Handle the trivial case
        if (value .eq. 0.) then
          sitlin(1:4) = ' 0.0'
          return

        else if (abs(value).ge.1.0e5 .or. abs(value) .lt. 0.001) then
! if really out of range write in scientific 1<= N <=10 format
          write(sitlin(:12),'(1pg12.5)') value

        else
! write out a 6 digit number
          frmat='(g12.6)'
          if(abs(value) .lt. 0.01) then
            frmat='(f11.8)'
          else if(abs(value) .lt. 0.1) then
            frmat='(f10.7)'
          else if(abs(value) .lt. 1.0) then
            frmat='(f9.6)'
          endif
          write(sitlin(:12),frmat) value
        endif


!         write(*,*) 'trimming  "',sitlin(:12),'"'
! Remove trailing zeros
!       remove trailing zeros in scientific notation
40      je = index(sitlin(:12),'0E') + index(sitlin(:12),'0D')
        if(je .gt. 0) then
          if(sitlin(je-1:je-1) .ne. '.') then
            sitlin(je:12) = sitlin(je+1:12)//'  '
            goto 40
          endif
        endif

!       remove trailing zeros in floating formats (no exponent)
        je = index(sitlin(:12),'E') + index(sitlin(:12),'D')
        jd = index(sitlin(:12),'.')
        if(je .eq. 0  .and.  jd .gt. 0) then
50        j  = index(sitlin(:12),'0 ')
!         trim a zero-space if before the exponent and after the decimal
!          if(j .gt. 0  .and.  j .gt. jd+1) then
          if(j .gt. 0) then
            sitlin(j:12) = ' '
            goto 50
          endif
        else

!         leading zero's exponent from the exponent
          j= index(sitlin(je+1:12),'+0') + index(sitlin(je+1:12),'-0')
          if(j .gt. 0) sitlin(je+j+1:12) = sitlin(je+j+2:12)//'  '
        endif

        return
        end
end MODULE wrtsitfil


!.... WRTSITE
      subroutine wrtsite(mode,i,sitlin)

      integer mode, i
      character sitlin*(*)

      integer, save            :: windx  = 0
      integer, save            :: oslen = 0, lbllen = 0
      integer                  :: j, k
      character (len=256),save :: outsit = " "
      character (len=80), save :: sitlbl = " "

      !include "jday.inc"

!   wrtsite
!      This routine writes the site file information to a disk file
!      input is loaded through repeated calls with mode set >0
!         mode =  0    write the record
!         mode =  1    save   the output file name and the write index
!         mode =  2    save   site record label
!         mode = 11    return the output file name and the write index
!         mode = 12    return site record label
!
!      The site files names MUST include the extension names so we can figure
!      out how to process the data.   NO extension then return without error
!
!       Modifications
!  Nov/2009  K. Killian
!    use the site nlaypgst for site input output rather than the dynamic nlaypg
!   Jan2007  K. Killian
!       corrected a bug in the rcelit calculation
!   9/2006  K. Killian
!       Added a fatal trap for r/w record =0
!       standrard open specifies convert='BIG_ENDIAN'
!   5/2006  K. Killian
!       removed aminrl from enhanced output list. It is recalculated every
!       month in simsom
!   4/2003  K. Killian
!       minor changes to ensure compatability with g77
!   6/2001  K. Killian
!       Added an output message that the ASCII file is being written

! NOTE:
!      The gfortran, Absoft, and ProFortran compilers all work with
!      convert='BIG_ENDIAN' in the binary open even if it is the native format.
!      The SUN compler f77 currently needs it removed.
!      The alternative is to read in native format and uncomment sitswend.

        ! write(*,*) 'wrtsite',mode,i," '",trim(sitlin),"'";

!       set the site file name and write index
        if(mode .eq. 1) then
           oslen = len_trim(sitlin)
           outsit = sitlin !  sitlin(:oslen)//' '
           windx = i

!       set the site label
        else if(mode .eq. 2) then
           sitlbl = sitlin
           lbllen = len(sitlin)

!       retrieve the site file name and index
        else if(mode .eq. 11) then
           sitlin = outsit
           i = windx

!       retrieve the site label
        else if(mode .eq. 12) then
           sitlin = sitlbl

!       write a binary site file
        else if(mode .eq. 0) then
          ! write(*,*) 'wrtsite outsit=',oslen, "  '",trim(outsit),"'"
          if(max(index(outsit,'.esa'),index(outsit,'.dsa')).ne.0) then

            ! preserve an important comment
            j = index(sitlbl, '!')
            if(j .ne. 0) then
              k = max(1, 40 - (len_trim(sitlbl) -j))
              sitlbl(k:) = sitlbl(j:)
            endif
            ! write(*,*) "sitlbl:'",trim(sitlbl),"'"
            call binsitrw(outsit(:oslen),0,windx, sitlbl, 17) ! daycent number >= 16

          ! write an ascii site file from stored info
          else if(mode .eq. 0 .and. oslen .gt. 0) then
            if(index(outsit,'.100') .eq. 0) then
              outsit(oslen+1:) = '.100'
              oslen = oslen+4
            endif
            call sitprint(outsit(:oslen),sitlbl,windx)
            call message("writing ASCII site file "//outsit(:oslen))
          endif

        else
           call message('Unknown wrtsite mode number')
        end if

      return
      end subroutine wrtsite


!.... sitprint
      subroutine sitprint(outsit,label,indx)
      USE wrtsitfil

      implicit none

!...Writes an ASCII site file in the enhanced format
!
!  Modifications
!  Jan 2009 K. Killian
!    allow printing to the standard output.

      include 'const.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'sitsoil.inc'       ! soils.in and sitpar.in variables
      include 'wth.inc'

      integer   i, j, k, indx, osf
      character label*(*), outsit*(*)
      character cindx*10
      real      ltc
      real      :: mrtfr, wmrtfr, arsum

      if(outsit .eq. 'STDOUT') then
        osf=6
      else
        osf=27
        open(osf, file=outsit, status="UNKNOWN", iostat=i)
        if(i .gt. 0) call abortrun('opening output site file '//trim(outsit))
      end if

!         call reratio(nelem,rcelit,rces1,rces2,rces3,som1ci,som1e,       &
!              som2ci,som2e,som3ci,som3e,metcis,metabe,strcis,struce)

         ! mature root fractions if we are doing a site file transfer
         if(bglcism(1) .eq. -1.0) then
           mrtfr = bglcism(2)
           bglcism = 0.
         else ! calculate the fraction from the live data.
           mrtfr  = sum(bglcism)
           arsum    = mrtfr + sum(bglcisj)
           if(arsum .gt. 0.0  .and.  mrtfr .gt. 0.0) then
             mrtfr  = mrtfr/arsum
           else
             mrtfr  = 0.0
           endif
         endif
         if(frtcism(1) .eq. -1.0) then
           wmrtfr = frtcism(2)
           frtcism = 0.
         else
           wmrtfr = sum(frtcism)
           arsum    = wmrtfr + sum(frtcisj)
           if(arsum .gt. 0.0  .and.  wmrtfr .gt. 0.0) then
             wmrtfr  = wmrtfr/arsum
           else
             wmrtfr  = 0.0
           endif
         endif

!**** write the first site 100 record
      write(cindx,'(i10)') indx
      write(osf,'(4a)') 'X    Archived site: ',trim(label),' #',trim(adjustl(cindx))
!      write(osf,'("X    Archived site: # ",2a)') adjustl(cindx),trim(label)
        k=2

!*** Climate parameters
      do i = 1, MONTHS
        call wrtlin(osf, k, precip(i))
      end do

      do i = 1, MONTHS
        call wrtlin(osf, k, prcstd(i))
      end do

      do i = 1, MONTHS
        call wrtlin(osf, k, prcskw(i))
      end do

      do i = 1, MONTHS
        call wrtlin(osf, k, tmn2m(i))
      end do

      do i = 1, MONTHS
        call wrtlin(osf, k, tmx2m(i))
      end do

!      do i = 1, MONTHS
!        call wrtlin(osf, k, cldcov(i+1))
!      end do

      do i = 1, MONTHS
        call wrtlin(osf, k, sradadj(i))
      end do

      call wrtlin(osf, k, hours_rain)

!*** Site and control parameters
      call wrtlin(osf, k, real(0.0)) ! IVAUTO = 0 ensures this data is read
      call wrtlin(osf, k, real(nelem))
      call wrtlin(osf, k, sitlat)
      call wrtlin(osf, k, sitlng)

      call wrtlin(osf, k, elevation)
      call wrtlin(osf, k, sitslp)
      call wrtlin(osf, k, aspect)
      call wrtlin(osf, k, ehoriz)
      call wrtlin(osf, k, whoriz)

!      call wrtlin(osf, k, sand)
!      call wrtlin(osf, k, silt)
!      call wrtlin(osf, k, clay)
!      call wrtlin(osf, k, rock)
!      call wrtlin(osf, k, bulkd)
!      call wrtlin(osf, k, real(nlayer))
      call wrtlin(osf, k, real(nlaypgst))
      call wrtlin(osf, k, drain)
      call wrtlin(osf, k, basef)
      call wrtlin(osf, k, stormf)
      call wrtlin(osf, k, precro)
      call wrtlin(osf, k, fracro)
      call wrtlin(osf, k, real(swflag))

!      do i = 1, CMXLYR
!        call wrtlin(osf, k, awilt(i))
!      end do
!      do i = 1, CMXLYR
!        call wrtlin(osf, k, afiel(i))
!      end do

!      call wrtlin(osf, k, ph)
      call wrtlin(osf, k, pslsrb)
      call wrtlin(osf, k, sorpmx)

      call wrtlin(osf, k, sublimscale)
      call wrtlin(osf, k, reflec)
      call wrtlin(osf, k, albedo)

      call wrtlin(osf, k, dmpflux)
      call wrtlin(osf, k, drainlag)
      call wrtlin(osf, k, hpotdeep)
      call wrtlin(osf, k, ksatdeep)

      call wrtlin(osf, k, dmp)
      call wrtlin(osf, k, Ncoeff)
      call wrtlin(osf, k, jdayStart)
      call wrtlin(osf, k, jdayEnd)
      call wrtlin(osf, k, N2Oadjust_fc)
      call wrtlin(osf, k, N2Oadjust_wp)
      call wrtlin(osf, k, MaxNitAmt)
      call wrtlin(osf, k, netmn_to_no3)
      call wrtlin(osf, k, wfpsdnitadj)
      call wrtlin(osf, k, N2N2Oadj)
!      call wrtlin(osf, k, flood_N2toN2O)
!      call wrtlin(osf, k, CO2_to_CH4)
!      call wrtlin(osf, k, frCH4emit) ! write the fraction methane emitted
!      call wrtlin(osf, k, frac_to_exudates)
!      call wrtlin(osf, k, Aeh)
!      call wrtlin(osf, k, Deh)
!      call wrtlin(osf, k, Beh_flood)
!      call wrtlin(osf, k, Beh_drain)
!      call wrtlin(osf, k, zero_root_frac)


!*** External nutrient input parameters
      call wrtlin(osf, k, epnfa(INTCPT))
      call wrtlin(osf, k, epnfa(SLOPE))
      call wrtlin(osf, k, epnfs(INTCPT))
      call wrtlin(osf, k, epnfs(SLOPE))
      call wrtlin(osf, k, satmos(INTCPT))
      call wrtlin(osf, k, satmos(SLOPE))
      call wrtlin(osf, k, sirri)
      call wrtlin(osf, k, afue)

!*** Organic matter initial values
      call wrtlin(osf, k, som1ci(SRFC,UNLABL))
      call wrtlin(osf, k, som1ci(SRFC,LABELD))
      call wrtlin(osf, k, som1ci(SOIL,UNLABL))
      call wrtlin(osf, k, som1ci(SOIL,LABELD))

      call wrtlin(osf, k, som2ci(SRFC,UNLABL))
      call wrtlin(osf, k, som2ci(SRFC,LABELD))
      call wrtlin(osf, k, som2ci(SOIL,UNLABL))
      call wrtlin(osf, k, som2ci(SOIL,LABELD))

      call wrtlin(osf, k, som3ci(UNLABL))
      call wrtlin(osf, k, som3ci(LABELD))

      do i = SRFC, SOIL
        do j = 1, MAXIEL
          if(j.le.nelem .and. som1e(i,j).ne.0) then
             rces1(i,j) = (som1ci(i,1)+som1ci(i,2)) / som1e(i,j)
          end if
          call wrtlin(osf, k, rces1(i,j))
        end do
      end do

!     rces = som2c/som2e
      do i = SRFC, SOIL
        do j = 1, MAXIEL
          if(j.le.nelem .and. som2e(i,j).ne.0) then
             rces2(i,j) = (som2ci(i,1)+som2ci(i,2)) / som2e(i,j)
          end if
          call wrtlin(osf, k, rces2(i,j))
        end do
      end do

      do j = 1, MAXIEL
        if(j.le.nelem .and. som3e(j).ne.0) then
           rces3(j) = (som3ci(1)+som3ci(2)) / som3e(j)
        end if
        call wrtlin(osf, k, rces3(j))
      end do

      call wrtlin(osf, k, clittr(SRFC,UNLABL))
      call wrtlin(osf, k, clittr(SRFC,LABELD))
      call wrtlin(osf, k, clittr(SOIL,UNLABL))
      call wrtlin(osf, k, clittr(SOIL,LABELD))

      do i = SRFC, SOIL
        ltc =  metcis(i,1)+metcis(i,2) + strcis(i,1)+strcis(i,2)
!        ltc =  clittr(i,1) + clittr(i,2)
!        ltc =  metabc(i) + strucc(i)
        do j = 1, MAXIEL
          if(j.le.nelem .and. metabe(i,j) + struce(i,j) .ne. 0) then
             rcelit(i,j) = ltc/(metabe(i,j) + struce(i,j))
          end if
          call wrtlin(osf, k, rcelit(i,j))
        end do
      end do

      call wrtlin(osf, k, aglcis(UNLABL))
      call wrtlin(osf, k, aglcis(LABELD))
      do i = 1, MAXIEL
        call wrtlin(osf, k, aglive(i))
      end do

      call wrtlin(osf, k, bglcisj(UNLABL)+bglcism(UNLABL))
      call wrtlin(osf, k, bglcisj(LABELD)+bglcism(LABELD))
      do i = 1, MAXIEL
        call wrtlin(osf, k, bglivej(i)+bglivem(i))
      end do

      call wrtlin(osf, k, stdcis(UNLABL))
      call wrtlin(osf, k, stdcis(LABELD))
      do i = 1, MAXIEL
        call wrtlin(osf, k, stdede(i))
      end do

      !*** Forest organic matter initial parameters
      call wrtlin(osf, k, rlvcis(UNLABL))
      call wrtlin(osf, k, rlvcis(LABELD))
      do i = 1, MAXIEL
        call wrtlin(osf, k, rleave(i))
      end do

      call wrtlin(osf, k, fbrcis(UNLABL))
      call wrtlin(osf, k, fbrcis(LABELD))
      do i = 1, MAXIEL
        call wrtlin(osf, k, fbrche(i))
      end do

      call wrtlin(osf, k, rlwcis(UNLABL))
      call wrtlin(osf, k, rlwcis(LABELD))
      do i = 1, MAXIEL
        call wrtlin(osf, k, rlwode(i))
      end do

      call wrtlin(osf, k, frtcisj(UNLABL)+frtcism(UNLABL))
      call wrtlin(osf, k, frtcisj(LABELD)+frtcism(LABELD))
      do i = 1, MAXIEL
        call wrtlin(osf, k, frootej(i)+frootem(i))
      end do

      call wrtlin(osf, k, crtcis(UNLABL))
      call wrtlin(osf, k, crtcis(LABELD))
      do i = 1, MAXIEL
        call wrtlin(osf, k, croote(i))
      end do

      call wrtlin(osf, k, wd1cis(UNLABL))
      call wrtlin(osf, k, wd1cis(LABELD))
      call wrtlin(osf, k, wd2cis(UNLABL))
      call wrtlin(osf, k, wd2cis(LABELD))
      call wrtlin(osf, k, wd3cis(UNLABL))
      call wrtlin(osf, k, wd3cis(LABELD))
      call wrtlin(osf, k, frstage)

      !*** Mineral initial parameters
      do i = 1, MAXIEL
        do j = 1, CMXLYR
          call wrtlin(osf, k, minerl(j,i))
        end do
      end do

      do i = 1, MAXIEL
        call wrtlin(osf, k, parent(i))
      end do

      do i = 1, MAXIEL
        call wrtlin(osf, k, secndy(i))
      end do

      call wrtlin(osf, k, occlud)

      do i = 1, CMXLYR
        call wrtlin(osf, k, rwcf(i))
      end do
      call wrtlin(osf, k, snlq)
      call wrtlin(osf, k, snow)
      call wrtlin(osf, k, SnowFlag)

      !*** soil.in soil water layers
!      dpthmx = 0;   lyblkd = 0;  fieldc = 0;   wiltpt = 0
!      ecoeff = 0;   tcoeff = 0;  sandfrac = 0;   clayfrac = 0
!      orgfrac = 0;   swclimit = 0;  satcond = 0;  lyrpH = 0

         do i = 1, MXSWLYR
           if(dpthmx(i) .eq. 0) then
             k=k+1
             cycle
           endif
           call wrtlin(osf, k, dpthmx(i))
           k = k+20
           call wrtlin(osf, k, lyblkd(i))
           k = k+20
           call wrtlin(osf, k, fieldc(i))
           k = k+20
           call wrtlin(osf, k, wiltpt(i))
           k = k+20
           call wrtlin(osf, k, ecoeff(i))
           k = k+20
           call wrtlin(osf, k, tcoeff(i))
           k = k+20
           call wrtlin(osf, k, sandfrac(i))
           k = k+20
           call wrtlin(osf, k, clayfrac(i))
           k = k+20
           call wrtlin(osf, k, orgfrac(i))
           k = k+20
           call wrtlin(osf, k, swclimit(i))
           k = k+20
           call wrtlin(osf, k, satcond(i))
           k = k+20
           call wrtlin(osf, k, lyrpH(i))
           k = k-231
         end do

      !*** Enhanced extend parameters
         k = k+231
         do j = SRFC, SOIL
           do i = UNLABL, LABELD
             if(clittr(j,i) .gt. 0.) then
               call wrtlin(osf, k, metcis(j,i)/clittr(j,i))
             else
               call wrtlin(osf, k, 0.0)
             endif
           end do
         end do

         do j = SRFC, SOIL
           do i = 1, MAXIEL
             if(metabe(j,i) .gt. 0.) then
               call wrtlin(osf, k, metabc(j)/metabe(j,i))
             else
               call wrtlin(osf, k, rces1(j,i)) ! 0.0)
             endif
           end do
         end do

         do i = SRFC, SOIL
           if(strucc(i) .gt. 0.) then
             call wrtlin(osf, k, strlig(i)/strucc(i))
           else
             call wrtlin(osf, k, 0.0)
           endif
         end do

         do i = 1, MAXIEL
           if(wood1e(i) .gt. 0.) then
             call wrtlin(osf, k, wood1c/wood1e(i))
           else
             call wrtlin(osf, k, 0.0)
           endif
         end do

         do i = 1, MAXIEL
           if(wood1e(i) .gt. 0.) then
             call wrtlin(osf, k, wood2c/wood2e(i))
           else
             call wrtlin(osf, k, 0.0)
           endif
         end do

         do i = 1, MAXIEL
           if(wood1e(i) .gt. 0.) then
             call wrtlin(osf, k, wood3c/wood3e(i))
           else
             call wrtlin(osf, k, 0.0)
           endif
         end do

         call wrtlin(osf, k, carbostg(CRPSYS,UNLABL))
         call wrtlin(osf, k, carbostg(CRPSYS,LABELD))
         call wrtlin(osf, k, carbostg(FORSYS,UNLABL))
         call wrtlin(osf, k, carbostg(FORSYS,LABELD))

         do i = 1, MAXIEL
           call wrtlin(osf, k, crpstg(i))
         end do

         do i = 1, MAXIEL
           call wrtlin(osf, k, forstg(i))
         end do


         ! use the actual mature/juvenile C fraction

         call wrtlin(osf, k,  mrtfr) ! wrtfrac
         call wrtlin(osf, k, wmrtfr) ! wmrtfrac

         call wrtlin(osf, k, real(ammonium))
         do i = 1, MXSWLYR
           call wrtlin(osf, k, real(nitrate(i)))
         end do

         do i = 1, MXSWLYR
           call wrtlin(osf, k, real(swc(i)))
         end do

      close(unit=osf)
      return
      end


!.... BINSITRW
      subroutine binsitrw(binam, rindx, windx, inlabel, vrsn)
      USE wrtsitfil

      implicit none

       character binam*(*), inlabel*(*)
       integer   rindx, windx, vrsn

!...This routine does the read/write access to the binary site file.
!     binam   name of file to open
!     rindx   index of record to read
!     windx   index of record to write
!     inlabel label stored in enhanced version
!     vrsn    version number
!
!  Modifications
!  July/2012 K. Killian
!    A negative vrsn makes the record checks non-fatal for sitarchive record scans
!    added a check on nelem and nlayer. If these are <=0  or greater than their
!    maxima then the record is assumed ot be bad (probably all zeros)
!  Feb/2008 K. Killian
!    initialize annual precipitation, prcann with climate data.
!      A value is needed for the initial crop lignin

!...includes that contain the site.100 variables
       include 'const.inc'
       include 'fertil.inc'
       include 'param.inc'
       include 'plot1.inc'
       include 'plot2.inc'
       include 'plot3.inc'
       include 'site.inc'
       include "sitsoil.inc"       ! soils.in and sitpar.in variables
       include 'wth.inc'

!...Local variables
       integer, parameter :: bin = 27
       integer      i, j, istat, lenam
       ! full daycent site input description, include  soils.in sitepar.in (627 lines)
       ! no usexdrvrs  that will be set/overridden by the weather file data
       ! no (-10) "pad variables for title lines, + version
       ! remove overwritten SAND, SILT, CLAY, BULKDEN, PH, AWILT and AFIEL (-25)
       ! 40 byte label
       integer, parameter  :: rlen2412 = 2412
       integer    :: rlen, filvrsn;
!       real      :: amnrl(MAXIEL) = 0
       character (len=40)  :: label
       character (len=8)   :: filstat
       character (len=256) :: mssg

       logical :: badsit
       logical :: fatal = .true.
       logical :: cheksitin, filexst
       real    :: mrtfr, wmrtfr

       real    :: null
       real    :: initfsw
       real    :: wtrtable(12) ! this will be ignored but its in the 16 format
           !  note cldcov was loaded as 1-12 in C leaving the 0 element open
           !       this translates to fortran elements 2:13
       real    :: cldcov(13)      ! this will be ignored but its in the 16 format
       real :: rmetcis(SOIL,LABELD), rstrlig(SOIL)
       real :: rcewod1(MAXIEL), rcewod2(MAXIEL), rcewod3(MAXIEL)

       ! write(*,*) "binsitrw '",binam,"'", rindx, windx,"'",inlabel,"'",vrsn

!      get the file name length for error processing
       lenam = len_trim(binam)
       label = " "
       rlen = rlen2412
       if(vrsn.lt.0) fatal = .false. ! make non fatal input checks
       filvrsn = abs(vrsn)


!****   set the file status to old on a read operation
       filstat='UNKNOWN'
       if(rindx.gt.0) filstat='OLD'
       if(rindx.le.0  .and.  windx.le.0) then
         mssg = 'please specify site record number for '//trim(binam)
         call abortrun(mssg)
       end if

       inquire(file=binam, exist=filexst)
       if(filexst) then
         filstat='OLD'
         ! call message('binsitrw open '//trim(binam)//' for version')
         open(unit=bin,file=binam,form='UNFORMATTED',status="OLD", &
            access='DIRECT',recl=44,iostat=istat,convert='BIG_ENDIAN')
         if(istat.ne.0) call abortrun('unable to open archive '//trim(binam))
         ! read the version from an existing file
         read(unit=bin, rec=1, iostat=istat, iomsg= mssg) label, filvrsn
         if(istat.ne.0) call abortrun('reading '//trim(mssg)//' file: '//trim(binam))
!DEF &    access='DIRECT',recl=rlen,iostat=istat,convert='BIG_ENDIAN')
!SUN &    access='DIRECT',recl=rlen,iostat=istat)
       else
         filstat='NEW'
         open(unit=bin,file=binam,form='UNFORMATTED',status="NEW", &
            access='DIRECT',recl=44,iostat=istat,convert='BIG_ENDIAN')
         write(unit=bin, rec=1, iostat=istat, iomsg= mssg)  &
               'New Archive                             ', filvrsn
         if(istat.ne.0) call abortrun('writing '//trim(mssg)//' file: '//trim(binam))
       endif
       close(unit=bin)
       if(filvrsn .eq. 16 .or. filvrsn .eq. 17) rlen = rlen2412

       open(unit=bin,file=binam,form='UNFORMATTED',status='OLD',  &
          access='DIRECT',recl=rlen,iostat=istat,convert='BIG_ENDIAN')
!DEF &    access='DIRECT',recl=rlen,iostat=istat,convert='BIG_ENDIAN')
!SUN &    access='DIRECT',recl=rlen,iostat=istat)

       if(istat.ne.0) then
          write(mssg,'("open failure ",i8,"  file: ",a)') istat, trim(binam)
          call abortrun(mssg)
       end if

       if(rindx.gt.0) then
         if(fatal) then
           write(mssg,'("reading v",I2," record ",i8,": ")') filvrsn,rindx
           mssg = trim(mssg)//trim(binam)
           call message (mssg)
         endif
         do i=1, len(label)
           if(iachar(label(i:i)) .lt. 32         .or.         &
              iachar(label(i:i)) .gt. 126) label(i:i) = ' '
         end do

         select case (filvrsn)
           case (16)
             call read16(bin, rindx, istat, mssg)

             N2N2Oadj     = 1.0       ! multiplier on N2/N2O ratio                         8/2013                        1.0
             netmn_to_no3 = 0.20      ! fraction of new net mineralization that goes to NO3 (0.0-1.0)
             wfpsdnitadj  = 1.        ! adjustment on inflection point for water filled pore space       1.0

           case (17)
             call read17(bin, rindx, istat, mssg)

           CASE DEFAULT
             write(mssg,*) filvrsn
             mssg = "unknown archive version v"//mssg
             call abortrun(trim(mssg)//': '//trim(binam))

         END SELECT

         if(istat.ne.0) then
           if(istat .lt. 0) then
             write(mssg,'(a,i8," from")') "EOF reading record ",rindx
           else
             i = len_trim(mssg)
             mssg = mssg(:i)//" record"
             write(mssg(i+8:),'(i8)') rindx
           end if
           mssg(len_trim(mssg)+2:) = binam
           call abortrun(mssg)
         end if

         ! convert from big-endian to little-endian
         ! call sitswend(vrsn,rivauto,rnelm,rnlyr,rnlypg,rswflg)

         ! remove any non printing characters from the label.
         ! not sure where they come from but clean existing archives
         do i=1, len(label)
           if(iachar(label(i:i)) .lt. 32         .or.         &
              iachar(label(i:i)) .gt. 126) label(i:i) = ' '
         end do
         inlabel = label//' '

         nlaypg   = nlaypgst

         badsit = .false.
         if (cheksitin(fatal) .or.      &
             nelem.lt.1 .or. nelem.gt.MAXIEL  .or.                        &
             nlaypg.le.0   .or. nlaypg.gt.CMXLYR) badsit = .true.

         numlyrs = maxloc(dpthmx,dim=1)  !  transfer(maxloc(dpthmx),1)

         if(.not. badsit) then
           do i=1, numlyrs
             if(sandfrac(i)<0 .or. sandfrac(i) >1. .or. &
                clayfrac(i)<0 .or. clayfrac(i) >1.) badsit = .true.
             ! not sure what condition this was looking for KLK 10May17
             !if(-0.3*sandfrac(i) + 15.7*clayfrac(i) + 3.1 .lt. 0.00001) badsit = .true.
           end do
         endif

         if (badsit) then
           write(mssg,*) "bad binary site record ",rindx
           if(fatal) then
             call abortrun(trim(mssg)//": "//trim(binam))
           else
             call message(mssg)
             return
           endif
         endif

!        added the code to set SOMnE from the C/E ratios
         som1c = sum(som1ci, dim=2)
         som2c = sum(som2ci, dim=2)
         som3c = sum(som3ci)

         do i=1,MAXIEL
           if(rcelit(SRFC,i).gt.0) struce(SRFC,i)=sum(clittr(SRFC,:))/rcelit(SRFC,i)-metabe(SRFC,i)
           if(rcelit(SOIL,i).gt.0) struce(SOIL,i)=sum(clittr(SOIL,:))/rcelit(SOIL,i)-metabe(SOIL,i)
           ! Compute N, P, and S for surface and soil som1, as well as for
           !  som2 and som3.   vek  08-91
           where (rces1(:,i).gt.0.) som1e(:,i) = som1c(:)/rces1(:,i)
           where (rces2(:,i).gt.0.) som2e(:,i) = som2c(:)/rces2(:,i)
           if(rces3(i).gt.0)        som3e(i) = som3c/rces3(i)
         end do

!------------------------------------------------------------------------
       else if(windx.gt.0) then
         label = inlabel//' ' !transfer 40 characters to the record label
         ! remove any non printing characters from the label.
         do i=1, len(label)
           if(iachar(label(i:i)) .lt. 32         .or.         &
              iachar(label(i:i)) .gt. 126) label(i:i) = ' '
         end do

         if(fatal) then
           write(mssg,'("writing v",i2," record ",i8)') filvrsn,windx
           mssg = trim(mssg)//"  label '"//trim(label)//"' : "//trim(binam)
           call message (mssg)
         endif

!****    calculate rcelit,rces1,rces2,rces3, the organic matter mineral ratios
         call reratio(rcelit,rces1,rces2,rces3,som1ci,som1e,       &
              som2ci,som2e,som3ci,som3e,metcis,metabe,strcis,struce)

!          swap to big-endian format
!           call sitswend(filvrsn,rivauto,rnelm,rnlyr,rnlypg,rswflg)

! SAND, SILT, CLAY, BULKDEN, PH, AWILT and AFIEL

         if(filvrsn.eq.16) then
           call write16(bin, windx, istat, mssg)
         else if(filvrsn.eq.17) then
           call write17(bin, windx, istat, mssg)
         end if

!        swap back to little-endian format
!         call sitswend(filvrsn,rivauto,rnelm,rnlyr,rnlypg,rswflg)
       else
         mssg= 'No site record index for binary read/write'
         call abortrun(mssg)
       end if

       if(istat.ne.0) then
         if(istat .lt.0) then
           mssg = 'EOF reading record '
           write(mssg(len_trim(mssg)+1:),'(i8," from")') windx
         else
           i = len_trim(mssg)
           write(mssg(i:),'(a,i8)') " record ",rindx
         end if
         mssg(len_trim(mssg)+2:) = binam
         call abortrun(mssg)
       end if
       ! write(*,*) "BINSITRW ",rindx, windx, strlig

      close(unit=bin)
      return

      CONTAINS
  !.... RERATIO
        subroutine reratio(rcelit,rces1,rces2,rces3,som1ci,som1e,    &
                           som2ci,som2e,som3ci,som3e,                      &
                           metcis,metabe,strcis,struce)

           real rcelit(2,3),rces1(2,3),rces2(2,3),rces3(3)
           real som1ci(2,2),som1e(2,3),som2ci(2,2),som2e(2,3)
           real som3ci(2),som3e(3)
           real metcis(2,2),metabe(2,3),strcis(2,2),struce(2,3)

          integer   i, j

  !****   calculate the slow and passive organic matter mineral ratios
          do j = 1, nelem
            do i= SRFC, SOIL
              ! litter C/E ratio
              if((metabe(i,j)+struce(i,j)) .ne. 0) rcelit(i,j) =          &
               (metcis(i,1)+metcis(i,2) + strcis(i,1)+strcis(i,2)) /      &
                          (metabe(i,j)+struce(i,j))
            end do
          end do

          ! som1 C/E ratio
          where(som1e(SRFC,1:nelem) .ne. 0.) rces1(SRFC,1:nelem) = sum(som1ci(SRFC,1:nelem))/som1e(SRFC,1:nelem)  ! som3 C/E ratio
          where(som1e(SOIL,1:nelem) .ne. 0.) rces1(SOIL,1:nelem) = sum(som1ci(SOIL,1:nelem))/som1e(SOIL,1:nelem)  ! som3 C/E ratio

          ! som2 C/E ratio
          where(som2e(SRFC,1:nelem) .ne. 0.) rces2(SRFC,1:nelem) = sum(som2ci(SRFC,1:nelem))/som2e(SRFC,1:nelem)  ! som3 C/E ratio
          where(som2e(SOIL,1:nelem) .ne. 0.) rces2(SOIL,1:nelem) = sum(som2ci(SOIL,1:nelem))/som2e(SOIL,1:nelem)  ! som3 C/E ratio

          where(som3e(1:nelem) .ne. 0.) rces3(1:nelem) = (som3ci(1)+som3ci(2))/som3e(1:nelem)  ! som3 C/E ratio
          return
        end


        subroutine write16(bin, indx, istat, mssg)
           integer :: bin, indx, istat
           character*(*) mssg
           real, parameter :: tbotmn=0.0, tbotmx=12.4, timlag=30.0 ! unused values

           real     N2Oadjust

           !store the average value and remove the old turnoverfrac
           ! This works for the old constant but you can't store the new
           ! limits in the old scaler
           N2Oadjust= (N2Oadjust_fc+N2Oadjust_wp)/2/0.02

           wtrtable = 0.
           !real rivauto, rnelm, rnlyr, rnlypg, rswflg
           write(unit=bin, rec=indx, iostat=istat, iomsg= mssg)           &
           label, filvrsn,                                                &
           !*** Climate parameters
           precip,prcstd,prcskw,tmn2m,tmx2m,cldcov(2:13),sradadj,hours_rain,     &
           !*** Site and control parameters
           0.0,real(nelem),sitlat,sitlng,elevation,sitslp,aspect,         &
           ehoriz, whoriz, rock, real(nlayer), real(nlaypgst),            &
           drain,basef,stormf,PRECRO,FRACRO,real(swflag),                 &
           pslsrb,sorpmx, sublimscale, reflec, albedo, tbotmn, tbotmx,    &
           dmp, timlag, ncoeff, real(jdayStart), real(jdayEnd), N2Oadjust,initfsw, &
           dmpflux,real(drainlag),hpotdeep,ksatdeep,wtrtable,             &
           !*** External nutrient input parameters
           epnfa,epnfs,satmos,sirri,                                      &
           !*** Organic matter initial values
           som1ci(SRFC,UNLABL),som1ci(SRFC,LABELD),                       &
           som1ci(SOIL,UNLABL),som1ci(SOIL,LABELD),                       &
           som2ci(SRFC,UNLABL),som2ci(SRFC,LABELD),                       &
           som2ci(SOIL,UNLABL),som2ci(SOIL,LABELD), som3ci,               &
           ((rces1(i,j), j= 1, MAXIEL), i= SRFC, SOIL),                   &
           ((rces2(i,j), j= 1, MAXIEL), i= SRFC, SOIL),rces3,             &
           clittr(SRFC,UNLABL),clittr(SRFC,LABELD),                       &
           clittr(SOIL,UNLABL),clittr(SOIL,LABELD),                       &
           ((rcelit(i,j), j= 1, MAXIEL), i= SRFC, SOIL),                  &
           aglcis,aglive,bglcisj,bglivej,stdcis,stdede,                   &
           !*** Forest organic matter initial parameters
           rlvcis,rleave,fbrcis,fbrche,rlwcis,rlwode,                     &
           frtcisj,frootej,crtcis,croote,wd1cis,wd2cis,wd3cis,            &
           !*** Mineral initial parameters
           minerl,parent,secndy,occlud,                                   &
           !*** Water initial parameters
           rwcf,snlq,snow,                                                &
           !*** soil parameters from soil.in
           dpthmx, lyblkd, fieldc, wiltpt, ecoeff, tcoeff,                &
           sandfrac, clayfrac, orgfrac, swclimit, satcond, lyrpH,         &
           !*** Enhanced extend parameters
           metcis(SRFC,UNLABL),metcis(SRFC,LABELD),                       &
           metcis(SOIL,UNLABL),metcis(SOIL,LABELD),                       &
           ((metabe(i,j), j= 1, MAXIEL), i= SRFC, SOIL),                  &
           strlig,wood1e,wood2e,wood3e,crpstg,forstg,                     &
           real(ammonium),real(nitrate),real(swc)
        end subroutine write16

        subroutine read16(bin, indx, istat, mssg)
           integer :: bin, indx, istat
           character*(*) mssg

           real rivauto, rnelm, rnlyr, rnlypg, rswflg
           real     ammonia             ! real ammonium
           real     lyrno3(MXSWLYR)     ! real nitrate
           real     swcinit(MXSWLYR)    ! real swc
           real     N2Oadjust, dnstart, dnend, drnlag
           real   :: tbotmn, tbotmx, timlag ! unused values

           read(unit=bin, rec=indx, iostat=istat, iomsg= mssg)            &
           label, filvrsn,                                                &
           !*** Climate parameters
           precip,prcstd,prcskw,tmn2m,tmx2m,cldcov(2:13),sradadj,hours_rain,    &
           !*** Site and control parameters
           rivauto,rnelm,sitlat,sitlng,elevation,sitslp,aspect,            &
           ehoriz, whoriz, rock, rnlyr, rnlypg,                           &
           drain,basef,stormf,PRECRO,FRACRO,rswflg,                       &
           pslsrb,sorpmx, sublimscale, reflec, albedo, tbotmn, tbotmx,    &
           dmp, timlag, ncoeff, dnstart, dnend, N2Oadjust,initfsw,        &
           dmpflux,drnlag,hpotdeep,ksatdeep,wtrtable,                   &
           !*** External nutrient input parameters
           epnfa,epnfs,satmos,sirri,                                      &
           !*** Organic matter initial values
           som1ci(SRFC,UNLABL),som1ci(SRFC,LABELD),                       &
           som1ci(SOIL,UNLABL),som1ci(SOIL,LABELD),                       &
           som2ci(SRFC,UNLABL),som2ci(SRFC,LABELD),                       &
           som2ci(SOIL,UNLABL),som2ci(SOIL,LABELD), som3ci,               &
           ((rces1(i,j), j= 1, MAXIEL), i= SRFC, SOIL),                   &
           ((rces2(i,j), j= 1, MAXIEL), i= SRFC, SOIL),rces3,             &
           clittr(SRFC,UNLABL),clittr(SRFC,LABELD),                       &
           clittr(SOIL,UNLABL),clittr(SOIL,LABELD),                       &
           ((rcelit(i,j), j= 1, MAXIEL), i= SRFC, SOIL),                  &
           aglcis,aglive,bglcisj,bglivej,stdcis,stdede,                   &
           !*** Forest organic matter initial parameters
           rlvcis,rleave,fbrcis,fbrche,rlwcis,rlwode,                     &
           frtcisj,frootej,crtcis,croote,wd1cis,wd2cis,wd3cis,            &
           !*** Mineral initial parameters
           minerl,parent,secndy,occlud,                                   &
           !*** Water initial parameters
           rwcf,snlq,snow,                                                &
           !*** soil parameters from soil.in
           dpthmx, lyblkd, fieldc, wiltpt, ecoeff, tcoeff,                &
           sandfrac, clayfrac, orgfrac, swclimit, satcond, lyrpH,         &
           !*** Enhanced extend parameters
           metcis(SRFC,UNLABL),metcis(SRFC,LABELD),                       &
           metcis(SOIL,UNLABL),metcis(SOIL,LABELD),                       &
           ((metabe(i,j), j= 1, MAXIEL), i= SRFC, SOIL),                  &
           strlig,wood1e,wood2e,wood3e,crpstg,forstg,                     &
           ammonia,lyrno3,swcinit

           if(istat .ne. 0) return      ! retrun and report i/O ERROR

           ivauto    = nint(rivauto)
           nelem     = nint(rnelm)
           nlayer    = nint(rnlyr)
           nlaypgst  = nint(rnlypg)
           swflag    = nint(rswflg)
           jdayStart = nint(dnstart)
           jdayEnd   = nint(dnend)
           drainlag  = nint(drnlag)
           ammonium  = ammonia
           nitrate   = lyrno3
           swc       = swcinit

           metabc = sum(metcis, dim=2) ! do here since sitein needs it to undo metabe
           strcis = clittr -metcis     !Set structural C as clittr - metcis
           strucc = sum(strcis, dim=2)

           ! set the new adjust parameters
           if(N2Oadjust_wp .eq. 0) then
             N2Oadjust_fc = N2Oadjust*0.02 ! maximum proportion of nitrified N lost as N2O @ field capacity0
             N2Oadjust_wp = N2Oadjust*0.02 ! minimum proportion of nitrified N lost as N2O @ wilting point
           endif

        end subroutine read16

        subroutine write17(bin, indx, istat, mssg)
           integer :: bin, indx, istat
           character*(*) mssg

           real :: mrtfr, wmrtfr, arsum
           real :: pad(3) = 0. ! words freed from v16  reserved for future use
           real :: null
           real, parameter :: tbotmn=0.0, tbotmx=12.4, timlag=30.0 ! unused values

           null = 0;
           rstrlig = 0
           strucc = sum(strcis, dim=2)
           where(strucc .gt. 0.) rstrlig = strlig/strucc

           rmetcis = 0
           where(clittr .gt. 0.) rmetcis = metcis/clittr

           where(metabe(SRFC,1:nelem) .gt. 0.) metabe(SRFC,1:nelem) = metabc(SRFC)/metabe(SRFC,1:nelem)
           where(metabe(SOIL,1:nelem) .gt. 0.) metabe(SOIL,1:nelem) = metabc(SOIL)/metabe(SOIL,1:nelem)

           where (wood1e .gt. 0.0)
             rcewod1 = sum(wd1cis) / wood1e
           else where
             rcewod1 = 0
           end where
           where (wood2e .gt. 0.0)
             rcewod2 = sum(wd2cis) / wood1e
           else where
             rcewod2 = 0
           end where
           where (wood3e .gt. 0.0)
             rcewod3 = sum(wd3cis) / wood1e
           else where
             rcewod3 = 0
           end where

           mrtfr    = sum(bglcism)
           arsum    = mrtfr + sum(bglcisj)
           if(arsum .gt. 0.0  .and.  mrtfr .gt. 0.0) then
             mrtfr  = mrtfr/arsum
           else
             mrtfr  = 0.0
           endif

         if(bglcism(1) .eq. -1.0) then
           mrtfr = bglcism(2)
           bglcism = 0.
         else ! calculate the fraction from the live data.
           mrtfr  = sum(bglcism)
           arsum    = mrtfr + sum(bglcisj)
           if(arsum .gt. 0.0  .and.  mrtfr .gt. 0.0) then
             mrtfr  = mrtfr/arsum
           else
             mrtfr  = 0.0
           endif
         endif
         if(frtcism(1) .eq. -1.0) then
           wmrtfr = frtcism(2)
           frtcism = 0.
         else
           wmrtfr = sum(frtcism)
           arsum    = wmrtfr + sum(frtcisj)
           if(arsum .gt. 0.0  .and.  wmrtfr .gt. 0.0) then
             wmrtfr  = wmrtfr/arsum
           else
             wmrtfr  = 0.0
           endif
         endif

           write(unit=bin, rec=indx, iostat=istat, iomsg= mssg)          &
           label, filvrsn,                                                &
           !*** Climate parameters
           precip,prcstd,prcskw,tmn2m,tmx2m,sradadj,hours_rain, &
           !*** Site and control parameters
           0.0,nelem,sitlat,sitlng,elevation,sitslp, &
           aspect,ehoriz,whoriz,0,nlaypgst,drain,  &
           basef,stormf,precro,fracro,swflag,pslsrb,sorpmx,  &
           sublimscale,reflec,albedo,dmpflux,drainlag,hpotdeep,ksatdeep,  &
           tbotmn,tbotmx,dmp,timlag,Ncoeff,jdayStart,jdayEnd,N2Oadjust_fc,  &
           N2Oadjust_wp,MaxNitAmt,SnowFlag,netmn_to_no3,wfpsdnitadj,N2N2Oadj,  &
! methane parameters in fix.100 save as place holders
!   flood_N2toN2O,CO2_to_CH4,C6H12O6_to_CH4,frac_to_exudates,Aeh,Deh,Beh_flood,Beh_drain,zero_root_frac,  &
           null,null,null,null,null,null,null,null,null, & !
           !*** External nutrient input parameters
           epnfa,epnfs,satmos,sirri,afue,  &
           !*** Organic matter initial values
           som1ci,som2ci,som3ci,rces1,rces2,rces3,clittr,rcelit,  &
           aglcis,aglive,bglcisj+bglcism,bglivej+bglivem,stdcis,stdede,  &
           !*** Forest organic matter initial parameters
           rlvcis,rleave,fbrcis,fbrche,rlwcis,rlwode,  &
           frtcisj+frtcism,frootej+frootem,crtcis,croote,wd1cis,wd2cis,wd3cis,  &
           !*** Mineral initial parameters
           minerl,parent,secndy,occlud,  &
           !*** Water initial parameters
           rwcf,snlq,snow,  &
           !*** soil parameters from soil.in
           dpthmx,lyblkd,fieldc,wiltpt,ecoeff,tcoeff,sandfrac,clayfrac,  &
           orgfrac,swclimit,satcond,lyrpH,  &
           !*** Enhanced extend parameters
           rmetcis,metabe,rstrlig,rcewod1,rcewod2,rcewod3,carbostg,crpstg,forstg,  &
           mrtfr, wmrtfr, real(ammonium),real(nitrate),real(swc), & ! real*4 versions
           frstage, pad ! write forest age and pad the original record length
        end subroutine write17


        subroutine read17(bin, indx, istat, mssg)
           integer :: bin, indx, istat
           character*(*) mssg

           real ::  ammonia          ! real ammonium
           real ::  lyrno3(MXSWLYR)  ! real nitrate
           real ::  swcinit(MXSWLYR) ! real swc
           ! unused variables; maintained for compatability
           integer :: nlayr     ! removed nlayer
           real :: pad(3)            ! words freed from v16  reserved for future use
           real :: methp(9)          ! moved methane parameters
           real :: tbotmn, tbotmx, timlag ! unused values
           integer :: ilyr
           real, parameter :: PARTDENS = 2.65 ! Particle Density (g/cm3)
           real            :: porosity        ! calculated porosity


           methp = 0

           mssg = ''
           read(unit=bin, rec=rindx, iostat=istat, iomsg= mssg)          &
           label, filvrsn,                                               &
           ! precip,prcstd,prcskw,tmn2m,tmx2m,cldcov,sradadj,hours_rain,
           precip,prcstd,prcskw,tmn2m,tmx2m,sradadj,hours_rain, &
           !*** Site and control parameters
           ivauto,nelem,sitlat,sitlng,elevation,sitslp, &
           aspect,ehoriz,whoriz,nlayr,nlaypgst,drain,  &
           basef,stormf,precro,fracro,swflag,pslsrb,sorpmx,  &
           sublimscale,reflec,albedo,dmpflux,drainlag,hpotdeep,ksatdeep,  &
           tbotmn,tbotmx,dmp,timlag,Ncoeff,jdayStart,jdayEnd,N2Oadjust_fc,  &
           N2Oadjust_wp,MaxNitAmt,SnowFlag,netmn_to_no3,wfpsdnitadj,N2N2Oadj,  &
! methane parameters in fix.100 save as place holders
!   flood_N2toN2O,CO2_to_CH4,C6H12O6_to_CH4,frac_to_exudates,Aeh,Deh,Beh_flood,Beh_drain,zero_root_frac,  &
           methp,  &
           !*** External nutrient input parameters
           epnfa,epnfs,satmos,sirri,afue,  &
           !*** Organic matter initial values
           som1ci,som2ci,som3ci,rces1,rces2,rces3,clittr,rcelit,  &
           aglcis,aglive,bglcisj,bglivej,stdcis,stdede,  &
           !*** Forest organic matter initial parameters
           rlvcis,rleave,fbrcis,fbrche,rlwcis,rlwode,  &
           frtcisj,frootej,crtcis,croote,wd1cis,wd2cis,wd3cis,  &
           !*** Mineral initial parameters
           minerl,parent,secndy,occlud,  &
           !*** Water initial parameters
           rwcf,snlq,snow,  &
           !*** soil parameters from soil.in
           dpthmx,lyblkd,fieldc,wiltpt,ecoeff,tcoeff,sandfrac,clayfrac,  &
           orgfrac,swclimit,satcond,lyrpH,  &
           !*** Enhanced extend parameters
           metcis,metabe,strlig,wood1e,wood2e,wood3e,carbostg,crpstg,forstg,  &
           mrtfr, wmrtfr, ammonia,lyrno3, swcinit, & ! store real*4 versions
           frstage, pad ! read new forest age and remaining extra words from the original record

           if(istat .ne. 0) return      ! return and report i/O ERROR

           metcis = metcis*clittr      ! metcis is really a ratio of metcis/clittr
           metabc = sum(metcis, dim=2) ! do here since we need it to undo metabe
           strcis = clittr - metcis
           strucc = sum(strcis, dim=2)

           ! we read a ratio. Convert the ratio to actual values
           where(metabe(SRFC,1:nelem) .gt. 0) metabe(SRFC,1:nelem) = metabc(SRFC)/metabe(SRFC,1:nelem)
           where(metabe(SOIL,1:nelem) .gt. 0) metabe(SOIL,1:nelem) = metabc(SOIL)/metabe(SOIL,1:nelem)

           strlig = strlig * strucc
!           where(rcelit(SRFC,1:nelem) >0) struce(SRFC,1:nelem) = sum(clittr(SRFC,1:nelem))/rcelit(SRFC,:) -metabe(SRFC,1:nelem)
!           where(rcelit(SRFC,1:nelem) >0) struce(SOIL,1:nelem) = sum(clittr(SRFC,1:nelem))/rcelit(SOIL,:) -metabe(SOIL,1:nelem)

           ! this version strlig is the ratio

           ! version 17 wood1e is stored as a ratio
           where (wood1e(:nelem) .gt. 0.0 ) wood1e(:nelem) = sum(wd1cis) / wood1e(:nelem)
           ! else where ; wood1e = 0 ; end where
           where (wood2e(:nelem) .gt. 0.0 ) wood2e(:nelem) = sum(wd2cis) / wood2e(:nelem)
           ! else where ; wood2e = 0 ; end where
           where (wood3e(:nelem) .gt. 0.0 ) wood3e(:nelem) = sum(wd3cis) / wood3e(:nelem)
           ! else where ; wood3e = 0 ; end where

           ! save the maure root fraction in the C array for DETIV

           ! transfer the fractions if we are doing a site file transfer
           bglcism = 0.
           frtcism = 0.
           if(mrtfr .gt. 0.0)  bglcism= [-1.0,  mrtfr]
           if(wmrtfr .gt. 0.0) frtcism= [-1.0,  wmrtfr]

           ammonium = ammonia
           nitrate = lyrno3
           swc = swcinit
           null = 0;

           ilyr = 1;
           do while (lyblkd(ilyr) > 0.)
             porosity = 0.99985 - lyblkd(ilyr) / PARTDENS;
             ! Calculation of thetas_bd added 9/18/00. -cindyk
             ! save as fraction since it is always used that way  KLK 28 Apr 2015
             thetas_bd(ilyr) = 0.95*porosity;

             ! check that field capacity is less than saturation. KLK 14 Nov 2016
             ! the inversion can occur in the Saxton equations with high (>70%) clay soils.
             if (fieldc(ilyr) > thetas_bd(ilyr)) then
!               if(ilyr==1 .or. lyblkd(ilyr) /= lyblkd(ilyr-1)) then
!write(*, '("Warning: lowering field capacity ",f9.7," TO ",f9.7," below ",f9.7," in layer ",i2," to match Bulk density",f9.7)') &
!                        fieldc(ilyr),porosity*(61.0/64.0),porosity,ilyr,lyblkd(ilyr);
!               endif
               fieldc(ilyr) = porosity*(121.0/128.0);
             endif
             ilyr = ilyr+1
           enddo

         if(sum(rlwcis).ne.0  .and.  frstage == 0) call message(&
            'restarting existing forest with juvenile growth');

        end subroutine read17
      end


!.... cheksitin
      logical function cheksitin(fatal)
      USE wrtsitfil

      include 'const.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'sitsoil.inc'
      include 'wth.inc'

      logical fatal

!...Checks the values input from a binary site archive for valid values
!    skip the values read from sitepar.in and soils.in. Although it is
!    possible for the transfer to corrupt them, or even the C structures
!    to be overwritten, in is MUCH more likely those errors will munge the
!    values than generate NaN and INF.
!


      integer   i, j, k
      integer   :: lbuf
      character :: erlst*132

!       count the variable we are looking at
        k=0
        erlst = ''
        lbuf  = 0
        cheksitin = .false.

      !*** Climate parameters
      do i = 1, MONTHS
        if(test754(precip(i))) call nambad("precip", erlst, lbuf)
      end do
      do i = 1, MONTHS
        if(test754(prcstd(i))) call nambad("prcstd", erlst, lbuf)
      end do
      do i = 1, MONTHS
        if(test754(prcskw(i))) call nambad("prcskw", erlst, lbuf)
      end do
      do i = 1, MONTHS
        if(test754(tmn2m(i))) call nambad("tmn2m", erlst, lbuf)
      end do
      do i = 1, MONTHS
        if(test754(tmx2m(i))) call nambad("tmx2m", erlst, lbuf)
      end do
      ! do i = 1, MONTHS
      !   if(test754(cldcov(i))) call nambad("cldcov", erlst, lbuf)
      ! end do
      ! do i = 1, MONTHS
      !   if(test754(sradadj(i))) call nambad("sradadj", erlst, lbuf)
      ! end do
      ! if(test754(hours_rain)) call nambad("hours_rain", erlst, lbuf)

      !*** Site and control parameters
      !if(test754(rivauto)) call nambad("rivauto", erlst, lbuf)
      !if(test754(rnelm)) call nambad("rnelm",  erlst, lbuf)
      !if(test754(rnlyr)) call nambad("rnlyr", erlst, lbuf)
      !if(test754(rnlypg)) call nambad("rnlypg", erlst, lbuf)
      !if(test754(rswflg)) call nambad("rswflg", erlst, lbuf)

      if(test754(sitlat)) call nambad("sitlat", erlst, lbuf)
      if(test754(sitlng)) call nambad("sitlng", erlst, lbuf)

       if(test754(elevation)) call nambad("elevation", erlst, lbuf)
       if(test754(sitslp)) call nambad("sitslp", erlst, lbuf)
       if(test754(aspect)) call nambad("aspect", erlst, lbuf)
       if(test754(ehoriz)) call nambad("ehoriz", erlst, lbuf)
       if(test754(whoriz)) call nambad("whoriz", erlst, lbuf)

      if(test754(drain)) call nambad("drain", erlst, lbuf)
      if(test754(basef)) call nambad("basef", erlst, lbuf)
      if(test754(stormf)) call nambad("stormf", erlst, lbuf)
      if(test754(precro)) call nambad("precro", erlst, lbuf)
      if(test754(fracro)) call nambad("fracro", erlst, lbuf)

      if(test754(pslsrb)) call nambad("pslsrb", erlst, lbuf)
      if(test754(sorpmx)) call nambad("sorpmx", erlst, lbuf)

      ! if(test754(sublimscale)) call nambad("sublimscale", erlst, lbuf)
      ! if(test754(reflec)) call nambad("reflec", erlst, lbuf)
      ! if(test754(albedo)) call nambad("albedo", erlst, lbuf)
      ! if(test754(tbotmn)) call nambad("tbotmn", erlst, lbuf)
      ! if(test754(tbotmx)) call nambad("tbotmx", erlst, lbuf)
      ! if(test754(dmp)) call nambad("dmp", erlst, lbuf)
      ! if(test754(timlag)) call nambad("timlag", erlst, lbuf)
      ! if(test754(Ncoeff)) call nambad("Ncoeff", erlst, lbuf)
      ! if(test754(jdayStart)) call nambad("jdayStart", erlst, lbuf)
      ! if(test754(jdayEnd)) call nambad("jdayEnd", erlst, lbuf)
      ! if(test754(N2Oadjust)) call nambad("N2Oadjust", erlst, lbuf)
      ! if(test754(fswcinit)) call nambad("fswcinit", erlst, lbuf)
      ! if(test754(dmpflux)) call nambad("dmpflux", erlst, lbuf)
      ! if(test754(hpotdeep)) call nambad("hpotdeep", erlst, lbuf)
      ! if(test754(ksatdeep)) call nambad("ksatdeep", erlst, lbuf)

      !*** External nutrient input parameters
      if(test754(epnfa(INTCPT))) call nambad("epnfa", erlst, lbuf)
      if(test754(epnfa(SLOPE))) call nambad("epnfa", erlst, lbuf)
      if(test754(epnfs(INTCPT))) call nambad("epnfs", erlst, lbuf)
      if(test754(epnfs(SLOPE))) call nambad("epnfs", erlst, lbuf)
      if(test754(satmos(INTCPT))) call nambad("satmos", erlst, lbuf)
      if(test754(satmos(SLOPE))) call nambad("satmos", erlst, lbuf)
      if(test754(sirri)) call nambad("sirri", erlst, lbuf)
      if(test754(afue)) call nambad("afue", erlst, lbuf)

      !*** Organic matter initial values
      if(test754(som1ci(SRFC,UNLABL))) call nambad("som1ci", erlst, lbuf)
      if(test754(som1ci(SRFC,LABELD))) call nambad("som1ci", erlst, lbuf)
      if(test754(som1ci(SOIL,UNLABL))) call nambad("som1ci", erlst, lbuf)
      if(test754(som1ci(SOIL,LABELD))) call nambad("som1ci", erlst, lbuf)
      if(test754(som2ci(SRFC,UNLABL))) call nambad("som2ci", erlst, lbuf)
      if(test754(som2ci(SRFC,LABELD))) call nambad("som2ci", erlst, lbuf)
      if(test754(som2ci(SOIL,UNLABL))) call nambad("som2ci", erlst, lbuf)
      if(test754(som2ci(SOIL,LABELD))) call nambad("som2ci", erlst, lbuf)
      if(test754(som3ci(UNLABL))) call nambad("som3ci", erlst, lbuf)
      if(test754(som3ci(LABELD))) call nambad("som3ci", erlst, lbuf)

      do i = SRFC, SOIL
        do j = 1, MAXIEL
          if(test754(rces1(i,j))) call nambad("rces1", erlst, lbuf)
        end do
      end do

!     rces = som2c/som2e
      do i = SRFC, SOIL
        do j = 1, MAXIEL
          if(test754(rces2(i,j))) call nambad("rces2", erlst, lbuf)
        end do
      end do

      do j = 1, MAXIEL
        if(test754(rces3(j))) call nambad("rces3", erlst, lbuf)
      end do

      if(test754(clittr(SRFC,UNLABL))) call nambad("clittr", erlst, lbuf)
      if(test754(clittr(SRFC,LABELD))) call nambad("clittr", erlst, lbuf)
      if(test754(clittr(SOIL,UNLABL))) call nambad("clittr", erlst, lbuf)
      if(test754(clittr(SOIL,LABELD))) call nambad("clittr", erlst, lbuf)

      do i = SRFC, SOIL
        do j = 1, MAXIEL
          if(test754(rcelit(i,j))) call nambad("rcelit", erlst, lbuf)
        end do
      end do

      if(test754(aglcis(UNLABL))) call nambad("aglcis", erlst, lbuf)
      if(test754(aglcis(LABELD))) call nambad("aglcis", erlst, lbuf)
      do i = 1, MAXIEL
        if(test754(aglive(i))) call nambad("aglive", erlst, lbuf)
      end do

      if(test754(bglcisj(UNLABL))) call nambad("bglcis", erlst, lbuf)
      if(test754(bglcisj(LABELD))) call nambad("bglcis", erlst, lbuf)
      do i = 1, MAXIEL
        if(test754(bglivej(i))) call nambad("bglive", erlst, lbuf)
      end do

      if(test754(stdcis(UNLABL))) call nambad("stdcis", erlst, lbuf)
      if(test754(stdcis(LABELD))) call nambad("stdcis", erlst, lbuf)
      do i = 1, MAXIEL
        if(test754(stdede(i))) call nambad("stdede", erlst, lbuf)
      end do

      !*** Forest organic matter initial parameters
      if(test754(rlvcis(UNLABL))) call nambad("rlvcis", erlst, lbuf)
      if(test754(rlvcis(LABELD))) call nambad("rlvcis", erlst, lbuf)
      do i = 1, MAXIEL
        if(test754(rleave(i))) call nambad("rleave", erlst, lbuf)
      end do

      if(test754(fbrcis(UNLABL))) call nambad("fbrcis", erlst, lbuf)
      if(test754(fbrcis(LABELD))) call nambad("fbrcis", erlst, lbuf)
      do i = 1, MAXIEL
        if(test754(fbrche(i))) call nambad("fbrche", erlst, lbuf)
      end do

      if(test754(rlwcis(UNLABL))) call nambad("rlwcis", erlst, lbuf)
      if(test754(rlwcis(LABELD))) call nambad("rlwcis", erlst, lbuf)
      do i = 1, MAXIEL
        if(test754(rlwode(i))) call nambad("rlwode", erlst, lbuf)
      end do

      if(test754(frtcisj(UNLABL))) call nambad("frtcis", erlst, lbuf)
      if(test754(frtcisj(LABELD))) call nambad("frtcis", erlst, lbuf)
      do i = 1, MAXIEL
        if(test754(frootej(i))) call nambad("froote", erlst, lbuf)
      end do

      if(test754(crtcis(UNLABL))) call nambad("crtcis", erlst, lbuf)
      if(test754(crtcis(LABELD))) call nambad("crtcis", erlst, lbuf)
      do i = 1, MAXIEL
        if(test754(croote(i))) call nambad("croote", erlst, lbuf)
      end do

      if(test754(wd1cis(UNLABL))) call nambad("wd1cis", erlst, lbuf)
      if(test754(wd1cis(LABELD))) call nambad("wd1cis", erlst, lbuf)
      if(test754(wd2cis(UNLABL))) call nambad("wd2cis", erlst, lbuf)
      if(test754(wd2cis(LABELD))) call nambad("wd2cis", erlst, lbuf)
      if(test754(wd3cis(UNLABL))) call nambad("wd3cis", erlst, lbuf)
      if(test754(wd3cis(LABELD))) call nambad("wd3cis", erlst, lbuf)

      !*** Mineral initial parameters
      do i = 1, MAXIEL
        do j = 1, CMXLYR
          if(test754(minerl(j,i))) call nambad("minerl", erlst, lbuf)
        end do
      end do

      do i = 1, MAXIEL
        if(test754(parent(i))) call nambad("parent", erlst, lbuf)
      end do

      do i = 1, MAXIEL
        if(test754(secndy(i))) call nambad("secndy", erlst, lbuf)
      end do

      if(test754(occlud)) call nambad("occlud", erlst, lbuf)

      !*** Water initial parameters
      do i = 1, CMXLYR
        if(test754(rwcf(i))) call nambad("rwcf", erlst, lbuf)
      end do

      if(test754(snlq)) call nambad("snlq", erlst, lbuf)
      if(test754(snow)) call nambad("snow", erlst, lbuf)

      !*** soil parameters from soil.in
         ! do i = 1, MAXIEL
         !   if(test754(dpthmx(i))) call nambad("dpthmx", erlst, lbuf)
         !   if(test754(lyblkd(i))) call nambad("lyblkd", erlst, lbuf)
         !   if(test754(fieldc(i))) call nambad("feldc", erlst, lbuf)
         !   if(test754(wiltpt(i))) call nambad("wilpt", erlst, lbuf)
         !   if(test754(ecoeff(i))) call nambad("ecoef", erlst, lbuf)
         !   if(test754(tcoeff(i))) call nambad("tcoef", erlst, lbuf)
         !   if(test754(sandfrac(i))) call nambad("sandf", erlst, lbuf)
         !   if(test754(clayfrac(i))) call nambad("clayf", erlst, lbuf)
         !   if(test754(orgfrac(i))) call nambad("orgf", erlst, lbuf)
         !   if(test754(swclimit(i))) call nambad("swclimit", erlst, lbuf)
         !   if(test754(satcond(i))) call nambad("satcnd", erlst, lbuf)
         !   if(test754(lyrpH(i))) call nambad("phlyr", erlst, lbuf)
         ! end do

      !*** Enhanced extend parameters
         if(test754(metcis(SRFC,UNLABL))) call nambad("metcis", erlst, lbuf)
         if(test754(metcis(SRFC,LABELD))) call nambad("metcis", erlst, lbuf)
         if(test754(metcis(SOIL,UNLABL))) call nambad("metcis", erlst, lbuf)
         if(test754(metcis(SOIL,LABELD))) call nambad("metcis", erlst, lbuf)

         do i = 1, MAXIEL
           if(test754(metabe(SRFC,i))) call nambad("metabe", erlst, lbuf)
         end do
         do i = 1, MAXIEL
           if(test754(metabe(SOIL,i))) call nambad("metabe", erlst, lbuf)
         end do

         if(test754(strlig(1))) call nambad("strlig", erlst, lbuf)
         if(test754(strlig(2))) call nambad("strlig", erlst, lbuf)
         do i = 1, MAXIEL
           if(test754(wood1e(i))) call nambad("wood1e", erlst, lbuf)
         end do

         do i = 1, MAXIEL
           if(test754(wood2e(i))) call nambad("wood2e", erlst, lbuf)
         end do

         do i = 1, MAXIEL
           if(test754(wood3e(i))) call nambad("wood3e", erlst, lbuf)
         end do

         do i = 1, MAXIEL
           if(test754(crpstg(i))) call nambad("crpstg", erlst, lbuf)
         end do

         do i = 1, MAXIEL
           if(test754(forstg(i))) call nambad("forstg", erlst, lbuf)
         end do

         if(test754(real(ammonium))) call nambad("ammoni", erlst, lbuf)
         do i = 1, MXSWLYR
           if(test754(real(nitrate(i)))) call nambad("nitrat", erlst, lbuf)
         end do

         if(erlst .ne. '') then
           cheksitin = .true.
           call message("IEEE floating exception (NaN or INF) in site variables ")
           if(fatal) then
              call abortrun(erlst)
           else
              call message('ERROR: '//erlst)
           endif
         endif

         return
      CONTAINS

        subroutine nambad(vnam, erlst, lbuf)
          character :: erlst*(*), vnam*(*)
          integer   :: lbuf

          ! add the variable name to the error list

          integer   :: l
          l = len(vnam)
          if(lbuf .eq. 0) then
            erlst = vnam
          elseif(lbuf+l .lt. len(erlst)  .and.          &
     &         index(erlst,vnam,back=.true.).eq.0) then
            erlst(lbuf:) = ', '//vnam
            lbuf = lbuf+l+2
          endif
          return
        end subroutine nambad


        logical function test754(value)
          ! IEEE 754 specifies that NaN's are not equil. Thus, unless the compiler
          !  optimizes away the IF statement, this is a quick test for NaN values
          !  use isnan f95 intrinsic
          ! Test for IEEE 754 NAN; (v .ne. v .or. 1.0e37/(1.+abs(v)) .eq. 0.0)
          test754 = (ISNAN(value)  .or.  abs(value) .gt. huge(value))
        end function test754

      end
