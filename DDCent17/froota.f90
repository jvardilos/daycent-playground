
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

      real function froota(a2drat,h2ogef,systype)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'

! ... Argument declarations
      integer systype
      real    a2drat(MAXIEL), h2ogef

! ... This function determines the fraction of production going to fine
! ... roots in crops, grasses, and woody plants based based on water
! ... and nutrient availability.
!
!    Modified great plains equation to avoid the pole in the equation KLK May 2014
!    This corrects problems with the equation in desert areas like the Southwest.
!
! ... CALLED FROM:  cropDynC
! ...               treeDynC
!
! ... a2drat(iel) - the ratio of available mineral to mineral demand by
! ...               the plant
! ... h2ogef      - the effect of water on root to shoot ratio,
! ...               computed in POTCRP or POTFOR
! ... h2oeff      - effect of water stress on fraction of root carbon
! ... ntreff      - effect of nutrient stress in fraction of root carbon
! ... rtsh        - root to shoot ratio
!
! ... For trees only:
! ...   tfrtcn(1) - maximum fraction of C allocated to fine roots under
! ...               maximum nutrient stress
! ...   tfrtcn(2) - minimum fraction of C allocated to fine roots with no
! ...               nutrient stress
! ...   tfrtcw(1) - maximum fraction of C allocated to fine roots under
! ...               maximum water stress
! ...   tfrtcw(2) - minimum fraction of C allocated to fine roots with no
! ...               water stress
!
! ... For grasses/crops only (crop.100):
! ...   frtcindx  - (0) Use Great Plains eqn,
! ...               (1) perennial plant,
! ...               (2) annual plant,
! ...               (3) perennial, use growing degree day implementation
! ...               (4) non-grain filling annual plant, growing degree day
! ...                   implementation, dynamic carbon allocation
! ...               (5) grain filling annual plant, growing degree day
! ...                   implementation, dynamic carbon allocation
! ...               (6) grain filling annual plant that requires a
! ...                   vernalization period (i.e. winter wheat), growing
! ...                   degree day implementation, dynamic carbon allocation
! ...   frtc(1)   - fraction of C allocated to roots at planting, with no
! ...               water or nutrient stress, used when FRTCINDX = 2, 4, 5, or
! ...               6
! ...   frtc(2)   - fraction of C allocated to roots at time FRTC(3), with no
! ...               water or nutrient stress, used when FRTCINDX = 2, 4, 5, or
! ...               6
! ...   frtc(3)   - time after planting (months with soil temperature greater
! ...               than RTDTMP) at which the FRTC(2) value is reached, used
! ...               when FRTCINDX = 2, 4, 5, or 6
! ...   frtc(4)   - maximum increase in the fraction of C going to the roots
! ...               due to water stress, used when FRTCINDX = 2, 4, 5, or 6
! ...   frtc(5)   - maximum increase in the fraction of C going to the roots
! ...               due to nutrient stress, used when FRTCINDX = 2, 4, 5, or 6
! ...   cfrtcn(1) - maximum fraction of C allocated to roots under maximum
! ...               water stress, used when FRTCINDX = 1 or 3
! ...   cfrtcn(2) - minimum fraction of C allocated to roots with no water
! ...               stress, used when FRTCINDX = 1 or 3
! ...   cfrtcw(1) - maximum fraction of C allocated to roots under maximum
! ...               nutrient stress, used when FRTCINDX = 1 or 3
! ...   cfrtcw(2) - minimum fraction of C allocated to roots with no nutrient
! ...               stress, used when FRTCINDX = 1 or 3

! ... Function declarations
      real ramp
      external ramp

! ... Local Variables
      integer iel
      real    h2oeff, ntreff
      real    rtsh
      real    effprc
      character  bffr*8

      do iel = 1, nelem
        if ((a2drat(iel) .lt. 0.0) .or. (a2drat(iel) .gt. 1.0)) then
          call abortrun('a2drat out of bounds in froota')
        endif
      end do

      froota = -1.0

      if (systype .eq. FORSYS) then
! ..... Effect of water limitation - inverse of growth curve
! ..... Compute water stress limitation on trees the same way
! ..... as we are computing water stress limitation for crops
! ..... and grasses cak - 09/12/03
!        h2oeff = frfrac(1) + ((frfrac(2) - frfrac(1)) * (1.0 - h2ogef))
        h2oeff = line(h2ogef, 0.0, tfrtcw(1), 1.0, tfrtcw(2))
! ..... Effect of nutrient limitation
        ntreff = 0.0
        do iel = 1, nelem
          ntreff = max(line(a2drat(iel), 0.0, tfrtcn(1), 1.0, tfrtcn(2)), ntreff)
        end do
! ..... Compute fraction of C to go to fine roots
        froota = max(h2oeff, ntreff)
! ..... Change fraction of C to go to fine roots by effect of CO2
        froota = froota * co2crs(FORSYS)
        froota = min(froota, 1.0)

      else if (systype .eq. CRPSYS) then

        if (frtcindx .eq. 0) then
! ....... Use Great Plains equation for root to shoot ratio
! ....... Done in cycle, no need to call this routine again here
!          call prcgrw(month)
          ! Limit the great plains equation          KLK 8/15/2002
          ! the current parameters for the great plains equation have a pole at
          ! about 5 cm annual precip. To avoid this force the yearly precipitation
          ! used to be greater than both the pole and zero
          effprc = max(grwprc, 1.0)
          if(agppb.ne.0.) effprc = max(effprc, -agppa/agppb+2.0)
          rtsh = (bgppa + effprc*bgppb)/(agppa + effprc*agppb)
          ! cap the root to shoot ratio to 5. or 83% of production
          if    (rtsh .gt. 5.) then
            rtsh = 5.
          elseif(rtsh .lt. 0.25) then
            rtsh = 0.25
          endif
          froota = rtsh / (rtsh +1.0) ! 1.0 / (1.0/rtsh + 1.0)

        elseif (frtcindx .eq. 1 .or. frtcindx .eq. 3) then
! ....... A perennial plant (grass)
! ....... Effect of water limitation
          h2oeff = line(h2ogef, 0.0, cfrtcw(1), 1.0, cfrtcw(2))
! ....... Effect of nutrient limitation
          ntreff = 0.0
          do iel = 1, nelem
            ntreff = max(line(a2drat(iel), 0.0, cfrtcn(1), 1.0, cfrtcn(2)), ntreff)
          end do
! ....... Compute fraction of C to go to roots
          froota = max(h2oeff, ntreff)
          froota = min(froota, 1.0)

        elseif (frtcindx .eq. 2 .or. frtcindx .ge. 4) then
! ....... An annual plant (crop)
! ....... Compute fraction of C to the roots when there is no
! ....... nutrient or water stress, cak - 09/12/03
          froota = ramp(real(plntcnt), 0.0, frtc(1), frtc(3), frtc(2))
! ....... Effect of water limitation
! ....... Increase fraction of C to roots due to water stress, cak - 09/12/03
!          h2oeff = line(h2ogef, 0.0, frtc(1), 1.0, frtc(2))
          h2oeff = line(h2ogef, 0.0, frtc(4), 1.0, 0.0)
! ....... Effect of nutrient limitation
! ....... Amount to increase fraction of C going to roots due to
! ....... nutrient stress, cak - 09/12/03
          ntreff = 0.0
          do iel = 1, nelem
            ntreff = max(line(a2drat(iel), 0.0, frtc(5), 1.0, 0.0), ntreff)
          end do
! ....... Compute fraction of C to go to roots
! ....... Fraction of C to go to roots is adjusted due to nutrient or
! ....... water stress, cak - 09/12/03
          froota = min(froota + max(h2oeff, ntreff), 1.0)

        else
          write(bffr,*) frtcindx
          call abortrun('Invalid value frtcindx in froota = '//bffr)
        endif

! ..... Change fraction of C to go to fine roots by effect of CO2
        froota = froota * co2crs(CRPSYS)

      else
      endif

      if (froota .lt. 0.0 .or. froota .gt. 1.0) then
        write(*,*) a2drat,h2ogef,systype, frtcindx
        write(bffr,*) froota
        call abortrun('froota: fine root fraction '//trim(bffr)//' out of bounds')
      endif

      return

      contains
        include 'line.f'
      end function froota
