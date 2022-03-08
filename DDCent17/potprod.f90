
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine potprod(trbasl, tavedly, petdly, tfrac, tavemth, curday, srad)

!    Modifications
!     Sept15: K. Killian
!              moved the tree basal area calculation to where it is used so it
!              responds daily and does not require the tree be defined at the
!              start of the month.

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'timvar.inc'

! ... Argument declarations
      real :: trbasl
      real :: tavedly
      real :: petdly
      real :: tfrac, tavemth
      integer :: curday
      double precision :: srad

! ... Local Variables
      real      :: pp, cancvr
      real      :: wdbmas

! ... Determine potential production if it's during the growth season.
! ... The variable, pp, is used below to recompute aglivb for use in h2olos.
      pp = 0.0

! ... Since we are outputing the water stress term add initialization
! ... for the h2ogef array so that it will output zero during periods
! ... of no growth, cak - 04/28/2006
      h2ogef = 0.0

      ! Added below for savanna model (rm)
      ! moved here from simsom so we can test to make sure the trees are growing
      cancvr = 0.0
      if (cursys .eq. SAVSYS  .and.  forgrw .eq. 1) then
        wdbmas = (fbrchc + rlwodc) * 2.0
        ! Can get a divide by zero error when there is no wood biomass, add a
        ! minimum wood biomass so that trees can grow from nothing, cak - 10/08/02
        if (wdbmas .le. 0.0) wdbmas = 50
        ! Change the way that tree basal area is calculated, cak 12/19/01
        trbasl = (wdbmas/(0.88 * ((wdbmas * 0.01)**0.635)))
        if (trbasl .lt. 250.0) trbasl = trbasl * basfct
        trbasl = wdbmas / trbasl
        cancvr = 1 - exp(-0.064 * trbasl)
        if (trbasl .le. 1.0E-6) then
          trbasl = 0.3
        endif
      endif

! ... For a Crop System...
      if (crpgrw .eq. 1) then
        if ((frtcindx .lt. 3) .or. ((frtcindx .ge. 3) .and. (.not. plntkill))) then
          call potcrp(cancvr, tavedly, petdly, tfrac, srad, curday, daylength(curday))
          pp = pcropc
        endif
      endif

! ... For a Forest System...
      if (forgrw .eq. 1) then
        call potfor(tavedly, petdly, tfrac, tavemth, srad, daylength(curday))
        pp = pforc
      endif

      if (cursys .eq. CRPSYS) then
        aglivb = aglivb + 0.25 * pp * 2.5
      endif

      return
      end
