
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine partit(cpart,recres,lyr,cdonor,edonor,frlign,friso)
      use calflow;

      implicit none
      include 'const.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'zztim.inc'

      ! Argument declarations
      integer :: lyr
      real :: cpart, frlign, friso
      real :: recres(MAXIEL), cdonor(ISOS), edonor(MAXIEL)

      ! Partition residue from cdonor and edonor into layer lyr
      ! of structural and metabolic.
      ! cpart is the amount of carbon in the residue.
      ! recres contains the n/c, p/c, and s/c ratios in the residue.
      ! frlign is the fraction of the incoming material which is lignin.
      ! friso is the fraction of cpart which is labeled.

      ! Local variables
      integer   iel, clyr
      real :: accum(ISOS), epart(MAXIEL), dirabs(MAXIEL)
      real :: caddm, cadds, eaddm, eadds
      real :: fligst, frmet, frn, rcetot, rlnres, delin
      real :: dellig, c13c12r, c13frac, c13lig, c13nlig
      real :: c13struc, c12struc, c13met, c12met
      real :: namt
      double precision  :: frac_nh4, frac_no3
      character(len=10) :: subname

      subname = 'partit    '
      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

      if (cpart .lt. 1.e-07) goto 999 ! skip the calculations if there is no C to apportion


      if (friso .lt. 0.0) then
        if (cdonor(2) .le. 0.0) then
          friso = 0.0
        else
          friso = 1.0
        endif
      endif

      ! Compute amount of element in residue.
      epart(1:nelem) = cpart * recres(1:nelem)

      ! Add residue input to the annual accumulators; For Cfarm KLK 06/2018
      clitad(lyr) = clitad(lyr) + cpart
      elitad(lyr,1:nelem) = elitad(lyr,1:nelem) + epart(1:nelem) ! yearly litter additions

      ! Direct absorption of mineral element by residue
      !  (mineral will be transferred to donor compartment then partitioned into
      !   structural and metabolic using flow routines.)

      ! If minerl(SRFC,iel) is negative then dirabs = zero.
      where (minerl(SRFC,1:nelem) .lt. 0.)
        dirabs(1:nelem) = 0.0
      elsewhere
        dirabs(1:nelem)=damr(lyr,1:nelem)*minerl(SRFC,1:nelem)* max(cpart/pabres,1.)
      endwhere


      ! For each mineral element...
      do 10 iel = 1, nelem

        ! If C/E ratio is too low, transfer just enough to make C/E of residue = damrmn
        if (epart(iel)+dirabs(iel) .le. 0.0) then
          rcetot = 0.0
        else
          rcetot = cpart/(epart(iel)+dirabs(iel))
        endif

        if (rcetot .lt. damrmn(iel)) then
          dirabs(iel) = max(cpart/damrmn(iel) - epart(iel), 0.)
        endif

        if (iel .eq. N) then
          namt = -1.0*dirabs(iel)
          clyr = 1
          call cmpnfrac(clyr,ammonium,nitrate,minerl,frac_nh4,frac_no3)
          call update_npool(clyr, namt, frac_nh4, frac_no3, ammonium, nitrate, subname)
        endif
        call flow(minerl(1,iel),edonor(iel),time,dirabs(iel))
10    continue

      ! Partition carbon into structural and metabolic fraction of
      !  residue (including direct absorption) which is nitrogen
      frn = (epart(1)+dirabs(1)) / (cpart*2.5)

      ! Lignin/nitrogen ratio of residue
      rlnres = frlign/frn

      ! Carbon added to metabolic
      !  Compute the fraction of carbon that goes to metabolic.
      frmet = spl(INTCPT)-spl(SLOPE)*rlnres

      ! Make sure the fraction of residue which is lignin isn't
      !  greater than the fraction which goes to structural.  -rm 12/91
      if (frlign .gt. (1.0 - frmet)) then
        frmet = (1.0 - frlign)
      endif

      ! Make sure at least 1% goes to metabolic
      if (frmet .lt. 0.20) then
        frmet = .20
      endif

      ! Compute amounts to flow
      caddm = cpart * frmet
      if (caddm .lt. 0) then
        caddm = 0.0
      endif
      cadds = cpart-caddm

      ! Adjust lignin content of structural.
      !  fligst is the fraction of incoming structural residue
      !  which is lignin; restricting it to a maximum of .8
      fligst = frlign/(cadds/cpart)

      ! Changed allowable maximum from .8 to .6 -rm 5/92
      !  Changed maximum fraction from .6 to 1.0  -lh 1/93
      if (fligst .gt. 1.0) then
        fligst = 1.0
      endif

      ! Determine what type of labeling is to be done
      if (labtyp .ne. 2) then

        ! Carbon added to metabolic
        accum(UNLABL) = 0.0
        accum(LABELD) = 0.0
        call csched(caddm,friso,1.0, &
                    cdonor(UNLABL),metcis(lyr,UNLABL), &
                    cdonor(LABELD),metcis(lyr,LABELD), &
                    1.0,accum)
        cinput = cinput + accum(UNLABL) + accum(LABELD)

        ! Carbon added to structural
        accum(UNLABL) = 0.0
        accum(LABELD) = 0.0
        call csched(cadds,friso,1.0, &
                    cdonor(UNLABL),strcis(lyr,UNLABL), &
                    cdonor(LABELD),strcis(lyr,LABELD), &
                    1.0,accum)
        cinput = cinput + accum(UNLABL) + accum(LABELD)

      else
        if (friso .gt. 0.0 .and. friso .lt. 1.0) then
          ! Calculate ratio of 13C to 12C in incoming material
          c13c12r = 1/(1/friso - 1)

          ! Calculate delta 13C of incoming material
          delin = ((c13c12r/PEEDEE) - 1) * 1000

          ! Calculate delta 13C of lignin
          dellig = delin + dligdf

          ! Calculate ratio of 13C to 12C in lignin
          c13c12r = dellig * PEEDEE * 1.0e-03 + PEEDEE

          ! Calculate fraction of 13C in lignin
          c13frac = 1/(1/c13c12r + 1)
        elseif (friso .eq. 0.0) then
          ! Set the fraction of 13C in lignin to 0
          c13frac = 0.0
        else
          ! (friso .eq. 1.0) so set the fraction of 13C in lignin to 1
          c13frac = 1.0
        endif

        ! Calculate 13C in lignin
        c13lig = c13frac * frlign * cpart

        ! Calculate 13C in non-lignin
        c13nlig = (cpart * friso) - c13lig

        ! Flow 13C to structural
        c13struc = cadds * fligst * c13frac + cadds * (1 - fligst) * &
                   c13nlig / ((1 - frlign) * cpart)
        call flow(cdonor(LABELD),strcis(lyr,LABELD),time,c13struc)
        cinput = cinput + c13struc

        ! Flow 12C to structural
        c12struc = cadds - c13struc
        call flow(cdonor(UNLABL),strcis(lyr,UNLABL),time,c12struc)
        cinput = cinput + c12struc

        ! Flow 13C to metabolic
        c13met = (cpart * friso) - c13struc
        call flow(cdonor(LABELD),metcis(lyr,LABELD),time,c13met)
        cinput = cinput + c13met

        ! Flow 12C to metabolic
        c12met = caddm - c13met
        call flow(cdonor(UNLABL),metcis(lyr,UNLABL),time,c12met)
        cinput = cinput + c12met
      endif

      call adjlig(strucc(lyr),fligst,cadds,strlig(lyr)) ! Adjust lignin

      ! Partition mineral elements into structural and metabolic
      do 20 iel = 1, nelem

        ! Flow into structural
        eadds = cadds/rcestr(iel)
        call flow(edonor(iel),struce(lyr,iel),time,eadds)
        ! Flow into metabolic
        eaddm = epart(iel)+dirabs(iel)-eadds
        call flow(edonor(iel),metabe(lyr,iel),time,eaddm)
20    continue

999   continue

      return
      end
