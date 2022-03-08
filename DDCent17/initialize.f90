
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine initialize(extndopt)

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

! ... Argument declarations
      integer extndopt

! ... This subroutine initializes several common block variables that were
!     being used in the code without initialization, - cak - 09/26/00
!     05/2014  K. Killian
!       initialize new fruit/nut harvest parameters

! ... Local variables

!...Initialize annual production accumulators
!   NOTE: month, dofrst, dofone and nelem are temporarily set so inprac and
!         annacc clear everything do variables are always reset in schedl
      dofone = .true.   ! reset in schedl
      dofrst = .true.   ! reset in schedl
      doplnt = .true.   ! reset in schedl
! ... Used in schedl subroutine, set to .true. so that crpgrw, msplt, and
! ... forgrw are initialized
      dolast = .true.   ! reset in schedl
      doflst = .true.   ! reset in schedl
      nelem = MAXIEL    ! reset in site file
      call inprac(CRPSYS)
      call inprac(FORSYS)

      !swold = 0                                                      ! parfs

!...Initialize monthly co2 accumlators
!...These plot1 variables need to initialized before calling SAVARP
        st1c2 = 0.0
        st2c2 = 0.0
        mt1c2 = 0.0
        mt2c2 = 0.0
        s11c2 = 0.0
        s12c2 = 0.0
        s21c2 = 0.0
        s22c2 = 0.0
        s3c2  = 0.0

! zero the plot3 forest retranslocation and harvest variables
        ! Fruit/Nut variables
        frnutc = 0.0  ! carbon
        frnute = 0.0  ! elements
        frnacc = 0.0  ! carbon growing season C production accumulator
        frncis = 0.0  ! unlabeled/labeled C
        afncis = 0.0  ! unlabeled/labeled C growing season accumulator
        fndday = 0.0  ! growing degree day accumulator
        fnftrm = 0.0  ! post growth fall time remaining
        fnhrvc = 0.0  ! harvested C
        fnhrve = 0.0  ! harvested E
        fnhrve = 0.0                                          ! plot3
        forstg = 0.0

! zero the universe variables in plot2
        csrsnk  = 0.0
        esrsnk  = 0.0

! ... The co2*(*) and cltfac(*) parameters from the plot1 and plot2
! ... common blocks should not be inititalized on an extend
      if (extndopt .eq. 0) then
! ..... These variables are set only if the proper conditions are met,
! ..... initialize with default values
        cltfac = 1.0  ! array

! ..... Initialize these variables for first print of plot commons
          co2crs = 1.0   ! array
          co2cpr = 1.0   ! array
          co2ctr = 1.0   ! array
          co2cce = 1.0   ! array
      endif

      !=======================================================
      call defsitpar()

      return
      end


   subroutine defsitpar()
     USE ISO_C_Binding
     include 'const.inc'
     include 'sitsoil.inc'
     !=======================================================
!     default values for new N2O and methane sitpar parameters
!     Make a stand alone subroutine so sitarchive can call it
!
!     these were not addressed in the 2013 inventory version
!     Need to maintain this until all the 2013 site archives are updated
!     CO2_to_CH4 is a local constant in methane_production and removed from input 10Mar2013

     N2Oadjust_fc     =   0.012  ! previous constant; use as flag
     N2Oadjust_wp     =   0.0    ! set to zero as a flag probably will change
     MaxNitAmt        =   1.0;   ! maximum daily nitrification amount (gN/m^2)
     SnowFlag         =   1;     ! snow insultes soil
     netmn_to_no3     =   0.2;   ! fraction of new net mineralization going to NO3
     wfpsdnitadj      =   1.0;   ! adjust inflection point for water filled pore space
     N2N2Oadj         =   1.0;   ! multiplier on N2/N2O ratio
     ammonium = 0.0
     nitrate = 0.0
   end subroutine defsitpar
