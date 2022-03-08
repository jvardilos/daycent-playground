
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

!  included the fruit/nut part in the forest system  4/2014  KLK
!  new names for number of layers in Century and soil water models that are
!  consistant with swconst.h  4/2014  KLK

! ... CRPSYS, FORSYS, SAVSYS are aliases for the system
      integer, parameter :: CRPSYS = 1, FORSYS = 2, SAVSYS = 3

      integer, parameter :: MAXIEL = 3  ! the total # of elements available
      integer, parameter :: N = 1,  P = 2,  S = 3  ! Nitrogen, Phosphorus, Sulphur

      integer, parameter :: CMXLYR = 10   ! maximum number of layers In Century soil
      integer, parameter :: MXSWLYR = 21  ! maximum number of layers in soil water model
      integer, parameter :: MONTHS = 12   ! the number of months

! ... SRFC, SOIL are enumerations for surface, soil
      integer, parameter :: SRFC = 1,  SOIL = 2

! ... CPARTS is the number of parts in the grassland/crop system:
      integer, parameter :: CPARTS = 3
      integer, parameter :: ABOVE = 1,   BELOW = 2 ! above and below ground
      integer, parameter :: BELOWJ = 2   !  BELOWJ  juvenile belowground (roots)
      integer, parameter :: BELOWM = 3   !  BELOWM  mature belowground (roots)

! ... FPARTS is the number of parts in the forest system:
      integer, parameter :: FPARTS = 7
      integer, parameter :: LEAF  = 1    !  LEAF    leaf forest part
      integer, parameter :: FROOTJ = 2   !  FROOTJ  juvenile fine root forest part
      integer, parameter :: FBRCH = 3    !  FBRCH   fine branch forest part
      integer, parameter :: LWOOD = 4    !  LWOOD   large wood forest part
      integer, parameter :: CROOT = 5    !  CROOT   coarse root forest part
      integer, parameter :: FRNUT = 6    !  FRNUT   tree fruit/nut part
      integer, parameter :: FROOTM = 7   !  FROOTM  mature fine root forest part
      integer, parameter :: FROOT = 2    ! for arrays that combine fine roots

! ... NEWFOR, OLDFOR are the new and old forests array pointers
      integer, parameter :: NEWFOR = 1, OLDFOR = 2

! ... IMIN, IMAX, IVAL are enumerations for indexing arrays
      integer, parameter :: IMIN = 1,  IMAX = 2,  IVAL = 3

! ... INTCPT, SLOPE are the intercept and slope array pointers
      integer, parameter :: INTCPT = 1,  SLOPE = 2

      integer, parameter :: ISOS = 2    ! number of isotopes (unlabeld, labeled)
      integer, parameter :: UNLABL = 1,  LABELD = 2 !  unlabeled, labeled array indices

! ... ESTOR, ESOIL, ENFIX, EFERT are enumerations used in
!      restrp.f, growth.f, trees.f
      integer, parameter :: ESTOR = 1, ESOIL = 2, ENFIX = 3, EFERT = 4

! ... Constant values
! ... Change value of PEEDEE constant used in del13c computations, cak - 12/19/01
      real, parameter :: PEEDEE=0.01119490 ! no longer 0.0112372

! ... Fraction of the total carbon that is C14 carbon
      real, parameter :: FRAC_C14=0.00000103 ! not 0.00103

! ... dconst.inc (Daily constants)
      integer, parameter :: NDAY=366, NMONTH=12
