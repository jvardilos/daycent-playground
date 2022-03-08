
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      common/pheno/accumdd, basetemp(2), cgrwdys, clsgres,               &
     &             curgdys, dayhrs, ddbase, decidgrow, fgrwdys,          &
     &             flsgres, furgdys, hrsinc, grnfill, grnfldys,          &
     &              grnhrvt, mnddhrv, mxddhrv, plntkill, soiltavewk,     &
     &             thermunits, tmpgerm, tmpkill, tmplff, tmplfs

! ... accumdd     - flag, set to true when germination occurs in
! ...               growing degree day implementation, is reset to
! ...               false when senensence and/or harvest occurs
! ... basetemp(1) - base temperature for crop growth, growing degree
! ...               days will accumulate only on days when the average
! ...               temperature is greater than the base temperature
! ...               for the crop (degrees C)
! ... basetemp(2) - ceiling on the maximum temperature used to
! ...               accumulate growing degree days (degrees C)
! ... cgrwdys     - number of days of grass/crop growth in the current
! ...               growing season
! ... clsgres     - late season crop/grass growth restriction factor (0.0 - 1.0)
! ... curgdys     - number of days after the start of growing season
! ...               of unrestricted growth of grass/crop components
! ... dayhrs      - number of hours in a day based on julian date
! ... ddbase      - thermal units (growing degree days) required to
! ...               trigger a senescence event for a perennial plant
! ...               (FRTCINDX 3) or a non-grain filling annual
! ...               (FRTCINDX = 4) or anthesis for a grain filling
! ...               annual (FRTCINDX 5 or 6)
! ... decidgrow   - true in deciduous growth period, between leaf out and leaf drop
! ... fgrwdys     - number of days that tree growth has occurred in the
! ...               current growing season
! ... flsgres     - late season forest growth restriction factor (0.0 - 1.0)
! ... furgdys     - number of days after leaf out of unrestricted growth
! ...               of woody tree components for a decidious tree type
! ... grnfill     - true if in grain fill period for a grain filling
! ...               annual (FRTCINDX 5 or 6)
! ... grnfldys    - number of days that grain filling has been occurring
! ... grnhrvt     - true if the grain is ready to harvest
! ... hrsinc      - indicates if the length of days in increasing or decreasing
! ... mnddhrv     - minimum number of degree days from anthesis to
! ...               harvest for grain filling annuals (FRTCINDX 5 or 6)
! ... mxddhrv     - maximum number of degree days from anthesis to
! ...               harvest for grain filling annuals (FRTCINDX 5 or 6)
! ... plntkill    - true if a killing frost has occurred
! ... soiltavewk  - average soil temperature in second DayCent soils.in
! ...               layer over the previous 7 days, degrees C
! ... thermunits  - growing degree day accumulator for thermal units implementation
! ... tmpgerm     - soil temperature required for seed germination, degrees C
! ... tmpkill     - air temperature at which plant growth stops, degrees C
! ... tmplff      - temperature at which deciduous leaf drop will occur, degrees C
! ... tmplfs      - temperature at which deciduous leaf out will occur, degrees C

      integer cgrwdys, curgdys, fgrwdys
      integer furgdys, grnfldys
      real    basetemp, clsgres, dayhrs, ddbase, flsgres, mnddhrv
      real    mxddhrv, soiltavewk, thermunits, tmpgerm
      real    tmpkill, tmplff, tmplfs
      logical accumdd, decidgrow, grnfill, grnhrvt, hrsinc, plntkill

      save /pheno/
