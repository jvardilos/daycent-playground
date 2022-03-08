
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

! Dec 2013
!    removed tafue in favor of a single sit level Automatic Fertilise use Efficiency
!    Added Century tree fruit/nut variables
!      fnftim      fruit/nut fall time in integer days
!      fngddl      fruit/nut growing degree day limit, tbase, tmax
!      tfhrvfr     Tree fruit harvest fraction


      common/parfs/basfct, btolai, ccefor, cerfor,                       &
     &    decid, decw1, decw2, decw3,                                    &
     &    fcfrac, forgrw,forrtf, fnftim, fngddl, tfhrvfr,                &
     &    klai, laitop, ldrmlt, leafdr,                                  &
     &    maxlai, maxldr, maxnp, sapk, swold,                            &
     &    tfrtcn, tfrtcw, tmxturn,                                       &
     &    wdlig, wmrtfrac, woodb, wooddr, wrdsrfc, twscoef

      integer decid, forgrw, fnftim

      real basfct,btolai, decw1,decw2,decw3
      real fngddl(3), tfhrvfr  ! fruit/nut control variables
      real forrtf(3), klai,laitop,ldrmlt,leafdr(12)
      real maxlai,maxldr,maxnp,sapk,swold,tmxturn
      real tfrtcn(2),tfrtcw(2),twscoef(2)
      real wmrtfrac,woodb,wrdsrfc
      ! real ccefor(2,5,3),cerfor(3,5,3), fcfrac(5,2) ! Why 1 less than FPARTS
      real ccefor(2,6,3), cerfor(3,6,3), fcfrac(6,2)  ! increased dimension for fruit/nut
      real wdlig(7),wooddr(7)

      save /parfs/


! These tree variables from PHENO should probably be here
!   dayhrs     - number of hours in a day based on julian date
!   decidgrow  - true in period between deciduous leaf out and leaf drop
!   hrsinc     - true if the day length is increasing
!   tmplff     - temperature below which deciduous leaf drop may occur, (C) C4.5
!   tmplfs     - temperature above which deciduous leaf out may occur,  (C) C4.5

! Century 4.5 variables added to Century but not used in DayCent
!   trecfrac   - C allocation fraction for new production  C4.5
!   leafout    - indicates deciduous leafout is occurring  C4.5
