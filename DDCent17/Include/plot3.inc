
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

!    added fruit/nut output variables

      common/plot3/ acrcis,afbcis,afrcisj,afrcism,alvcis,alwcis,         &
     & crootc,croote,crtacc,crtcis,crtprd,eupprt,                        &
     & fbrchc,fbrche,fbracc,fbrcis,fbrprd,                               &
     & fcacc,fcmth,fcprd,forstg,                                         &
     & frootcj,frootcm,frootej,frootem,                                  &
     & frtjacc,frtmacc,frtcisj,frtcism,frtjprd,frtmprd,                  &
     & frstc,frste,fsysc,fsyse,                                          &
     & rleavc,rleave,rlvacc,rlvcis,rlvprd,                               &
     & rlwodc,rlwode,rlwacc,rlwcis,rlwprd,                               &
     & sumrsp,tcrem,terem,                                               &
     & w1lig,w2lig,w3lig,w1mnr,w2mnr,w3mnr,wd1cis,                       &
     & wd2cis,wd3cis,wood1c,wood2c,wood3c,woodc,wood1e,                  &
     & wood2e,wood3e,woode,                                              &
     & frnutc, frnute, frnacc, frncis, afncis, fndday, fnftrm,           &
     & fnhrvc, fnhrve, frnprd, frstage

      real acrcis(2),afbcis(2),afrcisj(2),afrcism(2)
      real alvcis(2),alwcis(2)
      real crootc,croote(3),crtacc,crtcis(2),crtprd,eupprt(6,3)
      real fbrchc,fbrche(3),fbracc,fbrcis(2),fbrprd
      real fcacc,fcmth(12),fcprd,forstg(3)
      real frootcj,frootcm,frootej(3),frootem(3),frtjacc,frtmacc
      real frtcisj(2),frtcism(2),frtjprd,frtmprd,frstc,frste(3)
      real fsysc,fsyse(3),rleavc,rleave(3),rlvacc,rlvcis(2),rlvprd
      real rlwodc,rlwode(3),rlwacc,rlwcis(2),rlwprd
      real sumrsp,tcrem,terem(3)
      real w1lig,w2lig,w3lig, w1mnr(3),w2mnr(3),w3mnr(3)
      real wd1cis(2),wd2cis(2),wd3cis(2),wood1c,wood2c,wood3c,woodc
      real wood1e(3),wood2e(3),wood3e(3),woode(3)
      real frnutc, frnute(3), frnacc, frncis(2), afncis(2)             ! Fruit/nut
      real fndday, fnftrm, fnhrvc, fnhrve(3), frnprd, frstage           ! Fruit/nut

      save /plot3/

! Fruit/Nut variables
! frnutc    carbon
! frnute(3) elements
! frnacc    carbon growing season C production accumulator
! frnprd    fruit/Nut growing season production
! frncis(2) unlabeled/labeled C
! afncis(2) unlabeled/labeled C growing season accumulator
! fndday    growing degree day accumulator
! fnftrm    post growth fall time remaining
! fnhrvc    harvested C
! fnhrve(3) harvested E
