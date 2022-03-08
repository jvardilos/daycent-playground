
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine sumcar

      implicit none
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'timvar.inc'

      real          :: temp

! ... Sum unlabeled and labeled carbon to get totals.

      strucc(1)=sum(strcis(1,:))
      strucc(2)=sum(strcis(2,:))
      metabc(1)=sum(metcis(1,:))
      metabc(2)=sum(metcis(2,:))

! ... Sum som1 surface and soil isotopes separately.  vek 08-91
      som1c(1)=sum(som1ci(1,:))
      som1c(2)=sum(som1ci(2,:))
      som2c(1)=sum(som2ci(1,:))
      som2c(2)=sum(som2ci(2,:))
      som3c=   sum(som3ci)

      wood1c = sum(wd1cis)
      wood2c = sum(wd2cis)
      wood3c = sum(wd3cis)

      somtc=som1c(2) + som2c(2) + som3c + strucc(2) + metabc(2)
      aglivc=aglcis(1)+aglcis(2)
      stdedc=stdcis(1)+stdcis(2)
      bglivcj=bglcisj(1)+bglcisj(2)
      bglivcm=bglcism(1)+bglcism(2)

      temp = agcacc
      agcacc=agcisa(1)+agcisa(2)
      ! add current growth to monthly above ground crop production
      agcmth(month)  =  agcmth(month)  + (agcacc - temp)

      temp = bgcjacc
      bgcjacc=bgcisja(1)+bgcisja(2)
      ! add current growth to monthly below ground juvenile production
      bgcjmth(month) =  bgcjmth(month) + (bgcjacc - temp)
      ! bgcmth(month) =  bgcmth(month) + (bgcjacc - temp)

      temp = bgcmacc
      bgcmacc=bgcisma(1)+bgcisma(2)
      ! add current growth to monthly below ground mature production
      bgcmmth(month) =  bgcmmth(month) + (bgcmacc - temp)
      ! bgcmth(month) =  bgcmth(month) + (bgcjacc - temp)

      rleavc = sum(rlvcis)
      frootcj = sum(frtcisj)
      frootcm = sum(frtcism)
      fbrchc = sum(fbrcis)
      rlwodc = sum(rlwcis)
      crootc = sum(crtcis)
      frnutc = sum(frncis)

      rlvacc =  sum(alvcis)
      frtjacc = sum(afrcisj)
      frtmacc = sum(afrcism)
      fbracc =  sum(afbcis)
      rlwacc =  sum(alwcis)
      crtacc =  sum(acrcis)
      frnacc =  sum(afncis) ! fruit/nut

      temp   = fcacc ! add current growth to the forest production
      fcacc  = rlvacc + frtjacc + frtmacc + fbracc + rlwacc + crtacc + frnacc
      fcmth(month)   =  fcmth(month) + (fcacc - temp)

      return
      end
