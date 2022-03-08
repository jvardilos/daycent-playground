
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

!    removed cafue in favor of a single site level Automatic Fertilise use Efficiency


      common/parcp/aglivb,astrec,aglrem,astgc,astlbl,astlig,             &
     &    auirri,awhc,basfc2,bglrem,bioflg,biok5,biomax,                 &
     &    cfrtcn,cfrtcw,clteff,cmxturn,crpgrw,crprtf,cultra,             &
     &    efrgrn,eMax,fallrt,fawhc,fdfrem,fdgrem,fecf,feclig,            &
     &    flfrem,flghrv,flgrem,fligni,fnue,fret,frtc,                    &
     &    frtcindx,frtsh,fsdeth,fulcan,                                  &
     &    gfcret,gret,grwprc,grzeff,hibg,himax,himon,hiwsf,              &
     &    irramt,mrtfrac,omadtyp,                                        &
     &    pltmrf,pramn,pramx,prbmn,prbmx,                                &
     &    rdrj,rdrm,rdsrfc,rmvstr,rtdtmp,                                &
     &    sdethc,seedl,sfclit,stdead,vlossp,cwscoef

      integer auirri,bioflg,crpgrw,flghrv,grzeff,seedl,frtcindx
      integer himon(2)

      real aglivb,astrec(3),aglrem,astgc,astlbl
      real astlig,awhc,basfc2,bglrem,biok5,biomax
      real cfrtcn(2),cfrtcw(2),clteff(4),cmxturn
      real crprtf(3),cultra(7),cwscoef(2),efrgrn(3),eMax
      real fallrt,fawhc,fdfrem(4),fdgrem,fecf(3),feclig
      real flfrem,flgrem,fligni(2,3),fnue(2),fret(3,4),frtc(5),frtsh
      real fsdeth(4),fulcan,gfcret,gret(3),grwprc,hibg,himax,hiwsf
      real irramt, mrtfrac, omadtyp
      real pltmrf,pramn(3,2),pramx(3,2),prbmn(3,2),prbmx(3,2)
      real rdrj,rdrm,rdsrfc,rmvstr,rtdtmp
      real sdethc,sfclit,stdead,vlossp

      save /parcp/
