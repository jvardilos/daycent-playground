
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      common/plot1/ agcacc,agcmth,agcprd,agdefac,aglcis,aglivc,           &
     &     aglive,aminrl,amt1c2,amt2c2,anerb,annet,arspmth,as11c2,        &
     &     as12c2,as21c2,as22c2,as3c2,asmos,ast1c2,ast2c2,avh2o,          &
     &     bgcjacc,bgcmacc,bgcjmth,bgcmmth,bgcjprd,bgcmprd,               &
     &     bgdefac,bglcisj,bglcism,bglivcj,bglivcm,bglivej,bglivem,       &
     &     cautoresp,cgrain,cinput,clittr,cgrspflux,cmrspflux,            &
     &     co2cce,co2crs,co2cpr,co2ctr,cproda,cprodc,                     &
     &     cprodf,creta,crmvst,crpstg,crpval,dsomsc,                      &
     &     egrain,elimit,eprodc,eprodf,ermvst,                            &
     &     eupacc,eupaga,eupbga,eupprd,evap,                              &
     &     fautoresp,fertac,fertmth,fertot,fertprd,                       &
     &     fgrspflux,fmrspflux,grspann,grspflux,grspmth,                  &
     &     harmth,hi,irract,irrtot,metabc,metabe,metcis,minerl,           &
     &     mrspann,mrspflux,mrspmth,mt1c2,mt2c2,nfix,nfixac,occlud,       &
     &     parent,pet,petann,plabil,prcann,ptagc,ptbgc,pttr,              &
     &     rain,relyld,resp,respmth,runoff,rwcf,                          &
     &     s11c2,s12c2,s21c2,s22c2,s3c2,satmac,sclosa,scloss,sdrema,      &
     &     secndy,shrema,sirrac,snfxac,snlq,snow,soilnm,somsc,somse,      &
     &     somtc,som1c,som1ci,som1e,som2c,som2ci,som2e,som3c,som3ci,      &
     &     som3e,srspann,srspmth,stdcis,st1c2,st2c2,stemp,strcis,stream,  &
     &     strlig,strmac,strucc,struce,sumnrs,stdedc,stdede,              &
     &     tave,tlittr,tminrl,tnetmn,totc,tran,volgma,volexa,volpla,      &
     &     wd1c2,wd2c2,wd3c2,wdfxaa,wdfxas,                               &
     &     eproda, clitad, elitad ! yearly growth E and litter input

      real agcacc,agcmth(12),agcprd,agdefac,aglcis(2),aglivc
      real aglive(3),aminrl(3),amt1c2,amt2c2,anerb,annet,arspmth(2,2)
      real as11c2,as12c2,as21c2,as22c2,as3c2,asmos(10),ast1c2,ast2c2
      real avh2o(3)
      real bgcjacc,bgcmacc,bgcjmth(12),bgcmmth(12),bgcjprd,bgcmprd
      real bgdefac,bglcisj(2),bglcism(2)
      real bglivcj,bglivcm,bglivej(3),bglivem(3)
      real cautoresp(2),cgrain,cinput,clittr(2,2)
      real cgrspflux(3),cmrspflux(3)
      real co2cce(2,2,3),co2crs(2),co2cpr(2),co2ctr(2),cproda,cprodc
      real cprodf,creta,crmvst,crpstg(3),crpval,dsomsc
      real egrain(3),elimit,eprodc(3),eprodf(3),ermvst(3)
      real eupacc(3),eupaga(3),eupbga(3),eupprd(3),evap
      real fautoresp(2),fertac(2,3),fertmth(12,3),fertot(2,3)
      real fertprd(2,3),fgrspflux(7),fmrspflux(7)
      real grspann(2),grspflux(2), grspmth(2)
      real harmth,hi,irract,irrtot,metabc(2),metabe(2,3),metcis(2,2)
      real minerl(10,3),mrspann(2),mrspflux(2),mrspmth(2)
      real mt1c2(2),mt2c2(2),nfix,nfixac
      real occlud,parent(3),pet,petann,plabil,prcann,ptagc,ptbgc
      real pttr,rain,relyld,resp(2),respmth(2),runoff,rwcf(10)
      real s11c2(2),s12c2(2),s21c2(2),s22c2(2),s3c2(2),satmac,sclosa
      real scloss,sdrema,secndy(3),shrema,sirrac,snfxac(2),snlq,snow
      real soilnm(3),somsc,somse(3),somtc
      real som1c(2),som1ci(2,2),som1e(2,3)
      real som2c(2),som2ci(2,2),som2e(2,3),som3c,som3ci(2),som3e(3)
      real srspann(2),srspmth(2),stdcis(2),st1c2(2),st2c2(2),stemp
      real strcis(2,2),stream(8),strlig(2),strmac(8)
      real strucc(2),struce(2,3),sumnrs(3),stdedc,stdede(3)
      real tave,tlittr(2,2),tminrl(3),tnetmn(3),totc,tran
      real volgma,volexa,volpla
      real wd1c2(2),wd2c2(2),wd3c2(2),wdfxaa,wdfxas

      real eproda(3), clitad(2), elitad(2,3) ! yearly growth and litter input
      save /plot1/
