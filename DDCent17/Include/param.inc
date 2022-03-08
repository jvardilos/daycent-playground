
!               Copyright 1993 Colorado State University
!                       All Rights Reserved

!   afue        made the Automatic Fertilizer Use Efficiency a site level
!               parameter to help prevent orchard confusion.
!   undodran    undo values for DRAIN Added here since DRAIN
!               add another value (3) for the FLOD event
!   watrflag    0 = rain and irrigation events used to add water to the
!                   soil potentially bringing it to saturation
!               1 = water added to the system automatically to bring
!                   the soil water content to saturation

      common/param/afiel,afue,amov,autoresp1,autoresp2,awilt,basef,       &
     &    bulkd,co2ipr,co2ice,co2irs,co2itr,co2sys,co2tm,                 &
     &    drain,epnfa,epnfs,falprc,fracro,                                &
     &    hpttr,htran,ivauto,labtyp,labyr,                                &
     &    ckmrspmx,cmrspnpp,fkmrspmx,fmrsplai, cgresp,fgresp,             &
     &    maxtmp,mctemp,micosm,mintmp,                                    &
     &    nelem,newautoresp,Ninput,nlayer,nlaypg,nlaypgst,no3pref,        &
     &    npp2cs,Nstart,OMADinput,OMADstart,                              &
     &    ph,pHscalar,phstart,phsys,phtm,                                 &
     &    ppdf,prcskw,prcstd,prdx,precip,precipyr,precro,psloss,pslsrb,   &
     &    rcelit,rces1,rces2,rces3,remwsd,rock,                           &
     &    satmos,satmt,sirri,snfxmx,sorpmx,stamt,stormf,                  &
     &    strm5l,strm5u,ststart,stsys,swflag,tmxbio,trbasl,               &
     &    claypg,claypg_const,tlaypg, cmix,tmix, undodran, watrflag

      integer:: falprc,ivauto,labtyp,labyr,micosm,nelem
      integer:: Ninput,nlayer,nlaypg,nlaypgst,Nstart,OMADinput,OMADstart
      integer:: phsys,phtm,swflag,claypg,claypg_const,tlaypg,watrflag

      real   :: afiel(10),amov(10), awilt(10)
      real   :: autoresp1(2),autoresp2(2)
      real   :: afue, basef,bulkd
      real   :: co2ipr(2),co2ice(2,2,3),co2irs(2),co2itr(2),co2sys
      real   :: co2tm(2), drain,epnfa(2),epnfs(2),fracro
      real   :: hpttr(12),htran(12),ckmrspmx(3), cgresp(3)
      real   :: cmrspnpp(7), fkmrspmx(7),fmrsplai(7),fgresp(7) ! should be FPARTS (?)
      real   :: maxtmp(12),mctemp,mintmp(12)
      real   :: newautoresp(2),no3pref(2),npp2cs(2),ph,pHscalar(12)
      real   :: phstart,ppdf(4,2),prcskw(12),prcstd(12),prdx(2)
      real   :: precip(12),precipyr,precro,psloss,pslsrb
      real   :: rcelit(2,3),rces1(2,3),rces2(2,3),rces3(3),remwsd,rock
      real   :: satmos(2),satmt,sirri,snfxmx(2),sorpmx,stamt,stormf
      real   :: strm5l,strm5u,ststart,stsys,tmxbio,trbasl
      real   :: cmix,tmix
      real   :: undodran(3)

      save /param/
