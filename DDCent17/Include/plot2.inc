
!               Copyright 1993 Colorado State University
!                       All Rights Reserved



      common /plot2/ aagdefac, abgdefac, accrst, accrste, agcisa,         &
     &    aglcn, bgcisja, bgcisma, bglcnj, bglcnm,                        &
     &    carbostg, cgracc, cisgra, cltfac, csrsnk,                       &
     &    ch4mnox, ch4yrox, ch4mnem, ch4yrem, ch4mnpr, ch4yrpr,           &
     &    dautoresp, dbglivc, dbglivcj, dbglivcm, dblit,                  &
     &    dcarbostg, deloe, deloi, dfrootc, dfrootcj, dfrootcm,           &
     &    dmetc, dhetresp, dslit, dsoilresp, dsom1c, dsom2c,              &
     &    dsom3c,dsomtc,dstruc, egracc,ereta,esrsnk, gromin,              &
     &    lhzcac,lhzeac, metmnr, n2oacc, n2omth, n2oprd,                  &
     &    omadac, omadae, omadmte, omadmth,                               &
     &    omadprd, omadpre, omadtot, omaetot,                             &
     &    prcfal, rnpml1, sdrmae, sdrmai, shrmai, shrmae, somsci,         &
     &    somtci, somte, strmnr, s1mnr, s2mnr, s3mnr,                     &
     &    tcerat, tcnpro, tgzrte, tomres, totalc, totale, totsysc,        &
     &    totsyse, voleac, volex, volgac, volgm, volpac, volpl,           &
     &    wdfx, wdfxa, wdfxma, wdfxms, wdfxs



      real aagdefac, abgdefac, accrst, accrste(3), agcisa(2)
      real aglcn, bgcisja(2), bgcisma(2), bglcnj, bglcnm
      real carbostg(2,2), cgracc, cisgra(2), cltfac(4), csrsnk(2)
      real dautoresp(2), dbglivc, dbglivcj, dbglivcm, dblit
      real ch4mnox, ch4yrox, ch4mnem, ch4yrem, ch4mnpr, ch4yrpr
      real dcarbostg(2), deloe, deloi, dfrootc, dfrootcj, dfrootcm
      real dmetc(2), dhetresp, dslit, dsoilresp, dsom1c(2), dsom2c(2)
      real dsom3c, dsomtc, dstruc(2)
      real egracc(3), ereta(3),esrsnk(3), gromin(3)
      real lhzcac,lhzeac(3), metmnr(2,3), n2oacc,n2omth(12),n2oprd
      real omadac, omadae(3), omadmte(12,3), omadmth(12)
      real omadprd, omadpre(3), omadtot, omaetot(3), prcfal, rnpml1
      real sdrmae(3), sdrmai(2), shrmai(2)
      real shrmae(3), somsci(2), somtci(2), somte(3)
      real strmnr(2,3), s1mnr(2,3), s2mnr(2,3), s3mnr(3)
      real tcerat(3), tcnpro, tgzrte(3), tomres(2), totalc, totale(3)
      real totsysc, totsyse(3)
      real voleac, volex, volgac, volgm, volpac, volpl
      real wdfx, wdfxa, wdfxma, wdfxms, wdfxs

      save /plot2/
