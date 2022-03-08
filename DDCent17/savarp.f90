
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine savarp

      implicit none

! ... Compute variables for printing or plotting
      ! converted to F90 free form, uses F90 vector arithmetic, and sum
      ! converted del13out and del14out to a single routine removing the duplicate calls
      ! added the fruit/nut pool to output sums   KLK May 2014

      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

! ... Function declarations
      real      :: delLCout, fsfunc
      external  :: delLCout, fsfunc

! ... Local variables
      integer  :: iel, lyr, mm, nll
      real     :: fsol
      real, parameter :: toler = 1./(2.**52)

! ... Calculate total non-living C, minimum over the year
      totc = min(totc, sum(som1c) + sum(som2c) + som3c + sum(strucc) +   &
                 sum(metabc))

! ... Compute soil organic matter sums for plotting
      somsc = som1c(SOIL) + som2c(SOIL) + som3c
      somtc = somsc + strucc(SOIL) + metabc(SOIL)
      woodc = wood1c + wood2c + wood3c
      frstc = rleavc + frootcj + frootcm + fbrchc + rlwodc + crootc + frnutc
      fsysc = somtc + woodc + frstc + strucc(SRFC) + metabc(SRFC) + &
              som1c(SRFC) + som2c(SRFC)
      totsysc = fsysc + aglivc + bglivcj + bglivcm + stdedc

      somse(1:nelem) = som1e(SOIL,1:nelem) + som2e(SOIL,1:nelem) + som3e(1:nelem)
      somte(1:nelem) = somse(1:nelem) + struce(SOIL,1:nelem) + metabe(SOIL,1:nelem)
      woode(1:nelem) = wood1e(1:nelem) + wood2e(1:nelem) + wood3e(1:nelem)
      frste(1:nelem) = rleave(1:nelem) + frootej(1:nelem) + frootem(1:nelem) + &
           fbrche(1:nelem) + rlwode(1:nelem) + croote(1:nelem) + frnute(1:nelem)
      fsyse(1:nelem) = somte(1:nelem) + woode(1:nelem) + frste(1:nelem) + &
                   struce(SRFC,1:nelem) + metabe(SRFC,1:nelem) + &
                   som1e(SRFC,1:nelem) + som2e(SRFC,1:nelem)
      totsyse(1:nelem) = fsyse(1:nelem) + aglive(1:nelem) + bglivej(1:nelem) + &
                   bglivem(1:nelem) + stdede(1:nelem)

! ... Compute soil organic matter sums by isotope
      somsci = som1ci(SOIL,:) + som2ci(SOIL,:) + som3ci
      somtci = somsci + strcis(SOIL,:) + metcis(SOIL,:)

! ..... Add litter layer components, including surface som1 and surface
! ..... som2, to get total organic matter including residue.  vek 08-91
      tomres = somtci + strcis(SRFC,:) + metcis(SRFC,:) + &
               som1ci(SRFC,:) + som2ci(SRFC,:)

! ... Sum the co2 losses back into the appropriate source/sink, klk 1/05
      csrsnk = csrsnk + mt1c2 + mt2c2 + st1c2 + st2c2 + s11c2 + s12c2 + &
                      s21c2 + s22c2 + s3c2 + wd1c2 + wd2c2 + wd3c2

! ... Sum all state variables
! ... Include som1c(SRFC) and som2c(SRFC) since it is not included in
! ... somtc.   vek 08-91
! ... Remove strm5u and strm51 from this calculation as these variables
! ... are being added to the csrsnk(UNLABL) and csrsnk(LABELD)
! ... accumulators respectively in simsom. -mdh 7/19/01
      totalc = somtc + strucc(SRFC) + metabc(SRFC) + som1c(SRFC) + som2c(SRFC) + &
               aglivc + stdedc + bglivcj + bglivcm + sum(csrsnk) + woodc + frstc

! ... Calculate tminrl
      plabil = 0.0
      do 50 iel = 1, nelem
        tminrl(iel) = 0.0
        do lyr = 1, nlayer
          if (minerl(lyr,iel).gt.0.0) then
            tminrl(iel) = tminrl(iel) + minerl(lyr,iel)
            if (iel .eq. P) then
              fsol = fsfunc(minerl(lyr,P), pslsrb, sorpmx)
              plabil = plabil + minerl(lyr,P) * fsol
            endif
          endif
        end do

!...    To prevent underflow zero any layer with negligable nutrient
!****    (less than can be summed in an IEEE double precision)
        WHERE (abs(minerl(1:nlayer, iel)) .lt. tminrl(iel)*toler) minerl(1:nlayer,iel) = 0.0
50    continue

      nll = nlayer + 1
      do iel = 1, nelem

! ..... Include som1e(1,iel) and som2e(1,iel) since they are not
! ..... included in somte.  vek 08-91
! ..... Remove stream(iel+5) from this calculation as stream(iel+5) is
! ..... being added to the esrsnk(iel) accumulator in simsom. -mdh 7/19/01
        totale(iel) = tminrl(iel) + somte(iel) + struce(SRFC,iel) + &
                      metabe(SRFC,iel) + som1e(SRFC,iel) + &
                      som2e(SRFC,iel) + aglive(iel) + stdede(iel) + &
                      bglivej(iel) + bglivem(iel) + esrsnk(iel) + &
                      minerl(nll,iel) + parent(iel) + secndy(iel) + &
                      woode(iel) + frste(iel) + crpstg(iel) + forstg(iel)
        if (iel .eq. P) then
          totale(iel) = totale(iel) + occlud
        endif
      end do

! ... Above and below ground live C/N ratio
      if (aglive(N) .gt. 0.0 .and. bglivej(N) .gt. 0.0 .and. bglivem(N) .gt. 0.0) then
        aglcn = aglivc/aglive(N)
        bglcnj = bglivcj/bglivej(N)
        bglcnm = bglivcm/bglivem(N)
      else
        aglcn = -999.0
        bglcnj = -999.0
        bglcnm = -999.0
      endif

! ... Overall c/n, c/p, and c/s ratios in soil organic matter
      tcerat(1:nelem) = somtc/somte(1:nelem)

! ... Average annual value of agdefac and bgdefac, the decomposition
! ... factors which combine the effects of temperature and moisture
      aagdefac = 0.0
      do 80 mm = 1, MONTHS
! ..... A negative value indicates that agdefac has not yet been
! ..... calculated for month mm
        if (agdefacm(mm) .lt. 0.) then
          goto 90
        endif
        aagdefac = aagdefac + agdefacm(mm)
80    continue
      mm = 13
90    mm = mm-1
      if (mm .gt. 0) then
        aagdefac = aagdefac/float(mm)
      endif
      abgdefac = 0.0
      do 85 mm = 1, MONTHS
! ..... A negative value indicates that bgdefac has not yet been
! ..... calculated for month mm
        if (bgdefacm(mm) .lt. 0.) then
          goto 95
        endif
        abgdefac = abgdefac + bgdefacm(mm)
85    continue
      mm = 13
95    mm = mm-1
      if (mm .gt. 0) then
        abgdefac = abgdefac/float(mm)
      endif

! ... Litter output
      clittr = metcis + strcis
      tlittr = metcis + strcis + som1ci + som2ci

      if (labtyp .eq. 1 .or. labtyp .eq. 2) then
! ..... Delta 13C output (2)   or   Delta 14C output (1)
        dsomsc = delLCout(labtyp,somsci(LABELD), somsci(UNLABL), dsomsc)
        dsomtc = delLCout(labtyp,somtci(LABELD), somtci(UNLABL), dsomtc)
        dsom1c(SRFC) = delLCout(labtyp,som1ci(SRFC,LABELD), &
                                som1ci(SRFC,UNLABL), dsom1c(SRFC))
        dsom1c(SOIL) = delLCout(labtyp,som1ci(SOIL,LABELD), &
                                som1ci(SOIL,UNLABL), dsom1c(SOIL))
        dsom2c(SRFC) = delLCout(labtyp,som2ci(SRFC,LABELD), &
                                som2ci(SRFC,UNLABL), dsom2c(SRFC))
        dsom2c(SOIL) = delLCout(labtyp,som2ci(SOIL,LABELD), &
                                som2ci(SOIL,UNLABL), dsom2c(SOIL))
        dsom3c = delLCout(labtyp,som3ci(LABELD), som3ci(UNLABL), dsom3c)
        dslit = delLCout(labtyp,clittr(SRFC,LABELD) + som1ci(SRFC,LABELD), &
                         clittr(SRFC,UNLABL) + som1ci(SRFC,UNLABL), dslit)
        dblit = delLCout(labtyp,clittr(SOIL,LABELD) + som1ci(SOIL,LABELD), &
                         clittr(SOIL,UNLABL) + som1ci(SOIL,UNLABL), dblit)
        dstruc(SRFC) = delLCout(labtyp,strcis(SRFC,LABELD), &
                                strcis(SRFC,UNLABL), dstruc(SRFC))
        dstruc(SOIL) = delLCout(labtyp,strcis(SOIL,LABELD), &
                                strcis(SOIL,UNLABL), dstruc(SOIL))
        dmetc(SRFC) = delLCout(labtyp,metcis(SRFC,LABELD), &
                               metcis(SRFC,UNLABL), dmetc(SRFC))
        dmetc(SOIL) = delLCout(labtyp,metcis(SOIL,LABELD), &
                               metcis(SOIL,UNLABL), dmetc(SOIL))
        dcarbostg(CRPSYS) = delLCout(labtyp,carbostg(CRPSYS,LABELD), &
                                    carbostg(CRPSYS,UNLABL), dcarbostg(CRPSYS))
        dcarbostg(FORSYS) = delLCout(labtyp,carbostg(FORSYS,LABELD), &
                                    carbostg(FORSYS,UNLABL), dcarbostg(FORSYS))
        dautoresp(CRPSYS) = delLCout(labtyp,arspmth(CRPSYS,LABELD), &
                                     arspmth(CRPSYS,UNLABL), dautoresp(CRPSYS))
        dautoresp(FORSYS) = delLCout(labtyp,arspmth(FORSYS,LABELD), &
                                     arspmth(FORSYS,UNLABL), dautoresp(FORSYS))
        dhetresp = delLCout(labtyp,respmth(LABELD), respmth(UNLABL), dhetresp)
        dsoilresp = delLCout(labtyp,arspmth(CRPSYS,LABELD) + &
                            arspmth(FORSYS,LABELD) + respmth(LABELD), &
                            arspmth(CRPSYS,UNLABL) + arspmth(FORSYS,UNLABL) + &
                            respmth(UNLABL), dsoilresp)
        dbglivc = delLCout(labtyp,bglcisj(LABELD) + bglcism(LABELD), &
                           bglcisj(UNLABL) + bglcism(UNLABL), dbglivc)
        dbglivcj = delLCout(labtyp,bglcisj(LABELD), bglcisj(UNLABL), dbglivcj)
        dbglivcm = delLCout(labtyp,bglcism(LABELD), bglcism(UNLABL), dbglivcm)
        dfrootc = delLCout(labtyp,frtcisj(LABELD) + frtcism(LABELD), &
                           frtcisj(UNLABL) + frtcism(UNLABL), dfrootc)
        dfrootcj = delLCout(labtyp,frtcisj(LABELD), frtcisj(UNLABL), dfrootcj)
        dfrootcm = delLCout(labtyp,frtcism(LABELD), frtcism(UNLABL), dfrootcm)
        deloi = delLCout(labtyp,tlittr(SRFC,LABELD), tlittr(SRFC,UNLABL), deloi)
        deloe = delLCout(labtyp,tlittr(SOIL,LABELD) + som3ci(LABELD), &
                         tlittr(SOIL,UNLABL) + som3ci(UNLABL), deloe)
      endif

      return
      end


      real function delLCout(labtyp, labeled, unlabeled, oldval)

        implicit none
        include 'const.inc'

        ! labeled C output
        ! Combined del13out and del14out into a single function that automatically
        ! uses the correct constant rather than duplicating the calling code
        ! in savarp and simsom
        !this is a pure function.

        ! Argument declarations
        real,    intent(IN) :: labeled, unlabeled, oldval
        integer, intent(IN) :: labtyp

        delLCout = oldval
        if (unlabeled .gt. 0.0) then
          if(labtyp .eq. 1) then
!           delLCout = ((labeled * 1000.0) / (labeled + unlabeled) - 1) * 1000.0
            delLCout = ((labeled/unlabeled)/FRAC_C14 - 1.0) * 1000.0
          else !if(labtyp .eq. 2) then
            delLCout = ((labeled/unlabeled)/PEEDEE   - 1.0) * 1000.0
          endif
        endif

        return
      end
