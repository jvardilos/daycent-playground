
      !           Copyright 1993 Colorado State University
      !                   All Rights Reserved


    subroutine erosn(psloss,bulkd,edepth,enrich,lhzci,lhze,nelem,aceqcnt)
      use calflow;

      implicit none
      include 'const.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

      ! Argument declarations
      integer   nelem, aceqcnt
      real      psloss, bulkd, edepth, enrich, lhzci(3,2), lhze(3,3)

      ! Soil is removed from the system.

      ! psloss is the rate of soil loss (kg/m**2/m) per day -mdh 1/97
      ! bulkd is the bulk density of the soil (kg/liter)
      ! edepth is the depth of the soil (m)
      ! enrich is the enrichment factor for SOM losses
      ! scloss is the total carbon loss from soil organic
      !   matter and below ground litter for the day
      ! sclosa is the accumulator for scloss
      ! lhzci(pool,iso) is the carbon in the lower horizon pools
      ! lhze(pool,iel) is the organic N,P,S in the lower horizon pools
      ! lhzcac is the accumulator for carbon inputs from the
      !   lower horizon pools
      ! lhzeac(iel) is the accumulator for N,P,S inputs from the
      !   lower horizon pools

      ! Local variables
      integer :: iel, iso
      real    :: flost, input
      real    :: closs

      real, dimension(2,MAXIEL) :: eins2
      real, dimension(2)        :: cins2, couts2
      real, dimension(MAXIEL)   :: eins3
      real                      :: cins3, couts3

      cins3  = 0.0;    eins3  = 0.0;    couts3  = 0.0
      cins2  = 0.0;    eins2  = 0.0;    couts2  = 0.0

      ! Compute fraction of soil which is lost.
      flost = psloss/(1000.*bulkd*edepth) * enrich

      ! Total carbon loss
      scloss = somtc * flost
      sclosa = sclosa + scloss

      ! Soil losses from below ground som1, som2, som3, below ground
      ! metabolic and structural pools
      ! Argument after nelem tells how many layers of som1, som2, or
      ! som3 are being modeled.  Only the soil layer is used in soilos.
      ! vek 08/91
      call soilos(SOIL,flost,som1c(SOIL),som1ci,csrsnk,som1e,esrsnk, closs)
      call soilos(SOIL,flost,som2c(SOIL),som2ci,csrsnk,som2e,esrsnk, couts2(SOIL))
      call soilos(1,flost,som3c,som3ci,csrsnk,som3e,esrsnk, couts3)
      call soilos(SOIL,flost,metabc(SOIL),metcis,csrsnk,metabe,esrsnk, closs)
      call soilos(SOIL,flost,strucc(SOIL),strcis,csrsnk,struce,esrsnk, closs)

      ! This section commented out. It is assumed that losses from mineral
      ! pools are replaced by an equivalent amount from the next layer
      ! Calculate soil losses from mineral pools based on edepth
      !  real eloss, sum
      !  flost = psloss/(1000.*bulkd*edepth)
      !  do iel = 1, nelem
      !    eloss = parent(iel)*flost
      !    call flow(parent(iel),esrsnk(iel),time,eloss)
      !    eloss = secndy(iel)*flost
      !    call flow(secndy(iel),esrsnk(iel),time,eloss)
      !    if (iel .eq. P) then
      !      eloss = occlud*flost
      !      call flow(occlud,esrsnk(iel),time,eloss)
      !    endif
      !  enddo

      ! Calculate soil losses from mineral pools based on adep(1)
      ! Compute fraction of soil which is lost.
      !  flost = psloss/(1000.*bulkd*adep1)

      !  do iel = 1, nelem
      !    eloss = minerl(SRFC,iel)*flost
      !    call flow(minerl(SRFC,iel),esrsnk(iel),time,eloss)
      !  enddo

      ! Calculate input of organic matter from next soil horizon
      ! using an equivalent depth of soil to that eroded
      do iso = UNLABL, LABELD
        input = flost*lhzci(1,iso)
        lhzci(1,iso) = lhzci(1,iso) - input
        lhzcac = lhzcac + input
        call flow(csrsnk(iso),som1ci(SOIL,iso),time,input)
        input = flost*lhzci(2,iso)
        lhzci(2,iso) = lhzci(2,iso) - input
        lhzcac = lhzcac + input
        cins2(SOIL) = cins2(SOIL)+input
        call flow(csrsnk(iso),som2ci(SOIL,iso),time,input)
        input = flost*lhzci(3,iso)
        lhzci(3,iso) = lhzci(3,iso) - input
        lhzcac = lhzcac + input
        cins3 = cins3+input
        call flow(csrsnk(iso),som3ci(iso),time,input)
      end do

      do iel = 1, nelem
        input = flost*lhze(1,iel)
        lhze(1,iel) = lhze(1,iel) - input
        lhzeac(iel) = lhzeac(iel) + input
        call flow(esrsnk(iel),som1e(SOIL,iel),time,input)
        input = flost*lhze(2,iel)
        lhze(2,iel) = lhze(2,iel) - input
        lhzeac(iel) = lhzeac(iel) + input
        eins2(SOIL,iel) = eins2(SOIL,iel)+input
        call flow(esrsnk(iel),som2e(SOIL,iel),time,input)
        input = flost*lhze(3,iel)
        lhze(3,iel) = lhze(3,iel) - input
        lhzeac(iel) = lhzeac(iel) + input
        eins3(iel) = eins3(iel)   + input
        call flow(esrsnk(iel),som3e(iel),time,input)
      end do

      if(aceqcnt > 0) call EqAclratCE(cins3, eins3, couts3, cins2, eins2, couts2)
      return
    contains


      subroutine soilos(nlr,flost,somc,somci,csrsnk,some,esrsnk, closs)

  ! ... Argument declarations
        integer :: nlr
        real    :: flost, somc, somci(nlr,ISOS), csrsnk(ISOS)
        real    :: some(nlr,MAXIEL), esrsnk(MAXIEL)
        real    :: closs

  ! ... Compute soil loss for som1, som2, or som3.
  ! ... Only soil layer is considered here, so nlr, 2 for som1, som2, can be the index.
  !     som3 has no surface component so it MUST be 1 or memory corruption will occur!
  ! ...   vek  08-91
        ! klk 6 Apr 2016
        ! contain this routine in erosn, the only place it is called
        ! added the C and E loss output for equilibrium acceleration
        ! simplified the eloss. no need to multiply and divide by C amounts

! ... Local variables
        integer :: iel
        real    :: accum(ISOS)

        accum = 0.0

        if (somc .gt. 0.0001) then ! skip if there is no SOM to loose
          closs = somc * flost  ! Loss of carbon isotopes
          call csched(somc * flost, somci(nlr,LABELD),somc, &
                      somci(nlr,UNLABL),csrsnk(UNLABL),     &
                      somci(nlr,LABELD),csrsnk(LABELD),     &
                      1.0,accum)

          ! Loss for each other element is based on element/carbon ratio
          do iel = 1, nelem
            call flow(some(nlr,iel),esrsnk(iel),time, (flost * some(nlr,iel)))
          end do
        endif

        closs = sum(accum)

        return
      end subroutine soilos
    end subroutine erosn
