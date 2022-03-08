
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine grem()
      use calflow;

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'zztim.inc'

! ... Simulate removal of crop/grass by fire or grazing for the month.
      ! Fire events in forest and savanna systems will burn the litter layer
      ! as well as the dead fine branches and dead large wood, cak - 08/23/02

! ... Local variables
      integer   iel, lyr, clyr
      real      shremc, shreme(MAXIEL)
      real      litrme(MAXIEL), sdremc, sdreme(MAXIEL)
      real      eret(MAXIEL)
      real      recres(MAXIEL), ciso, friso
      real      accum(ISOS)
      real      esum
      double precision frac_nh4, frac_no3
      character        subname*10

      accum = 0.0

! ... NOTES:
      ! ciso tells what fraction of the C returned is labeled (C14).
      ! Initialize flrem, fdrem based on fire or grazing.
      ! Mod. fire routines, created 'litburn() subroutine  -mse 4-94.
      ! If dofire(cursys)=1 -> CRPSYS burns standing dead and litter.
        ! To eliminate potential uninitialized variable warnings and
        ! IMHO remove some obfuscation, deleted intermediate variables
        ! flrem, fdrem in favor of directly using the input parameters.



      ! Added for local initialization of variables which may not
      ! get initialized during a run. 8-31-90 -rm
      subname = 'grem      '
      ciso = 0.0
      shremc = 0.0
      sdremc = 0.0
        sdreme = 0.0 ! array
        shreme = 0.0 ! array
        litrme = 0.0 ! array

      ! shoot and standing dead removal based on fire/grazing parameters
      if (dofire(cursys)) then
        shremc = flfrem * aglivc
        sdremc = fdfrem(1) * stdedc
      else
        shremc = flgrem * aglivc
        sdremc = fdgrem * stdedc
      endif

      ! Shoots removed
      if (shremc .gt. 0.) then
        shrema = shrema + shremc ! carbon
        call csched(shremc,aglcis(LABELD),aglivc, &
                    aglcis(UNLABL),csrsnk(UNLABL), &
                    aglcis(LABELD),csrsnk(LABELD), &
                    1.0,shrmai)
        ciso = ciso + shremc*aglcis(LABELD)/aglivc

        ! elements
        do iel = 1, nelem
          shreme(iel) = shremc*aglive(iel)/aglivc
          shrmae(iel) = shrmae(iel) + shreme(iel)
          call flow(aglive(iel),esrsnk(iel),time,shreme(iel))
        end do
      endif

      ! Standing dead removed
      if (sdremc .gt. 0.) then
        sdrema = sdrema + sdremc ! carbon

        call csched(sdremc,stdcis(LABELD),stdedc, &
                    stdcis(UNLABL),csrsnk(UNLABL), &
                    stdcis(LABELD),csrsnk(LABELD), &
                    1.0,sdrmai)
        ciso = ciso + (sdremc*stdcis(LABELD)/stdedc)
        ! elements
        do iel = 1, nelem
          sdreme(iel) = sdremc*stdede(iel)/stdedc
          sdrmae(iel) = sdrmae(iel) + sdreme(iel)
          call flow(stdede(iel),esrsnk(iel),time,sdreme(iel))
        end do
      endif

      if (dofire(cursys)) then  ! ... FIRE
        call fire()
      else                      ! ... GRAZE
        call graze()
      endif

      return
      contains
        subroutine fire
          real ::  cgain, egain(MAXIEL)
          real ::  closs, eloss(MAXIEL)
          real ::  cpass, epass(MAXIEL), som2rmc

          cgain = 0.
          egain = 0
          cpass = 0.0
          epass = 0.0

          ! Residue (surface litter) removed by fire       vek 5/26/90
          call litburn(litrme, som2rmc)

          ciso= ciso + fdfrem(2)*(strcis(SRFC,LABELD)+metcis(SRFC,LABELD))

          ! Dead fine branches removed by fire, cak - 01/02
          if (wood1c .gt. 0.) then
            ! carbon
            closs = fdfrem(3) * wood1c
            tcrem = tcrem + closs
            call csched(closs,wd1cis(LABELD),wood1c, &
                        wd1cis(UNLABL),csrsnk(UNLABL), &
                        wd1cis(LABELD),csrsnk(LABELD), &
                        1.0,accum)
            ! elements
            do iel = 1, nelem
              eloss(iel) = closs * (wood1e(iel) / wood1c)
              terem(iel) = terem(iel) + eloss(iel)
              call flow(wood1e(iel),esrsnk(iel),time,eloss(iel))
            end do
          endif

          ! Dead large wood removed by fire, cak - 01/02
          if (wood2c .gt. 0.) then
            ! carbon
            closs = fdfrem(4) * wood2c
            tcrem = tcrem + closs
            call csched(closs,wd2cis(LABELD),wood2c, &
                        wd2cis(UNLABL),csrsnk(UNLABL), &
                        wd2cis(LABELD),csrsnk(LABELD), &
                        1.0,accum)
            ! elements
            do iel = 1, nelem
              eloss(iel) = closs * (wood2e(iel) / wood2c)
              terem(iel) = terem(iel) + eloss(iel)
              call flow(wood2e(iel),esrsnk(iel),time,eloss(iel))
            end do
          endif

          ! Carbon and nutrient return following removal by fire
          !   fret()    - fraction of element returned by fire
          ! The following variables have units g/m**2/month and are:
          !   sheret    - elemental return for shoots
          !   sderet    - elemental return for standing dead and litter
          !   eret(iel) - total elemental return for aboveground removal

          ! Return carbon from burning live shoots by the fire as charcoal
          ! to the passive SOM pool
          if (aglivc .gt. 0.) then
            cgain = flfrem * fret(1,1) * aglivc
            cpass = cpass + cgain
            call csched(cgain,aglcis(LABELD),aglivc, &
                        csrsnk(UNLABL),som3ci(UNLABL), &
                        csrsnk(LABELD),som3ci(LABELD), &
                        1.0,accum)
          endif
          ! Return carbon removed from standing dead by the fire as charcoal
          ! to the passive SOM pool
          if (stdedc .gt. 0.) then
            cgain = fdfrem(1) * fret(1,1) * stdedc
            cpass = cpass + cgain
            call csched(cgain,stdcis(LABELD),stdedc, &
                        csrsnk(UNLABL),som3ci(UNLABL), &
                        csrsnk(LABELD),som3ci(LABELD), &
                        1.0,accum)
          endif
          ! Return carbon from burning the structural component of surface
          ! litter by the fire as charcoal to the passive SOM pool
          if (strucc(SRFC) .gt. 0.) then
            cgain = fdfrem(2) * fret(1,1) * strucc(SRFC)
            cpass = cpass + cgain
            call csched(cgain,strcis(SRFC,LABELD),strucc(SRFC), &
                        csrsnk(UNLABL),som3ci(UNLABL), &
                        csrsnk(LABELD),som3ci(LABELD), &
                        1.0,accum)
          endif
          ! Return carbon from burning the metabolic component of surface
          ! litter by the fire as charcoal to the passive SOM pool
          if (metabc(SRFC) .gt. 0.) then
            cgain = fdfrem(2) * fret(1,1) * metabc(SRFC)
            cpass = cpass + cgain
            call csched(cgain,metcis(SRFC,LABELD),metabc(SRFC), &
                        csrsnk(UNLABL),som3ci(UNLABL), &
                        csrsnk(LABELD),som3ci(LABELD), &
                        1.0,accum)
          endif
          ! Return carbon from burning the surface component of active soil
          ! organic matter by the fire as charcoal to the passive SOM pool
          if (som1c(SRFC) .gt. 0.) then
            cgain = fdfrem(2) * fret(1,1) * som1c(SRFC)
            cpass = cpass + cgain
            call csched(cgain,som1ci(SRFC,LABELD),som1c(SRFC), &
                        csrsnk(UNLABL),som3ci(UNLABL), &
                        csrsnk(LABELD),som3ci(LABELD), &
                        1.0,accum)
          endif
          ! Return carbon from burning the surface component of intermediate
          ! soil organic matter by the fire as charcoal to the passive SOM pool
          if (som2c(SRFC) .gt. 0.) then
            cgain = fdfrem(2) * fret(1,1) * som2c(SRFC)
            cpass = cpass + cgain
            call csched(cgain,som2ci(SRFC,LABELD),som2c(SRFC), &
                        csrsnk(UNLABL),som3ci(UNLABL), &
                        csrsnk(LABELD),som3ci(LABELD), &
                        1.0,accum)
          endif
          ! Return carbon from burning the dead fine branches by the fire
          ! as charcoal to the passive SOM pool
          if (wood1c .gt. 0.) then
            cgain = fdfrem(3) * fret(2,1) * wood1c
            cpass = cpass + cgain
            call csched(cgain,wd1cis(LABELD),wood1c, &
                        csrsnk(UNLABL),som3ci(UNLABL), &
                        csrsnk(LABELD),som3ci(LABELD), &
                        1.0,accum)
          endif
          ! Return carbon from burning the dead large wood by the fire
          ! as charcoal to the passive SOM pool
          if (wood2c .gt. 0.) then
            cgain = fdfrem(4) * fret(3,1) * wood2c
            cpass = cpass + cgain
            call csched(cgain,wd2cis(LABELD),wood2c, &
                        csrsnk(UNLABL),som3ci(UNLABL), &
                        csrsnk(LABELD),som3ci(LABELD), &
                        1.0,accum)
          endif

          ! associated elements to passive pool based on max C/E ratio
          epass(1:nelem) = cpass/varat3(1,1:nelem)
          ! Add SOM3 return to erata accumulator
          ereta(1:nelem) = ereta(1:nelem) + epass(1:nelem)
          do iel = 1, nelem
            call flow(esrsnk(iel),som3e(iel),time,epass(iel))
          end do

          frac_nh4 = 0.5D0
          frac_no3 = 0.5D0
          ! Return nutrients
          do iel = 1, nelem
            ! burnt live shoots, standing dead, and surface litter
            ! sheret = fret(1,iel+1) * shreme(iel)
            sdreme(iel) = sdreme(iel) + litrme(iel)
            ! sderet = fret(1,iel+1) * sdreme(iel)
            ! eret(iel) = sheret + sderet
            eret(iel) = fret(1,iel+1) * (shreme(iel) + sdreme(iel))

            ! burnt dead fine branches and dead large wood, cak - 01/02
            egain(iel) = egain(iel) + &
                         fdfrem(3) * fret(2,iel+1) * wood1e(iel) + &
                         fdfrem(4) * fret(3,iel+1) * wood2e(iel)

            esum = eret(iel) + egain(iel) - epass(iel)
            if(esum .gt. 0.) then

              ! Add the mineral pool return to the accumulator
               ereta(iel)= ereta(iel) + esum

              if (iel .eq. N) then
                clyr = 1
                subname = 'grem1     '
                call update_npool(clyr, esum*eret(N)/(eret(N)+egain(N)), &
                     frac_nh4, frac_no3, ammonium, nitrate, subname)

                subname = 'grem2     '
                call update_npool(clyr, esum*egain(N)/(eret(N)+egain(N)), &
                     frac_nh4, frac_no3, ammonium, nitrate, subname)
              endif

              call flow(esrsnk(iel),minerl(1,iel),time,esum)
            endif
          end do

          ! END FIRE

          creta= creta + cpass  ! Accumulate amounts returned
          if(aceqcnt > 0) call EqAclratFire(cpass, epass, som2rmc)
          return

        end subroutine fire


        subroutine graze
        real      cret, feces, urine
          ! ... GRAZE NOTES:
          !   Carbon and nutrient return following removal by grazing.
          !   Grazing return with feces and urine explicitly separated.
          !   All carbon returned by grazing is in the form of feces.
          !     fcret     - fraction of carbon returned (input value gfcret directly KLK 17Aug11)
          !     gret(iel) - fraction of element returned by grazing
          !   The following variables have units g/m**2/month and are:
          !     cret      - amount of carbon returned to system
          !     sheret    - elemental return for shoots
          !     sderet    - elemental return for standing dead and litter
          !     eret(iel) - total elemental return for aboveground removal
          !     urine     - amount of urine returned
          !     feces     - amount of fecal material returned (N, P, S)
          !   To adjust for the changing lignin content of added material
          !   strucc(1) and strlig are recomputed.
          cret = gfcret * (shremc + sdremc)
          if (cret .le. 0.0) then
            cret = 0.0
            eret(1:nelem) = 0.0
          else
            frac_nh4 = 1.0D0
            frac_no3 = 0.0D0
            do iel = 1, nelem
              ! Fraction of N that is returned is a function of clay content, cak - 03/03/02
              if (iel .eq. N) then
                if (clay .lt. 0.0) then
                  gret(iel) = 0.7
                else if (clay .gt. 0.30) then
                  gret(iel) = 0.85
                else
                  gret(iel) = line(clay, 0.0, 0.7, 0.30, 0.85)
                endif
              endif
              ! sheret = gret(iel) * shreme(iel)
              ! sderet = gret(iel) * sdreme(iel)
              ! eret(iel) = sheret + sderet
              eret(iel) = gret(iel) * (shreme(iel) + sdreme(iel))
              tgzrte(iel) = tgzrte(iel) + eret(iel)
              ereta(iel) = ereta(iel) + eret(iel) ! update the return accumulator
              urine= (1.-fecf(iel)) * eret(iel)
              feces= fecf(iel) * eret(iel)
              recres(iel) = feces/cret
              if (iel .eq. N) then
                clyr = 1
                subname = 'grem3     '
                call update_npool(clyr, urine, frac_nh4, frac_no3, &
                                  ammonium, nitrate, subname)
              endif
              call flow(esrsnk(iel),minerl(1,iel),time,urine)

              ! Add the amount of N that is volatilized from excreted animal
              ! waste to the VOLPL and VOLPLA output variables, cak - 03/31/04
              if (iel .eq. N) then
                volpl = volpl + ((1.0 - gret(N)) * shreme(N)) + &
                                ((1.0 - gret(N)) * sdreme(N))
                volpla = volpla + ((1.0 - gret(N)) * shreme(N)) + &
                                  ((1.0 - gret(N)) * sdreme(N))
                volpac = volpac + ((1.0 - gret(N)) * shreme(N)) + &
                                  ((1.0 - gret(N)) * sdreme(N))
              endif
            end do

            ! Mod. to add structural & metabolic C into labeled (numerator)
            ! and total (denominator) C removed.  (vek  05/26/90)
            ! friso tells what fraction of the C returned is labeled
            lyr = 1
            friso = ciso / (shremc + sdremc)
            call partit(cret,recres,lyr,csrsnk,esrsnk,feclig,friso)

            creta= creta + cret  ! Accumulate amounts returned
          endif

        ! END GRAZE
        end subroutine graze

        include 'line.f'
      end subroutine grem
