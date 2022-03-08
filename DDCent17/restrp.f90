
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... RESTRP.F

      subroutine restrp(elimit, nelem, availm, favail, cerat, nparts, &
                        cfrac, potenc, rimpct, storage, snfxmx,       &
                        cprodl, eprodl, uptake, eup, plantNfix, relyld)

      implicit none
      include 'const.inc'
      include 'fertil.inc'
      include 'evivars.inc'

      ! Argument declarations
      integer   nelem, nparts
      real      favail(MAXIEL) ! Must be calles with the system specific vector!
      real      availm(MAXIEL), cerat(2,nparts,MAXIEL)
      real      cfrac(nparts), cprodl, elimit, eprodl(MAXIEL)
      real      plantNfix, potenc, relyld, rimpct, snfxmx
      real      storage(MAXIEL), uptake(4,MAXIEL)
      real ::   eup(nparts, MAXIEL) ! Added as a parameter   KLK May 2015

! ... Restrict the actual production based on C/E ratios.  Calculate
! ... minimum, and maximum whole plant nutrient concentrations.

!     Removed potent.inc and added eup to argument list KLK May 2015

! ... Local variables
! ... NOTE:  Local variables cannot have adjustable array size.  MINECI
! ...        and MAXECI are set to the largest array size which may occur.
      ! changed error messages to use message system KLK 28Nov12

      integer   iel, ipart
      real      afert(MAXIEL), eavail(MAXIEL)
      real      maxeci(FPARTS-1,MAXIEL), mineci(FPARTS-1,MAXIEL)
      character (len=40) :: string

! ... Definitions of Local variables
! ...   afert  - amount of automatic fertilizer to be added
! ...   eavail - E available
! ...   maxeci - maximum E/C ratio
! ...   mineci - minimum E/C ratio

      if ((nelem .le. 0) .or. (nelem .gt. 3)) then
        write(string,*) nelem
        call abortrun('nelem ('//trim(string)//') out of bounds in restrp'//        &
          'Check <site>.100')
      endif
      if (nparts .gt. FPARTS-1) then
        write(string,*) nparts
        call abortrun('nparts ('//trim(string)//') out of bounds in restrp')
      endif

! ... Reset variables to zero
      cprodl = 0.0
      eprodl(1:nelem) = 0.0
      afert(1:nelem) = 0.0
      eup(1:nparts, :) = 0.0

! ... There is no production if one of the mineral elements is not available.
      if(aufert .eq. 0) then ! aufert will provide something so don't test further
        ! quit growth if there is no available N and no symbiotic N fixation
        if (availm(N) + storage(N) .le. 1E-4 .and. snfxmx .eq. 0.0) goto 999

        ! check availability of other elements
        if(nelem .gt. N) then
          if(minval(availm(2:nelem) + storage(2:nelem)) .lt. 1.0E-4) goto 999
        end if
      end if

! ... Initialize cprodl
      cprodl = potenc

! ... Calculate soil available nutrients, based on a maximum fraction
! ... (favail) and the impact of root biomass (rimpct), adding storage.
      eavail(1:nelem) = max(0.0, (availm(1:nelem) * favail(1:nelem) *    &
                                 rimpct) + storage(1:nelem))

! ... Calculate E/C of plant parts (crop, grass, or tree)
      do iel = 1, nelem
        do ipart = 1, nparts
          if (cerat(IMAX,ipart,iel) .eq. 0.0) then
           write(string,'(i1,",",i1)') ipart,iel
           call abortrun('cerat(IMAX,'//trim(string)//') = 0.0 in restrp')
          endif
          if (cerat(IMIN,ipart,iel) .eq. 0.0) then
           write(string,'(i1,",",i1)') ipart,iel
           call abortrun('cerat(IMIN,'//trim(string)//') = 0.0 in restrp')
          endif
          mineci(ipart,iel) = 1.0 / cerat(IMAX,ipart,iel)
          maxeci(ipart,iel) = 1.0 / cerat(IMIN,ipart,iel)
        end do
      end do

! ... Compute the limitation
      call nutrlm(elimit, nelem, nparts, cfrac, eavail,                  &
                  maxeci, mineci, snfxmx, eup, cprodl,                   &
                  eprodl, plantNfix, afert, aufert, eviFlag)

! ... Calculate relative yield - reimplemented by AKM 18/8/2000
      relyld = cprodl/potenc

! ... Calculate uptakes from all sources: storage, soil, plantNfix, and autofert
      uptake(ENFIX,N) = plantNfix ! fixation is calculated in nutrlm
      uptake(EFERT,:) = afert     ! Uptake from automatic fertilization

      do iel = 1, nelem
! ..... If storage pool contains all needed for uptake
        if(storage(iel) .le. 0.0) then
          uptake(ESTOR,iel) = 0.0
          uptake(ESOIL,iel) = eprodl(iel) - uptake(ENFIX,iel) - uptake(EFERT,iel)
        elseif(storage(iel) .gt. eprodl(iel)) then
          uptake(ESTOR,iel) = eprodl(iel)
          uptake(ESOIL,iel) = 0.0
! ....... Otherwise, extra necessary from the soil pool
        else
          uptake(ESTOR,iel) = storage(iel)
! ......... subtract plantNfix -mdh 3/8/99
! ......... soil N uptake should not include monthly symbiotic N fixation amount
          uptake(ESOIL,iel) = eprodl(iel) - uptake(ENFIX,iel) -          &
                              uptake(EFERT,iel) - uptake(ESTOR,iel)
        endif
      end do

999   continue

      return
      end
