
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function leafa(rleavc_opt, rleavc, cprodfLeft, cprodf)

      implicit none
      include 'parfs.inc'

c ... Argument declarations
      real      rleavc, cprodfLeft, cprodf

c ... Compute the fraction of production going to leaves in crops
c ... and woody plants based on water and nutrient availability.
c ...   cprodf     - total potential C production for all tree parts
c ...                (gC/m^2)
c ...   cprodfLeft - amount of carbon still available for allocation (gC/m^2)
c ...   rleavc     - tree leaf carbon (gC/m^2)
!       rleavc_opt - optimum LAI C  calculated in treegrow so don't redo KLK 28 May

c ... Local Variables
      real rleavc_opt
      character bffr*8

      if (rleavc .lt. 0.0) then
        write(bffr,*) 'rleavc = ', rleavc
        call abortrun('leaf C rleavc = '//bffr//' < 0.0 in leafa')
      endif
      if(cprodfLeft .lt. 0.0) call abortrun('cprodfLeft < 0.0 in leafa')
      if (cprodf .le. 0.0) call abortrun('cprodf <= 0.0 in leafa')

c ... Calculate theoretical maximum for LAI based on large wood biomass,
c ... cak - 07/24/02

      leafa = 0.0 ! initialize leafa as no growth
      if (rleavc_opt .gt. rleavc) then
c ..... if possible, increase leaf biomass so that optimum is reached
        leafa = max(0.0,min(rleavc_opt-rleavc, cprodfLeft) / cprodf)

        if (leafa .gt. 1.0) then
          write(bffr,*) leafa
          call abortrun('leafa = '//bffr//' > 1.0 in leafa')
        endif
      endif

      return
      end
