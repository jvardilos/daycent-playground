
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... LACALC.F

      real function lacalc(fbrchc, rlwodc, maxlai, klai)

      implicit none

! ... Argument declarations
      real,intent(IN) :: fbrchc, rlwodc, maxlai, klai

! ... CALLED FROM:  potfor

! ... Calculate true LAI using leaf biomass and a biomass-to-LAI
!     conversion parameter which is the slope of a regression
!     line derived from LAI vs Foliar Mass for Slash Pine.

!     Calculate theoretical LAI as a function of large wood mass.
!     There is no strong consensus on the true nature of the relationship
!     between LAI and stemwood mass.  Version 3.0 used a negative exponential
!     relationship between leaf mass and large wood mass, which tended to
!     break down in very large forests.  Many sutdies have cited as "general"
!     an increase of LAI up to a maximum, then a decrease to a plateau value
!     (e.g. Switzer et al. 1968, Gholz and Fisher 1982).  However, this
!     response is not general, and seems to mostly be a feature of young
!     pine plantations.  Northern hardwoods have shown a monotonic increase
!     to a plateau  (e.g. Switzer et al. 1968).  Pacific Northwest conifers
!     have shown a steady increase in LAI with no plateau evident (e.g.
!     Gholz 1982).  In this version, we use a simple saturation fucntion in
!     which LAI increases linearly against large wood mass initially, then
!     approaches a plateau value.  The plateau value can be set very large to
!     give a response of steadily increasing LAI with stemwood.

! ... References:
!          1)  Switzer, G.L., L.E. Nelson and W.H. Smith 1968.
!              The mineral cycle in forest stands.  'Forest
!              Fertilization:  Theory and Practice'.  pp 1-9
!              Tenn. Valley Auth., Muscle Shoals, AL.
!
!          2)  Gholz, H.L., and F.R. Fisher 1982.  Organic matter
!              production and distribution in slash pine (Pinus
!              elliotii) plantations.  Ecology 63(6):  1827-1839.
!
!          3)  Gholz, H.L.  1982.  Environmental limits on aboveground
!              net primary production and biomass in vegetation zones of
!              the Pacific Northwest.  Ecology 63:469-481.

!  return a 0.0 instead of a divide by zero if C = 0 and klai is also 0 KLK 4/3/2017
!   This probably occurs before tree growth while the tree is not defined.

! ... Local variables
!        real rlai, tlai ! rlai - LAI from leaf biomass  tlai - LAI from large wood biomass
         real wdcarbon

! ... Due to the dynamic carbon allocation, return only the theoretical LAI, cak - 07/24/02
!        rlai = (rleavc * 2.5) * btolai

! ... wood carbon includes the fine branch in addition to the large wood, cak - 10/20/2006
        wdcarbon = fbrchc + rlwodc
        if(klai + wdcarbon .gt. 0.0) then
          lacalc = maxlai * wdcarbon / (klai + wdcarbon)
        else
          lacalc = 0.
        endif

      return
      end function lacalc
