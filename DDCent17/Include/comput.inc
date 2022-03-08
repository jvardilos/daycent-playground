
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... Variables used internally which must be shared among routines.
!     Change p1co2 to p1co2(2) to represent surface and soil layers.
!     Change ratnew to ratnew1 and ratnew2 for surface and soil decomposition MDH 9/24/2012
! ... baseNdep - base amount of atmospheric N deposition for the year
! ...            based on average annual precipitation only
! ... cercrp(mx/mn, nparts, nelem)
! ... fps1s3 - the fraction of som1 decomposing to som3
! ... fps2s3 - the fraction of som2 decomposing to som3
! ... lhzci(pool,iso) - the lower horizon carbon used as an input 
! ...                   to the surface soil in an erosion event
! ... lhze(pool,iel) - the lower horizon N,P,S used as an input 
! ...                  to the surface soil in an erosion event
! ... orglch - the fraction of organics that leach from soil som1 when
! ...          there is sufficient water flow.  Computed in predec.
! ... wc - the water capacity for the top layer afiel(1)-awilt(1)
! ... ratnew1(iel,1) - the C/E ratio for new material created when a
! ...                  lignin component decomposes to SRFC som1.
! ... ratnew2(iel,1) - the C/E ratio for new material created when a
! ...                  lignin component decomposes to SOIL som1.
! ... ratnew1(iel,2) - the C/E ratio for new material created when a
! ...                  lignin component decomposes to SRFC som2.
! ... ratnew2(iel,2) - the C/E ratio for new material created when a
! ...                  lignin component decomposes to SOIL som2.

      common/comput/agdefacm,bgdefacm,baseNdep,cercrp,eftext,            &
     &  fps1s3,fps2s3,lhzci,lhze,orglch,p1co2,h2ogef,wc,                 &
     &  ratnew1,ratnew2


      real agdefacm(12),  bgdefacm(12), baseNdep
      real cercrp(2,2,3), eftext,       fps1s3,       fps2s3
      real lhzci(3,2),   lhze(3,3),     ratnew1(3,2), ratnew2(3,2)
      real p1co2(2),      h2ogef(2),    orglch,       wc

      save /comput/
