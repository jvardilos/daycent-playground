
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


!      modified:
!      11 Oct 2014  K. Killian
!        made timary match fortran column-major storage format
!        This effects schedl.f90 readblk.f90 csa_detiv.f90 chrvar.inc
!      15Jan  K. Killian
!        combined schvar and chrvar since new Fortran specs don't care
!      May2011  K. Killian
!        F90 style declarations
!        increased the length of weather file variable
!        increased the option variable length
!        returned event array back to 4 characters.
!        added a parameter to set the length of these arrays
!        reordered the schvar common to put the arrays at the end.
!        added a pointer to the global events at the top end of the array
!        increased the length of weather file variable
!        removed fltary (foating argument) array
!        made wlabel a parameter array
!        included schvar in this file. schvar is called with chrvar and
!          this puts all the event and command variables in the same file
!          Put all the event, command, variables in common schvar.
!          Note 5.33 Fortran90 spec allows character and numeric in
!          the same common block.


!**** Scheduling option variables
!     evtptr	array index of next event in cmdary, timary, typary
!     mxevpt    array index for global clock events
!     fltary    array holding real (float) values for each event
!     rptyrs    number of years actually included in the block
!     timary	array holding year and step count (month, week) for each event
!     ttlind    total number of events read in a block

!     Character values needed in common block:
!     cmdary    array holding command for each event in the block
!     curcrp    the current crop
!     curtre    the current tree
!     initcp    the initial crop
!     initre    the initial tree
!     typary    array holding command type for each event in block
!     wlabel	string of weather labels "prec", "tmin", "tmax"
!     wthnam    if wthr = F, the weather file name
!     wthr      the weather choice

      integer, parameter :: lenary=2500

      integer            :: evtptr, glevpt, rptyrs, ttlind
      logical            :: dbgpsch
      integer            :: timary(2,lenary)

      character (len=4)  :: cmdary(lenary)
      character (len=28) :: curcrp,curtre, typary(lenary)
      character (len=4)  :: wthr
      character (len=132):: wthnam

      character (len=4), parameter :: wlabel(3)=(/'prec','tmin','tmax'/)

      common/schvar/evtptr, glevpt, rptyrs, ttlind, dbgpsch, wthr,       &
     &              wthnam, curcrp, curtre, timary, cmdary, typary


      save /schvar/

