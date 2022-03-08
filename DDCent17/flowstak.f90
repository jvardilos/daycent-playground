!               Copyright 1998 Colorado State University
!                       All Rights Reserved

! Flow system
!       schedule a 'flow' of material from one state variable to another;
!
!       Flow    stores these transferes in a stack.
!       flowl   is an alternative version which checks all scheduled flows
!       flowup  moves the stored quantity 'howmch' moved from 'source'
!           to 'sink' for those flows that were scheduled to occur before time.
!           Flows which are scheduled to occur after 'time' are pushed to the
!           top of the stack.
!
!       The variable 'when' tells the 'time' at which the flow is to
!           be completed; the first call to flowup with 'time' greater
!           then or equal to 'when' will cause the flow to complete.
!
!	flowup will complete the flows which
!

 module flowstak

   implicit none
!  private                      ! hide everything except the function calls
!  public flow, flowup, floclr, nflows
   public flowdat, nflows, stackhead, ierr


   type flowdat
!      private
!      integer                   :: when           !/* Time to flow */
      real                      :: when            !/* Time to flow */
!      character(LEN=16)         :: labl           !/* Flow frame label */
      real,   pointer           :: donr => null()  !/* Source */
      double precision, pointer :: ddonr => null() !/* Source */
      real,   pointer           :: rcvr => null()  !/* Destination */
      double precision, pointer :: drcvr => null() !/* Destination */
      double precision          :: amt             !/* Amount */
      type (flowdat), pointer   :: next => null()  !next entry
   end type flowdat

   integer                      :: ierr
   integer                      :: nflows = 0      !/* Number of flows */
   type (flowdat), POINTER      :: stackhead => null()
   save nflows, stackhead, ierr

 contains

   subroutine prntflw(stg,adv)
     character (len=*) :: stg
     character (len=1) :: adv

     if(adv .eq. "N") then
       call strout(trim(stg)//char(0)) ! write(*,'(a)', ADVANCE="NO") trim(stg);
     else
       call mssg(trim(stg)//char(0))   ! write(*,'(a)') trim(stg);
     endif
   end subroutine prntflw

   subroutine printflowdat(tmp_ptr)
     type (flowdat), pointer    :: tmp_ptr
     character (len=128) :: stg

     write(stg,'("at ",f9.3," move ",f15.7,":")') tmp_ptr%when,tmp_ptr%amt
     call strout(trim(stg)//char(0)) ! write(*,"(a)", ADVANCE='NO') stg
     if(associated(tmp_ptr%donr)) then
       write(stg,'(f15.7,z15,3x)') tmp_ptr%donr, loc(tmp_ptr%donr)
       call strout(trim(stg)//char(0)) ! write(*,"(a)", ADVANCE='NO') stg
     else
       write(stg,'(f15.7,z15,3x)') tmp_ptr%ddonr, loc(tmp_ptr%ddonr)
       call strout(trim(stg)//char(0)) ! write(*,"(a)", ADVANCE='NO') stg
     endif
write(stg,'("t: ",L3,L3,3x)') associated(tmp_ptr%rcvr),associated(tmp_ptr%drcvr)
call strout(trim(stg)//char(0)) ! write(*,"(a)", ADVANCE='NO') stg
     if(associated(tmp_ptr%rcvr)) then
       write(stg,'(" => ",f15.7,z15,3x)') tmp_ptr%rcvr, loc(tmp_ptr%rcvr)
       call mssg(trim(stg)//char(0)) ! write(*,"(a)") stg
     else
       write(stg,'(" => ",f15.7,z15,3x)') tmp_ptr%drcvr, loc(tmp_ptr%drcvr)
       call mssg(trim(stg)//char(0)) ! write(*,"(a)") stg
     endif

   end subroutine printflowdat

 end module flowstak


module calflow

!subroutine foo(a, b)
!real:: a
!integer::b
!  ...
!end subroutine foo

!subroutine bar(x, y)
!integer::x
!integer::y
!end subroutine bar

  interface flow
    module procedure flowff, flowdf,flowdfd, flowfd,flowfdd, flowdd,flowddd
  end interface

  interface flowl
    module procedure flowlff, flowldf, flowlfd, flowldd
  end interface

  interface flowbalanc
    module procedure balancd, balancf
  end interface

  contains


  subroutine flowff (source, sink, schtime, howmch)!, label)
     use flowstak

     implicit none
     real, target                  :: source, sink
     real, intent(IN)              :: schtime
!     integer, intent(IN)           :: schtime
     real, intent(IN)              :: howmch
!     character(LEN=16), intent(IN), optional   :: label

! Added the label string.  This is an optional tag the labels the flow
!       frame and allows you to put in a debugging text message
     type (flowdat), pointer    :: newflow

! ***  Don't bother with zero and insignificant, less than precision, flows
     if (abs(howmch) .gt. min(spacing(source),spacing(sink))) then
        allocate (newflow)          ! Create a new frame
        newflow%next => stackhead   ! it points to the current head
        stackhead => newflow        ! link it as the new head entry
        nflows = nflows + 1         ! increment flow count
!        nullify(stackhead%donr, stackhead%rcvr)
        nullify(stackhead%ddonr,stackhead%drcvr)

! Link in the appropriate data setting donor and receiver so that flow >0
        stackhead%when =  schtime
        if(howmch .gt. 0.) then
          stackhead%donr => source
          stackhead%rcvr => sink
          stackhead%amt  =  howmch
        else
          stackhead%donr => sink
          stackhead%rcvr => source
          stackhead%amt  =  abs(howmch)
        endif
!        stackhead%labl = ''; if(present(label)) stackhead%labl = label;
     end if
     return
  end subroutine flowff


  subroutine flowdd (source, sink, schtime, howmch)!, label)
     use flowstak
     implicit none
     double precision, target      :: source, sink
     real, intent(IN)              :: schtime
!     integer, intent(IN)           :: schtime
     real, intent(IN)  :: howmch

    call flowddd(source, sink, schtime, Dble(howmch))
  end subroutine flowdd

  subroutine flowddd (source, sink, schtime, howmch)!, label)
     use flowstak

     implicit none
     double precision, target      :: source, sink
     real, intent(IN)              :: schtime
!     integer, intent(IN)           :: schtime
     double precision, intent(IN)  :: howmch
!     character(LEN=16), intent(IN), optional   :: label

! Added the label string.  This is an optional tag the labels the flow
!       frame and allows you to put in a debugging text message
     type (flowdat), pointer    :: newflow

! ***  Don't bother with zero and insignificant, less than precision, flows
     if (abs(howmch) .gt. min(spacing(source),spacing(sink))) then
        allocate (newflow)          ! Create a new frame
        newflow%next => stackhead   ! it points to the current head
        stackhead => newflow        ! link it as the new head entry
        nflows = nflows + 1         ! increment flow count
        nullify(stackhead%donr, stackhead%rcvr)
!        nullify(stackhead%ddonr,stackhead%drcvr)

! Link in the appropriate data setting donor and receiver so that flow >0
        stackhead%when =  schtime
        if(howmch .ge. 0.) then
          stackhead%ddonr => source
          stackhead%drcvr => sink
          stackhead%amt  =  howmch
        else
          stackhead%ddonr => sink
          stackhead%drcvr => source
          stackhead%amt  =  abs(howmch)
        endif
!        stackhead%labl = ''; if(present(label)) stackhead%labl = label;
     end if
     return
  end subroutine flowddd


!====================================================================
  subroutine flowdf (dsource, sink, schtime, howmch)!, label)
     use flowstak
     implicit none
     double precision, target      :: dsource
     real, target                  :: sink
     real, intent(IN)              :: schtime
!     integer, intent(IN)           :: schtime
     real, intent(IN)              :: howmch

    call flowdfd(dsource, sink, schtime, Dble(howmch))
  end subroutine flowdf

  subroutine flowdfd (dsource, sink, schtime, howmch)!, label)
     use flowstak

     implicit none
     double precision, target      :: dsource
     real, target                  :: sink
     real, intent(IN)              :: schtime
!     integer, intent(IN)           :: schtime
     double precision, intent(IN)  :: howmch
!     character(LEN=16), intent(IN), optional   :: label

! Added the label string.  This is an optional tag the labels the flow
!       frame and allows you to put in a debugging text message
     type (flowdat), pointer    :: newflow

     if (howmch < 0.0) call flowfdd(sink, dsource, schtime, -howmch)

! ***  Don't bother with zero and insignificant, less than precision, flows
     if (abs(howmch) .gt. min(spacing(dsource),spacing(sink))) then
        allocate (newflow)          ! Create a new frame
        newflow%next => stackhead   ! it points to the current head
        stackhead => newflow        ! link it as the new head entry
        nflows = nflows + 1         ! increment flow count
        nullify(stackhead%donr, stackhead%drcvr)

! Link in the appropriate data setting donor and receiver so that flow >0
        stackhead%when =  schtime
        stackhead%ddonr => dsource
        stackhead%rcvr  => sink
        stackhead%amt   =  howmch
!        stackhead%labl = ''; if(present(label)) stackhead%labl = label;
     end if
     return
  end subroutine flowdfd


  subroutine flowfd (source, dsink, schtime, howmch)!, label)
     use flowstak
     implicit none
     real, target                  :: source
     double precision, target      :: dsink
     real, intent(IN)              :: schtime
!     integer, intent(IN)           :: schtime
     real, intent(IN)  :: howmch

    call flowfdd(source, dsink, schtime, Dble(howmch))
  end subroutine flowfd

  subroutine flowfdd (source, dsink, schtime, howmch)!, label)
     use flowstak

     implicit none
     real, target                  :: source
     double precision, target      :: dsink
     real, intent(IN)              :: schtime
!     integer, intent(IN)           :: schtime
     double precision, intent(IN)  :: howmch
!     character(LEN=16), intent(IN), optional   :: label

! Added the label string.  This is an optional tag the labels the flow
!       frame and allows you to put in a debugging text message
     type (flowdat), pointer    :: newflow

     if (howmch < 0.0) call flowdfd (dsink, source, schtime, -howmch)

! ***  Don't bother with zero and insignificant, less than precision, flows
     if (abs(howmch) .gt. min(spacing(source),spacing(dsink))) then
        allocate (newflow)          ! Create a new frame
        newflow%next => stackhead   ! it points to the current head
        stackhead => newflow        ! link it as the new head entry
        nflows = nflows + 1         ! increment flow count
        nullify(stackhead%ddonr, stackhead%rcvr)

! Link in the appropriate data setting donor and receiver so that flow >0
        stackhead%when =  schtime
        stackhead%donr  => source
        stackhead%drcvr => dsink
        stackhead%amt   =  howmch
!        stackhead%labl = ''; if(present(label)) stackhead%labl = label;
     end if
     return
  end subroutine flowfdd




  subroutine flowlff (var1, var2, schtime, howmch)!, label)
     use flowstak

     implicit none
     real, target, intent(IN)                  :: var1, var2
     real, intent(IN)                          :: schtime
!     integer, intent(IN)                       :: schtime
     real, intent(IN)                          :: howmch
!     character(LEN=16), intent(IN), optional   :: label

!
!       limited flow: schedule a 'flow' with a limit check on source;
!           records an amount so 'flowup' moves the quantity 'howmch'
!           from 'source' to 'sink'
!
!       This limited flowl scans all scheduled flows and guarantees this move
!           will not drive source negative; the amount moved is truncated as
!           necessary. howmch will not go negative, so this will not "fix"
!           it if previously scheduled events have driven source negative!
!
!       Calculates an err condition but it doesn't seem to return it.
!           This table shows the flow amount and err values for the final
!           source value, AFTER all outstanding flows are taken, and howmch;
!
!           source                 flow            ierr
!           source <= 0.            0.             -2
!           0 < source < howmch    source          -1
!           howmch <= source       howmch           0
!
!       N.B. if howmch is less than 0., the flow will be scheduled;
!           no guarantee is made regarding var2 not going negative!
!
!       The variable 'when' tells the 'time' when the flow is to complete.
!           flowup executes flows when 'time' is greater or equal to 'when'
!
!       KLK
!       Modified 31Dec10.
!         Changed the coding to track the original conditions more closely.
!           This corrects an apparent problem with the calls to flow.
!           (This coding is easier to analyze.)
!         If the flow is negative, change sense of the variables so the flow
!           is positive from source to sink. The check now works correctly.
!
!       Included for completeness. I have never seen it called.
!           Rewritten to do just the source balance then call flow
!

        real                         :: tmp, oput
        real, target                 :: source, sink

        if(howmch >= 0.) then
          source = var1; sink = var2;
        else
          source = var2; sink = var1;
        endif

! ***    Step through the stack
        call balancf (source, tmp, oput)
        tmp = source + tmp - oput

!       tmp is the var1 value after applying all currently outstanding flows
        if(tmp .le. 0.) then
!          if it is less than 0, don't schedule this flow, return -2
           ierr = -2
        else if(tmp .lt. abs(howmch)) then
!          requested removal exceeds value; limit to what is available, return -1
           call flow (source, sink, schtime, tmp)
           ierr = -1
        else
!          normal... schedule the flow and return 0
           call flow (source, sink, schtime, abs(howmch))
           ierr = 0
        endif

        return
  end subroutine flowlff

  subroutine flowldd (var1, var2, schtime, howmch)!, label)
     use flowstak

     implicit none
     double precision, target, intent(IN)      :: var1, var2
     real, intent(IN)                          :: schtime
!     integer, intent(IN)                       :: schtime
     double precision, intent(IN)              :: howmch
!     character(LEN=16), intent(IN), optional   :: label

!
!       limited flow: schedule a 'flow' with a limit check on source;
!           records an amount so 'flowup' moves the quantity 'howmch'
!           from 'source' to 'sink'
!
!       This limited flowl scans all scheduled flows and guarantees this move
!           will not drive source negative; the amount moved is truncated as
!           necessary. howmch will not go negative, so this will not "fix"
!           it if previously scheduled events have driven source negative!
!
!       Calculates an err condition but it doesn't seem to return it.
!           This table shows the flow amount and err values for the final
!           source value, AFTER all outstanding flows are taken, and howmch;
!
!           source                 flow            ierr
!           source <= 0.            0.             -2
!           0 < source < howmch    source          -1
!           howmch <= source       howmch           0
!
!       N.B. if howmch is less than 0., the flow will be scheduled;
!           no guarantee is made regarding var2 not going negative!
!
!       The variable 'when' tells the 'time' when the flow is to complete.
!           flowup executes flows when 'time' is greater or equal to 'when'
!
!       KLK
!       Modified 31Dec10.
!         Changed the coding to track the original conditions more closely.
!           This corrects an apparent problem with the calls to flow.
!           (This coding is easier to analyze.)
!         If the flow is negative, change sense of the variables so the flow
!           is positive from source to sink. The check now works correctly.
!
!       Included for completeness. I have never seen it called.
!           Rewritten to do just the source balance then call flow
!

        double precision             :: tmp, oput
        double precision, target     :: source, sink

        if(howmch >= 0.) then
          source = var1; sink = var2;
        else
          source = var2; sink = var1;
        endif

! ***    Step through the stack
        call balancd (source, tmp, oput)
        tmp = source + tmp - oput

!       tmp is the var1 value after applying all currently outstanding flows
        if(tmp .le. 0.) then
!          if it is less than 0, don't schedule this flow, return -2
           ierr = -2
        else if(tmp .lt. abs(howmch)) then
!          requested removal exceeds value; limit to what is available, return -1
           call flowddd (source, sink, schtime, tmp)
           ierr = -1
        else
!          normal... schedule the flow and return 0
           call flowddd (source, sink, schtime, abs(howmch))
           ierr = 0
        endif

        return
  end subroutine flowldd

  subroutine flowldf (source, sink, schtime, howmch)!, label)
     use flowstak

     implicit none
     double precision, target, intent(IN)      :: source
     real,   target, intent(IN)                :: sink
     real, intent(IN)                          :: schtime
!     integer, intent(IN)                       :: schtime
     double precision, intent(IN)              :: howmch
!     character(LEN=16), intent(IN), optional   :: label

!
!       limited flow: schedule a 'flow' with a limit check on source;
!           records an amount so 'flowup' moves the quantity 'howmch'
!           from 'source' to 'sink'
!
!       This limited flowl scans all scheduled flows and guarantees this move
!           will not drive source negative; the amount moved is truncated as
!           necessary. howmch will not go negative, so this will not "fix"
!           it if previously scheduled events have driven source negative!
!
!       Calculates an err condition but it doesn't seem to return it.
!           This table shows the flow amount and err values for the final
!           source value, AFTER all outstanding flows are taken, and howmch;
!
!           source                 flow            ierr
!           source <= 0.            0.             -2
!           0 < source < howmch    source          -1
!           howmch <= source       howmch           0
!
!       N.B. if howmch is less than 0., the flow will be scheduled;
!           no guarantee is made regarding var2 not going negative!
!
!       The variable 'when' tells the 'time' when the flow is to complete.
!           flowup executes flows when 'time' is greater or equal to 'when'
!
!       KLK
!       Modified 31Dec10.
!         Changed the coding to track the original conditions more closely.
!           This corrects an apparent problem with the calls to flow.
!           (This coding is easier to analyze.)
!         If the flow is negative, change sense of the variables so the flow
!           is positive from source to sink. The check now works correctly.
!
!       Included for completeness. I have never seen it called.
!           Rewritten to do just the source balance then call flow
!

        double precision             :: tmp, oput

        if(howmch < 0.) call flowlfd (sink, source, schtime, -howmch)

! ***    Step through the stack
        call balancd (source, tmp, oput)
        tmp = source + tmp - oput

!       tmp is the var1 value after applying all currently outstanding flows
        if(tmp .le. 0.) then
!          if it is less than 0, don't schedule this flow, return -2
           ierr = -2
        else if(tmp .lt. abs(howmch)) then
!          requested removal exceeds value; limit to what is available, return -1
           call flowdfd (source, sink, schtime, tmp)
           ierr = -1
        else
!          normal... schedule the flow and return 0
           call flowdfd (source, sink, schtime, abs(howmch))
           ierr = 0
        endif

        return
  end subroutine flowldf


  subroutine flowlfd (source, sink, schtime, howmch)!, label)
     use flowstak

     implicit none
     real,   target, intent(IN)                :: source
     double precision, target, intent(IN)      :: sink
     real, intent(IN)                          :: schtime
!     integer, intent(IN)                       :: schtime
     double precision, intent(IN)              :: howmch
!     character(LEN=16), intent(IN), optional   :: label

!
!       limited flow: schedule a 'flow' with a limit check on source;
!           records an amount so 'flowup' moves the quantity 'howmch'
!           from 'source' to 'sink'
!
!       This limited flowl scans all scheduled flows and guarantees this move
!           will not drive source negative; the amount moved is truncated as
!           necessary. howmch will not go negative, so this will not "fix"
!           it if previously scheduled events have driven source negative!
!
!       Calculates an err condition but it doesn't seem to return it.
!           This table shows the flow amount and err values for the final
!           source value, AFTER all outstanding flows are taken, and howmch;
!
!           source                 flow            ierr
!           source <= 0.            0.             -2
!           0 < source < howmch    source          -1
!           howmch <= source       howmch           0
!
!       N.B. if howmch is less than 0., the flow will be scheduled;
!           no guarantee is made regarding var2 not going negative!
!
!       The variable 'when' tells the 'time' when the flow is to complete.
!           flowup executes flows when 'time' is greater or equal to 'when'
!
!       KLK
!       Modified 31Dec10.
!         Changed the coding to track the original conditions more closely.
!           This corrects an apparent problem with the calls to flow.
!           (This coding is easier to analyze.)
!         If the flow is negative, change sense of the variables so the flow
!           is positive from source to sink. The check now works correctly.
!
!       Included for completeness. I have never seen it called.
!           Rewritten to do just the source balance then call flow
!

        real                         :: tmp, oput

!       scan the stack

        if(howmch < 0.) call flowldf (sink, source, schtime, -howmch)

! ***    Step through the stack
        call balancf (source, tmp, oput)
        tmp = source + tmp - oput

!       tmp is the var1 value after applying all currently outstanding flows
        if(tmp .le. 0.) then
!          if it is less than 0, don't schedule this flow, return -2
           ierr = -2
        else if(tmp .lt. abs(howmch)) then
!          requested removal exceeds value; limit to what is available, return -1
           call flowfdd (source, sink, schtime, dble(tmp))
           ierr = -1
        else
!          normal... schedule the flow and return 0
           call flowfdd (source, sink, schtime, abs(howmch))
           ierr = 0
        endif

        return
  end subroutine flowlfd


  subroutine balancf (var1, fiput, foput)
     use flowstak

     implicit none
     real, target, intent(IN)          :: var1
     real, intent(OUT)                 :: fiput, foput
!     character(LEN=16), intent(IN), optional   :: label


!       report the total flows in and out of a specified location. It
!       searches the flow stack summing the inputs into input and the
!       outputs into otput
!       based on the flowl routine
!
!       KLK
!       Written 24Jan11.
!         created to retrieve som2, som3, input for equilibrium calculations.

        integer i
        double precision          :: iput, oput

!       the start of the stack
        type (flowdat), pointer    :: tmp_ptr

        iput = 0; oput = 0. ! initialize the output variables

! ***    Step through the stack
        if (nflows .gt. 0) then
          tmp_ptr => stackhead
          do i= 1,nflows
            if(loc(var1) .eq. loc(tmp_ptr%donr)) then
              oput = oput + tmp_ptr%amt
! write(*,*) var1," out: ",tmp_ptr%amt
            else if(loc(var1) .eq. loc(tmp_ptr%rcvr)) then
              iput = iput + tmp_ptr%amt
! write(*,*) var1," In:  ",tmp_ptr%amt
            end if
            tmp_ptr => tmp_ptr%next
          end do
        end if

        fiput = real(iput); foput = real(oput)
        return
  end subroutine balancf


  subroutine balancd (var1, iput, oput)
     use flowstak

     implicit none
     double precision, target, intent(IN)      :: var1
     double precision, intent(OUT)             :: iput, oput
!     character(LEN=16), intent(IN), optional   :: label


!       report the total flows in and out of a specified location. It
!       searches the flow stack summing the inputs into iput and the
!       outputs into oput
!       based on the flowl routine
!
!       KLK
!       Written 24Jan11.
!         created to retrieve som2, som3, input for equilibrium calculations.


        integer i

!       the start of the stack
        type (flowdat), pointer    :: tmp_ptr

        iput = 0; oput = 0. ! initialize the output variables

! ***    Step through the stack
        if (nflows .gt. 0) then
          tmp_ptr => stackhead
          do i= 1,nflows
            if(loc(var1) .eq. loc(tmp_ptr%ddonr)) then
              oput = oput + tmp_ptr%amt
            else if(loc(var1) .eq. loc(tmp_ptr%drcvr)) then
              iput = iput + tmp_ptr%amt
            end if
            tmp_ptr => tmp_ptr%next
          end do
        end if

        return
  end subroutine balancd

end module calflow


  subroutine flowup (time)
     use flowstak

     implicit none
!     integer time, intent(IN) :: time
     real, intent(IN) :: time

     integer FlowsNotDone, i, istat
     type (flowdat), pointer    :: data_ptr
     type (flowdat), pointer    :: prevfram

     character (len=32) :: stg
     FlowsNotDone = 0

!
!       /* If there are any flows in the stack, determine which
!        * need to go now and do it.
!
!       This routine has a zero feature to help reduce problems from
!       precision (bit) errors. if the flow is approximately equal to
!       the donor, the donor is set to 0. The trigger for this is coded
!       as the absolute value of the residual being being one or two times
!       the donor spacing. (abs(donr - amt) < 3 *spacing(donr))
!       or alternately when loss precision(donr) exponent(donr) - exponent(donr - amt) <

     nullify(prevfram)
     data_ptr => stackhead
     if (nflows > 0) then
! *** set the active cell and clear the stack start
        do i=1,nflows
! *** step through the list
           if (time .ge. data_ptr%when) then

!        Error Check:  are the pointers associated
              if (associated(data_ptr%donr)) then
!                Remove the recorded flow from a single precision source
                 if(digits(data_ptr%donr)-(exponent(data_ptr%donr)-exponent(data_ptr%donr-data_ptr%amt)).lt. 3) then
!                 if(abs(data_ptr%donr - data_ptr%amt) .lt. 3*spacing(data_ptr%donr)) then
                   data_ptr%donr=0
                 else
                   data_ptr%donr = data_ptr%donr - data_ptr%amt
                 endif
              elseif (associated(data_ptr%ddonr)) then
!                Remove the recorded flow from a double precision source
                 if(digits(data_ptr%ddonr)-(exponent(data_ptr%ddonr)-exponent(data_ptr%ddonr-data_ptr%amt)).lt. 3) then
!                 if(abs(data_ptr%ddonr - data_ptr%amt) .lt. 3*spacing(data_ptr%ddonr)) then
                   data_ptr%ddonr=0
                 else
                   data_ptr%ddonr = data_ptr%ddonr - data_ptr%amt
                 endif
              else
                 stop 'Fatal Error: Unassociated "donr" pointer in flowup'
              endif

              if     (associated(data_ptr%rcvr)) then
                 data_ptr%rcvr = data_ptr%rcvr + data_ptr%amt
              elseif (associated(data_ptr%drcvr)) then
                 data_ptr%drcvr = data_ptr%drcvr + data_ptr%amt
              else
                 stop 'Fatal Error: Unassociated "rcvr" pointer in flowup'
              endif

              ! remove the completed flow from the stack
              if (associated(prevfram)) then
                 prevfram%next => data_ptr%next
                 deallocate (data_ptr, STAT = istat)     ! deallocate the flowdat
                 if(istat .ne. 0) then
                   write(stg,*) "istat = ",istat
                   call prntflw("flow deallocate error "// trim(stg),"")
                 endif
                 data_ptr => prevfram%next
              else
                 stackhead => data_ptr%next
                 deallocate (data_ptr, STAT = istat)     ! deallocate the flowdat
                 if(istat .ne. 0) then
                   write(stg,*) "istat = ",istat
                   call prntflw("flow deallocate error "// trim(stg),"")
                 endif
                 data_ptr => stackhead
              endif
           else
! *** This one doesn't need to be done yet
              FlowsNotDone = FlowsNotDone + 1

! *** set prevfram to the " last good" cell
              prevfram => data_ptr
              data_ptr => data_ptr%next
           endif
        end do
        nflows = FlowsNotDone
     endif
     return
  end subroutine flowup


  subroutine floclr
     use flowstak

     implicit none

     integer :: i
     type (flowdat), pointer    :: tmp_ptr

! *** Step through the stack dallocating all elements
     if (nflows .gt. 0) then
        do i= 1,nflows
           tmp_ptr => stackhead
           stackhead => tmp_ptr%next
           deallocate (tmp_ptr)
        end do
     end if
     nflows = 0
     return
  end subroutine floclr


  subroutine prstack
     use flowstak

     implicit none

     integer :: i
     type (flowdat), pointer    :: tmp_ptr
     character (len=32) :: stg

! *** Step through the stack dallocating all elements
     if (nflows .le. 0) then
       call prntflw("no flows",'')
     else
        tmp_ptr => stackhead
        i=0
        do while (associated(tmp_ptr))
          i=i+1
          write(stg,'(i5,3x)') i
          call prntflw(stg,'N')
          call printflowdat(tmp_ptr)
           tmp_ptr => tmp_ptr%next
        end do
     end if
     return
  end subroutine prstack

!  integer function flowcnt()
!     flowcnt=nflows
!  end function flowcnt

