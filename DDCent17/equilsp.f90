       MODULE eqAcc
         include 'const.inc'
         include 'cflows.inc'
         include 'dovars.inc'
         include 'parfx.inc'
         include 'param.inc'
         include 'plot1.inc'
         include 'plot2.inc'
         include 'timvar.inc'

         character*(256)         :: buffr ! output message buffer

!         real (c_float), save, bind(C, name="cNsoilt") :: cNsoilt(4);
         real (kind=8), save                      :: s3cin, s3cout
         real (kind=8), save, dimension(MAXIEL)   :: Etosom3;
         real, save                    :: som3st(2) = -1

         real, save                    :: convrat = 0   ! convergence ratio
         integer, save                 :: rptin   = 0   ! global self schedule repeat interval
         integer, save                 :: relaxtim = -2 ! post acceleration relaxation time

         integer, save :: ninpst = 0;
         real, parameter               :: tol = 1/(96.) ! default convergence ratio
       contains
         subroutine EqClear()
            s3cin   = 0.;
            s3cout  = 0.;

            som3st  = som3ci

            Etosom3   = 0;
         end subroutine EqClear
       END MODULE eqAcc


      integer function EqAcldat()
        use eqAcc
        implicit none

        EqAcldat = rptin
      end function EqAcldat

      subroutine EqAclratCE(cins3, eins3, couts3, cins2, eins2, couts2)
        use eqAcc
        implicit none

        real, dimension(2,MAXIEL) :: eins2
        real, dimension(2)        :: cins2, couts2
        real, dimension(MAXIEL)   :: eins3
        real                      :: cins3, couts3

          !character(len=9) :: date
          if(som3st(1) .eq. -1) call EqClear()

          s3cin  = s3cin  + cins3
          s3cout = s3cout + couts3
          Etosom3(1:nelem)  = Etosom3(1:nelem)  + eins3(1:nelem)
!  write(date,'(i3,"/",i5)') doy,cyear
!  write(*,*) 'EqAclratCE '//date, s3cin,Etosom3(1:nelem)," --- ",s3cout


          return
      end subroutine EqAclratCE


      subroutine EqAclrat(comand)
        use eqAcc
        implicit none

        ! Argument declarations
        character(len=*),           intent (IN) :: comand

        integer :: istat, lch, nv, sordr, m, ms, me
        character(len=9) :: date
        real (kind=8)    :: tmpval
        character (len=11),parameter :: keylist = "ARYI%?aryi%";
        character (len=4),parameter  :: evnt = "ACEQ";
        integer :: ReadInt
        real (kind=8) :: parsImedFld

        real (kind=8)                :: s3corRat = 0


        write(date,'(i3,"/",i5)') doy,cyear
        lch = index(trim(date)," ")
        date(lch:) = ADJUSTL(date(lch:))
        !write(*,*) date,comand,aceqcnt

        select case (comand(:3)) ! clear the equilibrium suns
          case ('AEQ')       ! update the SOM pools

           if(aceqcnt .eq. 0) return ! get out if we have shut this down

!      call message("EqAclrat(AEQ)  "//trim(date)//"   som3 =====")
!write(*,*) "EqAclrat(som3)  ",trim(date)," in:",s3cin,Etosom3(1:nelem),'   out:',s3cout
           if(s3cout>0) then
              call updatpool("EqAclrat(som3) "//trim(date),som3c, som3ci, som3st, s3cin, &
                                        s3cout, s3corRat, som3e, Etosom3, csrsnk, esrsnk)
             ! this is a reinitialization; accept the single precision values
             ! DPsom3e(1:nelem) = som3e(1:nelem)
             ! DPsom3ci = som3ci; DPsom3c = sum(DPsom3ci);
           endif

           ! removed the unused SOM2 acceleration code
           ! further development and evaluation of the impact of SOM2 to convergence is
           ! needed since current tests shows it to inhibit convergence.

           aceqcnt = max(aceqcnt -1, 0) !Decrement the iteration #
           if(convrat > 0  .and.  s3cout > 0) then
             if(abs(s3corRat -1.0) < convrat) then
               ! add the relaxation time to the current year
               if(relaxtim .le. 0) then
                 istat = cyear + rptin*abs(relaxtim)
               else
                 istat = cyear + relaxtim
               endif

               ! should we modify the end time? normally this is an early termination
               write(buffr,'(2i8)') cyear, istat
               if(istat<=tend) then
                 tend = istat ! shorten the run
                 ! log the
                 call message("ACEQ converged in "//buffr(:8)//": terminating at year "//buffr(9:))

               else
                 ! extend run limits for the cool down period.
                 tend = istat
                 if (blktnd < istat+1) blktnd = istat+1 ! the block must be long enough
                 ! Warn the user we extended the run
                 call message("ACEQ converged in "//buffr(:8)//": extending run to "//buffr(9:))
               endif

               aceqcnt = 0 ! terminate the acceleration
             end if
           end if

           call EqClear() ! clear the flow sums for the next iteration

!VERIFY(STRING, SET, back)  the position of the first character in STRING not in SET.
!SCAN(STRING, SET, back)    the position of the first character from SET in the string STRING.
          case ("SET")

            ! set the count from the **FIRST** command parsed
            if(aceqcnt .lt. 0) then
              istat = 0; lch   = 4
              call message('Parsing ACEQ comand: '//trim(comand(4:)))

              ! bare number; this is just the step count
              ! recall verify(st, set) returns the first element of st not in set
              if(VERIFY(trim(comand(4:)), '0123456789') .eq. 0) then
                aceqcnt = ReadInt(0,comand,lch,istat)

              else
                m  = 4
                me = len_trim(comand)

                sordr=0 ! set the sort order to 0
                do while (m .le. me) ! the last character is the terminator
                  !treat tab as a soft delimiter
                  ! remove spaces and parenthasis
                  if(index(' ()', comand(m:m)) .gt. 0) then; m=m+1; cycle; endif
                  ms = m ! save the field start for error reporting
                  ! keylist = "ARYI%?aryi%";
                  tmpval = parsImedFld(nv, comand, sordr,istat,m, ms,me,keylist,evnt)

                  select case (nv)
                    case (1) ! A automatic configuration set defaults input self scheduling
                      aceqcnt  = tmpval
                      convrat  = tol
                      relaxtim = -1

                    case (2) ! R input self scheduling
                      aceqcnt = tmpval

                    case (3) ! Y input relaxation time, years
                      relaxtim = tmpval

                    case (4) ! global self schedule repeat interval
                      rptin   = max(tmpval, 0.)

                    case (5) ! % convergence tolerance in percent
                      convrat = tmpval/100.

                    case DEFAULT
                      call abortrun('Unknown '//evnt//' key '//comand(3+ms:)//' in event '//comand(4:))
                  end select
                end do
                if (rptin == 0) then
                  call abortrun('ACEQ automatic repeat interval =0')
                else if(aceqcnt == 0) then
                  aceqcnt = int((tend-strtyr+1)/rptin)
                endif
              endif
            endif
            return

          case ("Cle")
            call EqClear()

          case default
            call message ("Unknown EqAclrat command "//comand)
        end select

        return
      contains

        subroutine updatpool (label,somC, somCi, somstrt, somCin, somCout, corRat, &
                              somE, somEin, csrsnk, esrsnk)

          character(len=*)                 :: label
          real (kind=4)                    :: somC, somCi(2)
          real (kind=8)                    :: somCin, somCout
          real (kind=8)                    :: corRat
          real (kind=8)                    :: tmp, dsomCi(2)
          real (kind=8), dimension(MAXIEL) :: somEin
          real (kind=4), dimension(MAXIEL) :: somE, esrsnk, delsome
          real (kind=4), dimension(2)      :: somstrt, csrsnk

        !==============================

           ! where (somEin(1:nelem).gt.0);  recin(1:nelem)  = somEin(1:nelem)/somCin;   end where

            tmp    = abs((dble(sum(somstrt))+(somCin-somCout))/dble(somC) - 1.d0)
        write(buffr,'(a,i5,3(a,f12.3),f12.3,a,1pE12.3)') label,aceqcnt,"  somstrt=",sum(somstrt)," somC",somC, &
              "  somCin,out=",somCin,-somCout," sum err:",tmp

            ! correct the C pool
            dsomCi    = somCi                  ! store the starting pt so we calculate change
            delsome(1:nelem) = somE(1:nelem);  ! store current E

            ! calculate the correction ratio
            if(somCout==0) then
              call message("EqAclratCE Warning: no C output from "//trim(label)//" skipping this adjustment")
              return
            else if (somCin == 0) then
              call message("EqAclratCE Warning: no C input to "//trim(label)//" skipping this adjustment")
              return
            endif

            corRat    = somCin/somCout  ! the correction ratio

            somCi  = somCi*corRat ! change somCi so flows match
            dsomCi = somCi - dsomCi       ! change in pool carbon
            somC   = sum(somCi)           ! sum the estimated pool
            csrsnk = csrsnk - dsomCi    ! update the C source/sink

            ! assume constant E/C; use the old total C, somC, to calculate somEout = somE * somCout/somC
            ! somE(1:nelem)  = somE * (somEin/(somE * (somCout/somC)))   ! change som E by C ratio
            !                = somEin/(somCout/somC)   ! change som E by C ratio
            somE(1:nelem)    = somC * (somEin(1:nelem)/somCin)   ! change som E by C ratio
            delsome(1:nelem) = somE(1:nelem) - delsome(1:nelem);    ! average input E/C ratio

        write(buffr(len_trim(buffr)+2:),'(4(a,f12.3))') " Ratio:",corRat,"  C/N ",somCin/somEin(1), &
        "  New",somC," adjust ", sum(dsomCi)

        call message(trim(buffr))

        if(sum(somCi) .eq. 0) call abortrun("EqAclratCE "//label//"updated C pool is zero;  no growth or error in accleration")
        !==============================

          return
        end subroutine updatpool
      end subroutine EqAclrat


      subroutine EqAclratFire(som3cin, som3ein, som2cout)
        use eqAcc

        real    :: som3cin, som3ein(MAXIEL), som2cout

         if(aceqcnt <=0) return

         s3cin            = s3cin   + som3cin
         Etosom3(1:nelem) = Etosom3(1:nelem) + som3ein(1:nelem)
! write(*,*) 'in EqAclratFire(',som3cin,som3ein
      end subroutine EqAclratFire


