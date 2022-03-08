
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


subroutine detiv(wrbin)
 USE ISO_C_Binding
 use, intrinsic :: iso_fortran_env ! get the precision parameters REAL32, REAL64, REAL128

 implicit none
 include 'chrvar.inc'
 include 'const.inc'
 include 'doubles.inc'
 include 'dovars.inc'
 include 'fertil.inc'
 include 'param.inc'
 include 'parcp.inc'
 include 'parfs.inc'
 include 'parfx.inc'
 include 'plot1.inc'
 include 'plot3.inc'
 include 'potent.inc'
 include 'seq.inc'
 include 'site.inc'
 include 'sitsoil.inc'
 include 'timvar.inc'
 include 'wth.inc'
 include 'wthdaily.inc'
 include 'zztim.inc'

! ... Determine name of schedule file, which contains the
! ... name of the site file, values of timing variables,
! ... and order of events

! Extend Options
!   Modified extend arguments to allow the extend record to be matched to
!     the starting time.  This requires passing strtyr to extend.
!   extndopt
!            0 no extend
!            1 Site file extend
!            2 copy data (normal)
!            3 initial condition
!            4 append to existing binary

! ... Fortran to C prototype
 INTERFACE

   SUBROUTINE initsw(sitlat, daylength, rwcf, adep)
     USE ISO_C_BINDING
     !MS$ATTRIBUTES ALIAS:'_initsw' :: initsw
     real(c_float)  :: sitlat
     real(c_float)  :: daylength(*), rwcf(*), adep(*)
   END SUBROUTINE initsw

   SUBROUTINE outfils()
   END SUBROUTINE outfils

   SUBROUTINE setasmos(asmos, nlayer, swc, numlyrs, avh2o, rwcf)
     USE ISO_C_BINDING
     !MS$ATTRIBUTES ALIAS:'_setasmos' :: setasmos
     INTEGER        :: nlayer, numlyrs
     real(c_float)  :: asmos(*), avh2o(*), rwcf(*)
     real(c_double) :: swc(*)
   END SUBROUTINE setasmos

   SUBROUTINE setlyrs(adep,nlayer,numlyrs, sand, silt, clay,       &
                      bulkd, ph, awilt, afiel, swflag,swdbgf)
     USE ISO_C_BINDING
     !MS$ATTRIBUTES ALIAS:'_setlyrs' :: setlyrs
     real(c_float)  :: adep(*)
     INTEGER        :: nlayer, numlyrs
     real(c_float)  :: sand, silt, clay, bulkd, ph
     real(c_float)  :: awilt(*), afiel(*)
     INTEGER        :: swflag
     integer        :: swdbgf(*)
   END SUBROUTINE setlyrs

   SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3,           &
              ammonium, nitrate, subname)
     USE ISO_C_BINDING
     !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
     INTEGER        :: clyr
     real(c_float)  :: amt
     real(c_double) :: frac_nh4, frac_no3, ammonium
     real(c_double) :: nitrate(*)
     CHARACTER      :: subname*10
   END SUBROUTINE update_npool

 END INTERFACE
 ! ... arguments
 logical     :: wrbin

 ! ... Function declarations
 integer ReadInt
 double precision ReadReal

 ! this common contains the library path
 character (len=256) filpath
 integer plen
 common /libpath/plen, filpath
 save /libpath/

 ! ... Local variables
 real        :: rt, mrtfr
 integer     :: oslen, rindx, soildepth, cld, sndx
 integer     :: nes, nee, nss, nse
 integer     :: clen, ii, istat, lch
 integer     :: extndopt, varcnt, outvar
 integer     :: swdbgf(2) = (/0, 0/) ! verbose, debug
 integer     :: lschnam(2), loc(2), srlen ! location of schedule file name root

 logical     :: goahead, ext_grid
 logical     :: sitparread, soilread, sitefirst

 character (len=84)  :: buffr
 character (len=256) :: inlin   ! input line to parse
 character (len=256) :: newbin, oldbin, schnam, ascnam, sitnam ! file names
 character (len=10)  :: subname = 'detiv'
 integer, parameter :: sext = 2
 character (len=4),  parameter :: sitext(sext) = (/'.100', '.dsa'/)
 character (len=1),  parameter :: FNSOIL(8) = 'soil.in'//char(0)
 character (len=1),  parameter :: FNSITE(11) = 'sitepar.in'//char(0)

 ! Initialize potential command line arguments
 extndopt = 0
 ext_grid = .false.
 ascnam = '{S}'
 schnam = '-'
 lschnam = 0
 newbin = ' '
 oldbin = ' '
 sitnam = ' '
 filpath = ' '
 inlin   = ' '
 ! ltrmeth = .false.
 sitparread = .FALSE.
 soilread   = .FALSE.
 sitefirst  = .FALSE.

 cursys = 0  ! clear current system
 varcnt =  0
 outvar =  0
 rindx  = 0  ! clear the binary site file read index
 sndx   = 0  ! clear site index
 oslen  =  0
 istat  = 0  ! clear the file I/O status
 rt     = -1

 ! parse the command line arguments
 call commandline()

! ... Check that minimal information was entered
 if (schnam .eq. '-' .or. schnam .eq. ' ') call abortrun('No schedule file name given.')

 !  Check that SOME output file is defined.
 !if (newbin .eq. ' ' .and. extndopt .gt. 0) then
   !   newbin = oldbin
 !end if
 if(newbin .eq. ' ' .and. varcnt .le. 0 .and. oslen .eq. 0) then
 ! No output file specified; define a binary file using the event file
   if(lschnam(1) .ne. 0) newbin = schnam(lschnam(1):lschnam(2))//'.bin'
   clen = clen+4
   call message('Writing output to binary file '//newbin(:clen))
 end if

!...  Make sure there is a list file if we specified output variables.
 if (varcnt .gt. 0 .and. ascnam(1:4) .eq. '{S}') then
!   No lis file specified; define a text file name using the event file
    if(lschnam(1) .ne. 0) ascnam = schnam(lschnam(1):lschnam(2))//'.lis'
    call message('The output list file is '//trim(ascnam))
 end if

!...open the ascii list file and write the header records
 call listvar(21,rt,schnam,ascnam)


! ... Add new routine to initialize common block variables to other than
! ... default values as necessary
 !=======================================================
 call initialize(extndopt)
 call defaultsitpar(dDO_fc, swc, swdbgf) ! set sitepar defaults; routine in initsit_tg


 !=======================================================
! ... Open the schedule file and read the header lines
 open(unit=15,file=schnam,status='OLD', IOSTAT=istat)
 ! is the schedule file really there
 if (istat.gt.0) call abortrun('opening schedule file "'//trim(schnam)//'"')

! ... Allow comments at the top of a schedule file for documentation.
 ! Converted to Kendrick's input routines. Allows both comments, anything
 ! after a '#' character, and blank lines

 lch = 0
 istat = 0
 strtyr = ReadInt(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('missing run start time in line:'//trim(inlin))
 ! all initialization is dated the day BEFORE the simulation starts
 ! Used in partit, flow, csched, calciv, and flowup without having been initialized
 cyear  = strtyr -1
 blktnd = strtyr
 month  = 12
 doy    = 365
 time   = strtyr ! = cyear + month*dtv

 ! read(15,*) tend
 lch = 0
 tend = ReadInt(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('missing run end time in line:'//trim(inlin))

 ! read(15,*) sitnam
 ! replaced variable directed read ('*') since it couldn't handle the special
 ! characters in a file path. Used READCLIN since it handles comments 30Dec11
 call READCLIN(15,inlin,lch,.TRUE.)
 if(sitnam .eq. ' ') sitnam = inlin
 do ii = 1,sext
   clen = index(sitnam,sitext(ii))
   if (clen .ne.0) then ! we have a known extension
     clen = clen+3
     inquire(file=sitnam(:clen), exist=goahead)
     if(goahead) then
       if(sitext(ii) .eq. '.dsa') then
         if(clen+1 < len(sitnam)) then
           lch = index(sitnam(clen+1:),'R')
           if(lch > 0) then
             lch   =lch+clen+1
             istat =0
             sndx = ReadInt(0,sitnam,lch,istat)
             if(istat <= 0  .or.  sndx .le. 0) sndx = 0
           endif
         endif
       endif
       sitnam(clen+1:) = ''
       exit;
     endif
     exit
   endif
 end do

! we haven't found a known extension
! search for a bare file name
 clen = len_trim(sitnam)
 do while (clen .ge. 1)   ! step through spaces (allows for files with embedded spaces)
   inquire(file=sitnam(:clen), exist=goahead);
   if(goahead) then
     if(clen+1 < len(sitnam)) then
       lch = index(sitnam(clen+1:),'R')
       if(lch > 0) then
         lch   =lch+clen+1
         istat =0
         sndx = ReadInt(0,sitnam,lch,istat)
         if(istat <= 0  .or.  sndx .le. 0) sndx = 0
       endif
     endif
     sitnam(clen+1:) = ''
     exit ! found so jump to the end of the outer loop
   else
     clen = index(sitnam(:clen),' ',.TRUE.)
     clen = len_trim(sitnam(:clen))
   endif
 end do
 if(.not. goahead) call abortrun("finding site file: "//trim(sitnam))
 call message('reading site: '//trim(sitnam))

 lch = 0
 labtyp = ReadInt(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('reading labeling type:'//trim(inlin))
 if(labtyp .lt. 0) labtyp = 0
 if(labtyp .gt. 2) then
   call abortrun('labtyp greater than 2 not defined:'//trim(inlin))
 endif
 lch = 0
 labyr = ReadInt(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('reading labeling year:'//trim(inlin))

 lch = 0
 mctemp = ReadReal(15,inlin,lch,istat)
 if (mctemp .ge. 0) then
   call abortrun('Microcosms are not implemented in this version'  &
               //' of Daily Century')
 endif
 micosm = 0

 lch = 0; istat = 0;
 co2sys = ReadReal(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('reading CO2 flag:'//trim(inlin))
 if (co2sys .gt. 0) then
   clen=101
   co2tm(1) = ReadReal(15, inlin, clen, istat)
   if(istat.le.0) call abortrun('reading CO2 time 1:'//trim(inlin))
   if (clen .lt. len(inlin))  then
     co2tm(2) = ReadReal(15, inlin, clen, istat)
     if(istat.le.0) call abortrun('reading CO2 time 2:'//trim(inlin))
   else
     co2tm(2) = co2tm(1) + 0.083
   endif
 endif

 ! Added error messages for header inputs  KLK 19 May 2015
 ! For CO2, PH, warming, N, OMAD weather input
 ! scalars, the starting time phtm, ststart, Nstart, OMADstart, wthstart,
 ! and the soil temperature amount, stamt can be placed on the same line
 ! as the flag allowing the number of header lines to remain constant.
 ! For compatibility with earlier schedules, the parser reads old format,
 ! 1 value per line, without error.

 ! New schedule header line for pH shift, cak - 08/02/02
 ! Change pH shift implementation to use scalar values, cak - 10/17/05
 lch = 0; istat = 0;
 phsys = ReadInt(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('reading PH effect flag: '//trim(inlin))
 if (phsys .gt. 0) then
   phtm = ReadInt(15,inlin,lch,istat)
   if(istat.le.0) call abortrun('reading PH effect start: '//trim(inlin))
 endif

 ! New schedule header lines for soil temperature warming experiments, cak - 07/02/03
 lch = 0; istat = 0;
 stsys = ReadReal(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('reading soil temperature warming flag: '//trim(inlin))
 if (stsys .gt. 0) then
   ststart = ReadReal(15,inlin,lch,istat)
   if(istat.le.0) call abortrun('reading soil temperature warming start: '//trim(inlin))
   stamt = ReadReal(15,inlin,lch,istat)
   if(istat.le.0) call abortrun('reading soil temperature warming amount: '//trim(inlin))
 endif

 ! New schedule header line for N input scalars, cak - 04/06/04
 lch = 0; istat = 0;
 Ninput = ReadInt(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('reading N input scalars flag: '//trim(inlin))
 if (Ninput .gt. 0) then
   Nstart = ReadInt(15,inlin,lch,istat)
   if(istat.le.0) call abortrun('reading N input scalars start: '//trim(inlin))
 endif

 ! New schedule header line for OMAD input scalars, cak - 04/06/04
 lch = 0; istat = 0;
 OMADinput = ReadInt(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('reading OMAD input scalars flag: '//trim(inlin))
 if (OMADinput .gt. 0) then
   OMADstart = ReadInt(15,inlin,lch,istat)
   if(istat.le.0) call abortrun('reading OMAD input scalars start: '//trim(inlin))
 endif

 ! New schedule header line for weather input scalars, cak - 10/18/05
 lch = 0; istat = 0;
 wthinput = ReadInt(15,inlin,lch,istat)
 if(istat.le.0) call abortrun('parsing weather input scalar flag: '//trim(inlin))
 if (wthinput .gt. 0) then
   wthstart = ReadInt(15,inlin,lch,istat)
   ! This is here
   if(istat.le.0) then
     write(*,'(a)') "Check for missing header lines in the event file:"
     write(*,'(a)') 'Schedule File Header Table'
     write(*,'(a,i6)')       'run start   strtyr      ', strtyr
     write(*,'(a,i6)')       'run end     tend        ', tend
     write(*,'(a,i6,a,i6)')  'label type  labtyp      ', labtyp,   '      labyr     ',labyr
     write(*,'(a,f10.3)')    'Microcosms  mctemp      ', mctemp
     write(*,'(a,f10.3,a,2f10.3)') 'CO2 flag    co2sys      ', co2sys,   '  CO2 times ',co2tm
     write(*,'(a,i6,a,i6)')  'PH flag     phsys       ', phsys,    '      phtm      ',phtm
     write(*,'(3(a,f10.3))') 't warm      stsys       ', stsys,   '  ststart   ',ststart, '  stamt ',stamt
     write(*,'(a,i6,a,i6)')  'N scaler    Ninput      ', Ninput,   '      Nstart    ',Nstart
     write(*,'(a,i6,a,i6)')  'OMAD scaler OMADinput   ', OMADinput,'      OMADstart ',OMADstart
     write(*,'(a,i6)')       'weather scaler wthinput ', wthinput

     call abortrun('reading weather input scalar flag: '//trim(inlin))
   endif
 endif

 read(15,*) decsys
 if (decsys .eq. SAVSYS) then
   decsys = FORSYS
 endif

 ! ... Obtain the initial values for the crop or forest system
 read(15,'(a)') curcrp
 curcrp = adjustl(curcrp)
 if(index(curcrp,' ') .gt. 0) curcrp(index(curcrp,' '):) = ' ' ! drop everything after space
 if (curcrp(:4) .eq. 'Init') curcrp = ' '

 read(15,'(a)') curtre
 curtre = adjustl(curtre)
 if(index(curtre,' ') .gt. 0) curtre(index(curtre,' '):) = ' '
 if (curtre(:4) .eq. 'Init') curtre = ' '


!==============================================
!     read and store any global events placed before the first block
!     sort the global events to prevent out of order entries
 rptyrs = tend   ! tend is now integer KLK 12/2013
 call readglobalevt()
!==============================================

! ... Read starting values from fixed parameter file
 call fixin

 ! Read starting site values from DayCent sitein and soil.in
 !Note: pass structures back into C routines with the first element address
 inquire(FILE = 'soils.in', EXIST=soilread)
 if(soilread) call inputlyrs(dDO_fc, swc, swdbgf);

 ! Read starting values from DayCent sitepar file now so site.100 can over ride
 ! moved this call from initsw
 inquire(FILE = "sitepar.in", EXIST=sitparread)
 if(sitparread) call initsite(dDO_fc, swc, swdbgf);

 ! Read starting values from site-specific file
 if(sndx>0  .and.  rindx<=0) rindx = sndx
 if(rindx.le.0   .or.  sitnam(:len_trim(sitnam)-4).eq.'.100') then
   call oldopen(7,sitnam) !open(unit=7,file=sitnam,status='OLD',err=1000)
   if (ext_grid) then
     call sitein_grid(extndopt)
   else
     call sitein(extndopt, swdbgf(1), sitparread, soilread)
   endif

   ! find the schedule and site file root names
   nes = index(schnam, '/', .TRUE.)+1
   nee = len_trim(schnam)
   nss = index(sitnam, '/', .TRUE.)+1
   nse = len_trim(sitnam)

   inlin = schnam(nes:nee)//"<-"//sitnam(nss:nse) ! initial site description
 else
   inlin = ' ' ! clear junk in inline
   ! Read starting values from a binary site archive
   call binsitrw(trim(sitnam), rindx, 0, inlin, 16) ! daycent files >= 16
   sitparread = .true. ! binary read sets sitepar variables
   soilread   = .true. ! binary read has soil layers

   ! find the schedule file root
   nes = index(schnam, '/', .TRUE.)+1
   nee = len_trim(schnam)
   ! save binary record tag
   if(index(inlin, schnam(nes:nee)) .eq. 0) then
     ! add a new schedule to the list if it wasn't there before
     inlin = schnam(nes:nee)//"<-"//trim(inlin)
   endif

  ! change the extend option if there is not a higher priority extend
  ! and some critical enhanced values are non zero. (added to prevent errors
  ! when a normal site is copied into a version 1 archive) 8/3/2001 KLK
    if(extndopt .le. 0) then
       if(.not. (som1e(1,2) .eq.0 .and. som1e(2,1) .eq.0 .and.     &
                 som1e(1,2) .eq.0 .and. som1e(2,1) .eq.0 .and.     &
                 som3e(1)   .eq.0)) then
          extndopt = 1
       else
          call message("incomplete enhanced archive")
       end if
    end if
 endif

 ! calculate nlayer
 soildepth = nint(maxval(dpthmx))
 if(numlyrs .eq. 0) numlyrs = maxloc(dpthmx,dim=1)  !  transfer(maxloc(dpthmx),1)
 cld = 0
 do ii=1, UBOUND(adep, DIM=1)
   cld = cld + nint(adep(ii))
   if(cld .eq. soildepth) then
     nlayer = ii
     exit
   endif
 end do

 call sitepostp() ! post processing of site data

 if(.not. sitparread) call message( &
   'WARNING: sitepar parameters unread - using defaults');
 if(.not. soilread) call abortrun('soil layer data missing');

 call wrtsite(2,1,trim(inlin)) ! store the updated site description

 nlaypgst = nlaypg ! retain the initial nlaypg for site file output

! ... Moved the read calls for the initial tree and crop to inside the extend
! ... if statement.  This is done to prevent a rather subtle bug that occurs
! ... when the initial crop/tree do not match the final values in the
! ... original schedule.  In that case, the derived output values, (crpval ...)
! ... do not match the current crop values.
! ... The crop/tree reads must occur before the calciv call on a normal run.
! ... 7/20/95  K. Killian

! ... Open binary file to write
 ! co-located the binary open with the extend, the first time the .bin file is used
 if (newbin .ne. ' ') then
   call message("new binary: "//trim(newbin))
   ! removed open status='OLD'; oldbin has already checked for existence
   ! binstat='UNKNOWN'; if (extndopt .gt. 1  .and.  newbin .eq. oldbin) binstat='OLD';
   open(unit=1,file=newbin,form='UNFORMATTED',IOSTAT=istat, action='READWRITE') ! status=binstat
   wrbin=.TRUE.
   if (istat.gt.0) then
     write(buffr,'(i4)') istat
     call abortrun(buffr(:4)//' opening Binary file "'//trim(newbin)//'"')
   endif
 endif

! ... Determine initial values
!      if      (extndopt .eq. 4) then
!        call extend(1,strtyr,0,wrbin)
!        ocrpval = crpval
!      else if (extndopt .eq. 3 .and. newbin .eq. oldbin) then
!        call extend(1,strtyr,2,wrbin)
!        ocrpval = crpval
!      else if (extndopt .gt. 1) then
 if (extndopt .gt. 1) then
   if (oldbin .ne. newbin) then
     ! open the source file for the extend/copy
     open(unit=3,file=oldbin,form='UNFORMATTED',status='OLD', iostat=istat)
     if (istat.gt.0) call abortrun('opening extend file "'//trim(oldbin)//'"')
     if(extndopt ==3) then; call extend(3,.false.); ! read to end; no copy
     else; call extend(3,wrbin);  ! read data; copy as required by output file
     endif;
     close(unit=3)        ! close the source
   else
     call extend(1,.FALSE.)     ! no need to copy data to the original file
     time = strtyr
     if(extndopt ==3) then
       rewind(UNIT= 1, iostat=istat); ! overwrite the file data
     else
       BACKSPACE(UNIT= 1, iostat=istat) ! clear the EOF flag
     endif
     if (istat.gt.0) call abortrun('resetting file "'//trim(oldbin)//'"')
   endif

! ..... Save the single precision secndy and occlud variables read into their
! ..... double precision counterparts, cak - 09/29/2005
   secndy_double(1) = secndy(1)
   secndy_double(2) = secndy(2)
   secndy_double(3) = secndy(3)
   occlud_double = occlud
 endif

 ! read the auxiliary output files list; moved from initsw 29 May 2015
 call outfils()

! =================================================================
! initialize the simulation

! zero some layer information below the soil depth
! clear input deeper than the soil profile; these inputs can break N balance and
!  generate SWC errors. This shouldn't happen but user input isn't always consistent
   minerl(NLAYER+1:,:) =0.
   rwcf(NLAYER+1:)     =0.
   nitrate(numlyrs+1:) =0.
   swc(numlyrs+1:)     =0.

! ... Added call to initsw for Daily water budget version of Century -mdh 9/94
 !moved call from initsw
 !Note: pass structures back into C routines with the first element address
 call initlyrs(dDO_fc, swc, swdbgf);

 ! Since SWC is available in sitsoil, eliminate the float copy swcinit
 call initsw(sitlat, daylength, rwcf, adep)

! ... Initialize soil properties based on structure of Daily Soil Water
! ... Model -mdh 10/96
 call setlyrs(adep,nlayer,numlyrs, sand, silt, clay, bulkd, ph,     &
              awilt, afiel, swflag, swdbgf)
 !=======================================================
 if(swdbgf(2) .eq. 1) call prntsitpar("set") !
 if(swdbgf(2) .eq. 1) then
   write(buffr,"('NOTICE: N2Oadjust_fc=',f6.3,'    N2Oadjust_wp=',f6.3)") N2Oadjust_fc,N2Oadjust_wp
   call message(buffr)
 endif

! ... Set the initial pH value based on the value returned from the setlyrs
! ... subroutine, cak - 08/02/02
 phstart = ph

! ... Initialize the layer beyond the last one the used for safety
! ... Zero out ALL layers below nlayer. -mdh 7/27/01
 minerl(nlayer+1:CMXLYR, :) = 0.0

   ! balance and/or initialize the minerl N and ammonium nitrate double accounting
   ! this has been recast to a contained subroutine
   call balinitN()

 frac_nh4_fert = 0.5
 frac_no3_fert = 0.5

! ... Initialize the fine root pools based on the information read from
! ... the site file and the crop/tree parameterization, cak - 05/24/2007
 if (curcrp .ne. ' ') then
   call cropin(curcrp)
   cursys = CRPSYS
 endif
 if (curtre .ne. ' ') then
   call treein(curtre)
   cursys = cursys + FORSYS
 endif

 if (extndopt .le. 1) then
   mrtfr = mrtfrac
   if(bglcism(1) .eq. -1) mrtfr = bglcism(2) ! use the listed ratio
   bglcism = bglcisj * mrtfr
   bglivem = bglivej * mrtfr
   bglcisj = bglcisj - bglcism
   bglivej = bglivej - bglivem

   mrtfr = wmrtfrac
   if(frtcism(1) .eq. -1) mrtfr = frtcism(2) ! use the listed ratio
   frtcism = frtcisj * mrtfr
   frootem = frootej * mrtfr
   frtcisj = frtcisj - frtcism
   frootej = frootej - frootem
 end if
 if(swdbgf(2) .eq. 1) call prntsitext('return:')

 if (extndopt .eq. 0) call calciv(swdbgf) ! calls sumcar
 ! if (extndopt .eq. 0) then; write(*,*) 'call calciv'; call calciv(swdbgf); endif! calls sumcar

! ... Sum up isotopes
 if (extndopt .le. 1) call sumcar ! do this on a non binary extend
 ! if (extndopt .le. 1) then; write(*,*) 'call sumcar'; call sumcar; endif ! do this on a non binary extend

! ... Do preliminary initializations and calculations
 call prelim

! ... Initialize soil moisture (asmos) based on Daily Soil Water Model -mdh 10/96
 ! 5/ 4/15  KLK  eliminate swctemp. pass setasmos swc directly
 call setasmos(asmos, nlayer, swc, numlyrs, avh2o, rwcf)

 ! initialize the temperature profile
 call tempinit (tmn2m, tmx2m, swc, daylength)

 return
 contains


   subroutine commandline ()
         integer     :: nargs
     character (len=8)   aswitch ! argument switch
     integer     :: oslen, windx
     logical     :: dpars = .false.

     windx  = 0
     aswitch = '  '

    ! ... Get command line arguments
     call retarg(nargs,'arg count')
     if (nargs .eq. 0) call manual()

    ! ... Process command line arguments
    ! ... Add a new argument to indicate reading an extended <site>.100 file
    ! ... which was created from a site.nc file from a Gridded Century run
    ! ... cak - 10/05/01
     do
       if(inlin .eq. ' ') call retarg(0, inlin)
       if(dpars)  call message('retarg: '//trim(inlin)//'  aswitch: '//trim(aswitch)) ! echo the input if requested

       ! SPECIAL CASES
       ! look for debug argument command
       if(inlin .eq. 'dpars'  .or.  inlin .eq. '--dpars') then
         call message("debugging parameter parsing")
         dpars = .true.              ! set the debug logical
         inlin = " "; cycle;         ! clear the input line and reread

       else if(inlin(:2) .eq. '-c') then ! Immediately read input from command file
         if(len_trim(inlin) .gt. 2) then
           inlin = inlin(3:)           ! Appended command file name
         else
           call retarg(0, inlin)       ! load the command file name
         endif
         inquire(file=inlin,exist=goahead) ! is there a file by that name?
         ! abort if there is nothing to open
         if (.not. goahead) call abortrun('command file not found: '//trim(inlin))
         call retarg(-21, inlin)     ! have the read package open the command file
         inlin = " ";    ! clear the input line
         cycle;          ! reread
       end if

       ! determine if this the next switch
       if(aswitch .eq. ' '  .or.  aswitch .eq. '-v') then
         if(inlin .eq. '?--?') then
           if(dpars)  call message("EOF")
           exit
         elseif(inlin(:2) .eq. '--') then
           aswitch = inlin
           inlin  = ' '
         elseif(inlin(:1) .eq. '-') then
           aswitch = inlin(:2)
           inlin   = inlin(3:)
         endif
       endif

       ! These switches don't require an argument
       select case (aswitch)
         case ('-h', '-?')
           call manual()

         case ('-q')
           call stdfil('century.log')
           call retarg(1,'echo')
           aswitch = ' '; inlin   = ' '
           cycle

         ! switch to increase sitepar and soils.in debug level
         case ('--sitdbg')
           if(swdbgf(1) .eq. 0) then
              swdbgf(1) = 1
           else
              swdbgf(2) = 1
           endif
           write(buffr,*) swdbgf; call message('debugging soil water swdbgf='//trim(buffr))
           aswitch = ' '; inlin   = ' '
           cycle

         ! switch to read soils.in and sitepar after site file
         case ('--update')
           sitefirst = .true.
           call message('update input order')
           aswitch = ' '; inlin   = ' '
           cycle

         ! Is this the delimiter flag
         ! Testing here requires the delimiter to be part of the -d argument; immediately
         ! following the -d switch. The advantage is that it never eats another field
         case ('-d')
           if(len_trim(inlin) .eq. 0) then
             call message("ignoring missing delimiter")
           else
             rt = 1
             ii = 1
             if(inlin(:1) .eq. "'") then
               ii = 2
             else if(len_trim(inlin) .gt. 1) then
               call message("Warning: unusual delimiter '"//trim(inlin)//"'")
             endif
             call listvar(0,rt,'delimiter'//inlin(ii:ii),' ')
           endif
           aswitch = ' '; inlin   = ' '
           cycle

         ! switch to change the litter pools used for methane production
         case ('--ltrm')
         !   ltrmeth = .true.
         !   write(*,'(a)') "Using only litter and active SOM CO2 for methane generation"
           aswitch = ' '; inlin   = ' '
           cycle
       end select

       if(inlin .eq. ' ') cycle ! remaining switches need to look at the arguments

       ! ---- BARE WORD,  No Pending flag; treat as a schedule file
       if(aswitch .eq. ' ') then
         if(schnam .eq. '-') then
           aswitch = '-s'
         else
           call message('Skipping unknown bare word argument: '//trim(inlin))
           inlin = ''
           cycle
         endif
       endif

       ! At this point everything should be defined
       clen = len_trim(inlin)
       if(dpars)  call message('** case aswitch: '//aswitch//'  inlin: '//inlin(:clen))

       select case (aswitch)

         ! arguments no longer supported
         case ('--ntrrd','--TPL','--PST','--PFT','--PDS','--TDP','--PSF')
           call thawpulse (aswitch)

         case ('-s')  ! set schnam and check for an allowed extension 'evt' or older 'sch'
           schnam = inlin(:clen)
           istat = max(index(schnam,'.evt'),index(schnam,'.sch'))

           ! No extension;  Append .evt and check for file existance
           if (istat.eq.0) then
             schnam(clen+1:) = '.evt'
             clen = clen+4
             inquire(file=schnam,exist=goahead)
             ! File not found: Replace the evt with .sch and check again
             if (.not. goahead) then
               schnam(clen-2:) = 'sch'
             endif
           endif
           call parsfilnam('root',schnam,lschnam,srlen)

         case ('-n', '-N') ! This is a new binary file
           if((inlin(:1) .eq. '-' .or. inlin .eq. '?--?')  .and.  lschnam(1) .ne. 0) then
             ! no file name given ... use default
             newbin = schnam(lschnam(1):lschnam(2))//'.bin'
           else
             newbin = inlin ! use the input name
             inlin = ' '    ! clear inlin since it has been used
             if (index(newbin,'.bin') .eq. 0) newbin = trim(newbin)//'.bin'

             ! look for schedule root symbol in binary file name
             ii= index(newbin,'{S}')
             if(ii > 0  .and.  lschnam(1) .ne. 0) newbin = newbin(:ii-1)// &
                     schnam(lschnam(1):lschnam(2))//newbin(ii+3:)
           endif

           if(aswitch .eq. '-n') then
             inquire(file=newbin,exist=goahead)
             if (goahead) call abortrun('The new binary file "'//trim(newbin)//'" exists')
           endif

         ! switch to increase sitepar and soils.in debug level
         case ('--site')
           sitnam = inlin(:clen)
           call message('replacement site file: '//trim(sitnam))

         case ('-t') ! The text output, .lis, file name
           if((inlin(:1) .eq. '-' .or. inlin .eq. '?--?')  .and.  lschnam(1) .ne. 0) then
             ! no file name given ... use default
             ascnam = schnam(lschnam(1):lschnam(2))//'.lis'
           else
             ascnam = inlin ! use the argument
             inlin = ' '    ! clear inlin since it has been used
             if (index(ascnam,'lis') .eq. 0) ascnam = trim(ascnam)//'.lis'

             ! look for schedule root symbol in list file name
             ii= index(ascnam,'{S}')
             if(ii > 0  .and.  lschnam(1) .ne. 0) ascnam = ascnam(:ii-1)//  &
                     schnam(lschnam(1):lschnam(2))//ascnam(ii+3:)
           endif
           ! test for the file name and print overwrite warning if exists
           inquire(file=ascnam,exist=goahead)
           if (goahead) call message('Overwriting: '//trim(ascnam))

         case ('-e', '-i')
           if (extndopt .ne. 0) call abortrun('Multiple extends specified')

           extndopt = index('#eia',aswitch(2:2))

           oldbin = inlin
           if (index(oldbin,'.bin').eq.0) oldbin(clen+1:clen+4) = '.bin'
           inquire(file=oldbin,exist=goahead)
           if (.not. goahead) call abortrun('extend file: '//trim(oldbin)//' can not be found.')

         case ('-l') ! store the library file path
           plen    = clen
           filpath = inlin(:clen)
           call parsfilnam('termdir',filpath,loc,plen)

         case ('-R') ! Handle the dual argument read record index (-R)
           ! parse a record number
           lch   =0
           istat =0
           rindx = ReadInt(0,inlin,lch,istat)
           if(istat <= 0  .or.  rindx .le. 0) then
             call abortrun('Missing Read index in'//inlin(:clen))
           endif

         case ('-W') ! write a Site file
           ! parse a record number
           if(windx .eq. 0) then
             lch   =0
             istat =0
             windx = ReadInt(0,inlin,lch,istat)
             if(istat > 0  .and.  lch > len_trim(inlin)) then
               inlin = ' '
               cycle
             endif
           endif

           ! look for the file name
           call wrtsite(1,windx,inlin(:clen))
           oslen = clen !retain length as proof that we saved a name

         ! parse the variable list (-v argument)
         case ('-v')

           ! parse the variable list
           call getlist(0,inlin,outvar,0,istat)

           if(outvar.eq.0 .and. istat.gt.0) call abortrun('bad Century output variable '//inlin(:clen))
           varcnt = varcnt + outvar
           inlin = ' '
           cycle

         case ('-g')
           ext_grid = .true.
           schnam = inlin(1:clen)
           if (index(schnam,'.sch').eq.0) schnam(clen+1:clen+4) = '.sch'
           inquire(file=schnam,exist=goahead)
           if (.not. goahead) call abortrun('schedule file not found; '//trim(schnam))

         case default
           call message('Skipping Unknown argument: '//aswitch)

       end select
       aswitch = ' '
       if(inlin(:1) .ne. '-') inlin   = ' '
     end do

     !... perform any close operations
     call retarg(-9,'')

   end subroutine commandline


   subroutine thawpulse (paswitch)
       ! change saturation fraction for saturated soils during freeze/thaw pulse
       ! switch to remove respiration restraint on denitrification during thaw pulse
     character (len=8) :: paswitch ! argument switch

     INTERFACE
       subroutine ssnowtrig(tnsad, tnstd) BIND(C) ! set original freeze/thaw pulse parameters
         USE ISO_C_BINDING
         real(c_float),   VALUE, INTENT(IN) :: tnsad, tnstd
       END subroutine ssnowtrig

       subroutine sfddtrig(plmfdd, plxfdd, ptrlyr) BIND(C) ! set pulse length frozen degree day limits
         USE ISO_C_BINDING
         real(c_float),   VALUE, INTENT(IN) :: plmfdd, plxfdd; ! pulse length frozen degree limits
         integer,         VALUE, INTENT(IN) :: ptrlyr          ! soil temperature trigger layer
       END subroutine sfddtrig

       subroutine stdpsf(tnpsif) BIND(C) ! set denitrification N2/N2O pore space inflection
         USE ISO_C_BINDING
         real(c_double),  VALUE, INTENT(IN) :: tnpsif
       END subroutine stdpsf

       subroutine stdpsfdd(dpsfdm, dpsfdx) BIND(C) ! set frozen degree day shift freezing degree day limits
         USE ISO_C_BINDING
         real(c_float),   VALUE, INTENT(IN) :: dpsfdm, dpsfdx; ! denitrification pore space frozen degree limits
       END subroutine stdpsfdd

       subroutine setdrrl(tndrrl) BIND(C) ! setting denitrification respiration restraint layer
         USE ISO_C_BINDING
         integer,       INTENT(IN) :: tndrrl
       END subroutine setdrrl

       subroutine sthwplsd(tplen) BIND(C) ! setting denitrification respiration restraint layer
         USE ISO_C_BINDING
         integer,       INTENT(IN) :: tplen
       END subroutine sthwplsd

       subroutine setsatfr(tnsatf) BIND(C) ! set N2N2O ratio saturated fraction
         USE ISO_C_BINDING
         real(c_float), INTENT(IN) :: tnsatf
       END subroutine setsatfr
     END INTERFACE


     ! local variables
     ! real :: decrate, decpool(8)
     integer        tni;    ! freeze thaw pulse integer parameter
     real(c_float)  tn(3);  ! freeze thaw pulse float parameters
     real(c_double) tndp;   ! freeze thaw pulse double parameter

     select case (paswitch)
       case ('--ntrrd')
         lch = 0; istat = 0;
         tni = ReadInt(0,inlin(:clen),lch,istat);
         if(istat .lt. 1 .or. istat .eq. 3) call abortrun( &
             'please indicate the layers to disable nitrification respiration restraint')
         call setdrrl(tni);

       ! input parameters to effect thaw N2O (water filled pore space)
       case ('--TPL')
         lch = 0; istat = 0;
         tni = ReadInt(0,inlin(:clen),lch,istat);
           if(istat .ne. 1) call abortrun('reading thaw pulse length')
         call sthwplsd(tni);

       case ('--PST') ! pulse Snow trigger
         lch = 0; istat = 0;
         tn(1) = readReal(0,inlin(:clen),lch,istat);
           if(istat .lt. 1) call abortrun('reading thaw pulse trigger depth')
         istat = 0;
         tn(2) = readReal(0,inlin(:clen),lch,istat);
           if(istat .lt. 1) call abortrun('reading thaw pulse start depth')
         call ssnowtrig(tn(1), tn(2)); ! set basic pulse description

       case ('--PFT') ! Pulse freezing degree day (temperature) trigger
         lch = 0; istat = 0;
         tn(1) = readReal(0,inlin(:clen),lch,istat);
           if(istat .lt. 0) call abortrun('reading thaw frozen degree day minimum')
         istat = 0;
         tn(2) = readReal(0,inlin(:clen),lch,istat);
           if(istat .lt. 0) call abortrun('reading thaw frozen degree day maximum')
         tni   = ReadInt(0,inlin(:clen),lch,istat)
         if(istat .lt. 0) tni = 0;
         call sfddtrig(tn(1),tn(2),tni);   ! set frozen degree day pulse limiting (tnsad, tnstd)

       case ('--PDS')
         lch = 0; istat = 0;
         tndp = readReal(0,inlin(:clen),lch,istat);
           if(istat .lt. 1) call abortrun('reading thaw pulse decomposition shift fraction')
         call stdpsf(tndp);           ! set denitrification N2/N2O pore space inflection

         if(lch .lt. clen-1) then
           istat = 0;
           tn(1) = readReal(0,inlin(:clen),lch,istat);
             if(istat .lt. 0) call abortrun('reading thaw frozen degree day minimum')
           istat = 0;
           tn(2) = readReal(0,inlin(:clen),lch,istat);
           if(istat .lt. 0) call abortrun('reading thaw frozen degree day maximum')
           call stdpsfdd(tn(1),tn(2));   ! set frozen degree day inflection modification
         endif

       ! ! change decomposition rate during freeze/thaw pulse
       ! case ('--TDP')
       !   decpool = 1.0
       !   lch = 0; istat = 0;
       !   decrate = readReal(0,inlin(:clen),lch,istat);
       !   if(istat<0) call abortrun("bad freeze/thaw pulse rate"//inlin(:clen))
       !   do while (lch <= clen)
       !     tni = ReadInt(0,inlin(:clen),lch,istat)
       !     decpool(tni) = decrate
       !    if(istat<0) call abortrun("bad freeze/thaw pulse rate"//inlin(:clen))
       !   end do
       !   call sethawdec(decpool) ! set the

       ! change saturation fraction for saturated soils during freeze/thaw pulse
       case ('--PSF')
         tn(1) = 0.0;
         lch = 0; istat = 0;
         tn(1) = readReal(0,inlin(:clen),lch,istat);
         if(istat<0) call abortrun("bad pulse N2N2O mixing rate"//inlin(:clen))
         call setsatfr(tn(1)) ! set the pulse n2n2o saturated mixing rate
     end select
     paswitch = ' '; inlin = ' ';
   end subroutine thawpulse


   subroutine prntsitext(tag)
   character*(*) tag
     write(*,*) tag," aglcis        ", aglcis
     write(*,*) tag," aglive        ", aglive

     write(*,*) tag," metcis        ", metcis
     write(*,*) tag," clittr        ", clittr
     write(*,*) tag," metabc        ", metabc
     write(*,*) tag," metabe        ", metabe
     write(*,*) tag," strlig        ", strlig
     write(*,*) tag," strucc        ", strucc
     write(*,*) tag," strcis        ", strcis
     write(*,*) tag," struce        ", struce
     write(*,*) tag," wood1e        ", wood1e
     write(*,*) tag," wood2e        ", wood2e
     write(*,*) tag," wood3e        ", wood3e
     write(*,*) tag," crpstg        ", crpstg
     write(*,*) tag," forstg        ", forstg
     write(*,*) tag," bglcism",bglcism,"   bglivem",bglivem
     write(*,*) tag," bglcisj",bglcisj,"   bglivej",bglivej
     write(*,*) tag," frtcisj",frtcisj,"   frootej",frootej
     write(*,*) tag," frtcism",frtcism,"   frootem",frootem
     write(*,*) tag," ammonium      ", ammonium
     write(*,*) tag," nitrate       ", nitrate
     write(*,*) tag," swc           ", swc
   end subroutine prntsitext

   subroutine prntsitpar(tag)
   character*(*) tag
     !write(*,'(g10.4,t14,a)') dDO_sat,"/ dDO_sat"
     write(*,'(g10.4,t14,a)') dDO_fc,"/ dDO_fc"
     write(*,'(g10.4,t14,a)') dDO_wp,"/ dDO_wp"
     write(*,'(i1,   t14,a)') usexdrvrs,"/ usexdrvrs"
     write(*,'(i1,   t14,a)') texture,"/ texture"
     write(*,'(2i4,  t14,a)') jdayStart,jdayEnd,"/ jdayStart jdayEnd"
     write(*,'(i2,   t14,a)') SnowFlag,"/ SnowFlag"
     write(*,'(g10.4,t14,a)') sublimscale,"/ sublimscale"
     write(*,'(g10.4,t14,a)') reflec,"/ reflec"
     write(*,'(g10.4,t14,a)') albedo,"/ albedo"
     write(*,'(g10.4,t14,a)') dmpflux,"/ dmpflux"
     write(*,'(g10.4,t14,a)') hours_rain,"/ hours_rain"
     write(*,'(g10.4,t14,a)') hpotdeep,"/ hpotdeep"
     write(*,'(g10.4,t14,a)') ksatdeep,"/ ksatdeep"
     write(*,'(g10.4,t14,a)') rlatitude,"/ rlatitude"
     write(*,'(a,12g10.4)')   "/ sradadj ",sradadj
     write(*,'(g10.4,t14,a)') dmp,"/ dmp"
     write(*,'(g10.4,t14,a)') Ncoeff,"/ Ncoeff"
     write(*,'(g10.4,t14,a)') N2Oadjust_fc,"/ N2Oadjust_fc"
     write(*,'(g10.4,t14,a)') N2Oadjust_wp,"/ N2Oadjust_wp"
     write(*,'(g10.4,t14,a)') MaxNitAmt,"/ MaxNitAmt"
     write(*,'(g10.4,t14,a)') netmn_to_no3,"/ netmn_to_no3"
     write(*,'(g10.4,t14,a)') wfpsdnitadj,"/ wfpsdnitadj"
     write(*,'(g10.4,t12,a)') N2N2Oadj,"/ N2N2Oadj"
     write(*,'(g10.4,t12,a)') elevation,"/ elevation"
     write(*,'(g10.4,t12,a)') sitslp,"/ sitslp"
     write(*,'(g10.4,t12,a)') aspect,"/ aspect"
     write(*,'(g10.4,t12,a)') ehoriz,"/ ehoriz"
     write(*,'(g10.4,t12,a)') whoriz,"/ whoriz"
     write(*,'(i1,   t12,a)') floodN2delay,"/ floodN2delay"
     write(*,'(g10.4,t12,a)') flood_N2toN2O,"/ flood_N2toN2O"
     write(*,'(g10.4,t12,a)') Aeh,"/ Aeh"
     write(*,'(g10.4,t12,a)') Deh,"/ Deh"
     write(*,'(g10.4,t12,a)') Beh_flood,"/ Beh_flood"
     write(*,'(g10.4,t12,a)') Beh_drain,"/ Beh_drain"
     write(*,'(g10.4,t12,a)') frCH4emit,"/ frCH4emit"
     write(*,'(g10.4,t12,a)') CO2_to_CH4,"/ CO2_to_CH4"
     write(*,'(g10.4,t12,a)') frac_to_exudates,"/ frac_to_exudates"
     write(*,'(g10.4,t12,a)') zero_root_frac,"/ zero_root_frac"
   end subroutine prntsitpar


   subroutine manual ()
     integer nst
     call message(' ')
     call message('DAILYDAYCENT SOIL ORGANIC MATTER and TRACEGAS MODEL')
     call message('         N2O/Methane/orchard  - May2018')
     call message(' ')
       call retarg(-5, inlin);
       nst = scan(inlin,'\/',BACK=.TRUE.)
       if(nst > 1) inlin = inlin(nst+1:)
     call message('usage:'//trim(inlin)//' [-s] schedule-file [options]')
     call message('')
     call message(' -s     <schedule file> required (-s flag optional)')
     call message('')
     call message(' -c     <command-file> option file')
     call message(' -d<x>  .lis field delimiter (no space): -dt tab delimited -dc comma delimited')
     call message(' -g     <gridded-schedule-file> (depricated)')
     call message(' -e     <old-binary-file> to extend')
     call message(' -i     <old-binary-file> providing initial conditions')
     call message(' -l     <directory> to search for library files')
     call message(' -n     <new-binary-output-file>(optional) new file')
     call message(' -N     <binary-output-file>(optional) allows overwrite')
     call message(' -t     <list-file>(.lis)(optional) to be written')
     call message('        -[n,N,t] use schedule name root if file name is omitted')
     call message(' -R     <read index> required for binary site file (.dsa)')
     call message(' -v     <output variable list> to include in the .lis file')
     call message('           shell requires parenthesis to be escaped or in quotes')
     call message(' -W     <updated-site-file> to be written')
     call message('        the file will be ASCII if the file name ends with .100')
     call message('        <record  binary-archive.dsa> record when given an index and a file name ending with .dsa')
     call message(' --site <site file> replacing file in schedule')
     call message(' --sitdbg  debug soil layer and sitepar input echo')
     call message('                 repeat add SW input echo')
     STOP
   end subroutine manual


   subroutine balinitN()
     real(kind=REAL64) :: nsum
     real(kind=REAL64), parameter :: nfrac = 0.8


     ! reconcile soil ammonium and nitrate with the mineral N
     ! The dual book keeping and input over specifies the input and both are
     ! rounded in the text site file.

     ! Make sure we don't throw an error in the first time step

     ! there are 3 possibilities for reconciling the minerl(N) layers with
     ! the ammonium and nitrate layers
     !  1 no nitrate an ammonium
     !    this is the old/default input and the old algorithm is used
     !  2 minerl(:,N) is zeroed but there is detailed nitrate and ammonium.
     !    Sum into the correct minerl layers
     !  3 both are read
     !    use the ammonium and nitrate to distribute the minerl pools to fine layers
     !    this is done to give the required input, minerl, priority.

     ! Coded using FORTAN vector processing. Tell the compiler to optimize these
     !  calculations for multiprocessing or vector capabilities.

     ! Remove unnecessary initial N balance errors by setting the level for small
     ! mineral inputs to 0. The Century 4 input limits values cause unnecessary
     ! N balance warnings and code works with 0 initial minerl. KLK 6/2018
     !  The old code was:
     !   if(nsum == 0) then
     !     ! no detailed N input, this is the old lower limit on minerl
     !     WHERE (minerl(1:nlayer,N) .lt. 0.05) minerl(1:nlayer,N) = 0.1
     !   else
     !     ! remove negative and set it to something non-zero
     !     WHERE (minerl(1:nlayer,N) .lt. 1.0e-6) minerl(1:nlayer,N) = 1.0e-6;
     !
     nsum = sum(nitrate)+ammonium
     ! the old lower limit on minerl was (minerl(1:nlayer,N) .lt. 0.05) minerl(1:nlayer,N) = 0.1
     WHERE (minerl(1:nlayer,N) .lt. 0.0) minerl(1:nlayer,N) = 0.
     if(nsum /= 0.  .and.  ammonium /= 0.) then
       ! make sure nitrate and ammonium not negative
       WHERE (nitrate .lt. 0.0)            nitrate = 0.0;
       ammonium = max(ammonium, 0.0);
     endif

     ! initialize unread Ammonium and nitrate from the minerl variables
     nsum = sum(nitrate(ubnd(1)+1:lbnd(1)+1))+ammonium
     if(nsum == 0.) then
       ! first layer, use any ammonium info we have
       call update_npool(1, minerl(1,N), nfrac, 1.0d0 - nfrac, ammonium, nitrate, subname)
     else
        nsum                         = nsum/minerl(1,N)
        ammonium                     = ammonium*nsum
        nitrate(ubnd(1)+1:lbnd(1)+1) = nitrate(ubnd(1)+1:lbnd(1)+1)*nsum
     endif

     do ii=2,nlayer
       nsum = sum(nitrate(ubnd(ii)+1:lbnd(ii)+1))
       if(nsum == 0) then
         ! lower layers, all of the initial mineral N is nitrate
         call update_npool(ii, minerl(ii,N), 0.0_REAL64, 1.0_REAL64, ammonium, nitrate, subname)
       else
         nitrate(ubnd(ii)+1:lbnd(ii)+1) = nitrate(ubnd(ii)+1:lbnd(ii)+1)*nsum/minerl(ii,N)
       endif
     end do

     return
   end subroutine balinitN

 end subroutine detiv
