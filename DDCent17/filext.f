c****	SUBROUTINE FILEXT
       subroutine filext(fnam,ext,appnd,lenam,status)
       character fnam*(*), ext*(*)
       integer lenam, status, appnd
C
       integer iext, maxext, lext, j
C -----------------------------------------------------------------------------
C this subroutine will process file names for extensions.
C (NOTE: since UNIX places no significance in a DOS/VMS like extension, the
C        routine searches backwards for a period.) The extension is compared
C        to a list of acceptable extensions.
C If there is a match, the file type is returned.
C The first (default) extension can be either appended to the name or
C      replace the existing extension
C
C Variables
C   fnam   character (40) for the file name
C   ext    character (20) the list of acceptable extensions. extensions must
C                         begin with a '.'
C   lenam   integer        length of the returned name
C   appnd  integer        flag to control the appending the default extension
C                          -1 do not append/replace any file extensions
C                           0 append/replace if no alternate extension is found
C                           1 append/replace the found extension 
C   status integer        variable containing the return condition
C                         >=1 number of the extension found
C                           0 no extension match found
C                          -1 error conditions
C
C
C  Coded by K. Killian 7/29/94
C  modified K. Killian 8/3/94 
C          made ext passed length to keep from including garbage on SUN
C
C -----------------------------------------------------------------------------
      status =-10
C
C determine if the extension list is specified
C
      if (ext(1:1) .ne. '.') then
        call message(' illegal list of extensions')
        status =-1
        return
      endif
C
C remove any leading blanks
C
      do 20 j= 1,len(fnam)
20      if (fnam(j:j) .ne. ' ') goto 25

        call message(' missing file name')
        status =-1
        return

 25   fnam = fnam(j:)

C
C determine length of the file name and any existing extensions
C
      iext  = 0
      lext  = 0
      lenam = 0
      do 30 j = len(fnam) ,1,-1
        if (lenam.eq.0  .and.  fnam(j:j).ne.' ') lenam = j
        if (fnam(j:j) .eq. '.') then
C
C found an extension. record and determine if it's in the list
C
          lext = j
          iext = index(ext,fnam(lext:lenam))
          goto 40
        endif
30    continue
      if (lext .eq. 0) lext = lenam+1

C
C identify the number of the extension found
C
40    status = 0
      maxext = index(ext(2:),'.')
      if (maxext .eq. 0) maxext = len(ext)
      if (iext .eq. 0) then
        iext = 1
      else 
        do 45 j = 1, iext
45        if (ext(j:j) .eq. '.') status = status +1
      endif

C
C process the append/replace operation
C
      if (appnd .gt. 0 .or. ((status .eq. 0) .and. (appnd .eq. 0))) then
        fnam(lext:) = ext(:maxext)
c note lext already counts the '.' so only add the 'extension' to the length
        lenam = lext + maxext -1
      endif
c     write (*,*) ' appnd, status ',appnd, status,' fnam ',fnam(:lenam)

      return
      end
