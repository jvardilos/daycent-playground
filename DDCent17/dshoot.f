
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine dshoot(bgwfunc, tfrac, curday)
      use calflow;

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      real    bgwfunc
      real    tfrac
      integer curday

c ... Simulate death of shoots for the month.
       ! moved dosene calculation to simsom;     KLK 15 Nov 2012

c ... Local variables
      integer   iel
      real      accum(ISOS), fdeth, sdethe, tostore
      real      death_volpl

      accum(UNLABL) = 0.0
      accum(LABELD) = 0.0


c ... Death of shoots.  Most shoots die during the month of senescence.
      if (aglivc .gt. 0) then
        if(dosene) then
          fdeth = fsdeth(2)
          thermunits = 0.0
          accumdd = .false.
c ....... For the growing degree implementation a senescence event also
c ....... triggers a LAST event
          if (frtcindx .ge. 3) then
            crpgrw = 0
            dolast = .true.
            plntcnt = 0
            lastday = curday
            cgrwdys = 0
            call inprac(CRPSYS)
          endif
        else
          ! During most months, only a small fraction of the shoots die.
          fdeth = fsdeth(1) * tfrac * (1. - bgwfunc)
          ! Increase the death rate of shoots to account for effect of shading.
          ! This is not done during senescence (when death rate is >= .4)
          ! Do shading here instead of later based on the death rate. KLK 4 Nov 13
          if(aglivc .gt. fsdeth(4)) fdeth = fdeth + fsdeth(3) * tfrac
        endif

c ..... Constrain the fraction
c ..... Allow senescence to kill 100% of the shoot, cak - 02/17/04
c        if (fdeth .gt. 0.95) then
        if (fdeth .gt. 1.0) then
          fdeth = 1.0
        endif

c ..... Calculate the amounts and flow
        sdethc = aglivc * fdeth
        call csched(sdethc,aglcis(LABELD),aglivc,
     &              aglcis(UNLABL),stdcis(UNLABL),
     &              aglcis(LABELD),stdcis(LABELD),
     &              1.0,accum)
        do 10 iel = 1, nelem
          sdethe = fdeth * aglive(iel)
c ....... Use the local variable death_volpl so that volatilization that
c ....... occurs at harvest and senescence can both be tracked, see harvst,
c ....... cak - 01/02
          if (iel .eq. N) then
            death_volpl = vlossp * sdethe
            call flow(aglive(iel),esrsnk(iel),time,death_volpl)
c ......... volpl added here -mdh 8/1/00
            volpl = volpl + death_volpl
            volpla = volpla + death_volpl
            volpac = volpac + death_volpl
            sdethe = sdethe - death_volpl
          endif
          tostore = sdethe * crprtf(iel)
          call flow(aglive(iel),crpstg(iel),time,tostore)
          sdethe = sdethe * (1 - crprtf(iel))
          call flow(aglive(iel),stdede(iel),time,sdethe)
10      continue
      endif

      return
      end
