! File VersionID:
!   $Id: readswap.f90 325 2017-03-16 15:59:06Z kroes006 $
! ----------------------------------------------------------------------
      subroutine readswap
! ----------------------------------------------------------------------
!     Date               : April 2014   
!     Purpose            : read main input file .SWP
! ----------------------------------------------------------------------
      use variables
      implicit none

      integer posarg, numchar
      integer swp,i,datea(6),getun,getun2,runf,swrunon
      integer swirgfil,datefix(2),ifnd,swyrvar,swmonth
      integer swerror,swbbcfile
      integer lay,swdc,ini,tss
      integer isublay(macp),bbc,il, sol, j

      real(4) fsec
      real(8) outdate1,rainflux(30),raintime(30)
      real(8) ores(maho),osat(maho),alfa(maho),npar(maho)
      real(8) maxdepth,lexp(maho),alfaw(maho)
      real(8) dates(mabbc),gwlevel(mabbc),cseeparr(mabbc)
      real(8) haquif(mabbc),hbottom(mabbc),ttop(mabbc)
      real(8) hhtab(mabbc),qhtab(mabbc),qboti(mabbc),tbot(mabbc)
      real(8) hthr,help,term1
      real(8) headtab(matab),thetatab(matab),conductab(matab)
      real(8) dydx(matab),pondmxtb(mairg),datepmx(mairg)
      real(8) dum
      character(len=200) filenamesophy(maho)
      character(len=200) filnam,rufil
      character(len=80)  filtext
      character(len=400) swpfilnam,logfilnam
      character(len=32)  bbcfil,irgfil,swpfil,tsoilfile
      character(len=11)  tmp
      character(len=400) messag

      logical   toscr, flsat, rdinqr

      real(8) nihil,vsmall
      parameter (vsmall = 1.d-16)
      parameter (nihil = 1.d-24)

!     pressure head where interpolation is used to calculate K from VG and ksatexm
      data hthr  /-2.0d0/
! ----------------------------------------------------------------------

! --- delete existing file Swap.ok
      call delfil ('Swap.ok',.false.)

! --- write message running to screen
      write (*,'(/,a)') '  running swap ....'

! --- path and filename of executable through argument command line
      PosArg = 1
      Call Get_Command_Argument (PosArg,swpfil,NumChar)
      if (NumChar.lt.1) swpfil = 'swap'
      project = swpfil

      if(NumChar.gt.3) then
         if(swpfil(NumChar-3:NumChar).eq.'.swp' .or.                    &
     &                       swpfil(NumChar-3:NumChar).eq.'.SWP') then
            swpfilnam = trim(swpfil)
         else
            swpfilnam = trim(swpfil)//'.swp'
         end if
      else 
        swpfilnam = trim(swpfil)//'.swp'
      end if
      Numchar = len(trim(swpfilnam))
      logfilnam = swpfilnam(1:NumChar-4)//'_swap.log'

! --- open log file
      logf = getun (20,99)
      call fopens(logf,logfilnam,'new','del')

! --  write head to logfile
      filtext = 'Main input variables'
      call writehead(logf,1,swpfilnam,filtext,project)
      write (logf,*)
      
! --- open swp file
      swp = getun2 (10,90,2)
      call rdinit(swp,logf,swpfilnam)

! -   environment
      call rdscha ('project',project)
      call rdscha ('pathwork',pathwork)
      call rdscha ('pathatm',pathatm)
      call rdscha ('pathcrop',pathcrop)
      call rdscha ('pathdrain',pathdrain)
      call rdsinr ('swscre',0,2,swscre)
      call rdsinr ('swerror',0,1,swerror)
      if (swerror .eq. 1) then
        toscr = .true.
      else
        toscr = .false.
      endif
      call messini(toscr,.true.,logf)

! --  Shared simulation
      flSwapShared = .false.
      if(rdinqr('flSwapShared')) then
        call rdslog ('flSwapShared',flSwapShared)
        if(flSwapShared) then
          messag = 'Simulation with shared data exchange (flSwapShared)'
          call warn ('Readswap',messag,logf,swscre)
        endif
      endif

! -   simulation period
      call rdstim('tstart',tstart)
      call dtdpar (tstart+0.1d0, datea, fsec)
      iyear = datea(1)
      imonth = datea(2)
      call dtdpar (tstart, datea, fsec)
      call rdstim('tend',tend)

! -   check begin and end date of simulation
      if ((tend - tstart) .lt. 0.0d0) then
        messag = 'The end date of the simulation'//                     &
     &                    ' should be larger than the begin date!'
        call fatalerr ('readswap',messag)
      endif

! -   Output dates for balances
      call rdsinr ('swyrvar',0,1,swyrvar)
      if (swyrvar .eq. 0) then 
        call rdfinr('datefix',1,31,datefix,2,2)
        datea(1) = iyear
        datea(2) = datefix(2)
        datea(3) = datefix(1)
        fsec = 0.0
        call dtardp (datea, fsec, outdate1)
        if (outdate1 .lt. tstart) then
          datea(1) = datea(1) + 1
          call dtardp (datea, fsec, outdate1)
        endif
        i = 1
        outdat(i) = outdate1
        datea(1) = datea(1) + 1
        call dtardp (datea, fsec, outdate1)
        do while (outdate1 .lt. (tend + 0.1d0))
          i = i + 1
          outdat(i) = outdate1
          datea(1) = datea(1) + 1
          call dtardp (datea, fsec, outdate1)
        end do
      else
        call rdatim ('outdat',outdat,maout,ifnd)
      endif

! -   Intermediate output dates
      call rdsinr ('nprintday',1,1000,nprintday)

! --- output each time interval dt
      flprintdt = .false.
      if(rdinqr('flprintdt')) then
        call rdslog ('flprintdt',flprintdt)
      endif

      call rdsinr ('swmonth',0,1,swmonth)
      if (swmonth .eq. 1) then 
        datea(1) = iyear
        datea(2) = imonth
        if (datea(2) .lt. 12) then
          datea(2) = datea(2) + 1          
        else
          datea(1) = datea(1) + 1          
          datea(2) = 1          
        endif
        datea(3) = 1
        fsec = 0.0
        call dtardp (datea, fsec, outdate1)
        i = 0
        do while ((outdate1 - 1.d0) .lt. (tend + 0.1d0))
          i = i + 1
          outdatint(i) = outdate1 - 1.d0
          if (datea(2) .lt. 12) then
            datea(2) = datea(2) + 1          
          else
            datea(1) = datea(1) + 1          
            datea(2) = 1          
          endif
          call dtardp (datea, fsec, outdate1)
        enddo
        period = 0
        swres = 0
      else
        call rdsinr ('period',0,366,period)
        call rdsinr ('swres',0,1,swres)
      endif

      call rdsinr ('swodat',0,1,swodat)
      if (swodat.eq.1 .and. swmonth .eq.0)                              &
     &   call rdatim ('outdatint',outdatint,maout,ifnd)

! -   output files
      call rdscha ('outfil',outfil)
      call rdsinr ('swheader',0,1,swheader)

      call rdsinr ('swafo',0,3,swafo)
      call rdsinr ('swaun',0,2,swaun)
      if (swaun.ge.1 .or. swafo.ge.1) then
        call rdsinr ('swdiscrvert',0,1,swdiscrvert)
        if (swdiscrvert.eq.1) then

! -       read values for numnodNew and dzNew
          call rdsinr ('numnodNew',1,macp,numnodNew)
          call rdfdor ('dzNew',1.0d-6,5.0d2,dzNew,macp,numnodNew)
        endif
        call rdsdor ('CritDevMasBal',                                   &
     &                1.0d-30, 1.0d0,CritDevMasBal)
      endif
      call rdsinr ('swvap',0,1,swvap)
      call rdsinr ('swate',0,1,swtem)
      call rdsinr ('swblc',0,1,swblc)
      call rdsinr ('swbma',0,1,swbma) 

! --- output of recharge/storage info for ModFlow
      swoutputmodflow = 0
      if(rdinqr('swoutputmodflow')) then
        call rdsinr ('swoutputmodflow',0,1,swoutputmodflow)
        if(swoutputmodflow.eq.1) then
          messag = 'simulation with addtional output for ModFlow'
          call warn ('Readswap',messag,logf,swscre)
        endif
      endif


! -   meteo
      call rdscha ('metfil',metfil)
      call rdsdor ('lat',-90.0d0,90.0d0,lat)
      call rdsdor ('alt',-400.0d0,3000.0d0,alt)
!     optional input of angstrom coefficients a and b
      angstroma = 0.25d0
      angstromb = 0.50d0
      if(rdinqr('angstroma')) then
        call rdsdor ('angstroma',0.0d0,1.0d0,angstroma)
      endif
      if(rdinqr('angstromb')) then
        call rdsdor ('angstromb',0.0d0,1.0d0,angstromb)
      endif
      if(rdinqr('swdivide')) then
        call rdsinr ('swdivide',0,1,swdivide)
      endif

!     detailed meteo input as option
      call rdsinr ('swmetdetail',0,1,swmetdetail)
      if (swmetdetail.eq.1) then
        call rdsinr ('nmetdetail',1,96,nmetdetail)
        swrain = 0
        swetr = 0
      else
        call rdsinr ('swrain',0,3,swrain)
        call rdsinr ('swetsine',0,1,swetsine)
        call rdsinr ('swetr',0,1,swetr)
      endif

!     data only required for Penman
      if (swetr.eq.0) then
        call rdsdor ('altw',0.0d0,99.0d0,altw)
      endif

      if (swrain .eq. 1) then
        call rdador ('time',0.d0,366.d0,raintime,30,ifnd)
        call rdfdor ('rainflux',0.d0,1000.d0,rainflux,30,ifnd)
        do i = 1, ifnd
          raintab(i*2) = rainflux(i)
          raintab(i*2-1) = raintime(i)
        end do
      endif

      if (swrain .eq. 3) then
        call rdscha ('rainfil',rainfil)
      endif


! -   crop rotation scheme


! -   switch crop simulation
      swcrop = 1   ! default: with crop simulation
      if(rdinqr('swcrop')) then
        call rdsinr ('swcrop',0,1,swcrop)
      endif

      if (swcrop .eq. 1) then

!       crop calendar
        ifnd = 0
        
        call rdatim ('cropstart',cropstart,macrop,ifnd)
        call rdftim ('cropend',cropend,macrop,ifnd)
        call rdfcha ('cropname',cropname,macrop,ifnd)
        call rdfcha ('cropfil',cropfil,macrop,ifnd)
        call rdfinr ('croptype',1,3,croptype,macrop,ifnd)
!       Germination option as crop specific switch
        initcrp = 1   ! default starting date: emergence
        if(rdinqr('initcrp')) then
          call rdfinr ('initcrp',1,2,initcrp,macrop,ifnd)
        endif
        
        do i = 1,ifnd
          write(tmp,'(i11)') i
          tmp = adjustl(tmp)
          if ((cropstart(i+1)-cropend(i)).lt.0.5d0                      &
     &                                           .and. i.lt.ifnd) then
            write(tmp,'(i11)') i+1
            tmp = adjustl(tmp)
            messag = 'The begin date of crop number '//trim(tmp)//      &
     &      ' should be larger than the end date of the former crop!'
            call fatalerr ('readswap',messag)
          endif
! -       check combination detailed crop with LAT>60.
          if(croptype(i).ne.1 .and. lat.gt.60.0d0) then
            messag = 'Fatal combination for crop number '//trim(tmp)//  &
     &      'detailed crop module within polar circle (LAT>60)!'
            call fatalerr ('readswap',messag)
          endif
        end do

! ---   Nitrogen in Crop
! -     flag for nutrients in crop and soil
        flCropNut = .false.
        if(rdinqr('flCropNut')) then
          call rdslog ('flCropNut',flCropNut)
        endif

        flCropReadFile = .true.
        flCropOpenFile = .true.

      endif
!

! --- fixed irrigation events
      swirgfil = 0
      irgfil = ' ' 
      call rdsinr ('swirfix',0,1,swirfix)
      if (swirfix .eq. 1) then
        do i = 1,mairg
          irdate(i) = 0.d0
        end do
        call rdsinr ('swirgfil',0,1,swirgfil)
        if (swirgfil .eq. 0) then
          call rdatim ('irdate',irdate,mairg,ifnd)
          call rdfdor ('irdepth',0.d0,1000.d0,irdepth,mairg,ifnd)
          call rdfdor ('irconc',0.d0,1000.d0,irconc,mairg,ifnd)
          call rdfinr ('irtype',0,1,irtype,mairg,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,irdate,tend,tstart,'irdate',              &
     &                   'readswap/swirfix=1      ')
! ---   convert mm irrigation to cm
          do i = 1, ifnd
            irdepth(i) = irdepth(i) / 10.d0
          end do
        else if (swirgfil .eq. 1) then
          call rdscha ('irgfil',irgfil)
        endif
      endif

! =================================================================
!     Section soil profile

! --- initial water conditions
      call rdsinr ('swinco',1,3,swinco)
      if (swinco.eq.1) then
        call rdador ('zi',-1.d5,0.d0,zi,macp,ifnd)
        call rdfdor ('h',-1.d10,1.d4,h,macp,ifnd)
        nhead = ifnd
      elseif (swinco.eq.2) then
        call rdsdor ('gwli',-10000.0d0,1000.0d0,gwli)
      elseif (swinco.eq.3) then
!       read data from file with results from previous simulation
        call rdscha ('inifil',inifil)
      endif

! --  drainage resistance of surface runoff
      call rdsdor ('rsro',0.001d0,1.0d0,rsro)
      call rdsdor ('rsroexp',0.01d0,10.0d0,rsroexp)

! --- ponding
      swpondmx = 0
      if(rdinqr('swpondmx')) then
        call rdsinr ('swpondmx',0,1,swpondmx)
        if(swpondmx.eq.1) then
          messag = 'Simulation with time dependent ponding-threshold'
          call warn ('Readswap',messag,logf,swscre)
          call rdatim ('datepmx',datepmx,mairg,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,datepmx,tend,tstart,'datepmx ',           &
     &                      'readswap/swpondmx=1')
          call rdfdor('pondmxtb',0.0d0,1000.d0,pondmxtb,mairg,ifnd)
          do i = 1, ifnd
            pondmxtab(i*2-1) = datepmx(i)
            pondmxtab(i*2)   = pondmxtb(i)
          enddo
          pondmx = pondmxtb(1)
        endif
      endif
      if(swpondmx.eq.0) then
        call rdsdor ('pondmx',0.0d0,1000.0d0,pondmx)
      endif

! --  drainage resistance of surface runoff
      call rdsdor ('rsro',0.001d0,1.0d0,rsro)
      call rdsdor ('rsroexp',0.01d0,10.0d0,rsroexp)
! --  check combination of time-dependent ponding and runoff-resistance
      if(swpondmx.eq.1 .and. rsro.lt.1.0d-02) then
        messag = 'Fatal error: time-dependent threshold for ponding '// &
     &       ' (SWPONDMX=1) requires resistance (RSRO) > 0.01 d-1!'
        call fatalerr ('readswap',messag)
      endif

! --- soil evaporation
      call rdsinr ('swcfbs',0,1,swcfbs)
      if (swcfbs.eq.1) then
        call rdsdor ('cfbs',0.0d0,1.5d0,cfbs)
      else
        cfbs = 1.0d0
      endif
      if (swdivide.eq.1) call rdsdor ('rsoil',0.0d0,1.0d3,rsoil)
      call rdsinr ('swredu',0,2,swredu)
      if (swredu .gt.0) then
        call rdsdor ('cofred',0.0d0,1.0d0,cofred)
        call rdsdor ('rsigni',0.0d0,1.0d0,rsigni)  
      endif
      cfevappond = 1.25d0
      if(rdinqr('cfevappond')) then
        call rdsdor ('cfevappond',0.0d0,3.0d0,cfevappond)
      endif
 
! --- vertical discretization of soil profile
      call rdainr ('isoillay',1,maho,isoillay,macp,ifnd)
      call rdfinr ('isublay',1,macp,isublay,macp,ifnd)
      call rdfdor ('hsublay',0.d0,10000.d0,hsublay,macp,ifnd)
      call rdfdor ('hcomp',0.d0,1000.d0,hcomp,macp,ifnd)
      call rdfinr ('ncomp',0,macp,ncomp,macp,ifnd)
      nsublay = ifnd
      numlay = isoillay(ifnd)

! --- check input vertical discretization of soil profile
      do i = 1, nsublay
        if (abs(hsublay(i) - hcomp(i)*ncomp(i)) .gt. 1.d-6) lay = 1
        if (i .gt. 1)then
          if( (isoillay(i) .ne. isoillay(i-1) .or.                      &
     &        isoillay(i) .ne. isoillay(i-1) + 1)) lay = 1
        end if
      end do

! --- Soil Hydraulic relation: as MVG-functions or as Tables
      swsophy = 0
! -   Tables for each soil layer
      if(rdinqr('swsophy')) then
        call rdsinr ('swsophy',0,1,swsophy)
        if(swsophy.eq.1) then
!         not allowed: table-option combined with output with 
!         adjusted vertical discrectization (swdiscrvert=1)
          if(swdiscrvert.eq.1) then
            messag = 'combi of sophys-tables (swsophy=1) '//            &
     &                'and swdiscrvert=1 not allowed'
            call fatalerr ('Readswap',messag)
          endif
          if (numlay .gt. 1) then
            call rdacha ('filenamesophy',filenamesophy,maho,numlay)
          else
            call rdscha ('filenamesophy',filenamesophy(1))
          endif
          messag = 'simulation with tables for soil hydraulic relations'
          call warn ('Readswap',messag,logf,swscre)
        endif
      endif

! --- Switch/flag to reduce capillary rise in compartment below rootzone
      swcaprise = .false.
      if(rdinqr('swcaprise')) then
        call rdslog ('swcaprise',swcaprise)
      endif
      swcapriseoutput = .false.
      if(rdinqr('swcapriseoutput')) then
        call rdslog ('swcapriseoutput',swcapriseoutput)
      endif

! ----  bulk density required for oxygen stress and solute transport
      if(rdinqr('bdens')) then
        call rdfdor ('bdens',100.0d0,10000.0d0,bdens,maho,numlay)
      endif

! -   MVG-functions: parameters of functions of each soil layer
      if(swsophy.eq.0) then
        call rdfdor ('osat',0.d0,1.0d0,osat,maho,numlay)
        call rdfdor ('ores',0.d0,1.0d0,ores,maho,numlay)
        call rdfdor ('alfa',1.d-4,100.d0,alfa,maho,numlay)
        call rdfdor ('alfaw',1.d-4,100.d0,alfaw,maho,numlay)
        call rdfdor ('npar',1.001d0,9.d0,npar,maho,numlay)
        call rdfdor ('lexp',-25.d0,25d0,lexp,maho,numlay)
        call rdfdor ('h_enpr',-40.d0,0.d0,h_enpr,maho,numlay)
!       to allow downward compatibility from Swap3.2.23
        if(rdinqr('ksatfit')) then
          call rdfdor ('ksatfit',1.d-5,1.d5,ksatfit,maho,numlay)
        else
          call rdfdor ('ksat',1.d-5,1.d5,ksatfit,maho,numlay)
        endif
        flksatexm  =.false.
        if(rdinqr('ksatexm')) then
          call rdfdor ('ksatexm',1.d-5,1.d5,ksatexm,maho,numlay)
          flksatexm  =.true.
        endif

!       assign sophy-values to paramvg (input of cofgen)
        paramvg = 0.0d0
        do lay = 1,numlay
          paramvg(1,lay) = ores(lay)
          paramvg(2,lay) = osat(lay)
          paramvg(3,lay) = ksatfit(lay)
          paramvg(4,lay) = alfa(lay)
          paramvg(5,lay) = lexp(lay)
          paramvg(6,lay) = npar(lay)
          paramvg(7,lay) = 1.d0 - (1.d0 / paramvg(6,lay))
          paramvg(8,lay) = alfaw(lay)
          paramvg(9,lay) = h_enpr(lay)
          paramvg(10,lay) = -999.d0
          if (flksatexm) then
            paramvg(10,lay) = ksatexm(lay)
            help = abs(hthr * paramvg(4,lay))**paramvg(6,lay)
            help = (1.0d0 + help) ** paramvg(7,lay)
            relsatthr(lay) = 1.0d0/help
            term1 = ( 1.0d0-relsatthr(lay)**(1.0/paramvg(7,lay)) )      &
     &                                                  **paramvg(7,lay)
            ksatthr(lay)=paramvg(3,lay)*(relsatthr(lay)**paramvg(5,lay))&
     &                                     * (1.0d0-term1)*(1.0d0-term1)
            if(ksatthr(lay).ge.ksatexm(lay)) then
              write(messag,'(a,i2)') 'ksatexm < ksatthr for layer ',lay
              call fatalerr ('readswap',messag)
            endif
          endif
        end do
      endif

! --- hysteresis
      call rdsinr ('swhyst',0,2,swhyst)
      if (swhyst.gt.0) then
        call rdsdor ('tau',0.0d0,1.0d0,tau)
      endif

!     combination of soilphysical tables and hysteresis is not possible
      if (swhyst.eq.1 .and. swsophy.eq.1) then
         messag = 'Combination of hysteresis (swhyst=1) and tabulated ' &
     &//'soil physics (swsophy=1) is not possible !'
          call fatalerr('Readswap',messag)
      end if

! --- rooting depth limitation
      call rdsdor ('rds',1.d0,5000.0d0,rdmax)

! --- preferential flow due to macropores
      flInitDraBas = .false.
      call rdsinr ('swmacro',0,1,swmacro)
      if (swmacro.eq.1) then
!- a. PARAMETERS FOR GEOMETRY MACROPORES
         MaxDepth= 0.d0
         do il = 1, nsublay
            MaxDepth= MaxDepth - hsublay(il)
         enddo
         call rdsdor('Z_Ah',maxdepth,0.d0,Z_Ah)
         call rdsdor('Z_Ic',maxdepth,0.d0,Z_Ic)
         call rdsdor('Z_St',maxdepth,0.d0,Z_St)
         call rdsdor('VlMpStSs',0.0d0,0.5d0,VlMpStSs)
         call rdsdor('PpIcSs',0.0d0,0.99d0,PpIcSs)
         call rdsinr('NumSbDm',0,(MaDm-2),NumSbDm)
         if (NumSbDm.eq.0 .and. PpIcSs.gt.0.d0) then
            messag = ' NumSbDm .eq.0 .and. PpIcSs.gt.0.d0'
            call fatalerr('MacroRead',messag)
         endif
         call rdsdor('DiPoMi',0.1d0,1000.0d0,DiPoMi)
         call rdsdor('DiPoMa',0.1d0,1000.0d0,DiPoMa)
         call rdsdor('PndmxMp',0.0d0,10.0d0,PndmxMp)    
!   - optional: 
         call rdsdor('PowM',0.0d0,100.0d0,PowM)    ! default 1.0
         call rdsdor('Rzah',0.0d0,1.0d0,Rzah)      ! default 0.0
         call rdsdor('Spoint',0.0d0,1.0d0,Spoint)  ! default 1.0
         call rdsinr('SwPowM',0,1,SwPowM)          ! default 0
         call rdsdor('PndmxMp',0.0d0,10.0d0,PndmxMp)    
         Z_MB50 = 0.5d0*(Z_Ic+Z_St)
         if(rdinqr('Z_MB50')) then
           call rdsdor('Z_MB50',Z_St,Z_Ic,Z_MB50)
         endif
!- b. PARAMETERS FOR SHRINKAGE CHARACTERISTICS
         call rdfinr('SwSoilShr',0,2,SwSoilShr,MaHo,numlay)
         call rdfinr('SwShrInp',1,3,SwShrInp,MaHo,numlay)
         call rdfdor('ThetCrMP',0.0d0,1.0d0,ThetCrMp,MaHo,numlay)
         call rdfdor('GeomFac',0.0d0,10.0d0,GeomFac,MaHo,numlay)
         call rdfdor('ShrParA',0.0d0,10.0d0,ShrParA,MaHo,numlay)
         call rdfdor('ShrParB',-10.0d0,100.0d0,ShrParB,MaHo,numlay)
         call rdfdor('ShrParC',-10.0d0,100.0d0,ShrParC,MaHo,numlay)
         call rdfdor('ShrParD',0.0d0,100.0d0,ShrParD,MaHo,numlay)
         call rdfdor('ShrParE',-10.0d0,10.0d0,ShrParE,MaHo,numlay)
         call rdsdor('ZnCrAr',Z_Ah,0.0d0,ZnCrAr)
!- c. PARAMETERS FOR SORPTIVITY
         call rdfinr('SwSorp',1,2,SwSorp,MaHo,numlay)
         call rdfdor('SorpFacParl',0.0d0,1.0d2,SorpFacParl,MaHo,numlay)
         call rdfdor('SorpMax',0.0d0,100.0d0,SorpMax,MaHo,numlay)
         call rdfdor('SorpAlfa',-10.0d0,10.0d0,SorpAlfa,MaHo,numlay)
!- d. SHAPE FACTOR for saturated exchange between macropores and matrix
         call rdsdor('ShapeFacMp',0.0d0,100.0d0,ShapeFacMp)
!- e. CRITICAL value for undersaturation volume
         call rdsdor('CritUndSatVol',0.0d0,10.0d0,CritUndSatVol)
         call rdsinr('SwDarcy',0,1,SwDarcy)
!- f. PARAMETERS FOR RAPID DRAINAGE
!       only possible when at least one drainege level!!!
         call rdsinr('SwDrRap',0,1,SwDrRap)
         if (SwDrRap.eq.1) then
!          at this moment for 1 system only, may be extended to other systems 
            call rdsdor ('RapDraResRef',0.d0,1.d10,RapDraResRef(1))
            do il= 2,madr
               RapDraResRef(il) = RapDraResRef(1)
            enddo
            call rdsdor('RapDrareaExp',0.0d0,100.0d0,RapDraReaExp)
            call rdsinr('NumLevRapDra',1,5,NumLevRapDra)
!
!-    set flag for initialisation of drainage basis rapid drainage
            flInitDraBas = .true.
         endif
      endif

! -   snow and frost conditions 
      call rdsinr ('swsnow',0,1,swsnow)
      if (swsnow.eq.1) then
        call rdsdor ('snowinco',0.0d0,1000.0d0,snowinco)
        call rdsdor ('teprrain',0.0d0,10.0d0,teprrain)
        call rdsdor ('teprsnow',-10.0d0,0.0d0,teprsnow)
        call rdsdor ('snowcoef',0.0d0,10.0d0,snowcoef)
        swsublim = 0
        if(rdinqr('swsublim')) then
          call rdsinr ('swsublim',0,1,swsublim)
        endif
      endif
      call rdsinr ('swfrost',0,1,swfrost)
      if (swfrost.eq.1) then
        call rdsdor ('tfroststa',-10.0d0,5.0d0,tfroststa)
        call rdsdor ('tfrostend',-10.0d0,5.0d0,tfrostend)
      endif
!     combination of frost and macropore-flow is not possible (yet)
      if (swfrost.eq.1 .and. swmacro.eq.1) then
         messag = 'Combination of frost (swfrost=1) and macropore-flow '&
     &//'(swmacro=1) is not operational !'
          call fatalerr('Readswap',messag)
      end if

!     combination of snow and Et variation during the day (SwEtSine)
      if (swmetdetail.eq.1 .and. swsnow.eq.1) then
        messag = 'In case of snow the sublimation is not simulated'//   &
     & 'correctly if short time meteorological records are used'//      &
     & '(swmetdetail=1 and swsnow=1), please adapt input !'
        call warn ('Readswap',messag,logf,swscre)
      end if

! --- parameters numerical scheme
      call rdsdor ('dtmin', 1.0d-7,0.1d0,dtmin)
      call rdsdor ('dtmax', dtmin, 1.0d0,dtmax)
      call rdsdor ('gwlconv',1.0d-5,1000.0d0,gwlconv)
      call rdsdor ('CritDevPondDt',1.0d-6,1.0d-01,CritDevPondDt)
! Convergence criteria (optional input)
      CritDevh1Cp = 1.0d-2
      if(rdinqr('CritDevh1Cp')) then
        call rdsdor ('CritDevh1Cp',1.0d-10,1.0d3,CritDevh1Cp)
      endif
      CritDevh2Cp = 1.0d-1
      if(rdinqr('CritDevh2Cp')) then
        call rdsdor ('CritDevh2Cp',1.0d-10,1.0d3,CritDevh2Cp)
      endif
! flag to generate additional output about convergence-warnings from subr Headcalc
      fldumpconvcrit = .false.
      if(rdinqr('fldumpconvcrit')) then
        call rdslog ('fldumpconvcrit',fldumpconvcrit)
      endif
! Maximum number of iterations [5,100]
      call rdsinr ('MaxIt',5,100,MaxIt)
! Maximum number of back track cycles within an iteration cycle [1,10]
      call rdsinr ('MaxBackTr',1,10,MaxBackTr)
! Maximum number of iterations: no input
      msteps = 100000000
! Switch for mean of hydraulic conductivity 
!  SwkMean=1,2:unweighted arithmic mean,weighted arithmic mean
!  SwkMean=3,4:unweighted geometric mean, weighted geometric mean
      call rdsinr ('SWkmean',1,4,SWkmean)
! Switch for implicit solution with hydraulic conductivity: 0 = explicit, 1 = implicit
      call rdsinr ('SwkImpl',0,1,SwkImpl)
! Maximum cputime, introduced to be able to interrupt (near) endless iterations
      flMaxIterTime = .false.
      if(rdinqr('flMaxIterTime')) then
        call rdslog ('flMaxIterTime',flMaxIterTime)
        MaxIterTime = 2419200                ! 4 weeks = 60*60*24*7*4=2419200 secs)
        if(flMaxIterTime) then
          call rdsinr ('MaxIterTime',1,2419200,MaxIterTime)
        endif
      endif

! =================================================================
!     Lateral drainage section
!     extended or basic drainage
      call rdsinr ('swdra',0,2,swdra)
      if (swdra .ne. 0) call rdscha ('drfil',drfil)
      if (SwDrRap.eq.1 .and. SwDra.eq.0) then
          messag = ' There are no drainage levels, so rapid drainage'// &
     &'is not possible !'
          call fatalerr('MacroRead',messag)
      endif

!     runon from external source (field)
      call rdsinr ('swrunon',0,1,swrunon)
      if (swrunon .eq. 1) then
        call rdscha ('rufil',rufil)
      endif

!     output-options drainage fluxes, surface reservoir 
      if (swdra.eq.2) then
        call rdsinr ('swdrf',0,1,swdrf)
        call rdsinr ('swswb',0,1,swswb)
      endif


! =================================================================
!     Bottom boundary section

! --- Initialise
      do i = 1,2*mabbc
        gwltab(i) = 0.0D0
        hbotab(i) = 0.0D0
        qbotab(i) = 0.0D0
        haqtab(i) = 0.0D0
        cseeptab(i) = 0.0D0
      end do

! --- option for input of bottom boundary condition
      call rdsinr ('swbbcfile',0,1,swbbcfile)
      if (swbbcfile .eq. 1) then
        call rdscha ('bbcfil',bbcfil)
      endif

! =================================================================
!     Section heat flow

! --- Switch whether simulation includes heat simulation or not
      call rdsinr ('swhea',0,1,swhea)

      if (swhea .eq. 1) then
! ---   analytical or numerical method 
        call rdsinr ('swcalt',1,2,swcalt)

        if (swcalt.eq.1) then
          call rdsdor ('tampli',0.0d0,50.0d0,tampli)
          call rdsdor ('tmean',-10.0d0,30.0d0,tmean)
          call rdsdor ('timref',0.0d0,366.0d0,timref)
          call rdsdor ('ddamp',1.0d0,500.0d0,ddamp)
        else
          call rdador ('zh',-1.0d5,0.0d0,zh,macp,ifnd)
          call rdfdor ('tsoil',-50.0d0,50.0d0,tsoil,macp,ifnd)
          call rdfdor ('psand',0.0d0,1.0d0,psand,maho,numlay)
          call rdfdor ('psilt',0.0d0,1.0d0,psilt,maho,numlay)
          call rdfdor ('pclay',0.0d0,1.0d0,pclay,maho,numlay)
          call rdfdor ('orgmat',0.0d0,1.0d0,orgmat,maho,numlay)
          nheat = ifnd

!   -     top boundary temperature
          call rdsinr ('SwTopbHea',1,2,SwTopbHea)
          Tsoilfile = ' '
          if (swtopbhea .eq. 2) then
            call rdscha ('TSoilFile',tsoilfile)
          endif

!   -     bottom boundary temperature
          call rdsinr ('SwBotbHea',1,2,SwBotbHea)
          if (SwBotbHea.eq.2) then
             call rdatim ('datet',dates,mabbc,ifnd)
! -     at least one date must be within simulation period
             call checkdate(ifnd,dates,tend,tstart,'datet ',            &
     &                      'readswap/swbotbhea=2    ')
             call rdfdor('tbot',-50.0d0,50.d0,tbot,mabbc,ifnd)
!
             do i = 1, ifnd
               tembtab(i*2)   = tbot(i)
               tembtab(i*2-1) = dates(i)
             enddo
!
           endif
!
        endif
      endif

! -   fatal error when frost or snow is simulated without heat flow
      if ((swfrost.eq.1 .or. swsnow.eq.1) .and. swhea.eq.0) then
        messag = 'In case of snow or frost the soil heat flow should '//&
     &  'be simulated! Adapt .swp input file.'
        call fatalerr ('readswap',messag)
      endif

! =================================================================
!     Section solute transport

! --- Switch whether simulation includes solute transport or not
      call rdsinr ('swsolu',0,1,swsolu)

      if (swsolu .eq. 1) then

! ---   top boundary and initial condition
        call rdsdor ('cpre',0.0d0,100.0d0,cpre)
        call rdsdor ('cdrain',0.0d0,100.0d0,cdrain)
!       bottom condition with 3 options:
!       SWBOTBC = 0: lateral (cdrain as input), bottom set to lateral (cseep=cdrain)
!       SWBOTBC = 1: separate input of bottom conc as cseep
!       SWBOTBC = 2: separate input of dynamic bottom conc as (datec,cseeparr)
        swbotbc = 0
        if (rdinqr('swbotbc')) then
          call rdsinr ('swbotbc',0,2,swbotbc)
        endif
        if (swbotbc.eq.0) then
          cseep = cdrain
        elseif (swbotbc.eq.1) then
          call rdsdor ('cseep',0.0d0,100.0d0,cseep)
        elseif (swbotbc.eq.2) then
          call rdatim ('datec',dates,mabbc,ifnd)
          call rdfdor ('cseeparr',0.0d0,100.0d0,cseeparr,mabbc,ifnd)
          do i = 1, ifnd
            cseeptab(i*2) = cseeparr(i)
            cseeptab(i*2-1) = dates(i)
          end do
        endif

        call rdador ('zc',-1.0d5,0.0d0,zc,macp,ifnd)
        call rdfdor ('cml',0.0d0,1000.0d0,cml,macp,ifnd)
        nconc = ifnd

! ---   diffusion, dispersion and solute uptake by roots
        call rdsdor ('ddif',0.0d0,10.0d0,ddif)
        call rdfdor ('ldis',0.0d0,100.0d0,ldis,maho,numlay) 
        call rdsdor ('tscf',0.0d0,10.0d0,tscf)

! ---   sorption
        call rdsinr ('swsp',0,1,swsp)
        if (swsp .eq. 1) then
          call rdsdor ('frexp',0.0d0,   10.0d0,frexp)  
          call rdsdor ('cref', 0.0d0, 1000.0d0,cref)
          call rdfdor ('kf',   0.0d0,10000.0d0,kf,maho,numlay)
          call rdfdor ('bdens',100.0d0,10000.0d0,bdens,maho,numlay)
        else
          cref  = 1.0d0
          frexp = 1.0d0
          do lay = 1,numlay
            kf(lay) = 0.0d0
          end do
        endif

! ---   decomposition
        call rdsinr ('swdc',0,1,swdc)
        if (swdc .eq. 1) then
          call rdfdor ('decpot',0.0d0,10.0d0,decpot,maho,numlay)
          call rdsdor ('gampar',0.0d0,0.5d0,gampar)
          call rdsdor ('rtheta',0.0d0,0.4d0,rtheta)
          call rdsdor ('bexp',0.0d0,2.0d0,bexp)
          call rdfdor ('fdepth',0.0d0,1.0d0,fdepth,maho,numlay)
        else 
          gampar = 0.0d0
          rtheta = 0.5d0
          bexp = 0.0d0
          do lay = 1,numlay
            decpot(lay) = 0.0d0
            fdepth(lay) = 0.0d0
          end do
        endif

! ---   breakthrough
        call rdsinr ('swbr',0,1,swbr)
        if (swbr .eq. 0) then
          daquif = 100.0d0
          poros = 1.0d0
          kfsat = 0.0d0
          decsat = 0.0d0
        else
          call rdsdor ('daquif',0.0d0,10000.0d0,daquif)
          call rdsdor ('poros',0.0d0,0.6d0,poros)
          call rdsdor ('kfsat',0.0d0,100.0d0,kfsat)
          call rdsdor ('decsat',0.0d0,10.0d0,decsat)
          call rdsdor ('cdraini',0.0d0,100.0d0,cdrain)
          cseep = cdrain
        endif

      endif

! =================================================================
!     Section Ageing according to Goode (1996): 
!      "Direct simulation of groundwater age, WRR vol.32, p 289-296"
      flAgeTracer = .false.
      if (rdinqr('flAgeTracer')) then
        call rdslog ('flAgeTracer',flAgeTracer)
        if (flAgeTracer) then

! ---     top and bottom boundary and initial condition
          cpre = 0.0d0
          cirr = 0.0d0
          cdrain = 0.0d0
          cseep = 0.0d0
          call rdador ('zc',-1.0d5,0.0d0,zc,macp,ifnd)
          call rdfdor ('cml',0.0d0,1000.0d0,cml,macp,ifnd)
          nconc = ifnd

! ---     diffusion, dispersion and uniform solute uptake by roots
          call rdsdor ('ddif',0.0d0,10.0d0,ddif)
          call rdfdor ('ldis',0.0d0,100.0d0,ldis,maho,numlay) 
          tscf = 1.0d0

!         warnings and errors
          messag = 'simulation with Age Tracer option (flAgeTracer)'
          call warn ('Readswap',messag,logf,swscre)
          if(swSolu.eq.1) then
            write(messag,'(3a)')                                        &
     &      'Combination of solute transport and Ageing is not allowed',&
     &      '(flAgeTracer=.true.  AND   swSolu=1 is not allowed !' 
            call fatalerr ('Readswap',messag)
          endif
          if(swSnow.eq.1) then
            write(messag,'(3a)')                                        &
     &       'Combination of snow and Ageing is not allowed',           &
     &       '(flAgeTracer=.true.  AND   swSnow=1 is not allowed !'
            call fatalerr ('Readswap',messag)
          endif
        endif
      endif


! =================================================================
!     Section bottom boundary condition

      if (swbbcfile .eq. 1) then
        close (swp)
        bbc = getun2 (10,90,2)
        filnam = trim(pathwork)//trim(bbcfil)//'.bbc'
        call rdinit(bbc,logf,filnam)
      endif

! --- option for bottom boundary condition
      call rdsinr ('swbotb',1,8,swbotb)

! --- given groundwaterlevel
      if (swbotb.eq.1) then
        call rdatim ('date1',dates,mabbc,ifnd)
        call rdfdor ('gwlevel',-10000.0d0,1000.d0,gwlevel,mabbc,ifnd)

! -     check: at least one date must be within simulation period
        call checkdate(ifnd,dates,tend,tstart,'date1 ',                 &
     &                      'readswap//swbotb=1      ')
! -     check if gwlevel is above soil surface, eliminate certain combinations
        flsat = .false.
        do i = 1,ifnd
          if (gwlevel(i) .gt. -0.5d0*(hcomp(1)))  flsat = .true.
        enddo
        if(flsat) then
           if (swKimpl.eq.1) then
              write(messag,'(3a)')                                      &
     &         'Groundwaterlevel above soil surface, combined with ',   &
     &         'Implicit Conductivity (swbotb=1 AND swKimpl=1 AND  ',   &
     &         'gwl>z(1))  is not allowed !'
              call fatalerr ('Readswap',messag)
           endif
           if (swmacro.eq.1) then
              write(messag,'(3a)')                                      &
     &         'Groundwaterlevel above soil surface, combined with ',   &
     &         'MacroPore flow (swbotb=1 AND swMacro=1 AND gwl>z(1)) ', &
     &         'is not allowed !' 
              call fatalerr ('Readswap',messag)
           endif
           if (swfrost.eq.1) then
              write(messag,'(3a)')                                      &
     &         'Groundwaterlevel above soil surface, combined with ',   &
     &         'frost conditions (swbotb=1 AND swforst=1 AND gwl>z(1))',&
     &         ' not well tested, be aware of balance errors !' 
              call warn ('Readswap',messag,logf,swscre)
           endif
        endif

! -     store values in gwltab
        do i = 1,ifnd
          gwltab(i*2) = gwlevel(i) 
          gwltab(i*2-1) = dates(i)
        enddo
      endif

! --- regional bottom flux is given                             
      if (swbotb.eq.2) then
        call rdsinr ('sw2',1,2,sw2)
        if (sw2 .eq. 1) then
          call rdsdor ('sinave',-10.0d0,10.0d0,sinave)
          call rdsdor ('sinamp',-10.0d0,10.0d0,sinamp)
          call rdsdor ('sinmax',0.d0,366.d0,sinmax)
        else
! -       read tabular data
          call rdatim ('date2',dates,mabbc,ifnd)
          call rdfdor ('qbot2',-100.0d0,100.0d0,qboti,mabbc,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,dates,tend,tstart,'date2 ',               &
     &                      'readswap//swbotb=2      ')
! -       fill qbotab table
          do i = 1,ifnd
            qbotab(i*2) = qboti(i) 
            qbotab(i*2-1) = dates(i)
          enddo
        endif
      endif

! --- calculated flux through the bottom of the profile
      if (swbotb.eq.3) then

!       Switch for implicit solution with lower boundary option 3 (Cauchy): 0 = explicit, 1 = implicit
        call rdsinr ('swbotb3Impl',0,1,swbotb3Impl)

        call rdsdor ('shape',0.0d0,1.0d0,shape)
        if (swbotb3Impl.eq.1 .and. abs(shape-1.0d0).gt.1.0d-7) then
          write(messag,'(a)')                                           &
     &    ' Possible lower boundary inconsistency using SwBotb=3: ',    &
     &    '  Combination of swbotb3Impl AND shape not equal 1.0',       &
     &    '  This is not recommended. Suggestion is: swbotb3Impl=0'
          call warn ('readswap',messag,logf,swscre)
        endif

        call rdsdor ('hdrain',-1.0d4,0.0d0,hdrain)
        call rdsdor ('rimlay',0.0d0,1.0d5,rimlay)

!       Switch to suppress addition of vertical resistance 
!                      between bottom of model and groundwater level
        call rdsinr ('SwBotb3ResVert ',0,1,SwBotb3ResVert)

        call rdsinr ('sw3',1,2,sw3)
        if (sw3 .eq. 1) then
          call rdsdor ('aqave',-10000.0d0,1000.0d0,aqave)
          call rdsdor ('aqamp',0.0d0,1000.0d0, aqamp)
          call rdsdor ('aqtmax',0.0d0,366.d0,aqtmax)
          call rdsdor ('aqper',0.0d0,366.0d0,aqper)
        else
! -       read tabular data
          call rdatim ('date3',dates,mabbc,ifnd)
          call rdfdor ('haquif',-10000.0d0,1000.d0,haquif,mabbc,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,dates,tend,tstart,'date3 ',               &
     &                   'readswap//swbotb=3      ')
! -       fill haqtab table
          do i = 1,ifnd
            haqtab(i*2) = haquif(i) 
            haqtab(i*2-1) = dates(i)
          enddo
        endif
        call rdsinr ('sw4',0,1,sw4)
        if (sw4 .eq. 1) then
          if (swbotb3Impl.eq.1) then
            write(messag,'(3a)')                                        &
     &       ' Implicit solution of Cauchy, combined with fluxes ',     &
     &       'is active !(swbotb3Impl=1 AND sw4=1)'
            call warn ('readswap',messag,logf,swscre)
          endif
! -       read tabular data
          call rdatim ('date4',dates,mabbc,ifnd)
          call rdfdor ('qbot4',-100.0d0,100.d0,qboti,mabbc,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,dates,tend,tstart,'date4 ',               &
     &                   'readswap//swbotb=3/sw4=1')
! -       fill qbotab table
          do i = 1,ifnd
            qbotab(i*2) = qboti(i) 
            qbotab(i*2-1) = dates(i)
          enddo
        endif
      endif

! --- flux-groundwater level relationship 
      if (swbotb.eq.4) then
        call rdsinr ('swqhbot',1,2,swqhbot)
        if (swqhbot.eq.1) then
          call rdsdor ('cofqha',-100.0d0,100.0d0,cofqha)
          call rdsdor ('cofqhb',-1.0d0,1.0d0,cofqhb)
          swcofqhc = 0
          cofqhc = 0.0d0
          if (rdinqr('cofqhc')) then
            swcofqhc = 1
            call rdsdor ('cofqhc',-10.0d0,10.0d0,cofqhc)
          endif
        else if (swqhbot.eq.2) then
          call rdador ('qtab',-100.0d0,100.d0,qhtab,mabbc,ifnd)
          call rdfdor ('htab', -1.0d4, 0.0d0, hhtab,mabbc,ifnd)
          do i = 1,ifnd
            qbotab(i*2) = qhtab(i) 
            qbotab(i*2-1) = abs(hhtab(i))
          enddo
        endif
      endif

! --- pressure head of lowest compartment is given 
      if (swbotb.eq.5) then
        call rdatim ('date5',dates,mabbc,ifnd)
        call rdfdor('hbot5',-1.0d10,1000.d0,hbottom,mabbc,ifnd)
! -     at least one date must be within simulation period
        call checkdate(ifnd,dates,tend,tstart,'date5 ',                 &
     &                 'readswap//swbotb=5      ')
! -     store pressure head values in hbotab
        do i = 1, ifnd
          hbotab(i*2) = hbottom(i) 
          hbotab(i*2-1) = dates(i)
        end do
      endif

      if (swbbcfile .eq. 0) then
        close (swp)
      elseif (swbbcfile .eq. 1) then
        close (bbc)
      endif

! -   Tables for each soil layer
      if(swsophy.eq.1) then
          ientrytablay = 0
          do lay = 1,numlay
            filnam = filenamesophy(lay)
            sol = getun2 (50,90,2)
            call rdinit(sol,logf,filnam)
            call rdador ('headtab',-1.0d09,1.0d09,headtab,matab,ifnd)
            call rdfdor ('thetatab',0.0d0,1.0d0,thetatab,matab,ifnd)
            call rdfdor                                                 &
     &                ('conductab',0.0d0,1000.0d0,conductab,matab,ifnd)
            close (sol)
            numtablay(lay) = ifnd

!           verify incremental sequence of values
            do i = 2,numtablay(lay)
              if((headtab(i)-headtab(i-1)) .le. nihil) then
                messag = ' No incremental values for head soil physic'  &
     &          //'in input file '//trim(filnam)
              call fatalerr('Readswap',messag)
              endif
              if((thetatab(i)-thetatab(i-1)) .le. vsmall) then
                messag = ' No incremental values for theta soil physic' &
     &          //'in input file '//trim(filnam)
              call fatalerr('Readswap',messag)
              endif
              if((conductab(i)-conductab(i-1)) .le. nihil) then
                messag = ' No incremental values for conduc soil physic'&
     &          //'in input file '//trim(filnam)
              call fatalerr('Readswap',messag)
              endif
            enddo

!           preprocess table
!---- sorting of arrays to an ascending sequence

            Do i=1,numtablay(lay)-1
               Do j=i+1,numtablay(lay)
                  If(headtab(i) .gt. headtab(j))Then
                     Dum          = headtab(j)
                     headtab(j)   = headtab(i)
                     headtab(i)   = Dum
                     Dum          = thetatab(j)
                     thetatab(j)  = thetatab(i)
                     thetatab(i)  = Dum
                     Dum          = conductab(j)
                     conductab(j) = conductab(i)
                     conductab(i) = Dum
                  End If
               End Do
            End Do
            ientrytablay(lay,0)=numtablay(lay)
            do i = 1,numtablay(lay)
              if(headtab(i).gt.-1.0d-1)then
                 j=0
              else
                 j = int(1000*(log10(-headtab(i))+1.d0))+1
              end if
              ientrytablay(lay,j) = i
              sptablay(1,lay,i) = headtab(i)
              sptablay(2,lay,i) = thetatab(i)
              sptablay(3,lay,i) = conductab(i)
            enddo
            ientrytablay(lay,1) = 0
            do j=matabentries-1,1,-1
               if(ientrytablay(lay,j) .eq. 0)                           &
     &            ientrytablay(lay,j) = ientrytablay(lay,j+1)
            end do
            call PreProcTabulatedFunction(1,                            &
     &                             numtablay(lay),headtab,thetatab,dydx)
            do i = 1,numtablay(lay)
              sptablay(4,lay,i) = dydx(i)
            enddo
            call PreProcTabulatedFunction(2,                            &
     &                           numtablay(lay),headtab,conductab,dydx)
            do i = 1,numtablay(lay)
              sptablay(5,lay,i) = dydx(i)
            enddo
!           tables must have values for a head=0
            if (sptablay(1,lay,numtablay(lay)) .lt. -1.d-20) then
               messag = ' No values for head=0 in tabulated soil physic'&
     &          //'in input file '//trim(filnam)
              call fatalerr('Readswap',messag)
            end if
          enddo
      endif


! =================================================================
!     Read data of drainage input file
      if (swdra.eq.1)                                                   &
     &    call rddrb (swallo,swdtyp,swmacro,dramet,dra,ipos,nrlevs,     &
     &    logf,swdivd,l,zbotdr,owltab,drares,infres,qdrtab,basegw,      &
     &    wetper,khtop,khbot,kvtop,kvbot,zintf,entres,geofac,drfil,     &
     &    pathdrain,numlay,cofani,swnrsrf,swscre,cofintfl,expintfl,     &
     &    swdislay,swtopdislay,ztopdislay,ftopdislay,shape,SwTopnrsrf,  &
     &    swdivdinf,FacDpthInf)

! =================================================================
! -   read data of fixed irrigation from separate file
      if (swirgfil .eq. 1 .and. swirfix .eq. 1) then
        filnam = trim(pathwork)//trim(irgfil)//'.irg'
        irg = getun2 (10,90,2)
        call rdinit(irg,logf,filnam)
        call rdatim ('irdate',irdate,mairg,ifnd)
        call rdfdor ('irdepth',0.d0,1000.d0,irdepth,mairg,ifnd)
        call rdfdor ('irconc',0.d0,1000.d0,irconc,mairg,ifnd)
        call rdfinr ('irtype',0,1,irtype,mairg,ifnd)
        close (irg)
! -     at least one date must be within simulation period
        call checkdate(ifnd,irdate,tend,tstart,'irdate',                &
     &                 'readswap//swirgfil=1     ')
! ---   convert mm irrigation to cm
        do i = 1, ifnd
          irdepth(i) = irdepth(i) / 10.d0
        end do
      endif

! --- read data from file with results from previous simulation
! --- set initial ponding conditions
      if (swinco.eq.3) then
        ini = getun2 (10,90,2)
        call rdinit(ini,logf,inifil)
!        write (messag,107)                                              &
!     &          '*  I/O of variables from file (SWINCO=3)',             &
!     &          '   should be considered carefully',                    &
!     &          '   this option is incomplete and poorly tested!'
 107    format (3a)
!        call warn ('readswap',messag,logf,swscre)

        call rdsdor ('ssnow',0.0d0,1000.0d0,ssnow)
        if (swsnow.ne.1) then
          write (messag,107)                                            &
     &         'No Simulation of snow, therefore : Initial',            &
     &         '  storage of snow set to 0.0',                          &
     &         '  neglecting ssnow-value from file '//trim(inifil)
          call warn ('readswap',messag,logf,swscre)
          ssnow = 0.0d0
        endif
        call rdsdor ('slw',0.0d0,1000.0d0,slw)
        call rdsdor ('pond',0.0d0,100.0d0,pond)
        call rdador ('z_h',-1.0d5,0.0d0,zi,macp,ifnd)
        call rdfdor ('h',-1.0d10,1.0d4,h,macp,ifnd)
        nhead = ifnd
        if (swhea.eq.1 .and. swcalt.eq.2) then
          call rdador ('z_Tsoil',-1.0d5,0.0d0,zh,macp,ifnd)
          call rdfdor ('Tsoil',-50.0d0,50.0d0,Tsoil,macp,ifnd)
        endif
        if (swsolu.eq.1 .or. flAgeTracer) then
          call rdador ('z_Cml',-1.0d5,0.0d0,zc,macp,ifnd)
          call rdfdor ('Cml',0.0d0,1000000.0d0,cml,macp,ifnd)
        endif
!       surface water level
        if(swdra.eq.2) then
          if(rdinqr('wls')) then
            call rdsdor ('wls',-1000.0d0,1000.0d0,wls)
            messag = 'Initial wls read from file (SWINCO=3)'//          &
     &      ' not implemented yet'
            call warn ('readswap',messag,logf,swscre)
          endif
        endif
!       write soil evaporation reservoirs
        if(swredu.eq.1) then
          messag = 'Initial ldwet read from file (SWINCO=3)'//          &
     &     ' if absent, then default of 0.0 is assumed'
          call warn ('reduceva',messag,logf,swscre)
          if(rdinqr('ldwet')) then
!           Time after significant rainfall (d)  (Black)
            call rdsdor ('ldwet',0.0d0,300.0d0,ldwet)
          else
            ldwet = 1.0d0
          endif
        endif
        if(swredu.eq.2) then
          messag = 'Initial spev read from file (SWINCO=3)'//           &
     &     ' if absent, then default of 0.0 is assumed'
          call warn ('reduceva',messag,logf,swscre)
          if(rdinqr('spev')) then
!           Rainfall excess (cm) (Boesten/Stroosnijder)
            call rdsdor ('spev',0.0d0,1000.0d0,spev)
          else
            spev = 0.0d0
          endif
        endif
!       length of final timestep (d)
        if(rdinqr('dt')) then
          !call rdsdor ('dt',dtmin,dtmax,dt) ! Aanpassing H.M. Mulder 2016-10-31 probleem bij TTUTIL tijdens inlezen dt (strenge controle niet nodig vanwege ini-file geproduceerd met SWAP)
          call rdsdou ('dt',dt)
          if (dt .lt. dtmin) dt = dtmin
          if (dt .gt. dtmax) dt = dtmax
        endif

!       close file
        close(ini)

!         macropore variables
        if(swmacro.eq.1) then
          write (messag,107)                                            &
     &          '*  I/O of macropore variables from file ',             &
     &          '   (SWINCO=3) not implemented yet !'
          call warn ('readswap',messag,logf,swscre) 
        endif

      endif

! --- runon from external file
      if (swrunon .eq. 1) then
        runf = getun2 (10,90,2)
        call rdinit(runf,logf,trim(rufil))
        call rdador ('runoff',0.d0,1000.d0,runonarr,maday,ifnd)
        flrunon = .true.
        close (runf)
      else
        flrunon = .false.
      endif

! --- read soil surface temperatures
      if (swtopbhea .eq. 2) then
        filnam = trim(pathwork)//trim(tsoilfile)//'.tss'
        tss = getun2 (10,90,2)
        call rdinit(tss,logf,filnam)
        call rdatim ('datet',dates,mabbc,ifnd)
        call rdfdor('ttop',-50.0d0,50.d0,ttop,mabbc,ifnd)
        do i = 1, ifnd
           temtoptab(i*2)   = ttop(i)
           temtoptab(i*2-1) = dates(i)
        enddo
        close (tss)
      endif

! --- copy content of key-file to log-file
      write (logf,14)  
 14   format('*',70('-'),'*',/,' Echo of input file:',/)
      call copfl2(swp,swpfilnam,logf,.true.)

      return
      end

! ----------------------------------------------------------------------
      subroutine rddrb (swallo,swdtyp,swmacro,dramet,dra,ipos,nrlevs,   &
     &    logf,swdivd,l,zbotdr,owltab,drares,infres,qdrtab,basegw,      &
     &    wetper,khtop,khbot,kvtop,kvbot,zintf,entres,geofac,drfil,     &
     &    pathdrain,numlay,cofani,swnrsrf,swscre,cofintfl,expintfl,     &
     &    swdislay,swtopdislay,ztopdislay,ftopdislay,shape,SwTopnrsrf,  &
     &    swdivdinf,FacDpthInf)
! ----------------------------------------------------------------------
!     Date               : July 2002
!     Purpose            : read input data for basic drainage    
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

      integer   swallo(madr),swdtyp(madr),dramet,dra,ipos,nrlevs,logf
      integer   swdivd,swdivdinf,numlay,swnrsrf,swscre,swmacro

      real(8)   l(madr),zbotdr(madr),owltab(madr,2*maowl)
      real(8)   drares(madr),infres(madr),shape
      real(8)   qdrtab(50),basegw,wetper(madr),khtop,khbot,kvtop,kvbot
      real(8)   zintf,entres,geofac,cofani(maho)
      real(8)   cofintfl,expintfl
      integer   swdislay,swtopdislay(madr),SwTopnrsrf
      real(8)   ztopdislay(Madr),ftopdislay(Madr),FacDpthInf

      character(len=16) drfil
      character(len=80) pathdrain
      logical   rdinqr
! ----------------------------------------------------------------------
! --- local
      integer   ifnd,i,getun2, swintfl
      real(8)   datowl(maowl),level(maowl),qdrain(25),gwl(25)
      character(len=80)  filnam
      character(len=200) messag

! ----------------------------------------------------------------------

! --- open file with drainage data
      filnam = trim(pathdrain)//trim(drfil)//'.dra'
      dra = getun2 (10,90,2)
      call rdinit(dra,logf,filnam)

! --- method to establish drainage/infiltration fluxes
      call rdsinr ('dramet',1,3,dramet) 

      if (dramet.ne.3) nrlevs = 1

! --- division of drainage fluxes
      call rdsinr ('swdivd',0,1,swdivd)
      if (swdivd.eq.0) then
          write(messag,'(3a)')                                          &
     &    ' Variabel SWDIVD=0 in input file : ',trim(filnam),           &
     &    ' this is not recommended and may cause numerical instability'
          call warn ('rddrb',messag,logf,swscre)
      endif
      if (swdivd.eq.1 .and. dramet.eq.1) then
          write(messag,'(3a,i3)')                                       &
     &    ' Variabel SWDIVD=1 and DRAMET=1 in inputfile :',trim(filnam),&
     &    ' this is not allowed for drainage method (dramet)= ',dramet
          call fatalerr ('Rddrb',messag)
      endif
      if (swdivd .eq. 1) then
        if(rdinqr('swdivdinf')) then
          call rdsinr ('swdivdinf',0,1,swdivdinf)
            if (swdivdinf.eq.1 .and. dramet.ne.3) then
              write(messag,'(3a,i3)')                                   &
     &        ' Variabel SWDIVDINF=1 and DRAMET not 1 in inputfile :',  &
     &        trim(filnam),' this option is not allowed for drainage ', &
     &        ' method (dramet)= ',dramet
              call fatalerr ('Rddrb',messag)
            endif
          call rdsdor ('FacDpthInf',0.d0,1.d0,FacDpthInf)
        else
          swdivdinf = 0
        endif
        call rdfdor ('cofani',1.d-4,1000.d0,cofani,maho,numlay)
      endif


! --- input table of drainage flux as function of groundwater level
      if (dramet.eq.1) then  
        call rdsdor ('lm1',1.0d0,1000.d0,l(1))
        call rdador ('gwl',-10000.0d0,10.0d0,gwl,25,ifnd)
        call rdfdor ('qdrain',-100.d0,1000.0d0,qdrain,25,ifnd)

        l(1) = 100.0d0*l(1)
        do  i = 1,50
          qdrtab(i) = 0.0d0
        end do
        do i = 1,ifnd
          qdrtab(i*2-1) = abs(gwl(i))
          qdrtab(i*2) = qdrain(i)
        end do

!   - in case of rapid drainage due to macropore flow: read zbotdr
        if (swmacro.eq.1) then
           call rdsinr ('swdtyp' ,1,2,swdtyp(1))
           call rdsdor ('zdrabas',-1000.0d0,0.0d0,zbotdr(1))
        endif

! --- input drainage formula of Hooghoudt or Ernst
      elseif (dramet.eq.2) then
   
! ---   read drain characteristics
        call rdsdor ('lm2',1.0d0,1000.0d0,l(1))
        l(1) = 100.0d0*l(1)
        call rdsdor ('shape',0.0d0,1.0d0,shape)
        call rdsdor ('wetper',0.0d0,1000.0d0,wetper(1))
        call rdsdor ('zbotdr',-1000.0d0,0.0d0,zbotdr(1))
        call rdsdor ('entres',0.0d0,1000.0d0,entres)

! ---   read profile characteristics
        call rdsinr ('ipos',1,5,ipos)
        call rdsdor ('basegw',-1.0d4,0.0d0,basegw)
        call rdsdor ('khtop',0.0d0,1000.0d0,khtop)

        if (ipos.ge.3) then
          call rdsdor ('khbot',0.0d0,1000.0d0,khbot)
          call rdsdor ('zintf',-1.0d4,0.0d0,zintf)
        endif
        if (ipos.ge.4) then
          call rdsdor ('kvtop',0.0d0,1000.0d0,kvtop)
          call rdsdor ('kvbot',0.0d0,1000.0d0,kvbot)  
        endif
        if (ipos.eq.5) then
          call rdsdor ('geofac',0.0d0,100.0d0,geofac)
        endif

! --- drainage and infiltration resistance
      elseif (dramet.eq.3) then
   
        call rdsinr ('nrlevs',1,Madr,nrlevs)
        if (nrlevs.gt.5) then
          write(messag,'(3a)')                                          &
     &    ' Number of Drainage levels >5 in inputfile : ',trim(filnam), &
     &    '   part of the output is limited to 5 levels'
          call warn ('rddrb',messag,logf,swscre)
        endif

! --- type of highest drainage level
        call rdsinr ('swintfl',0,1,swintfl)
!   - 'swintfl' is only used as input variabele, its function is taken over in the entire code by 'swintfl'    
        swnrsrf = swintfl
        if (swnrsrf.eq.1) then
          call rdsdor ('cofintflb',0.01d0,10.0d0,cofintfl)     
          call rdsdor ('expintflb',0.1d0,1.0d0,expintfl) 
        endif

! kroes 20080707 : allow disabling of option
        if (swdivd .eq. 1) then
          if (swnrsrf.gt.0)  call rdsinr ('SwTopnrsrf',0,1,SwTopnrsrf)
        endif

! ---   drainage level 1
        if (nrlevs.ge.1) then
          call rdsdor ('drares1',1.0d0,1.0d5,drares(1))
          call rdsdor ('infres1',0.0d0,1.0d5,infres(1))
          call rdsinr ('swallo1',1,3,swallo(1))
          if (swdivd .eq. 1) then
            call rdsdor ('l1',1.0d0,100000.0d0,l(1))
            l(1) = 100.0d0*l(1)
          endif
          call rdsdor ('zbotdr1',-10000.0d0,0.0d0,zbotdr(1))
          call rdsinr ('swdtyp1',1,2,swdtyp(1))
          if (swdtyp(1).eq.2) then
            call rdatim ('datowl1',datowl,maowl,ifnd)
            call rdfdor ('level1',-10000.0d0,200.0d0,level,maowl,ifnd)
! -         store values in table
            do i = 1,ifnd
              owltab(1,i*2) = level(i)
              owltab(1,i*2-1) = datowl(i)
            end do
          endif
        endif 

        if (nrlevs.ge.2) then
          call rdsdor ('drares2',1.0d0,1.0d5,drares(2))
          call rdsdor ('infres2',0.0d0,1.0d5,infres(2))
          call rdsinr ('swallo2',1,3,swallo(2))
          if (swdivd .eq. 1) then
            call rdsdor ('l2',1.0d0,100000.0d0,l(2))
            l(2) = 100.0d0*l(2)
          endif
          call rdsdor ('zbotdr2',-10000.0d0,0.0d0,zbotdr(2))
          call rdsinr ('swdtyp2',1,2,swdtyp(2))
          if (swdtyp(2).eq.2) then
            call rdatim ('datowl2',datowl,maowl,ifnd)
            call rdfdor ('level2',-10000.0d0,200.0d0,level,maowl,ifnd)
! -         store values in table
            do i = 1,ifnd
              owltab(2,i*2) = level(i)
              owltab(2,i*2-1) = datowl(i)
            end do
          endif
        endif 

        if (nrlevs.ge.3) then
          call rdsdor ('drares3',1.0d0,1.0d5,drares(3))
          call rdsdor ('infres3',0.0d0,1.0d5,infres(3))
          call rdsinr ('swallo3',1,3,swallo(3))
          if (swdivd .eq. 1) then
            call rdsdor ('l3',1.0d0,100000.0d0,l(3))
            l(3) = 100.0d0*l(3)
          endif
          call rdsdor ('zbotdr3',-10000.0d0,0.0d0,zbotdr(3))
          call rdsinr ('swdtyp3',1,2,swdtyp(3))
          if (swdtyp(3).eq.2) then
            call rdatim ('datowl3',datowl,maowl,ifnd)
            call rdfdor ('level3',-10000.0d0,200.0d0,level,maowl,ifnd)
! -         store values in table
            do i = 1,ifnd
              owltab(3,i*2) = level(i)
              owltab(3,i*2-1) = datowl(i)
            end do
          endif
        endif 

        if (nrlevs.ge.4) then
          call rdsdor ('drares4',1.0d0,1.0d5,drares(4))
          call rdsdor ('infres4',0.0d0,1.0d5,infres(4))
          call rdsinr ('swallo4',1,3,swallo(4))
          if (swdivd .eq. 1) then
            call rdsdor ('l4',1.0d0,100000.0d0,l(4))
            l(4) = 100.0d0*l(4)
          endif
          call rdsdor ('zbotdr4',-10000.0d0,0.0d0,zbotdr(4))
          call rdsinr ('swdtyp4',1,2,swdtyp(4))
          if (swdtyp(4).eq.2) then
            call rdatim ('datowl4',datowl,maowl,ifnd)
            call rdfdor ('level4',-10000.0d0,200.0d0,level,maowl,ifnd)
! -         store values in table
            do i = 1,ifnd
              owltab(4,i*2) = level(i)
              owltab(4,i*2-1) = datowl(i)
            end do
          endif
        endif 

        if (nrlevs.ge.5) then
          call rdsdor ('drares5',1.0d0,1.0d5,drares(5))
          call rdsdor ('infres5',0.0d0,1.0d5,infres(5))
          call rdsinr ('swallo5',1,3,swallo(5))
          if (swdivd .eq. 1) then
            call rdsdor ('l5',1.0d0,100000.0d0,l(5))
            l(5) = 100.0d0*l(5)
          endif
          call rdsdor ('zbotdr5',-10000.0d0,0.0d0,zbotdr(5))
          call rdsinr ('swdtyp5',1,2,swdtyp(5))
          if (swdtyp(5).eq.2) then
            call rdatim ('datowl1',datowl,maowl,ifnd)
            call rdfdor ('level1',-10000.0d0,200.0d0,level,maowl,ifnd)
! -         store values in table
            do i = 1,ifnd
              owltab(5,i*2) = level(i)
              owltab(5,i*2-1) = datowl(i)
            end do
          endif
        endif

!       kro 20141230: allow limitation of infiltration head to the waterdepth in the channel
        swnrsrf = 0
        if(rdinqr('swnrsrf')) then
          call rdsinr ('swnrsrf',0,1,swnrsrf)
        endif 

      endif 

!     top of model dicharge layer, determined by factor or direct input
      if (swdivd .eq. 1) then
        call rdsinr ('swdislay',0,2,swdislay)
        if (swdislay .eq. 1) then
          call rdfinr ('swtopdislay',0,1,swtopdislay,madr,nrlevs)
          call rdfdor('ztopdislay',-1.0d4,0.0d0,ztopdislay,madr,nrlevs)
        elseif (swdislay .eq. 2) then
          call rdfinr ('swtopdislay',0,1,swtopdislay,madr,nrlevs)
          call rdfdor('ftopdislay',0.0d0,1.0d0,ftopdislay,madr,nrlevs)
        endif
      endif


! --- close file with drainage data
      close (dra)         

      return
      end
      
! ----------------------------------------------------------------------
      subroutine readcropfixed (crpfil,pathcrop,idev,lcc,tsumea,        &
     & tsumam,tbase,kdif,kdir,gctb,swgc,cftb,                           &
     & cfeictb,fimin,siccaptb,swcf,                                     &
     & rdtb,rdctb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,          &
     & adcrh,adcrl,kytb,cofab,logf,schedule,swinter,pfreetb,pstemtb,    &
     & scanopytb,avprectb,avevaptb,cumdens,chtb,albedo,swetr,flsolute,  &
     & alphacrit,swdrought,wiltpoint,rootradius,rootcoefa,rsw,          &
     & q10_root,q10_microbial,specific_resp_humus,                      &
     & c_mroot,srl,dry_mat_cont_roots,air_filled_root_por,              &
     & spec_weight_root_tissue,var_a,f_senes,                           &
     & swhea,swcalt,swoxygen,swtopsub,nrstaring,oxygenslope,            &
     & oxygenintercept,swoxygentype,rooteff,swsalinity,                 &
     & swhydrlift,criterhr,stephr,kroot,rxylem,taccur,kstem,            &
     & wrtb,mrftb,swrootradius,root_radiusO2,bdens,inifil,              &
     & nofd,daycrop,dvs,tsum,t1900,tstart,swinco,                       &
     & cropstart,saltmax,saltslope,salthead)
! ----------------------------------------------------------------------
!     Update             : October 2012
!     Update             : July 2009
!     date               : July 2002             
!     purpose            : get crop parameters from cropfile
! ----------------------------------------------------------------------
      implicit  none
      include  'arrays.fi'
      
      integer   crp,idev,i,lcc,swgc,swcf,logf,ifnd,getun2,schedule
      integer   swinter,swetr,swdrought,swhea,swcalt,swoxygen
      integer   swtopsub,nrstaring,swoxygentype,swsalinity,swhydrlift
      integer   swrootradius,nofd,ini,daycrop
      logical   flsolute 
      real(8)   gctb(2*magrs),cftb(2*magrs),cfeictb(2*magrs)
      real(8)   fimin,siccaptb(2*magrs),mrftb(2*magrs)
      real(8)   rdtb(2*magrs),kytb(2*magrs),wrtb(2*magrs)
      real(8)   adcrh,adcrl,tbase,tsumam,tsumea,rdctb(22)
      real(8)   hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc
      real(8)   sum,afgen,kdif,kdir,cofab,chtb(2*magrs)
      real(8)   tinter(magrs),pfree(magrs),pstem(magrs)
      real(8)   scanopy(magrs),avprec(magrs),avevap(magrs)
      real(8)   pfreetb(2*magrs),pstemtb(2*magrs),scanopytb(2*magrs)
      real(8)   avprectb(2*magrs),avevaptb(2*magrs)
      real(8)   depth,rootdis(202),cumdens(202),rooteff
      real(8)   dvsinput(magrs),cfinput(magrs),chinput(magrs)
      real(8)   oxygenslope(6)
      real(8)   albedo,alphacrit,dvs,tsum
      real(8)   saltmax,saltslope,salthead
      real(8)   wiltpoint,rootradius,rootcoefa,rsw,oxygenintercept(6)
      real(8)   q10_root,q10_microbial,kstem,root_radiusO2
      real(8)   c_mroot,srl,dry_mat_cont_roots,bdens(maho)
      real(8)   spec_weight_root_tissue,var_a,f_senes
      real(8)   specific_resp_humus,air_filled_root_por
      real(8)   criterhr,stephr,kroot,rxylem,taccur
      logical   rdinqr
      character(len=*)   crpfil,pathcrop
      character(len=200) inifil

! locals
      character(len=200) message,filnam
      real(8)   siccap(magrs)
      real(8)   t1900,tstart,cropstart
      integer   swinco
! ----------------------------------------------------------------------

! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)

! --- phenology
      call rdsinr ('idev',1,2,idev)
      if (idev.eq.1) then
        call rdsinr ('lcc',1,366,lcc)
      elseif (idev.eq.2) then
        call rdsdor ('tsumea',0.0d0,10000.0d0,tsumea)
        call rdsdor ('tsumam',0.0d0,10000.0d0,tsumam)
        call rdsdor ('tbase',-10.0d0, 30.0d0,tbase)
      endif

! --- assimilation                        
      call rdsdor ('kdif',0.0d0,2.0d0,kdif)
      call rdsdor ('kdir',0.0d0,2.0d0,kdir)
     
! --- LAI or soil cover fraction 
      call rdsinr ('swgc',1,2,swgc)
      if (swgc.eq.1) then
        call rdador ('gctb',0.0d0,12.0d0,gctb,(2*magrs),ifnd)
      elseif (swgc.eq.2) then
        call rdador ('gctb',0.0d0,2.0d0,gctb,(2*magrs),ifnd)
        do i = 2,ifnd,2
          if (gctb(i).gt.1.0d0) then
            message = 'Fatal Error in crop-(.crp)-file: SoilCover > 1 '
            call fatalerr ('ReadCropFixed',message)
          endif
        enddo
      endif

! --- Crop factor or crop height
      call rdsinr ('swcf',1,3,swcf)

! --- check use of crop factors in case of ETref
      if (swetr.eq.1 .and. swcf.eq.2) then
        message = 'If ETref is used (SWETR = 1), always define crop '// &
     &           'factors (SWCF = 1 or 3)' 
        call fatalerr ('ReadCropFixed',message)
      endif

      if (swcf.eq.1 .or. swcf.eq.3) then
! ---   crop factor is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---   store values in cftb
        do i = 1,ifnd
          cftb(i*2) = cfinput(i) 
          cftb(i*2-1) = dvsinput(i)
        enddo
        if (swcf.eq.3) then
! ---     crop factor of wet crop is input
          call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
          call rdfdor ('cfw',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---     store values in cfwtb with crop factor of wet crop is input
          do i = 1,ifnd
            cfeictb(i*2) = cfinput(i) 
            cfeictb(i*2-1) = dvsinput(i)
          enddo
        endif
        chtb = -99.99d0
      else
! ---   crop height is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('ch',0.0d0,1.0d4,chinput,(magrs),ifnd)
! ---   store values in chtb
        do i = 1,ifnd
          chtb(i*2) = chinput(i) 
          chtb(i*2-1) = dvsinput(i)
        enddo
        cftb = -99.99d0
      endif

! --- reflection coefficient and crop resistance
      if (swcf.eq.1 .or. swcf.eq.3) then
! ---   use standard values for ETref
        albedo = 0.23d0
        rsc = 70.0d0
        rsw = 0.0d0
      else
! ---   use crop specific values
        call rdsdor ('albedo',0.0d0,1.0d0,albedo)
        call rdsdor ('rsc',0.0d0,1.0d6,rsc)
        call rdsdor ('rsw',0.0d0,1.0d6,rsw)
      endif

! --- rooting depth
      call rdador ('rdtb',0.0d0,1000.0d0,rdtb,(2*magrs),ifnd)

! --- yield response
      call rdador ('kytb',0.0d0,5.0d0,kytb,(2*magrs),ifnd)

! --- oxygen stress
      swoxygen = 1
      if(rdinqr('swoxygen')) then
        call rdsinr ('swoxygen',0,2,swoxygen)
! -     fatal error when physical oxygen stress is simulated without numerical heat flow
        if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1)) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), soil heat flow should be'//&
     &   ' numerically simulated: SwCalT=2! Adapt .swp input file.'
          call fatalerr ('readcropfixed',message)
        endif
! -     fatal error when physical oxygen stress is simulated without realistic bdens value
        if (swoxygen.eq.2 .and. bdens(1).lt.100.d0) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), bulk density must have'//  &
     &   ' realistic values; adjust BDENS-value(s) in .swp input file.'
          call fatalerr ('readcropfixed',message)
        endif
      endif

      if (swoxygen.eq.1) then                                          !
! ---   oxygen stress according to Feddes
        call rdsdor ('hlim1' ,-100.0d0,100.0d0,hlim1)                   !
        call rdsdor ('hlim2u',-1000.0d0,100.0d0,hlim2u)                 !
        call rdsdor ('hlim2l',-1000.0d0,100.0d0,hlim2l)                 !
      endif

      if (swoxygen.eq.2) then                                          
! ---   oxygen stress according to Bartholomeus

        swoxygentype = 1
        if(rdinqr('swoxygentype')) then
          call rdsinr ('swoxygentype',1,2,swoxygentype)
        endif

        if (swoxygentype.eq.1) then
! ---      use physical processes
           call rdsdor ('q10_root',1.0d0,4.0d0,q10_root)
           call rdsdor ('q10_microbial',1.0d0,4.0d0,q10_microbial)
           call rdsdor ('specific_resp_humus',0.0d0,1.0d0,              &
     &                   specific_resp_humus)
           call rdsdor ('c_mroot',0.0d0,1.0d0,c_mroot)
           call rdsdor ('srl',0.0d0,1.0d10,srl)
           call rdsdor ('f_senes',0.0d0,1.0d0,f_senes)         
           call rdador ('wrtb',0.0d0,10.0d0,wrtb,(2*magrs),ifnd)
           call rdador ('mrftb',0.0d0,100.0d0,mrftb,(2*magrs),ifnd)
           call rdsinr ('swrootradius',1,2,swrootradius)
           if (swrootradius .eq. 1) then
             call rdsdor ('dry_mat_cont_roots',0.0d0,1.0d0,             &
     &                   dry_mat_cont_roots)
             call rdsdor ('air_filled_root_por',0.0d0,1.0d0,            &
     &                   air_filled_root_por)
             call rdsdor ('spec_weight_root_tissue',0.0d0,1.0d5,        &
     &                   spec_weight_root_tissue)
             call rdsdor ('var_a',0.0d0,1.0d0,var_a)
           else
             call rdsdor ('root_radiusO2',1.0d-6,0.1d0,root_radiusO2)
           endif
        else
! ---      use reproduction functions
           call rdsinr ('SwTopSub',1,2,SwTopSub)
           call rdsinr ('NrStaring',1,18,NrStaring)
           call oxygen_dat (SwTopSub,NrStaring,OxygenSlope,             &
     &                      OxygenIntercept)
        endif
      endif

! --- drought stress
      swdrought = 1
      if(rdinqr('swdrought')) then
        call rdsinr ('swdrought',1,2,swdrought)                         
      endif

      if (swdrought.eq.1) then                                          
! ---    drought stress according to Feddes
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
         call rdsdor ('adcrh',0.0d0,5.0d0,adcrh)                        
         call rdsdor ('adcrl',0.0d0,5.0d0,adcrl)                        
! --     Criticial stress index for compensation of root water uptake (-)
         alphacrit = 1.0d0
         if(rdinqr('alphacrit')) then
           call rdsdor ('alphacrit',0.1d0,1.0d0,alphacrit)
         endif
      else
! ---    drought stress according to De Jong van Lier
         call rdsdor ('wiltpoint',-1.0d8,-1.0d2,wiltpoint)
         call rdsdor ('kstem',1.0d-10,1.0d1,kstem)               
         call rdsdor ('rxylem',1.0d-4,1.d0,rxylem)               
         call rdsdor ('rootradius',1.0d-4,1.0d0,rootradius)            
         call rdsdor ('kroot',1.0d-10,1.0d10,kroot)            
         call rdsdor ('rootcoefa',0.0d0,1.0d0,rootcoefa)
         call rdsinr ('swhydrlift',0,1,swhydrlift)
         call rdsdor ('rooteff',0.0d0,1.0d0,rooteff)                 
         call rdsdor ('stephr',0.0d0,10.d0,stephr)                 
         call rdsdor ('criterhr',0.0d0,10.d0,criterhr)                 
         call rdsdor ('taccur',0.0d-5,10.d-2,taccur)                 
      endif                                                             


! --- salt stress
      if (flsolute) then

        if(rdinqr('swsalinity')) then
          call rdsinr ('swsalinity',0,2,swsalinity)
        endif
 
        if (swsalinity .eq. 1) then
! ---     input for Maas and Hoffman salt reduction function
          call rdsdor ('saltmax',0.0d0,100.0d0,saltmax)
          call rdsdor ('saltslope',0.0d0,1.0d0,saltslope)
        elseif (swsalinity .eq. 2) then
! ---     osmotic head salinity stress concept
          call rdsdor ('salthead',0.0d0,1.0d3,salthead)
          if (swdrought.eq.1) then
           message = 'In case salinity stress is calculated with'//     &
     &     ' osmotic head (SwSalinity = 2), the drought stress'//       &
     &     ' should be calculated according to De Jong van Lier:'//     &
     &     ' SwDrought = 2! Adapt .crp input file.'
           call fatalerr ('readcropfixed',message)
          endif
        endif
      endif


! --- interception
      call rdsinr ('swinter',0,3,swinter)
      if (swinter .eq. 1) then
        call rdsdor ('cofab',0.0d0,2.0d0,cofab)
      else if (swinter .eq. 2) then
        call rdador ('t',0.d0,366.d0,tinter,(magrs),ifnd)
        call rdfdor ('pfree',0.d0,1.d0,pfree,(magrs),ifnd)
        call rdfdor ('pstem',0.d0,1.d0,pstem,(magrs),ifnd)
        call rdfdor ('scanopy',0.d0,10.d0,scanopy,(magrs),ifnd)
        call rdfdor ('avprec',0.d0,100.d0,avprec,(magrs),ifnd)
        call rdfdor ('avevap',0.d0,10.d0,avevap,(magrs),ifnd)
        do i = 1, ifnd
          pfreetb(i*2) = pfree(i)
          pfreetb(i*2-1) = tinter(i)
          pstemtb(i*2) = pstem(i)
          pstemtb(i*2-1) = tinter(i)
          scanopytb(i*2) = scanopy(i)
          scanopytb(i*2-1) = tinter(i)
          avprectb(i*2) = avprec(i)
          avprectb(i*2-1) = tinter(i)
          avevaptb(i*2) = avevap(i)
          avevaptb(i*2-1) = tinter(i)
        end do
      else if (swinter .eq. 3) then
        call rdsdor ('fimin',0.0d0,1.0d0,fimin)
        call rdador ('t',0.d0,366.d0,tinter,(magrs),ifnd)
        call rdfdor ('siccap',0.d0,1.d0,siccap,(magrs),ifnd)
        do i = 1, ifnd
          siccaptb(i*2)   = siccap(i)
          siccaptb(i*2-1) = tinter(i)
        enddo      
      endif

! --- read table with root distribution coefficients
      call rdador ('rdctb',0.0d0,100.0d0,rdctb,22,ifnd)


! --- determine whether irrigation scheduling is applied
      call rdsinr ('schedule',0,1,schedule)

      if (schedule.eq.1 .and. swdrought.eq.2) then
! ---    read limiting pressure heads for irrigation scheduling
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
      endif

! --- close file with crop data
      close (crp)


! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION
      if (swdrought.eq.1) then 
! ---   root water extraction according to Feddes function  


! ---   specify array ROOTDIS with root density distribution
        do i = 0,100
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
        enddo

! ---   calculate cumulative root density function
        do i = 1,202,2
! ---     relative depths
          cumdens(i) = rootdis(i)
        enddo
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0                 &
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
        enddo

! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo
      endif

! --- read crop data of former day from *.END file
      if (t1900 - tstart .lt. 1.d-3 .and. swinco .eq. 3 .and. .not.     &
     &   (t1900 - cropstart .gt. -1.d-3 .and.                           &
     &                t1900 - cropstart .lt. 1.d-3)) then

!       open inifile
        ini = getun2 (10,90,2)
        call rdinit(ini,logf,inifil)

        call rdsinr ('nofd',0,7,nofd)
        call rdsinr ('daycrop',0,366,daycrop)
        call rdsdor ('dvs',0.0d0,2.0d0,dvs)
        call rdsdor ('tsum',0.0d0,1.0d4,tsum)

!       close inifile
        close(ini)
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine readwofost (crpfil,pathcrop,swcf,cftb,idsl,dlo,dlc,    &
     &  initcrp,tsumemeopt,tbasem,teffmx,                               &
     &  hdrygerm,hwetgerm,agerm,cgerm,bgerm,                            &
     &  tsumea,tsumam,dtsmtb,dvsend,tdwi,laiem,rgrlai,slatb,spa,        &
     &  ssa,span,tbase,kdif,kdir,eff,amaxtb,tmpftb,tmnftb,cvl,cvo,cvr,  &
     &  cvs,q10,rml,rmo,rmr,rms,rfsetb,frtb,fltb,fstb,fotb,perdl,rdrrtb,&
     &  rdrstb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,       &
     &  adcrl,cofab,rdi,rri,rdc,rdctb,logf,schedule,cumdens,chtb,albedo,&
     &  swetr,flsolute,relmf,flpotrelmf,cfeictb,siccaplai,alphacrit,    &
     &  swdrought,wiltpoint,rootradius,rootcoefa,rsw,                   &
     &  q10_microbial,specific_resp_humus,                              &
     &  srl,dry_mat_cont_roots,air_filled_root_por,                     &
     &  spec_weight_root_tissue,var_a,                                  &
     &  swhea,swcalt,swoxygen,swtopsub,nrstaring,oxygenslope,           &
     &  oxygenintercept,swoxygentype,rooteff,swsalinity,                &
     &  swhydrlift,criterhr,stephr,kroot,rxylem,taccur,kstem,           &
     &  swrootradius,root_radiusO2,swdmi2rd,bdens,                      &
     &  rd,rdpot,dvs,flanthesis,tsum,ilvold,ilvoldpot,wrt,wrtpot,tadw,  &
     &  tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,lai,laipot,dwrt,&
     &  dwrtpot,dwlv,dwlvpot,dwst,dwstpot,gasst,gasstpot,mrest,mrestpot,&
     &  cwdm,cwdmpot,sla,slapot,lvage,lvagepot,lv,                      &
     &  lvpot,inifil,daycrop,tsumgerm,nofd,atmin7,swsoybean,mg,         &
     &  dvsi,dvrmax1,dvrmax2,flrfphotoveg,tmaxdvr,tmindvr,toptdvr,      &
     &  popt,pcrt,flphenodayl,                                          &
     &  flCropCalendar,flCropEmergence,flCropHarvest,swinco,            &
     &  t1900,tstart,cropstart,verndvs,vernsat,vernbase,verntb,         &
     &  saltmax,saltslope,salthead)
! ----------------------------------------------------------------------
!     Update:  March 2017
!     Update:  December 2009
!     purpose: read parameters for wofost crop growth routine
! ----------------------------------------------------------------------
      implicit none
      include  'arrays.fi'

      integer crp,i,idsl,logf,ifnd,swcf,getun2,schedule,swetr
      integer swdrought,swhea,swcalt,swoxygen,swtopsub,nrstaring
      integer swsalinity,swoxygentype,swhydrlift,ini,swanthesis,count
      integer swrootradius,initcrp,swdmi2rd,ilvold,ilvoldpot
      integer daycrop,nofd
      logical flsolute,rdinqr,flanthesis
      real(8) cftb(2*magrs),cfeictb(2*magrs),siccaplai
      real(8) dtsmtb(30),slatb(30),amaxtb(30),tmpftb(30)
      real(8) tmnftb(30),rfsetb(30),frtb(30),fltb(30),fstb(30),fotb(30)
      real(8) rdrrtb(30),rdrstb(30),kdif,kdir,cofab,laiem
      real(8) adcrh,adcrl,cvl,cvo,cvr,cvs,dlc,dlo,dvsend,eff
      real(8) perdl,q10,rdc,rdi,rgrlai,rml,rmo,rmr,rms,rri,spa,span
      real(8) ssa,tbase,tdwi,tsumam,tsumea,rdctb(22),chtb(2*magrs)
      real(8) hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc
      real(8) sum,afgen,depth,rootdis(202),cumdens(202),albedo
      real(8) oxygenslope(6)
      real(8) relmf,alphacrit,atmin7(7),tsumgerm
      real(8) wiltpoint,rootradius,rootcoefa,rsw,oxygenintercept(6)
      real(8) q10_microbial,root_radiusO2
      real(8) srl,dry_mat_cont_roots,saltmax,saltslope,salthead
      real(8) spec_weight_root_tissue,var_a,bdens(maho)
      real(8) specific_resp_humus,air_filled_root_por,rooteff
      real(8) criterhr,stephr,kroot,rxylem,taccur,kstem
      real(8) tsumemeopt,TBASEM,TEFFMX
      real(8) hdrygerm,hwetgerm,agerm,cgerm,bgerm
      real(8) rd,rdpot,dvs,tsum,wrt,wrtpot,tadw,dwrt,mrestpot
      real(8) tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,lai,laipot
      real(8) dwrtpot,dwlv,dwlvpot,dwst,dwstpot,gasst,gasstpot,mrest
      real(8) cwdm,cwdmpot,lv(366)
      real(8) lvpot(366),sla(366),slapot(366),lvage(366),lvagepot(366)
      character(len=*)   crpfil,pathcrop
      character(len=200) inifil
! --- only for soybean
      logical   flpotrelmf         ! Flag indicating calculation of limited attainable yield instead of potential yield (due to management factor)
      real(8) mg,dvsi,dvrmax1,dvrmax2,tmaxdvr,tmindvr,toptdvr
      integer swsoybean
      logical flrfphotoveg
      logical flphenodayl     ! Flag to allow input of POPT and PCRT or using empirical relation from Setiyono et al
      real(8) popt            ! optimal daylength for phenological development (hr)
      real(8) pcrt            ! critical daylength for phenological development (hr)
!
      integer swCropEmergence
      logical flCropCalendar,flCropEmergence
      integer swcropharvest
      logical flCropHarvest
      integer swinco
! --- only for vernalisation
      real(8)   verndvs            ! critical development stage after which the effect of vernalisation is halted [-]
      real(8)   vernsat            ! saturated vernalisation requirement [d]
      real(8)   vernbase           ! base vernalisation requirement [d]
      real(8)   verntb(30)         ! table with rate of vernalisation as function of tav [days/degrees]
      
      real(8) t1900,tstart,cropstart

! locals
      real(8)   dvsinput(magrs),cfinput(magrs),chinput(magrs)
      real(8)   laiinput(magrs),cfeicinput(magrs)
      character(len=200) message,filnam

! ----------------------------------------------------------------------

! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)


      ! --- crop factor or crop height
      call rdsinr ('swcf',1,3,swcf)

! --- check use of crop factors in case of ETref
      if (swetr.eq.1 .and. swcf.eq.2) then
        message = 'If ETref is used (SWETR = 1), always define crop '// &
     &           'factors (SWCF = 1 or 3)' 
        call fatalerr ('ReadWofost',message)
      endif

      if (swcf.eq.1) then
! ---   crop factor is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---   store values in cftb
        do i = 1,ifnd
          cftb(i*2) = cfinput(i) 
          cftb(i*2-1) = dvsinput(i)
        enddo
        chtb = -99.99d0
      elseif (swcf.eq.2)then
! ---   crop height is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('ch',0.0d0,1.0d4,chinput,(magrs),ifnd)
! ---   store values in chtb
        do i = 1,ifnd
          chtb(i*2) = chinput(i) 
          chtb(i*2-1) = dvsinput(i)
        enddo
        cftb = -99.99d0
      elseif (swcf.eq.3)then
! ---   LAI-dependent crop factors for dual crop coefficient
        call rdador ('lai',0.0d0,10.0d0,laiinput,36,ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,36,ifnd)
        call rdfdor ('cfeic',0.0d0,2.0d0,cfeicinput,36,ifnd)

        do i = 1,ifnd
          cftb(i*2)        = cfinput(i) 
          cftb(i*2-1)      = laiinput(i)
          cfeictb(i*2)     = cfeicinput(i) 
          cfeictb(i*2-1)   = laiinput(i)
        enddo

! ---   LAI-dependent crop height
        call rdfdor ('ch',0.0d0,1.0d2,chinput,36,ifnd)

        do i = 1,ifnd
          chtb(i*2)     = chinput(i)
          chtb(i*2-1)   = laiinput(i)
        enddo
!
! ---   Interception capacity per unit of LAI
        call rdsdor ('vxiclai',0.0d0,1.0d-2,siccaplai)
        siccaplai = siccaplai          
      endif

! --- reflection coefficient and crop resistance
      if (swcf.eq.1 .or. swcf.eq.3) then
! ---   use standard values for ETref
        albedo = 0.23d0
        rsc = 70.0d0
        rsw = 0.0d0
      else
! ---   use crop specific values
        call rdsdor ('albedo',0.0d0,1.0d0,albedo)
        call rdsdor ('rsc',0.0d0,1.0d6,rsc)
        call rdsdor ('rsw',0.0d0,1.0d6,rsw)
      endif

! --- phenology
! --- only for soybean:  maturity group parameters
      swsoybean = 0
      if(rdinqr('swsoybean')) then
        call rdsinr ('swsoybean',0,1,swsoybean)
      endif
      if (swsoybean.eq.1) then
        call rdsdor ('mg',0.1d0,9.0d0,mg)
        call rdsdor ('dvsi',0.0d0,2.0d0,dvsi)
        call rdsdor ('dvrmax1',0.0d0,1.0d0,dvrmax1)
        call rdsdor ('dvrmax2',0.0d0,1.0d0,dvrmax2)
        call rdsdor ('tmaxdvr',0.0d0,45.0d0,tmaxdvr)
        call rdsdor ('tmindvr',0.0d0,tmaxdvr,tmindvr)
        call rdsdor ('toptdvr',tmindvr,tmaxdvr,toptdvr)
        flrfphotoveg = .true.
        if(rdinqr('flrfphotoveg')) then
           call rdslog('flrfphotoveg',flrfphotoveg)
        endif
        flphenodayl = .false.
        if(rdinqr('flphenodayl')) then
           call rdslog('flphenodayl',flphenodayl)
           if(flphenodayl) then
              call rdsdor ('POPT',0.0d0,24.0d0,POPT)
              call rdsdor ('PCRT',0.0d0,24.0d0,PCRT)
           endif
        endif
      endif
! --- for non-soybean crops:
!     idsl =  Switch for crop development before anthesis: 
!             0 depends on temperature; 
!             1 depends on temperature and day length; 
!             2 depends on temperature, day length and vernalisation factor
      if (swsoybean.eq.0) then
        call rdsinr ('idsl',0,2,idsl)
        if (idsl.eq.1.or.idsl.eq.2) then
          call rdsdor ('dlo',0.0d0,24.0d0,dlo)
          call rdsdor ('dlc',0.0d0,24.0d0,dlc)
        endif
        if (idsl.eq.0.or.idsl.eq.2) then
          call rdsdor ('tsumea',0.0d0,10000.0d0,tsumea)
          call rdsdor ('tsumam',0.0d0,10000.0d0,tsumam)
        endif
        call rdador ('dtsmtb',0.0d0,100.0d0,dtsmtb,30,ifnd)
      endif

      call rdsdor ('dvsend',0.0d0,3.0d0,dvsend)
!     vernalisation
      if (idsl.eq.2) then
        call rdsdor ('verndvs',0.0d0,0.4d0,verndvs)            ! critical development stage after which the effect of vernalisation is halted [-]
        call rdsdor ('vernsat',0.0d0,200.0d0,vernsat)          ! saturated vernalisation requirement [d]
        call rdsdor ('vernbase',0.0d0,200.0d0,vernbase)        ! base vernalisation requirement [d]
        call rdador ('verntb',-100.0d0,100.0d0,verntb,30,ifnd) ! rate of vernalisation as function of tav [days/degrees]
      endif

!---  Germination option : timing parameters
      if (initcrp .eq. 2) then
        call rdsdor ('tsumemeopt',0.0d0,1.0d3,tsumemeopt)
        call rdsdor ('TBASEM',-20.0d0,4.0d1,TBASEM)
        call rdsdor ('TEFFMX',0.0d0,4.0d1,TEFFMX)
        call rdsdor ('hdrygerm',-1000.0d0,-1.0d-2,hdrygerm)
        call rdsdor ('hwetgerm',-100.0d0,-1.0d-2,hwetgerm)
        call rdsdor ('agerm',1.0d0,1.0d3,agerm)
        call rdsdor ('cgerm',-1.0d3,-1.0d0,cgerm)
        call rdsdor ('bgerm',1.0d0,1.0d3,bgerm)
      endif
! --- initial
      call rdsdor ('tdwi',0.0d0,10000.0d0,tdwi)
      call rdsdor ('laiem',0.0d0,10.0d0,laiem)
      call rdsdor ('rgrlai',0.0d0,1.0d0,rgrlai)

! --- green area
      call rdador ('slatb',0.0d0,2.0d0,slatb,30,ifnd)
      call rdsdor ('spa',0.0d0,1.0d0,spa)
      call rdsdor ('ssa',0.0d0,1.0d0,ssa)
      call rdsdor ('span',0.0d0,366.0d0,span)
      call rdsdor ('tbase',-10.0d0,30.0d0,tbase)

! --- assimilation
      call rdsdor ('kdif',0.0d0,2.0d0,kdif)
      call rdsdor ('kdir',0.0d0,2.0d0,kdir)
      call rdsdor ('eff',0.0d0,10.0d0,eff)
      call rdador ('amaxtb',0.0d0,100.0d0,amaxtb,30,ifnd)
      call rdador ('tmpftb',-10.0d0,50.0d0,tmpftb,30,ifnd)
      call rdador ('tmnftb',-10.0d0,50.0d0,tmnftb,30,ifnd)

! --- conversion of assimilates into biomass
      call rdsdor ('cvl',0.0d0,1.0d0,cvl)
      call rdsdor ('cvo',0.0d0,1.0d0,cvo)
      call rdsdor ('cvr',0.0d0,1.0d0,cvr)
      call rdsdor ('cvs',0.0d0,1.0d0,cvs)

! --- maintenance respiration
      call rdsdor ('q10',0.0d0,5.0d0,q10)
      call rdsdor ('rml',0.0d0,1.0d0,rml)
      call rdsdor ('rmo',0.0d0,1.0d0,rmo)
      call rdsdor ('rmr',0.0d0,1.0d0,rmr)
      call rdsdor ('rms',0.0d0,1.0d0,rms)
      call rdador ('rfsetb',0.0d0,3.0d0,rfsetb,30,ifnd)

! --- partitioning
      call rdador ('frtb',0.0d0,3.0d0,frtb,30,ifnd)
      call rdador ('fltb',0.0d0,3.0d0,fltb,30,ifnd)
      call rdador ('fstb',0.0d0,3.0d0,fstb,30,ifnd)
      call rdador ('fotb',0.0d0,3.0d0,fotb,30,ifnd)

! --- death rates
      call rdsdor ('perdl',0.0d0,3.0d0,perdl)
      call rdador ('rdrrtb',0.0d0,3.0d0,rdrrtb,30,ifnd)
      call rdador ('rdrstb',0.0d0,3.0d0,rdrstb,30,ifnd)

! --- oxygen stress
      swoxygen = 1
      if(rdinqr('swoxygen')) then
        call rdsinr ('swoxygen',0,2,swoxygen)
! -     fatal error when physical oxygen stress is simulated without numerical heat flow
        if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1)) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), soil heat flow should be'//&
     &   ' numerically simulated: SwCalT=2! Adapt .swp input file.'
          call fatalerr ('readwofost',message)
        endif
! -     fatal error when physical oxygen stress is simulated without realistic bdens value
        if (swoxygen.eq.2 .and. bdens(1).lt.100.0d0) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), bulk density must have'//  &
     &   ' realistic values; adjust BDENS-value(s) in .swp input file.'
          call fatalerr ('readwofost',message)
        endif
      endif

      if (swoxygen.eq.1) then
! ---   oxygen stress according to Feddes
        call rdsdor ('hlim1' ,-100.0d0,100.0d0,hlim1)
        call rdsdor ('hlim2u',-1000.0d0,100.0d0,hlim2u)
        call rdsdor ('hlim2l',-1000.0d0,100.0d0,hlim2l)
      endif

      if (swoxygen.eq.2) then                                          
! ---   oxygen stress according to Bartholomeus

        swoxygentype = 1
        if(rdinqr('swoxygentype')) then
          call rdsinr ('swoxygentype',1,2,swoxygentype)
        endif

        if (swoxygentype.eq.1) then
! ---      use physical processes
           call rdsdor ('q10_microbial',1.0d0,4.0d0,q10_microbial)
           call rdsdor ('specific_resp_humus',0.0d0,1.0d0,              &
     &                   specific_resp_humus)                        
           call rdsdor ('srl',0.0d0,1.0d10,srl)                        
           call rdsinr ('swrootradius',1,2,swrootradius)
           if (swrootradius.eq.1) then
              call rdsdor ('dry_mat_cont_roots',0.0d0,1.0d0,            &
     &                   dry_mat_cont_roots)                        
              call rdsdor ('air_filled_root_por',0.0d0,1.0d0,           &
     &                   air_filled_root_por)                        
              call rdsdor ('spec_weight_root_tissue',0.0d0,1.0d5,       &
     &                   spec_weight_root_tissue)
              call rdsdor ('var_a',0.0d0,1.0d0,var_a)
           else 
              call rdsdor ('root_radiusO2',0.0d0,1.0d0,root_radiusO2)
           endif  
        else
! ---      use reproduction functions
           call rdsinr ('SwTopSub',1,2,SwTopSub)
           call rdsinr ('NrStaring',1,18,NrStaring)
           call oxygen_dat (SwTopSub,NrStaring,OxygenSlope,             &
     &                      OxygenIntercept)
        endif
      endif

! --- drought stress
      swdrought = 1
      if(rdinqr('swdrought')) then
        call rdsinr ('swdrought',1,2,swdrought)                         
      endif

      if (swdrought.eq.1) then                                          
! ---    drought stress according to Feddes
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
         call rdsdor ('adcrh',0.0d0,5.0d0,adcrh)                        
         call rdsdor ('adcrl',0.0d0,5.0d0,adcrl)                        
! --     Criticial stress index for compensation of root water uptake (-)
         alphacrit = 1.0d0
         if(rdinqr('alphacrit')) then
           call rdsdor ('alphacrit',0.2d0,1.0d0,alphacrit)
         endif
      else 
! ---    drought stress according to De Jong van Lier
         call rdsdor ('wiltpoint',-1.0d8,-1.0d2,wiltpoint)
         call rdsdor ('kstem',1.0d-10,1.0d1,kstem)               
         call rdsdor ('rxylem',1.0d-4,1.d0,rxylem)               
         call rdsdor ('rootradius',1.0d-4,1.0d0,rootradius)            
         call rdsdor ('kroot',1.0d-10,1.0d10,kroot)            
         call rdsdor ('rootcoefa',0.0d0,1.0d0,rootcoefa)
         call rdsinr ('swhydrlift',0,1,swhydrlift)
         call rdsdor ('rooteff',0.0d0,1.0d0,rooteff)                 
         call rdsdor ('stephr',0.0d0,10.d0,stephr)                 
         call rdsdor ('criterhr',0.0d0,10.d0,criterhr)                 
         call rdsdor ('taccur',0.0d-5,10.d-2,taccur)                 
      endif

! --- salt stress
      if (flsolute) then

        if(rdinqr('swsalinity')) then
          call rdsinr ('swsalinity',0,2,swsalinity)
        endif
 
        if (swsalinity .eq. 1) then
! ---     input for Maas and Hoffman salt reduction function
          call rdsdor ('saltmax',0.0d0,100.0d0,saltmax)
          call rdsdor ('saltslope',0.0d0,1.0d0,saltslope)
        elseif (swsalinity .eq. 2) then
! ---     osmotic head salinity stress concept
          call rdsdor ('salthead',0.0d0,1.0d3,salthead)
          if (swdrought.eq.1) then
           message = 'In case salinity stress is calculated with'//     &
     &     ' osmotic head (SwSalinity = 2), the drought stress'//       &
     &     ' should be calculated according to De Jong van Lier:'//     &
     &     ' SwDrought = 2! Adapt .crp input file.'
           call fatalerr ('readcropfixed',message)
          endif
        endif
      endif

! --- management factor to account for other forms of stress
      relmf = 1.0d0
      if(rdinqr('relmf')) then
        call rdsdor ('relmf',0.0d0,1.0d0,relmf) 
      endif

! --- factor to account for difference between theoretical potential and attainable yield
      flpotrelmf = .false.
      if(rdinqr('flpotrelmf')) then
        call rdslog ('flpotrelmf',flpotrelmf)
      endif
      
! --- interception
      if(swcf.ne.3) then
        call rdsdor ('cofab',0.0d0,1.0d0,cofab)
      endif

! --- rooting
      call rdador ('rdctb',0.0d0,100.0d0,rdctb,22,ifnd)
      call rdsdor ('rdi',0.0d0,1000.0d0,rdi)
      call rdsdor ('rri',0.0d0,100.0d0,rri)
      call rdsdor ('rdc',0.0d0,1000.0d0,rdc)

! --- rooting depth influenced by dry matter increase
      swdmi2rd = 0
      if(rdinqr('swdmi2rd')) then
        call rdsinr ('swdmi2rd',0,1,swdmi2rd)
      endif

! --- determine whether irrigation scheduling is applied
      call rdsinr ('schedule',0,1,schedule)

      if (schedule.eq.1 .and. swdrought.eq.2) then
! ---    read limiting pressure heads for irrigation scheduling
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
      endif

      close (crp)

! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION

      if (swdrought.eq.1) then

! ---   specify array ROOTDIS with root density distribution
        do i = 0,100
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
        enddo

! ---   calculate cumulative root density function
        do i = 1,202,2
! ---     relative depths
          cumdens(i) = rootdis(i)
        enddo
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0                 &
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
        enddo

! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo
      endif 

! --- read crop data of former day from *.END file
      if (t1900 - tstart .lt. 1.d-3 .and. swinco .eq. 3 .and. .not.     &
     &   (t1900 - cropstart .gt. -1.d-3 .and.                           &
     &                t1900 - cropstart .lt. 1.d-3)) then
         
        if (flCropCalendar) then
        
          ini = getun2 (10,90,2)
          call rdinit(ini,logf,inifil)
          
          if (rdinqr('swcropharvest')) then
          
            call rdsinr ('swcropharvest',0,1,swcropharvest)
            if (swcropharvest .eq. 0) then
              flCropHarvest = .false.
            else
              flCropHarvest = .true.
            endif
            
            call rdsinr ('swCropEmergence',0,1,swCropEmergence)
            if (swCropEmergence .eq. 0) then
              flCropEmergence = .false.
            else
              flCropEmergence = .true.
            endif
          
            if (.not. flCropHarvest) then
              call rdsdor ('rd',0.0d0,1.0d4,rd)               
              call rdsdor ('rdpot',0.0d0,1.0d4,rdpot)               
              call rdsdor ('dvs',0.0d0,3.0d0,dvs)               
              call rdsinr ('daycrop',0,366,daycrop)               
              call rdsinr ('swanthesis',0,1,swanthesis)               
              call rdsdor ('tsum',0.0d0,1.0d4,tsum)               
              call rdsinr ('ilvold',0,366,ilvold)               
              call rdsinr ('ilvoldpot',0,366,ilvoldpot)               
              call rdsdor ('wrt',0.0d0,1.0d8,wrt)               
              call rdsdor ('wrtpot',0.0d0,1.0d8,wrtpot)               
              call rdsdor ('tadw',0.0d0,1.0d8,tadw)               
              call rdsdor ('tadwpot',0.0d0,1.0d8,tadwpot)               
              call rdsdor ('wst',0.0d0,1.0d8,wst)               
              call rdsdor ('wstpot',0.0d0,1.0d8,wstpot)               
              call rdsdor ('wso',0.0d0,1.0d8,wso)               
              call rdsdor ('wsopot',0.0d0,1.0d8,wsopot)               
              call rdsdor ('wlv',0.0d0,1.0d8,wlv)               
              call rdsdor ('wlvpot',0.0d0,1.0d8,wlvpot)               
              call rdsdor ('laiexp',0.0d0,1.0d8,laiexp)               
              call rdsdor ('lai',0.0d0,1.0d8,lai)
              call rdsdor ('laipot',0.0d0,1.0d8,laipot)
              call rdsdor ('dwrt',0.0d0,1.0d8,dwrt)
              call rdsdor ('dwrtpot',0.0d0,1.0d8,dwrtpot)
              call rdsdor ('dwlv',0.0d0,1.0d8,dwlv)
              call rdsdor ('dwlvpot',0.0d0,1.0d8,dwlvpot)
              call rdsdor ('dwst',0.0d0,1.0d8,dwst)
              call rdsdor ('dwstpot',0.0d0,1.0d8,dwstpot)
              call rdsdor ('gasst',0.0d0,1.0d10,gasst)
              call rdsdor ('gasstpot',0.0d0,1.0d10,gasstpot)
              call rdsdor ('mrest',0.0d0,1.0d10,mrest)
              call rdsdor ('mrestpot',0.0d0,1.0d10,mrestpot)
              call rdsdor ('cwdm',0.0d0,1.0d10,cwdm)
              call rdsdor ('cwdmpot',0.0d0,1.0d10,cwdmpot)
              call rdsdor ('tsumgerm',0.0d0,1.0d4,tsumgerm)
              call rdsinr ('nofd',0,366,nofd)               
              call rdfdor ('atmin7',-50.d0,60.d0,atmin7,7,7)
              
              count = max(ilvold,ilvoldpot)
              call rdfdor ('sla',0.d0,1.0d3,sla,366,count)
              call rdfdor ('slapot',0.d0,1.0d3,slapot,366,count)
              call rdfdor ('lvage',0.d0,1.0d3,lvage,366,count)
              call rdfdor ('lvagepot',0.d0,1.0d3,lvagepot,366,count)
              call rdfdor ('lv',0.d0,1.0d3,lv,366,count)
              call rdfdor ('lvpot',0.d0,1.0d3,lvpot,366,count)
              
              if (swanthesis .eq. 0) then
                flanthesis = .false.
              else
                flanthesis = .true.
              endif
            else

              daycrop = 0
              flCropEmergence = .false.

            endif
          else

            daycrop = 0
            flCropEmergence = .false.

          endif
        endif
!       close inifile
        close(ini)
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine readgrass (crpfil,pathcrop,tdwi,laiem,rgrlai,slatb,    &
     &  ssa,span,tbase,kdif,kdir,eff,amaxtb,tmpftb,tmnftb,cvl,cvr,cvs,  &
     &  q10,rml,rmr,rms,rfsetb,frtb,fltb,fstb,perdl,rdrrtb,             &
     &  rdrstb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl, &
     &  cofab,rdi,rdc,rdctb,rlwtb,logf,schedule,cumdens,flsolute,       &
     &  SeqGrazMow,swharvest,dmharvest,daylastharvest,dmlastharvest,    &
     &  dateharvest,flmowingtb,dmmowtb,maxdaymow,daymow,daymowpot,      &
     &  DelayRegrowthTab,dmgrazing,LSDb,DaysGrazing,UptGrazing,lsda,    &
     &  DaysGrazingtab,UptGrazingtab,tagprest,LossGrazing,              &
     &  wrtmax,relmf,flpotrelmf,LossGrazingtab,                         &
     &  swcf,swetr,cftb,chtb,cfeictb,siccaplai,alphacrit,               &
     &  swdrought,wiltpoint,rootradius,rootcoefa,rsw,                   &
     &  q10_microbial,specific_resp_humus,                              &
     &  srl,dry_mat_cont_roots,air_filled_root_por,                     &
     &  spec_weight_root_tissue,var_a,                                  &
     &  swhea,swcalt,swoxygen,swtopsub,nrstaring,oxygenslope,           &
     &  oxygenintercept,swoxygentype,rooteff,swsalinity,                &
     &  swhydrlift,criterhr,stephr,kroot,rxylem,taccur,kstem,           &
     &  swrootradius,root_radiusO2,bdens,fltsum200,                     &
     &  rd,rdpot,dvs,flanthesis,tsum,ilvold,ilvoldpot,wrt,wrtpot,tadw,  &
     &  tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,lai,laipot,dwrt,&
     &  dwrtpot,dwlv,dwlvpot,dwst,dwstpot,gasst,gasstpot,mrest,mrestpot,&
     &  cwdm,cwdmpot,sla,slapot,lvage,lvagepot,lv,                      &
     &  lvpot,inifil,daycrop,tsumgerm,nofd,atmin7,rid,flgrazingpot,     &
     &  idregr,idregrpot,laiexppot,laimax,tagp,tagppot,mowrest,         &
     &  tagpt,tagptpot,tsum200,iharvest,idaysgraz,                      &
     &  iseqgm,iseqgmpot,flgrazing,idaysgrazpot,                        &
     &  flCropCalendar,t1900,tstart,tend,swinco,cropstart,              &
     &  fltsumttd,tsumtime,tsumtemp,tsumdepth,cuptgraz,cuptgrazpot,     &
     &  saltmax,saltslope,salthead)

! ----------------------------------------------------------------------
!     Update             : August 2014   
!     date               : november 2004   
!     purpose            : read parameters for grass growth routine
! ----------------------------------------------------------------------
      implicit none
      include  'arrays.fi'

      integer crp, logf,ifnd,getun2,schedule,swdrought
      integer swcf,swetr,swhea,swcalt,swoxygen,swtopsub,nrstaring
      integer daylastharvest,swharvest,SeqGrazMow(366)
      integer maxdaymow,daymow,daymowpot,daydelay(100)
      integer swsalinity,swhydrlift,swrootradius,swoxygentype
      integer ilvold,ilvoldpot,daycrop,nofd,idaysgrazpot
      integer idregr,idregrpot,iseqgm,iseqgmpot,iharvest,idaysgraz
      real(8) slatb(30),amaxtb(30),tmpftb(30),depth,rootdis(202)
      real(8) tmnftb(30),rfsetb(30),frtb(30),fltb(30),rdrrtb(30)
      real(8) rdrstb(30),kdif,kdir,laiem,cofab,fstb(30)
      real(8) hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rdctb(22),        &
     &        rlwtb(22)
      real(8) sum,adcrh,adcrl,afgen,cvl,cvr,cvs,eff,rmr
      real(8) perdl,q10,rdc,rdi,rgrlai,rml,rms,span,ssa,tbase,tdwi
      real(8) cumdens(202),oxygenslope(6),oxygenintercept(6)
      real(8) cftb(2*magrs),chtb(2*magrs),cfeictb(2*magrs),siccaplai
      real(8) root_radiusO2,tagp,tagppot,tagprest
      real(8) albedo,rsc,rsw,alphacrit,rooteff,LossGrazingtab(2*100)
      real(8) wiltpoint,rootradius,rootcoefa,LossGrazing(100)
      real(8) dmharvest,dmlastharvest,dateharvest(999),dmgrazing
      real(8) dmmowtb(20),mowrest
      real(8) dmmowdelay(100),DelayRegrowthTab(2*100)
      real(8) lsdb(100),DaysGrazing(100),UptGrazing(100),lsda(366)
      real(8) DaysGrazingtab(2*100),UptGrazingtab(2*100)
      real(8) relmf,wrtmax,saltmax,saltslope,salthead
      real(8) q10_microbial,srl,dry_mat_cont_roots
      real(8) spec_weight_root_tissue,var_a,bdens(maho)
      real(8) specific_resp_humus,air_filled_root_por
      real(8) criterhr,stephr,kroot,rxylem,taccur,kstem
      real(8) rd,rdpot,dvs,tsum,wrt,wrtpot,tadw,rid
      real(8) tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,lai,laipot
      real(8) dwrtpot,dwlv,dwlvpot,dwst,dwstpot,gasst,gasstpot,mrest
      real(8) cwdm,cwdmpot,dwrt,mrestpot
      real(8) laiexppot,laimax,tagpt,tsumgerm,tagptpot,tsum200
      real(8) sla(366),slapot(366),lvage(366),lvagepot(366),lv(366)
      real(8) lvpot(366),atmin7(7),cuptgraz,cuptgrazpot
      logical flsolute,rdinqr,fltsum200,flanthesis
      logical flgrazing,flCropCalendar,flgrazingpot
      character(len=*)   crpfil,pathcrop
      character(len=200) inifil
      logical   flmowingtb
!      logical   flGrassGrowth      ! flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]
      logical   fltsumttd          ! Flag indicating suppress grass growth until criteria for temperature, time and depth are met
      logical   flpotrelmf         ! Flag indicating calculation of limited attainable yield instead of potential yield (due to management factor)
      integer   tsumtime           ! time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
      real(8)   tsumtemp           ! temperature limit to initiate grass growth  [0.0..20.0 grC, R]
      real(8)   tsumdepth          ! depth at which temp above tsumtemp for grass growth [0.0..100.0 cm below soil surface, R]

! locals
      integer   i,nrofSeqGM,ini,count
      integer   swcropharvest,swgrazing,swgrazingpot,swanthesis
      real(8)   dnrinput(magrs),cfinput(magrs),chinput(magrs)
      real(8)   laiinput(magrs),cfeicinput(magrs)
      real(8)   t1900,tstart,tend,cropstart
      integer   swinco
      logical   swSeqGrazMow

      character(len=200) message,filnam
! ----------------------------------------------------------------------


! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)

! --- ET related params  ---------

! --- crop factor or crop height
      call rdsinr ('swcf',1,3,swcf)

! --- check use of crop factors in case of ETref
      if (swetr.eq.1 .and. swcf.eq.2) then
        message = 'If ETref is used (SWETR = 1), always define crop '// &
     &           'factors (SWCF = 1 or 3)' 
        call fatalerr ('ReadGrass',message)
      endif

      if (swcf.eq.1) then
! ---   crop factor is input
        call rdador ('dnr',0.0d0,366.0d0,dnrinput,(magrs),ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---   store values in cftb
        do i = 1,ifnd
          cftb(i*2) = cfinput(i) 
          cftb(i*2-1) = dnrinput(i)
        enddo
        chtb = -99.99d0
      elseif (swcf.eq.2) then
! ---   crop height is input
        call rdador ('dnr',0.0d0,366.0d0,dnrinput,(magrs),ifnd)
        call rdfdor ('ch',0.0d0,1.0d4,chinput,(magrs),ifnd)
! ---   store values in chtb
        do i = 1,ifnd
          chtb(i*2) = chinput(i) 
          chtb(i*2-1) = dnrinput(i)
        enddo
        cftb = -99.99d0
      elseif (swcf.eq.3) then
! ---   LAI-dependent crop factors for dual crop coefficient
        call rdador ('lai',0.0d0,10.0d0,laiinput,36,ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,36,ifnd)
        call rdfdor ('cfeic',0.0d0,2.0d0,cfeicinput,36,ifnd)

        do i = 1,ifnd
          cftb(i*2)        = cfinput(i) 
          cftb(i*2-1)      = laiinput(i)
          cfeictb(i*2)     = cfeicinput(i) 
          cfeictb(i*2-1)   = laiinput(i)
        enddo

! ---   LAI-dependent crop height
        call rdfdor ('ch',0.0d0,1.0d2,chinput,36,ifnd)

        do i = 1,ifnd
          chtb(i*2)     = chinput(i)
          chtb(i*2-1)   = laiinput(i)
        enddo
!
! ---   Interception capacity per unit of LAI
        call rdsdor ('vxiclai',0.0d0,1.0d-2,siccaplai)
        siccaplai = siccaplai
      endif

! --- reflection coefficient and crop resistance
      if (swcf.eq.1 .or. swcf.eq.3) then
! ---   use standard values for ETref
        albedo = 0.23d0
        rsc = 70.0d0
        rsw = 0.0d0
      else
! ---   use crop specific values
        call rdsdor ('albedo',0.0d0,1.0d0,albedo)
        call rdsdor ('rsc',0.0d0,1.0d6,rsc)
        call rdsdor ('rsw',0.0d0,1.0d6,rsw)
      endif


! --- crop growth related params ---------

! --- initial
      call rdsdor ('tdwi',0.0d0,10000.0d0,tdwi)
      call rdsdor ('laiem',0.0d0,10.0d0,laiem)
      call rdsdor ('rgrlai',0.0d0,1.0d0,rgrlai)

! -   yearly start of growth 
! -   using tsum200 
      fltsum200 = .true.
      if(rdinqr('fltsum200')) then
        call rdslog ('fltsum200',fltsum200)
      endif
! -   grass growth initiated by temperature, time and depth
! -   using fltsumttd implies automatic suppress of tsum200
      fltsumttd = .false.
      if(rdinqr('fltsumttd')) then
        call rdslog ('fltsumttd',fltsumttd)
        if (fltsumttd) then                                 
           fltsum200 = .false.
           call rdsdor ('tsumtemp',0.0d0,20.0d0,tsumtemp)    ! temperature limit to initiate grass growth (grC)
           call rdsinr ('tsumtime',1,20,tsumtime)            ! time (nrs of sequential days) with temp above tsumtemp (days)
           call rdsdor ('tsumdepth',0.0d0,100.0d0,tsumdepth) ! depth at which temp above tsumtemp (cm below soil surface)
        endif
      endif

! --- green area
      call rdador ('slatb',0.0d0,366.0d0,slatb,30,ifnd)
      call rdsdor ('ssa',0.0d0,1.0d0,ssa)
      call rdsdor ('span',0.0d0,366.0d0,span)
      call rdsdor ('tbase',-10.0d0,30.0d0,tbase)

! --- assimilation
      call rdsdor ('kdif',0.0d0,2.0d0,kdif)
      call rdsdor ('kdir',0.0d0,2.0d0,kdir)
      call rdsdor ('eff',0.0d0,10.0d0,eff)
      call rdador ('amaxtb',0.0d0,366.0d0,amaxtb,30,ifnd)
      call rdador ('tmpftb',-10.0d0,50.0d0,tmpftb,30,ifnd)
      call rdador ('tmnftb',-10.0d0,50.0d0,tmnftb,30,ifnd)

! --- conversion of assimilates into biomass
      call rdsdor ('cvl',0.0d0,1.0d0,cvl)
      call rdsdor ('cvr',0.0d0,1.0d0,cvr)
      call rdsdor ('cvs',0.0d0,1.0d0,cvs)

! --- maintenance respiration
      call rdsdor ('q10',0.0d0,5.0d0,q10)
      call rdsdor ('rml',0.0d0,1.0d0,rml)
      call rdsdor ('rmr',0.0d0,1.0d0,rmr)
      call rdsdor ('rms',0.0d0,1.0d0,rms)
      call rdador ('rfsetb',0.0d0,366.0d0,rfsetb,30,ifnd)

! --- partitioning
      call rdador ('frtb',0.0d0,366.0d0,frtb,30,ifnd)
      call rdador ('fltb',0.0d0,366.0d0,fltb,30,ifnd)
      call rdador ('fstb',0.0d0,366.0d0,fstb,30,ifnd)

! --- death rates
      call rdsdor ('perdl',0.0d0,3.0d0,perdl)
      call rdador ('rdrrtb',0.0d0,366.0d0,rdrrtb,30,ifnd)
      call rdador ('rdrstb',0.0d0,366.0d0,rdrstb,30,ifnd)

! --- oxygen stress
      swoxygen = 1
      if(rdinqr('swoxygen')) then
        call rdsinr ('swoxygen',0,2,swoxygen)
! -     fatal error when physical oxygen stress is simulated without numerical heat flow
        if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1)) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), soil heat flow should be'//&
     &   ' numerically simulated: SwCalT=2! Adapt .swp input file.'
          call fatalerr ('readgrass',message)
        endif
! -     fatal error when physical oxygen stress is simulated without realistic bdens value
        if (swoxygen.eq.2 .and. bdens(1).lt.100.d0) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), bulk density must have'//  &
     &   ' realistic values; adjust BDENS-value(s) in .swp input file.'
          call fatalerr ('readgrass',message)
        endif
      endif

      if (swoxygen.eq.1) then                                          !
! ---   oxygen stress according to Feddes
        call rdsdor ('hlim1' ,-100.0d0,100.0d0,hlim1)                   !
        call rdsdor ('hlim2u',-1000.0d0,100.0d0,hlim2u)                 !
        call rdsdor ('hlim2l',-1000.0d0,100.0d0,hlim2l)                 !
      endif

      if (swoxygen.eq.2) then                                          
! ---   oxygen stress according to Bartholomeus

        swoxygentype = 1
        if(rdinqr('swoxygentype')) then
          call rdsinr ('swoxygentype',1,2,swoxygentype)
        endif

        if (swoxygentype.eq.1) then                                 
! ---      use physical processes
           call rdsdor ('q10_microbial',1.0d0,4.0d0,q10_microbial)  
           call rdsdor ('specific_resp_humus',0.0d0,1.0d0,              &
     &                   specific_resp_humus)                        
           call rdsdor ('srl',0.0d0,1.0d10,srl)                        
           call rdsinr ('swrootradius',1,2,swrootradius)
           if (swrootradius.eq.1) then
              call rdsdor ('dry_mat_cont_roots',0.0d0,1.0d0,            &
     &                   dry_mat_cont_roots)                        
              call rdsdor ('air_filled_root_por',0.0d0,1.0d0,           &
     &                   air_filled_root_por)                        
              call rdsdor ('spec_weight_root_tissue',0.0d0,1.0d5,       &
     &                   spec_weight_root_tissue)                   
              call rdsdor ('var_a',0.0d0,1.0d0,var_a)
           else
              call rdsdor ('root_radiusO2',0.0d0,1.0d0,root_radiusO2)
           endif  
        else
! ---      use reproduction functions
           call rdsinr ('SwTopSub',1,2,SwTopSub)
           call rdsinr ('NrStaring',1,18,NrStaring)
           call oxygen_dat (SwTopSub,NrStaring,OxygenSlope,             &
     &                      OxygenIntercept)
        endif
      endif

! --- drought stress
      swdrought = 1
      if(rdinqr('swdrought')) then
        call rdsinr ('swdrought',1,2,swdrought)                         
      endif

      if (swdrought.eq.1) then                                          
! ---    drought stress according to Feddes
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
         call rdsdor ('adcrh',0.0d0,5.0d0,adcrh)                        
         call rdsdor ('adcrl',0.0d0,5.0d0,adcrl)                        
! --     Criticial stress index for compensation of root water uptake (-)
         alphacrit = 1.0d0
         if(rdinqr('alphacrit')) then
           call rdsdor ('alphacrit',0.2d0,1.0d0,alphacrit)
         endif
      else                                                              !
! ---    drought stress according to De Jong van Lier
         call rdsdor ('wiltpoint',-1.0d8,-1.0d2,wiltpoint)
         call rdsdor ('kstem',1.0d-10,1.0d1,kstem)               
         call rdsdor ('rxylem',1.0d-4,1.d0,rxylem)               
         call rdsdor ('rootradius',1.0d-4,1.0d0,rootradius)            
         call rdsdor ('kroot',1.0d-10,1.0d10,kroot)            
         call rdsdor ('rootcoefa',0.0d0,1.0d0,rootcoefa)
         call rdsinr ('swhydrlift',0,1,swhydrlift)
         call rdsdor ('rooteff',0.0d0,1.0d0,rooteff)                 
         call rdsdor ('stephr',0.0d0,10.d0,stephr)                 
         call rdsdor ('criterhr',0.0d0,10.d0,criterhr)                 
         call rdsdor ('taccur',0.0d-5,10.d-2,taccur)                 
      endif

! --- salt stress
      if (flsolute) then

        if(rdinqr('swsalinity')) then
          call rdsinr ('swsalinity',0,2,swsalinity)
        endif
 
        if (swsalinity .eq. 1) then
! ---     input for Maas and Hoffman salt reduction function
          call rdsdor ('saltmax',0.0d0,100.0d0,saltmax)
          call rdsdor ('saltslope',0.0d0,1.0d0,saltslope)
        elseif (swsalinity .eq. 2) then
! ---     osmotic head salinity stress concept
          call rdsdor ('salthead',0.0d0,1.0d3,salthead)
          if (swdrought.eq.1) then
           message = 'In case salinity stress is calculated with'//     &
     &     ' osmotic head (SwSalinity = 2), the drought stress'//       &
     &     ' should be calculated according to De Jong van Lier:'//     &
     &     ' SwDrought = 2! Adapt .crp input file.'
           call fatalerr ('readcropfixed',message)
          endif
        endif
      endif

! --- interception
      if(swcf.ne.3) then
        call rdsdor ('cofab',0.0d0,1.0d0,cofab)
      endif

! --- rooting
      call rdador ('rdctb',0.0d0,100.0d0,rdctb,22,ifnd)
      call rdador ('rlwtb',0.0d0,5000.0d0,rlwtb,22,ifnd)
      call rdsdor ('rdi',0.0d0,1000.0d0,rdi)
      call rdsdor ('rdc',0.0d0,1000.0d0,rdc)


! --- management factor to account for other forms of stress (pests, diseases, nutrients, etc ..)
      relmf = 1.0d0
      if(rdinqr('relmf')) then
        call rdsdor ('relmf',0.0d0,1.0d0,relmf) 
      endif

! --- factor to account for difference between theoretical potential and 
!     attainable yield according to ASG field-experiments
      flpotrelmf = .false.
      if(rdinqr('flpotrelmf')) then
        call rdslog ('flpotrelmf',flpotrelmf)
      endif

! --- Management of mowing and grazing
! -   0. Define periods with mowing or grazing 
!     read (yearly) sequence of periods with mowing-grazing (1=grazing, 2=mowing)
      call rdainr ('SeqGrazMow',1,2,SeqGrazMow,366,nrofSeqGM)
      
! -   1. Mowing settings
      swSeqGrazMow = .false.
      do i = 1,nrofSeqGM
         if (SeqGrazMow(i) .eq. 2) swSeqGrazMow = .true.
      end do

      if (swSeqGrazMow) then
        call rdsdor ('mowrest',0.0d0,5000.0d0,mowrest)
        call rdsinr ('swharvest',1,2,swharvest)
        if(swharvest.eq.1) then
          flmowingtb = .false.
          if (rdinqr('dmmowtb')) then
            flmowingtb = .true.
            call rdador ('dmmowtb',0.0d0,100000.0d0,dmmowtb,20,ifnd)
            call rdsinr ('maxdaymow',1,366,maxdaymow)
          else
            call rdsdor ('dmharvest',0.0d0,100000.0d0,dmharvest)
            call rdsinr ('daylastharvest',1,366,daylastharvest)
            call rdsdor ('dmlastharvest',0.0d0,100000.0d0,dmlastharvest)
          endif
        elseif(swharvest.eq.2) then
          call rdatim ('dateharvest',dateharvest,999,ifnd)
          dateharvest(ifnd + 1) = tend + 1.d0
        endif
        
        ! days delay of regrowth after mowing
        call rdainr ('daydelay',0,366,daydelay,100,ifnd)
        call rdador ('dmmowdelay',0.0d0,100000.0d0,dmmowdelay,100,ifnd)
        
        ! store values in DelayRegrowthTab
        do i = 1,ifnd
          DelayRegrowthTab(i*2)   = dble(daydelay(i))
          DelayRegrowthTab(i*2-1) = dmmowdelay(i)
        enddo
      
      endif

! -   2. Grazing settings
      swSeqGrazMow = .false.
      do i = 1,nrofSeqGM
         if (SeqGrazMow(i) .eq. 1) swSeqGrazMow = .true.
      end do

      if (swSeqGrazMow) then
         call rdsinr ('swharvest',1,2,swharvest)
         if(swharvest.eq.1) then
            call rdsdor ('dmgrazing',0.0d0,100000.0d0,dmgrazing)
         elseif(swharvest.eq.2) then
           call rdatim ('dateharvest',dateharvest,999,ifnd)
           dateharvest(ifnd + 1) = tend + 1.d0
         endif

        call rdsdor ('tagprest',0.0d0,100000.0d0,tagprest)
!       LiveStock Density basal with days and uptake during grazing
        call rdador ('LSDb',0.0d0,1000.0d0,LSDb,100,ifnd)
        call rdfdor ('daysgrazing',0.0d0,366.0d0,daysgrazing,100,ifnd)
        call rdfdor ('uptgrazing',0.0d0,1000.0d0,uptgrazing,100,ifnd)
        call rdfdor ('lossgrazing',0.0d0,1000.0d0,lossgrazing,100,ifnd)
!       store values in daysgrazingtab, uptgrazingtab and lossgrazingtab 
        do i = 1,ifnd
          daysgrazingtab(i*2)   = daysgrazing(i)
          daysgrazingtab(i*2-1) = LSDb(i)
        enddo
        do i = 1,ifnd
          uptgrazingtab(i*2)   = uptgrazing(i)
          uptgrazingtab(i*2-1) = LSDb(i)
        enddo
        do i = 1,ifnd
          lossgrazingtab(i*2)   = lossgrazing(i)
          lossgrazingtab(i*2-1) = LSDb(i)
        enddo
!       LiveStock Density actual
        call rdador ('LSDa',0.0d0,1000.0d0,LSDa,366,ifnd)
!       verify length of array, should correspond to periods with grazing
        if (ifnd.ne.nrofSeqGM) then
          message = 'Dynamic grassland, input for LSDa'//               &
     &   ' (LiveStock Density actual); nr of values must be equal'//    &
     &   ' to the nr of values in SEQGRAZMOW! Adapt .crp input file.'
        endif
      endif

! --- maximum weight of roots (kg/ha dm)
      call rdsdor ('wrtmax',0.0d0,100000.0d0,wrtmax)

! --- determine whether irrigation scheduling is applied
      call rdsinr ('schedule',0,1,schedule)

      if (schedule.eq.1 .and. swdrought.eq.2) then
! ---    read limiting pressure heads for irrigation scheduling
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
      endif

      close (crp)

! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION

      if (swdrought.eq.1) then 

! ---   specify array ROOTDIS with root density distribution
        do i = 0,100
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
        enddo

! ---   calculate cumulative root density function
        do i = 1,202,2
! ---     relative depths
          cumdens(i) = rootdis(i)
        enddo
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0                 &
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
        enddo

! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo

      endif

! --- read crop data of former day from *.END file

      if (t1900 - tstart .lt. 1.d-3 .and. swinco .eq. 3 .and. .not.     &
     &   (t1900 - cropstart .gt. -1.d-3 .and.                           &
     &                t1900 - cropstart .lt. 1.d-3)) then

        if (flCropCalendar) then  
          ini = getun2 (10,90,2)
          call rdinit(ini,logf,inifil)

          call rdsdor ('rd',0.0d0,1.0d4,rd)               
          call rdsdor ('rdpot',0.0d0,1.0d4,rdpot)               
          call rdsdor ('dvs',0.0d0,3.0d0,dvs)               
          call rdsinr ('daycrop',0,366,daycrop)               
          call rdsinr ('swanthesis',0,1,swanthesis)               
          call rdsdor ('tsum',0.0d0,1.0d4,tsum)               
          call rdsinr ('ilvold',0,366,ilvold)               
          call rdsinr ('ilvoldpot',0,366,ilvoldpot)               
          call rdsdor ('wrt',0.0d0,1.0d8,wrt)               
          call rdsdor ('wrtpot',0.0d0,1.0d8,wrtpot)               
          call rdsdor ('tadw',0.0d0,1.0d8,tadw)               
          call rdsdor ('tadwpot',0.0d0,1.0d8,tadwpot)               
          call rdsdor ('wst',0.0d0,1.0d8,wst)               
          call rdsdor ('wstpot',0.0d0,1.0d8,wstpot)               
          call rdsdor ('wso',0.0d0,1.0d8,wso)               
          call rdsdor ('wsopot',0.0d0,1.0d8,wsopot)               
          call rdsdor ('wlv',0.0d0,1.0d8,wlv)               
          call rdsdor ('wlvpot',0.0d0,1.0d8,wlvpot)               
          call rdsdor ('laiexp',0.0d0,1.0d8,laiexp)               
          call rdsdor ('lai',0.0d0,1.0d8,lai)
          call rdsdor ('laipot',0.0d0,1.0d8,laipot)
          call rdsdor ('dwrt',0.0d0,1.0d8,dwrt)
          call rdsdor ('dwrtpot',0.0d0,1.0d8,dwrtpot)
          call rdsdor ('dwlv',0.0d0,1.0d8,dwlv)
          call rdsdor ('dwlvpot',0.0d0,1.0d8,dwlvpot)
          call rdsdor ('dwst',0.0d0,1.0d8,dwst)
          call rdsdor ('dwstpot',0.0d0,1.0d8,dwstpot)
          call rdsdor ('gasst',0.0d0,1.0d10,gasst)
          call rdsdor ('gasstpot',0.0d0,1.0d10,gasstpot)
          call rdsdor ('mrest',0.0d0,1.0d10,mrest)
          call rdsdor ('mrestpot',0.0d0,1.0d10,mrestpot)
          call rdsdor ('cwdm',0.0d0,1.0d10,cwdm)
          call rdsdor ('cwdmpot',0.0d0,1.0d10,cwdmpot)
          call rdsdor ('tsumgerm',0.0d0,1.0d4,tsumgerm)
          call rdsinr ('nofd',0,366,nofd)               
          call rdfdor ('atmin7',-50.d0,60.d0,atmin7,7,7)
          call rdsdor ('rid',0.0d0,1.0d4,rid)
          call rdsinr ('idregr',0,366,idregr)               
          call rdsinr ('idregrpot',0,366,idregrpot)               
          call rdsdor ('laiexppot',0.0d0,1.0d2,laiexppot)
          call rdsdor ('laimax',0.0d0,1.0d2,laimax)
          call rdsdor ('tagp',0.0d0,1.0d5,tagp)
          call rdsdor ('tagppot',0.0d0,1.0d5,tagppot)
          call rdsdor ('tagpt',0.0d0,1.0d5,tagpt)
          call rdsdor ('tagptpot',0.0d0,1.0d5,tagptpot)
          call rdsinr ('daymow',0,366,daymow)
          call rdsinr ('daymowpot',0,366,daymowpot)
          call rdsdor ('cuptgraz',0.0d0,1.0d5,cuptgraz)
          call rdsdor ('cuptgrazpot',0.0d0,1.0d5,cuptgrazpot)
          call rdsdor ('tsum200',0.0d0,1.0d4,tsum200)
          call rdsinr ('swcropharvest',0,1,swcropharvest)               
          call rdsinr ('iseqgm',0,366,iseqgm)               
          call rdsinr ('iseqgmpot',0,366,iseqgmpot)               
          call rdsinr ('swgrazing',0,1,swgrazing)
          call rdsinr ('swgrazingpot',0,1,swgrazingpot)
          call rdsinr ('iharvest',0,366,iharvest)
          call rdsinr ('idaysgraz',0,366,idaysgraz)
          call rdsinr ('idaysgrazpot',0,366,idaysgrazpot)
          
          count = max(ilvold,ilvoldpot)
          call rdfdor ('sla',0.d0,1.0d3,sla,366,count)
          call rdfdor ('slapot',0.d0,1.0d3,slapot,366,count)
          call rdfdor ('lvage',0.d0,1.0d3,lvage,366,count)
          call rdfdor ('lvagepot',0.d0,1.0d3,lvagepot,366,count)
          call rdfdor ('lv',0.d0,1.0d5,lv,366,count)
          call rdfdor ('lvpot',0.d0,1.0d5,lvpot,366,count)
          
          if (swanthesis .eq. 0) then
            flanthesis = .false.
          else
            flanthesis = .true.
          endif
          if (swgrazing .eq. 0) then
            flgrazing = .false.
          else
            flgrazing = .true.
          endif
          if (swgrazingpot .eq. 0) then
            flgrazingpot = .false.
          else
            flgrazingpot = .true.
          endif
        
        else
          
          daycrop = 0
          
        endif  
!       close inifile
        close(ini)
      
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine rddre (drfil,pathdrain,nrlevs,nrpri,nrsec,l,           &
     & zbotdr,widthr,taludr,rdrain,rinfi,rentry,rexit,gwlinf,swdtyp,    &
     & wlptab,swsec,wls1,osswlm,nmper,wlstar,impend,swman,wlsman,       &
     & wscap,swqhr,hbweir,alphaw,betaw,nqh,hqhtab,qqhtab,dropr,         &
     & gwlcrit,wlstab,sttab,swstini,swst,wlsbak,swsrf,nphase,hcrit,     &
     & hdepth,vcrit,wlp1,nodhd,numnod,dz,wldip,numadj,logf,intwl,       &
     & swnrsrf,rsurfdeep,rsurfshallow,t1900,cofintfl,expintfl,swscre,   &
     & swdivd,cofani,numlay,tstart,tend,swdislay,swtopdislay,ztopdislay,&
     & ftopdislay,SwTopnrsrf,                                           & 
     & swdivdinf,FacDpthInf)

! ----------------------------------------------------------------------
!     UpDate             : 20080109
!     Date               : 20010605                   
!     Purpose            : reading multilevel drainage characteristics 
!                          and specification of the surface water system
!                          for a period up to one year;  
!
! --- 1 Reading extended drainage input from .dra-file
! --- 2 Initializations
!
! --- Initializations:
! -1- wlp1: water level in primary system   (SWSRF = 3)
! -2- wls1: water level in secondary system (SWSRF = 2 or 3, SWSEC = 1 or 2) 
! -3- HBWEIR(IMPER) in case of table discharge relation (SWQHR = 2)
! -4- NUMADJ: number of target level adjustments
! -5- sttab(22,2) table with storage as a function of water level
! ---    sttab(i,1) contains levels:
! ---    i=1: +100 cm; i=2: 0 cm; i=22: bottom of deepest dr. medium
! ---    sttab(i,2) contains storage expressed as surface layer 
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer   nrpri,nrsec,nrlevs,swdtyp(Madr),nmper,swman(mamp)
      integer   swqhr,nqh(mamp),swsec,swsrf,nphase(mamp)
      integer   nodhd(mamp),numnod,numadj,logf,intwl(mamp),swnrsrf
      integer   swdivd,numlay,swscre, swdivdinf
      real(8)   l(Madr),zbotdr(Madr),widthr(Madr),taludr(Madr)
      real(8)   rdrain(Madr),rinfi(Madr),rentry(Madr),rexit(Madr)
      real(8)   gwlinf(Madr),wlptab(2*mawlp)
      real(8)   rsurfdeep,rsurfshallow,impend(mamp),t1900
      real(8)   wls1,osswlm,wlstar,wscap(mamp),hbweir(mamp),alphaw(mamp)
      real(8)   hqhtab(mamp,mamte),qqhtab(mamp,mamte)
      real(8)   betaw(mamp),dropr(mamp*mamte),hdepth(mamp*mamte)
      real(8)   wlsman(mamp,mamte),wlstab(2*mawls),sttab(22,2),wlsbak(4)
      real(8)   swstini,swst,hcrit(mamp,mamte),vcrit(mamp,mamte)
      real(8)   afgen,wlp1,dz(macp),gwlcrit(mamp,mamte),wldip(mamp)
      real(8)   cofintfl,expintfl,cofani(maho)
      real(8)   tend,tstart,FacDpthInf
      integer   swdislay,swtopdislay(madr),SwTopnrsrf
      real(8)   ztopdislay(Madr),ftopdislay(madr)
      character(len=16) drfil
      character(len=80) pathdrain

! --- local
      integer   dra,level(Madr),itab,i,ilev,getun2
      integer   iph,nrman1,nrman2,node,imper,imperb,imperi,ifnd
      integer   imper_4b(mamp),imper_4c(mamp),imper_4d(mamp)
      integer   imper_4e1(mamp*mamte),imper_4e2(mamp*mamte)
      integer   imptab(mamte),impphase(mamp*mamte), nmper2
      real(8)   wlp(mawlp),wls(mawlp),zb, hhtab(mamp),qhtab(mamp)
      real(8)   wlsman_t(mamp*mamte),gwlcrit_t(mamp*mamte)
      real(8)   hcrit_t(mamp*mamte),vcrit_t(mamp*mamte)
      real(8)   dates(mabbc)
      real(8)   wdepth,wbreadth,wvolum,dep,swstlev,sofcu,altcu
      logical   flweir(mamp),exists(mamp),flzero(mamp), rdinqr
      character(len=80)  filnam
      character(len=200) messag
      real(8)   small
      data      small     /0.0001d0/

! ----------------------------------------------------------------------

! --- open file with extended drainage data
      filnam = trim(pathdrain)//trim(drfil)//'.dra'
      dra = getun2 (10,90,2)
      call rdinit(dra,logf,filnam)

! --- division of drainage fluxes
      call rdsinr ('swdivd',0,1,swdivd)
      if (swdivd.eq.0) then
          write(messag,'(3a)')                                          &
     &    ' Variabel SWDIVD=0 in input file : ',trim(filnam),           &
     &    ' this is not recommended and may cause numerical instability'
          call warn ('rddre',messag,logf,swscre)
      endif
      if (swdivd .eq. 1) then
        if(rdinqr('swdivdinf')) then
          call rdsinr ('swdivdinf',0,1,swdivdinf)
          call rdsdor ('FacDpthInf',0.d0,1.d0,FacDpthInf)
        else
          swdivdinf = 0
        endif

        call rdfdor ('cofani',0.d0,1000.d0,cofani,maho,numlay)
      endif

! --- altitude of control unit (relative to reference level)
      call rdsdor ('altcu',-3.0d5,3.0d5,altcu)

! --- part 1

! --  number of drainage levels
      call rdsinr ('nrsrf',1,Madr,nrlevs)
      if (nrlevs.gt.5) then
          write(messag,'(3a)')                                          &
     &    ' Number of Drainage levels >5 in inputfile : ',trim(filnam), &
     &    '   part of the output is limited to 5 levels'
          call warn ('rddre',messag,logf,swscre)
      endif

! --  top of model dicharge layer, determined by factor or direct input
      if (swdivd .eq. 1) then
        call rdsinr ('swdislay',0,2,swdislay)
        if (swdislay .eq. 1) then
          call rdfinr ('swtopdislay',0,1,swtopdislay,madr,nrlevs)
          call rdfdor('ztopdislay',-1.0d4,0.0d0,ztopdislay,madr,nrlevs)
        elseif (swdislay .eq. 2) then
          call rdfinr ('swtopdislay',0,1,swtopdislay,madr,nrlevs)
          call rdfdor('ftopdislay',0.0d0,1.0d0,ftopdislay,madr,nrlevs)
        endif
      endif


! --- characteristics of each drainage level 
      call rdfinr ('lev',1,Madr,level,Madr,nrlevs)
      call rdfinr ('swdtyp',0,1,swdtyp,Madr,nrlevs)
      call rdfdor ('l',1.0d0,100000.0d0,l,Madr,nrlevs)
      call rdfdor ('zbotdre',(altcu-1.0d3),(altcu-1.0d-2),              &
     &              zbotdr,Madr,nrlevs)
      call rdfdor ('gwlinf',-10000.0d0,0.0d0,gwlinf,Madr,nrlevs)
      call rdfdor ('rdrain',1.0d0,1.0d5,rdrain,Madr,nrlevs)
      call rdfdor ('rinfi',1.0d0,1.0d5,rinfi,Madr,nrlevs)
      call rdfdor ('rentry',0.0d0,10.0d0,rentry,Madr,nrlevs)
      call rdfdor ('rexit',0.0d0,10.0d0,rexit,Madr,nrlevs)
      call rdfdor ('widthr',0.0d0,10000.0d0,widthr,Madr,nrlevs)
      call rdfdor ('taludr',1.0d-2,5.0d0,taludr,Madr,nrlevs)

! -   conversions and security checks....
      do i = 1, nrlevs
!       conversions
        l(i) = l(i)*100.0d0
        zbotdr(i) = zbotdr(i) -altcu
        if (swdivd.eq.1 .and. swdislay.eq.1) then
          ztopdislay(i) = ztopdislay(i) - altcu
          if(ztopdislay(i).lt.zbotdr(i)) then
            messag = 'Ztopdislay cannot be below Zbotdre, verify input!' 
            call fatalerr ('Rddra',messag)
          endif
        endif
! -     security checks....
! -     0 - verify levels-index
! -     1 - levels must be ordered, starting with the deepest
! -     2 - zbotdr always below surface level
! -     3 - gwlinf must be below bottom of deepest drainage medium
! -     4 - widthr must be greater than zero for non-tube drain systems
        if (level(i).ne.i) then
          messag = 'Drainage level index is not consistent' 
          call fatalerr ('Rddra',messag)
        endif
        if (i .gt. 1) then
          if (zbotdr(i) .lt. zbotdr(i-1)) then
            messag = 'Levels must be ordered, starting with deepest' 
            call fatalerr ('Rddra',messag)
          endif
        endif
        if (gwlinf(i).gt.zbotdr(i)) then
            messag = 'Gwlinf should be lower than Zbotdr !' 
            call fatalerr ('Rddra',messag)
        endif
        if (swdtyp(i).eq.0 .and.widthr(i).lt.small) then
            write(messag,'(a,f10.5)')                                   &
     &              'widthr must be greater than',small
            call fatalerr ('Rddra',messag)
        endif

      enddo

! --- type of highest drainage level
      call rdsinr ('swnrsrf',0,2,swnrsrf)
      if (swnrsrf.eq.1) then
        call rdsdor ('rsurfdeep',0.001d0,1000.0d0,rsurfdeep)
        call rdsdor ('rsurfshallow',0.001d0,1000.0d0,rsurfshallow)
      else if (swnrsrf.eq.2) then
        call rdsdor ('cofintfl',0.01d0,10.0d0,cofintfl)     
        call rdsdor ('expintfl',0.1d0,1.0d0,expintfl) 
      endif


! kroes 20080707 : allow disabling of option
      if (swdivd .eq. 1) then
        if (swnrsrf.gt.0)  call rdsinr ('SwTopnrsrf',0,1,SwTopnrsrf)
      endif


! --- part 2a
      call rdsinr ('swsrf',1,3,swsrf)
      if (swsrf.eq.3) then
        nrpri = 1
        nrsec = nrlevs-nrpri
      elseif (swsrf.eq.2) then
        nrpri = 0
        nrsec = nrlevs-nrpri
      elseif (swsrf.eq.1) then
! ---   no surface water system...ready
        close (dra)
        return
      endif

      if (swdtyp(1+nrpri).ne.0) then
        messag = 'Deepest sec. level must be open.' 
        call fatalerr ('Rddra',messag)
      endif

! --- part 2b

! --- read table with water levels in the primary system (SWSRF=3)
      if (swsrf.eq.3) then
! ---   init table
        do 8 i = 1,2*mawlp
          wlptab(i) = 0.0d0
    8   continue
! -     read and store
        call rdatim ('date1',dates,mawlp,ifnd)
! -     at least one date must be within simulation period
        call checkdate(ifnd,dates,tend,tstart,'date1 ',                 &
     &                 'readswap//swsrf=3        ')
        call rdfdor                                                     &
     &        ('wlp',(altcu-1000.0d0),(altcu+200.0d0),wlp,mawlp,ifnd)
        do i = 1, ifnd
          wlptab(i*2) = wlp(i) - altcu
          wlptab(i*2-1) = dates(i)
        enddo
! ---   set initial value
        wlp1 = afgen (wlptab,2*mawlp,t1900-1.d0)

! ---   Ready in case of only a primary system ??
!       if (NRSEC.LT.1) then
!          CLOSE(DRE)
!          Return
!       endif 

      endif

! --- part 2c

      call rdsinr ('swsec',1,2,swsec)

! --- part 3

      if (swsec.eq.1) then
! ---   surface water level of secondary system is input
! ---   position file pointer

! ---   init table
        do 60 i = 1,2*mawls
          wlstab(i) = 0.0d0
   60   continue
! -     read and store
        call rdatim ('date2',dates,mawlp,ifnd)
! -     at least one date must be within simulation period
        call checkdate(ifnd,dates,tend,tstart,'date2 ',                 &
     &                 'readswap//swsec=1        ')

        call rdfdor                                                     &
     &         ('wls',(altcu-1000.0d0),(altcu+200.0d0),wls,mawlp,ifnd)
        do i = 1, ifnd
          wlstab(i*2) = wls(i) -altcu
          wlstab(i*2-1) = dates(i)
        enddo
! ---   set initial value
        wls1 = afgen (wlstab,2*mawls,t1900-1.0d0)

! ---   part 4a

! --- surface water level of secondary system is simulated
      elseif (swsec.eq.2) then
        call rdsdor ('wlact',dble(zbotdr(1+nrpri)+altcu)                &
     &       ,dble(altcu),wls1)
        wls1 = wls1 - altcu
        call rdsdor ('osswlm',0.0d0,10.0d0,osswlm)

! ---   part 4b

        call rdsinr ('nmper',1,mamp,nmper)
        wlstar = wls1

        call rdfinr ('imper_4b',1,nmper,imper_4b,mamp,nmper)
        call rdftim ('impend',impend,mamp,nmper)
        call rdfinr ('swman',1,2,swman,mamp,nmper)
        call rdfdor ('wscap',0.0d0,10.0d0,wscap,mamp,nmper)
        call rdfdor ('wldip',0.0d0,100.0d0,wldip,mamp,nmper)
        call rdfinr ('intwl',1,31,intwl,mamp,nmper)

! -     for each type of management: count number of periods 
        nrman1 = 0
        nrman2 = 0
        do imper = 1, nmper
          if (swman(imper).eq.2 .and. intwl(imper).lt.1) then
            messag = 'intwl (management interval) must be > = 1d'
            call fatalerr ('Rddra',messag)
          endif
          wldip(imper) = abs(wldip(imper))
          if (swman(imper).eq.1) then 
            nrman1 = nrman1+1
          elseif (swman(imper).eq.2) then
            nrman2 = nrman2+1
          else
            messag = 'SWMAN is out of range.'
            call fatalerr ('Rddra',messag)
          endif
        enddo
        if ((nrman1+nrman2).ne.nmper) then
          messag = 'nrman1+nrman2 does not match nmper'
          call fatalerr ('Rddra',messag)
        endif

! ---   type of discharge relationship
        call rdsinr ('swqhr',1,2,swqhr)

        IF (SWQHR.EQ.1) THEN

! ---     part 4c

          call rdsdor ('sofcu',0.1d0 ,100000.0d0,sofcu)

          call rdfinr ('imper_4c',1,nmper,imper_4c,mamp,nmper)
          zb = min(zbotdr(1),zbotdr(2))
          call rdfdor                                                   &
     &           ('hbweir',(altcu+zb),(altcu+100.0d0),hbweir,mamp,nmper)
          call rdfdor ('alphaw',0.1d0,50.0d0,alphaw,mamp,nmper)
          call rdfdor ('betaw',0.5d0,3.0d0,betaw,mamp,nmper)

! ---     convert and check values
!         initialise
          imperb = 0
          imperi = 0
          do imper = 1, nmper
            hbweir(imper) = hbweir(imper)-ALTCU
! ---       correction for units
            alphaw(imper) = alphaw(imper)* (8.64d0* 100.0d0**           &
     &                      (1.0d0-betaw(imper))/SOFCU)     
            if (hbweir(imper).lt.zbotdr(1+NRPRI)) then
              messag = 'Weir crest level is below bottom of'            &
     &                 //' deepest channel of secondary system.'
              call fatalerr ('Rddra',messag)
            endif
! --- check for target level above channel bottom when supply is 
! --- attempted (system may never become dry in this case)
            if (swman(imper).eq.1 .and. wscap(imper).gt.1.d-7 .and.     &
     &         (hbweir(imper)-wldip(imper)).lt.                         &
     &         (zbotdr(1+NRPRI)+1.d-4)) then
              messag = 'HBWEIR/WLDIP !'                                 &
     &                 //' Supply not possible =< zbotdr !'
              call fatalerr ('Rddra',messag)
            endif

            if (imper.ne.imperb) then
              imperb = imper
              imperi = imperi + 1
            else
              messag = '4c imper not unique'
              call fatalerr ('Rddra',messag)
            endif
          enddo

! ---     check number of records.. 
          if (imperi.ne.nmper) then
              messag = 'part 4c - not enough records'
              call fatalerr ('Rddra',messag)
          endif

        elseif (swqhr.eq.2) then

! ---     part 4d

! --      initialise
          imperb = 0
          imperi = 0
          do 34 imper = 1,nmper
            nqh (imper) = 0
            flweir (imper) = .false.
            exists (imper) = .false.
            flzero (imper) = .false.
34        continue

! --      read and store
          call rdainr ('imper_4d',1,nmper,imper_4d,mamp,ifnd)
          call rdfinr ('imptab',1,mamte,imptab,mamte,ifnd)
          call rdfdor                                                   &
     &         ('htab',(altcu-1000.0d0),(altcu+100.0d0),hhtab,mamp,ifnd)
          call rdfdor ('qtab',0.0d0,500.0d0,qhtab,mamp,ifnd)

! ---     convert
          do i = 1, ifnd
            hqhtab(imper_4d(i),imptab(i)) = hhtab(i) - Altcu
            qqhtab(imper_4d(i),imptab(i)) = qhtab(i)
          enddo

! ---     check values
          do i = 1, ifnd
            imper = imper_4d(i)
            itab = imptab(i)
            nqh(imper) = nqh(imper)+1
            exists(imper) = .true.
            if (qqhtab(imper,itab).lt.0.000001d0) flzero(imper)=.true.
            if (hqhtab(imper,itab).lt.zbotdr(1+nrpri)) then
              messag = 'Level in q-h table below bottom of'             &
     &                 //' deepest channel of secondary system.'
              call fatalerr ('Rddra',messag)
            endif
            if (imper.ne.imperb) then
              imperb = imper
              imperi = imperi + 1
            endif

! -         establish hbweir (at level where qtab = 0) 
            if (qqhtab(imper,itab).lt.0.0001d0.and..not.flweir(imper))  &
     &      then
              hbweir(imper) = hqhtab(imper,itab)
              flweir(imper) = .true.
            endif 

! -         consistency checks
            if (nqh(imper).ne.itab) then
              messag = 'qh-table / imper - itab mismatch'
              call fatalerr ('Rddra',messag)
            endif
            if (itab.eq.1) then
              if (abs(hqhtab(imper,itab)-100.0d0).gt.0.00001d0) then
              messag = 'First value in htab should be altcu+100.0'
              call fatalerr ('Rddra',messag)
              endif
            endif
            if (itab.gt.1) then
              if ((hqhtab(imper,itab).ge.hqhtab(imper,itab-1)).or.      &
     &           (qqhtab(imper,itab).gt.qqhtab(imper,itab-1))) then
                 messag = 'qh-table - no descending values'
                 call fatalerr ('Rddra',messag)
              endif
            endif

          enddo

! ---     check number of periods.. 
          if (imperi.ne.nmper) then
            messag = 'part 4d - number of periods incorrect.'
            call fatalerr ('Rddra',messag)
          endif

! ---     check that QQHTAB goes down to zero
          do imper = 1,nmper
            if (exists(imper).and..not.flzero(imper)) then
              messag = 'qqhtab not going down to zero.'
              call fatalerr ('Rddra',messag)
            endif
          enddo

        endif

! ---   part 4e1

        if (nrman2.gt.0) then

          imperb = 0
          imperi = 0

! ---     read table with drop rates (length of table must be equal 
!            to the nr of management periods referred to in tabel 4e2
          call rdainr('imper_4e1',1,nmper,imper_4e1,(mamp*mamte),nmper2)
          call rdfdor('dropr',0.0d0,100.0d0,dropr,(mamp*mamte),nmper2)
          call rdfdor('hdepth',-100d0,0.0d0,hdepth,(mamp*mamte),nmper2)

          do i = 1, nmper2
            imper = imper_4e1(i)
            hdepth(i) = -abs(hdepth(i))

! ---       determine compartment number related to hdepth
            dep = 0.0d0
            node = 1
            do while (hdepth(i).lt.(dep-1.0d-6).and.node.le.numnod) 
              nodhd(imper) = node
              dep = dep-dz(node)
              node = node + 1
            enddo

            if (swman(imper) .ne. 2) then
              messag = '#4e1 swman - imper mismatch.'
              call fatalerr ('Rddra',messag)
            endif
            if (imper.ne.imperb) then
              imperb = imper
              imperi = imperi + 1
            else
              messag = 'Two drop rates at the same period.'
              call fatalerr ('Rddra',messag)
            endif
          enddo
          if (imperi .ne. nrman2) then
            messag = '#4e1 number of periods for drop rate incorrect.'
            call fatalerr ('Rddra',messag)
          endif

! ---   part 4e2

          imperb = 0
          imperi = 0

          do 54 imper = 1,nmper
            nphase(imper) = 0
   54     continue

! --      read and store
          call rdainr('imper_4e2',1,nmper,imper_4e2,(mamte*mamp),ifnd)
          call rdfinr('impphase',1,mamte,impphase,(mamte*mamp),ifnd)
          call rdfdor                                                   &
     &       ('wlsman',(altcu-500.0d0),altcu,wlsman_t,(mamp*mamte),ifnd)
          call rdfdor                                                   &
     &       ('gwlcrit',-500.0d0,0.0d0,gwlcrit_t,(mamp*mamte),ifnd)
          call rdfdor('hcrit',-1000.0d0,0.0d0,hcrit_t,(mamp*mamte),ifnd)
          call rdfdor('vcrit',0.0d0,20.0d0,vcrit_t,(mamp*mamte),ifnd)

! ---     convert
          do i = 1, ifnd
            wlsman(imper_4e2(i),impphase(i)) = wlsman_t(i) - Altcu
            gwlcrit(imper_4e2(i),impphase(i)) = gwlcrit_t(i)
            hcrit(imper_4e2(i),impphase(i)) = hcrit_t(i)
            vcrit(imper_4e2(i),impphase(i)) = vcrit_t(i)
          enddo

! ---     check values
          do i = 1, ifnd
            imper = imper_4e2(i)
            iph = impphase(i)

            nphase(imper) = nphase(imper) + 1
            if (wlsman(imper,iph).lt.zbotdr(1+nrpri)) then
              messag = 'Level of automatic weir below bottom'           &
     &                 //' deepest channel of secondary system.'
              call fatalerr ('Rddra',messag)
            endif
            if (swman(imper) .ne. 2) then
              messag = '#4e2 swman - imper mismatch'
              call fatalerr ('Rddra',messag)
            endif
            if (imper.ne.imperb) then
              imperb = imper
              imperi = imperi + 1
            endif
          enddo

          if (imperi .ne. nrman2) then
            messag = '#4e2 inconsistency between the nr of periods'     &
     &            //'with automatic weir (tabel 4e (IMPER4e2)) and'     &
     &            //' swman in tabel4b (IMPER_4b)'
            call fatalerr ('Rddra',messag)
          endif

! ---     consistency checks WLSMAN, GWLCRIT, HCRIT and VCRIT
          do 446 imper = 1,nmper
            if (swman(imper).eq.2) then
              if (abs(gwlcrit(imper,1)) .gt. 0.01d0) then
                messag = '#4e2 - gwlcrit(1) must be 0.'
                call fatalerr ('rddra',messag)
              endif
              if (abs(vcrit(imper,1)) .gt. 0.01d0) then
                messag = '#4e2 - vcrit(1) must be 0.'
                call fatalerr ('rddra',messag)
              endif
              if (abs(hcrit(imper,1)) .gt. 0.01d0) then
                messag = '#4e2 - hcrit(1) must be 0.'
                call fatalerr ('rddra',messag)
              endif
              if (hbweir(imper).gt. (wlsman(imper,1)-0.99999d0)) then
                messag = '#4e2 - HBWEIR within 1 cm of wlsman(1)'
                call fatalerr ('Rddra',messag)
              endif

              do 448 iph = 2,nphase(imper)
                if (wlsman(imper,iph).lt.wlsman(imper,iph-1)) then
                  messag = '#4e2 - wlsman inconsistent'
                  call fatalerr ('rddra',messag)
                endif
                if (gwlcrit(imper,iph).gt.gwlcrit(imper,iph-1)) then
                  messag = '#4e2 - gwlcrit inconsistent'
                  call fatalerr ('rddra',messag)
                endif
                if (hcrit(imper,iph).gt.hcrit(imper,iph-1)) then
                  messag = '#4e2 - hcrit inconsistent'
                  call fatalerr ('rddra',messag)
                endif
                if (vcrit(imper,iph).lt.vcrit(imper,iph-1)) then
                  messag = '#4e2 - vcrit inconsistent'
                  call fatalerr ('rddra',messag)
                endif
448           continue
            endif 
446       continue
        endif
      endif
! ----------------------------------------------------------------------
! --- initialize counter for number of target level adjustments
      numadj = 0

! --- sttab(i,1) contains depths

! --- sw-levels, to be used in piece-wise linear functions (first level)
! --- is 100 cm above soil surface to allow for situations with ponding)
      sttab(1,1) = 100.0d0
      sttab(2,1) =   0.0d0
      do 100 i = 3,22
! ---   layer between surface level and deepest
! ---   drain/channel bottom is divided into 20 compartments
        sttab(i,1) = zbotdr(1+nrpri)*(i-2)/20.0d0
100   continue

! --- sttab(i,2) contains storage expressed in cm

! --- calculation of surface water storage (in cm, i.e. volume
! --- per unit area), as a function of sw-level, only for open channels:
      do 104 i = 1,22
        sttab(i,2) = 0.0d0
        do 108 ilev = 1+nrpri,nrlevs
          if (swdtyp(ilev).eq.0.and.sttab(i,1).gt.zbotdr(ilev)) then
! ---       for levels above soil surface the volume-increment is
! ---       computed for a rectangle, and not a trapezium
            if (sttab(i,1) .le. 0.0d0) then
              wdepth = sttab(i,1)-zbotdr(ilev)
              wvolum = wdepth*(widthr(ilev)+wdepth/taludr(ilev))
            else
              wdepth = -zbotdr(ilev)
              wvolum = wdepth*(widthr(ilev)+wdepth/taludr(ilev))
              wbreadth = widthr(ilev) + 2*wdepth/taludr(ilev)
              wdepth = sttab(i,1)
              wvolum = wvolum + wbreadth*wdepth
            endif
            sttab(i,2) = sttab(i,2)+wvolum/l(ilev)
          endif
108     continue
104   continue

! --- initial storage swstini
      swstini = swstlev (wls1,sttab)
      swst = swstini

! --- initialize memorization of wls for most recent 4 timesteps 
      do i=1,4
        wlsbak(i) = 0.0d0
      enddo

! --- close input file with lateral boundary conditions
      CLOSE (DRA)         

      RETURN
      END

! ----------------------------------------------------------------------
      subroutine checkdate(ifnd,dates,tend,tstart,namedat,topic) 
! ----------------------------------------------------------------------
!     Date               : April 2006   
!     Purpose            : check range of input dates to range of simulation period
! ----------------------------------------------------------------------
      implicit none

! --- global
      integer   ifnd
      real(8)   dates(ifnd), tend, tstart
      character(len=200) messag
      character(len=5)   namedat
      character(len=*)   topic

! ----------------------------------------------------------------------
! --- local
      integer i
      logical fldaterr

! ----------------------------------------------------------------------

!   - at least one input date must be within simulation period or 
!     simulation period should be completely within range of input dates
      fldaterr = .true.
      i = 1
      do while (fldaterr .and. i.le.ifnd)
         if(dates(i).gt.tstart-1.d-6.and.dates(i).lt.tend+1.d-6)        &
     &      fldaterr=.false.
         i = i + 1
      enddo
      if(fldaterr) then
         if(dates(1).lt.tstart+1.d-6 .and. dates(ifnd).gt.tend-1.d-6)   &
     &      fldaterr=.false.
      endif
      if(fldaterr) then
        messag = 'Fatal '//namedat//                                    &
     &           ', no input date within simulation period'
        call fatalerr(topic,messag)
      endif

      return
      end
