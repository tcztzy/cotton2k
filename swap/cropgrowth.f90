! File VersionID:
!   $Id: cropgrowth.f90 328 2017-06-03 21:05:50Z kroes006 $
! ----------------------------------------------------------------------
      subroutine CropGrowth(task) 
! ----------------------------------------------------------------------
!     UpDate             : May 2014   
!     Date               : Aug 2004   
!     Purpose            : Call proper crop routines for initialization,
!                          calculation of rate/state variables and output
! ----------------------------------------------------------------------

      use variables
      implicit none

      integer task,node !,oldcrop
!      character(len=200) messag
!      include   'params.fi'

      select case (task)
      case (1)

! === initialization =========================================================

! --- set crop conditions 
      
! --- check growing season

      icrop = 1
      flCropCalendar = .false.
      do while (.not. flCropCalendar) 
        
        if (cropstart(icrop) .lt. 1.d0) exit
        
        if (t1900 - cropstart(icrop) .gt. -1.d-3                        &
     &                 .and. t1900 - cropend(icrop) .lt. 1.d-3) then
          flCropCalendar = .true.
        else
          icrop = icrop + 1
        endif
      enddo
      
!---  reset if new crop
      if (flCropCalendar) then
        if (t1900 - cropstart(icrop) .gt. -1.d-3 .and.                  &
     &                t1900 - cropstart(icrop) .lt. 1.d-3) then
          daycrop  = 0
          tsumgerm = 0.d0
          flCropEmergence = .false.
          flCropHarvest   = .false.
          flHarvestDay    = .false.
          flCropReadFile  = .true.
        
        endif
      endif

!---  initialize crop conditions
      if (flCropCalendar .and. .not. flCropHarvest) then
        if (flCropReadFile) then
! ---     fixed crop development -----------------------------------------------  
          if (croptype(icrop) .eq. 1) call CropFixed(1)
! ---     detailed crop growth -------------------------------------------------
          if (croptype(icrop) .eq. 2) call Wofost(1)
! ---     detailed grass growth  -----------------------------------------------
          if (croptype(icrop) .eq. 3) call Grass(1)
          flCropReadFile = .false.
        
          call CropOutput(1)
        endif
      endif

! --- set bare soil condition  -----------------------------------------------
      if (.not. flCropCalendar .or. flCropHarvest) then
        call nocrop (rd,lai,cf,ch,albedo,rsc)
      endif
      
      return

      case (2)

! === calculation of potential crop rate and state variables =================

      if (flCropHarvest) return

! --- detailed crop growth -------------------------------------------------
      if (croptype(icrop) .eq. 2) then

!       Init Wofost with startdate = emergence time 
        if (initcrp(icrop) .eq. 1) then
          flCropEmergence = .true.

!       Init Wofost with startdate = sowing time
!          implies: delay of growth and Simulate germination time 
        elseif (initcrp(icrop) .eq. 2) then
          if (.not.flCropEmergence .and. .not. flHarvest) then 
!             pressure head of rootzone determines timing of emergence
              node = 1
              hrz1 = 0.0d0
              do while ( abs(z(node)-0.5d0*dz(node)) .lt. 10.00001d0)
                hrz1 = hrz1 + h(node)
                node = node + 1
              enddo
              hrz1 = hrz1 / (node-1)
!             simulate germination time
              call wofostarablelandgerm(2,hrz1)
!             delay growth until tsumgerm is reached
              if (tsumgerm .lt. tsumemeopt) then
                ch  = 0.0D0     ! pvw
                cf  = 0.0D0     ! pvw
                rd = 0.0d0
                rdpot = 0.0D0  ! pvw
                lai = 0.0d0
                laipot = 0.0D0   ! pvw
              else
                flCropEmergence = .true.
!               initialise for next crop/germination period
                tsumgerm = 0.0d0
!               initial rooting depth and LAI equal to input values for emergence
                rd = rdi
                rdpot = rdi     ! pvw
                lai = laiem
                laipot = laiem  ! pvw  
              endif
          endif
        endif

!       start growth of arable crop with Wofost
        if(flCropEmergence) call Wofost(2)

      else if (croptype(icrop) .eq. 3) then
        call Grass(2)
      endif

      return

      case (3)

      if (flCropHarvest) return
      
! === calculation of actual crop rate and state variables =================
      
! --- fixed crop development -----------------------------------------------
      if (croptype(icrop).eq.1) call CropFixed(3)
! --- detailed crop growth, establish N-demand from soil -------------------
      if (croptype(icrop).eq.2)then
         if(flCropEmergence .or. flHarvestDay) then
!            if (laipot .gt. (laiem+nihil)) call Wofost(3)    ! pvw  20160225
            call Wofost(3)                                    ! jkro 20160426
         endif
      end if
! --- detailed grass growth  -----------------------------------------------
      if (croptype(icrop).eq.3) call Grass(3)

      return

      case (4)

      if (flCropHarvest) return
      
! --- detailed crop growth -------------------------------------------------
      if (croptype(icrop).eq.2)then
         if (flCropEmergence .or. flHarvestDay) then
           call Wofost(4)  
         endif
      end if

      case default
         call fatalerr ('CropGrowth', 'Illegal value for TASK')
      end select

      return

      end
!
! ----------------------------------------------------------------------
      subroutine cropfixed (task)
! ----------------------------------------------------------------------
!     date               : august 2004                           
!     purpose            : simple crop growth routine for swap 
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local variables
      integer i,task,lcc,swhydrlift
      
      real(8) gctb(2*magrs),rdtb(2*magrs),kytb(2*magrs),dummy,          &
     &        afgen,dtsum,dvr,watcon,dvssend
      parameter (dvssend=2.d0)
      include 'params.fi'
      save
! ----------------------------------------------------------------------

      select case (task)
      case (1)

! === initialization ===================================================
      call InitializeCrop
      cftb = 0.0d0
      chtb = 0.0d0
      gctb = 0.0d0
      kytb = 0.0d0
      rdtb = 0.0d0

! --- read crop data
      call readcropfixed (cropfil(icrop),pathcrop,idev,lcc,tsumea,      &
     & tsumam,tbase,kdif,kdir,gctb,swgc,cftb,cfeictb,fimin,siccaptb,    &
     & swcf,rdtb,rdctb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,     &
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
     & nofd,daycrop,dvs,tsum,t1900,tstart,swinco,cropstart(icrop),      &
     & saltmax,saltslope,salthead)
     
! --- development stage and tsum
      if (t1900 - tstart .gt. 1.d-3 .or. swinco .ne. 3 .or.             &
     &   (t1900 - cropstart(icrop) .gt. -1.d-3 .and.                    &
     &                t1900 - cropstart(icrop) .lt. 1.d-3)) then
        dvs = 0.0d0
        tsum = 0.0d0
      endif

! --- initial lai or sc
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif

! --- initial crop factor or crop height
      cf = afgen (cftb,(2*magrs),dvs)
      ch = afgen (chtb,(2*magrs),dvs)
      if (swcf.eq.3) then
        cfeic = afgen (cfeictb,(2*magrs),dvs)
      endif

! --- initial storage on canopy, for NHI option (if used).ge.span
      if (swinter.eq.3) then
        siccapact = afgen(siccaptb,(2*magrs),t)
        sicact    = 0.
      endif

! --- actual rooting depth [cm]
      rd = min (rdmax,afgen (rdtb,(2*magrs),dvs))

! --- initial dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),dvs)

! --- initial ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),dvs)

! --- initialize matric flux potential and hleaf
      if (swdrought .eq. 2) then                                        
        call MatricFlux(1,h(1),1,dummy)
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,numnod
         twilt(i) = watcon(i,wiltpoint,cofgen,swsophy,numtab,sptab,     &
     &                     ientrytab)
         hroot(i) = h(i)
        enddo
        hleaf = -2000.d0
      endif                                                             

      return          

      case (2)
      continue

! === calculate potential rate and state variables ======================
      case (3)

! === calculate actual rate and state variables ======================

! --- increase in temperature sum
      dtsum = max (0.0d0,tav-tbase)

! --- development rate
      if (idev.eq.1) then
        dvr = 2.0/lcc
      elseif (idev.eq.2) then
        if (dvs.lt.1.0d0) then
          dvr = dtsum/tsumea
        else
          dvr = dtsum/tsumam
        endif
      endif

! --- water stress
      if(abs(ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(min(tra/ptra,1.0d0),0.0d0)
      endif

! ----integrals of the crop --------------------------------------------

! --- phenological development stage
      dvs = min(dvs+dvr,dvssend)
      tsum = tsum + dtsum

! --- leaf area index or soil cover fraction    
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif

! --- crop factor or crop height
      cf        = afgen (cftb,(2*magrs),dvs)
      ch        = afgen (chtb,(2*magrs),dvs)
      if (swcf.eq.3) then
        cfeic     = afgen (cfeictb,(2*magrs),dvs)
      endif
      if (swinter.eq.3) then
        siccapact = afgen(siccaptb,(2*magrs),t)
      endif

! --- rooting depth [cm]
      rd = min (rdmax,afgen (rdtb,(2*magrs),dvs))

! --- dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),dvs)

! --- ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),dvs)


      case default
         call fatalerr ('CropFixed', 'Illegal value for TASK')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine cropoutput(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : open and write crop output files 
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- local variables ------------------
      integer task,getun   !,numcrop
      character(len=200) messag
      character(len=160) filnam,filtext

      select case (task)
      case (1)

! === open crop output file and write headers =====================

! --- open crop output file
      if (flCropOpenFile) then
! ---   open crop output file and write general header
        if (trim(outfil).eq.trim(cropfil(1))) then
          Messag = 'The name of the input crop-file (''//trim(cropfil'//&
     &   '(icrop))//'') cannot be equal to the name of'                 &
     &   //'the output crop-file '//trim(outfil)//' Adjust a filename !'
          call fatalerr ('crops',messag)
        endif
        filnam = trim(pathwork)//trim(outfil)//'.crp'
        crp = getun (20,90)
        call fopens(crp,filnam,'new','del')
        filtext = 'output data of simple or detailed crop growth model'
        call writehead (crp,1,filnam,filtext,project)

! ---   write header fixed crop growth
        if (croptype(icrop) .eq. 1)                                     &
     &    call OutCropFixed(1,date,t,daycrop,dvs,tsum,lai,cf,rd,crp,ch)           
! ---   write header detailed crop 
        if (croptype(icrop) .eq. 2)                                     &
     &    call OutWofost(1,date,daycrop,crp,t,dvs,tsum,laipot,lai,      &
     &            rdpot,rd,ch,cf,cwdmpot,cwdm,wsopot,wso,               &
     &            wlvpot,wlv,wstpot,wst,wrtpot,wrt,                     &
     &            dwlvCrop,dwlvSoil,dwst,dwrt,dwso,HarLosOrm_tot)
! ---   write header detailed grass growth
        if (croptype(icrop) .eq. 3)                                     &
     &    call OutGrass(1,date,daycrop,crp,t,dvs,tsum200,laipot,lai,    &
     &            rdpot,rd,ch,cf,tagppot,tagp,tagptpot,tagpt,           &
     &            wlvpot,wlv,wstpot,wst,wrtpot,wrt,cuptgraz,cuptgrazpot)

        flCropOpenFile = .false.

      else
! ---   header for second and subsequent crops

! ---   write header fixed crop growth
        if (croptype(icrop).eq.1 .and. swheader.eq.1)                   &
     &    call OutCropFixed(1,date,t,daycrop,dvs,tsum,lai,cf,rd,crp,ch)
! ---   write header detailed crop growth 
        if (croptype(icrop).eq.2 .and. swheader.eq.1)                   &
     &    call OutWofost(1,date,daycrop,crp,t,dvs,tsum,laipot,lai,      &
     &            rdpot,rd,ch,cf,cwdmpot,cwdm,wsopot,wso,               &
     &            wlvpot,wlv,wstpot,wst,wrtpot,wrt,                     &
     &            dwlvCrop,dwlvSoil,dwst,dwrt,dwso,HarLosOrm_tot)
! ---   write header detailed grass growth
        if (croptype(icrop).eq.3 .and. swheader.eq.1)                   &
     &    call OutGrass(1,date,daycrop,crp,t,dvs,tsum200,laipot,lai,    &
     &            rdpot,rd,ch,cf,tagppot,tagp,tagptpot,tagpt,           &
     &            wlvpot,wlv,wstpot,wst,wrtpot,wrt,cuptgraz,cuptgrazpot)

      endif

      return

      case (2)

! --- write actual data ----------------------------------------------------

! --- fixed crop file
      if (croptype(icrop) .eq. 1)                                       &
     &  call OutCropFixed(2,date,t,daycrop,dvs,tsum,lai,cf,rd,crp,ch)

! --- detailed crop growth 
      if (croptype(icrop) .eq. 2)                                       &
     &  call OutWofost(2,date,daycrop,crp,t,dvs,tsum,laipot,lai,        &
     &          rdpot,rd,ch,cf,cwdmpot,cwdm,wsopot,wso,                 &
     &          wlvpot,wlv,wstpot,wst,wrtpot,wrt,                       &
     &          dwlvCrop,dwlvSoil,dwst,dwrt,dwso,HarLosOrm_tot)
 
! --- detailed grass growth
      if (croptype(icrop) .eq. 3)                                       &
     &  call OutGrass(2,date,daycrop,crp,t,dvs,tsum200,laipot,lai,      &
     &          rdpot,rd,ch,cf,tagppot,tagp,tagptpot,tagpt,             &
     &          wlvpot,wlv,wstpot,wst,wrtpot,wrt,cuptgraz,cuptgrazpot)
      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (crp)

      case default
         call fatalerr ('CropOutput', 'Illegal value for TASK')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine nocrop (rd,lai,cf,ch,albedo,rsc)
! ----------------------------------------------------------------------
      real(8) rd,lai,cf,ch,albedo,rsc
! ----------------------------------------------------------------------
      rd = 0.0d0
      lai = 0.0d0
      cf = 0.0d0
      ch = 12.d0
      albedo = 0.23d0
      rsc = 70.d0

      return
      end

! ----------------------------------------------------------------------
      subroutine wofost(task)
! ----------------------------------------------------------------------
!     update             : march 2015
!     date               : october 2004
!     purpose            : detailed crop growth routine
! ----------------------------------------------------------------------
      use variables
      use wofost_soil_interface

      implicit none
 
      integer   i1,task,swhydrlift,i,uco2,ifnd

      real(8)   afgen,amax,asrc,ccheck,cosld,cvf   !,admi
      real(8)   laicr,lasum,mres,rdm
      real(8)   dalv,dayl,delt,dmi
      real(8)   drrt,drst,dslv,dteff,dtga,dtsum,dvr   !,dslv1,dslv2,dslvt
      real(8)   dvred,fl,fo,fr,fs,drlv !,fcheck
      real(8)   fysdel,gass,gla,glaiex,grlv,grrt,grst,grso  !,glasol
      real(8)   gwrt,gwso,gwst,pgass,rmres,rr   !,rest
      real(8)   sinld,slat,teff,twlv,twst
      real(8)   lasumpot,watcon,dtgapot
      real(8)   pgasspot,gasspot,rmrespot,mrespot,asrcpot,dmipot,rrpot
      real(8)   admipot,grrtpot,drrtpot,gwrtpot,grlvpot
      real(8)   dslvpot,restpot,dalvpot,drlvpot,gwsopot
      real(8)   glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real(8)   dslvtpot,twlvpot,twstpot
      real(8)   dummy
      logical   flpotrelmf    ! Flag indicating calculation of limited attainable yield instead of potential yield (due to management factor)
      real(8)   relmf         ! relative Management factor that reduces crop growth (-)
!      character(len=15) tmp
      character(len=200) messag

!     n-p-k use
      real(8) LRNR, LSNR, NLAI   !, LRPR, LSPR, LRKR, LSKR
      real(8) NLUE, NMAXSO, NPART, NFIXF
      real(8) NSLA, RNFLV, RNFRT, RNFST, TCNT      !, TCPT
      real(8) DVSNLT, DVSNT, RDRNS, FNTRT, FRNX    !, TCKT
      integer getun2
      real(8) INSW
      real(8) NMAXLV,NMAXST,NMAXRT,NNI,FSTR
      real(8) NMXLV(30)
      integer ILNMXL
      real(8) Fstress
!      real(8) FERNTAB(90), NRFTAB(30)
      integer nut    !ILFERN, ILNRFT
!      real(8) DAY
      real(8) NDEMTO   !,NMINT
      real(8) ANLV,ANST,ANRT,ANSO,ANLVI,ANSTI,ANRTI,ANSOI
      real(8) RNLDRT,RNLDST,RNLDLV,NDEML,NDEMS,NDEMR,NDEMSO
      real(8) NSUPSO,ATN,RNSO,NLIMIT,NUPTR
      real(8) NFIXTR,ATNLV,ATNST,ATNRT,RNTLV,RNTST,RNTRT
      real(8) RNULV,RNUST,RNURT,NUPTT,NFIXTT
      real(8) RNLV,RNST,RNRT  !,FERTN,NRF,FERTNS
      real(8) NLOSSL,NLOSSR,NLOSSS  
      real(8) NBALAN !,npki
!!!      real(8) NLOSST,NROOT,NLIVT
      real(8) NdemandBioFix
      real(8) ombalan,wlvt0,wstt0,wsot0,wrtt0,storagediff,drso
      real(8) FraHarLosOrm_lv,FraHarLosOrm_st,FraHarLosOrm_so
      real(8) FraDeceasedLvToSoil
      real(8) HarLosOrm_rt, HarLosOrm_lv, HarLosOrm_st, HarLosOrm_so
      real(8) HarLosNit_rt, HarLosNit_lv, HarLosNit_st, HarLosNit_so
      real(8) HarLosOrm_dwrt, HarLosOrm_dwlv
      real(8) HarLosOrm_dwst, HarLosOrm_dwso
      real(8) HarLosNit_dwrt, HarLosNit_dwlv 
      real(8) HarLosNit_dwst, HarLosNit_dwso
      real(8) idwlvCrop,idwlvSoil, NLOSSLDeceasedLvToSoil
      character(len=200) filnam
      logical rdinqr
      real(8) CO2, effc
      integer ifindi, indexyr

! --- only for soybean
      integer swsoybean
      logical flrfphotoveg
      real(8) mg,dvsi,dvrmax1,dvrmax2,tmaxdvr,tmindvr,toptdvr
      real(8) rfmgphotop, rfmgtemp
      logical flphenodayl     ! Flag to allow input of POPT and PCRT or using empirical relation from Setiyono et al
      real(8) popt            ! optimal daylength for phenological development (hr)
      real(8) pcrt            ! critical daylength for phenological development (hr)

      parameter (delt=1.0d0)
      include 'params.fi'

! --- only for vernalisation
      logical flvernalised
      real(8) r,vern,vernfac,vernrate,interpol

!      logical RunRoutine 
      
      save
! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initialization ====================================================
      call InitializeCrop
      
! --- read general crop data
      call readwofost (cropfil(icrop),pathcrop,swcf,cftb,idsl,dlo,dlc,  &
     &  initcrp(icrop),tsumemeopt,tbasem,teffmx,                        &
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
     &  t1900,tstart,cropstart(icrop),verndvs,vernsat,vernbase,verntb,  &
     &  saltmax,saltslope,salthead)


! --- if crop based on calendar is still active, but already harvested
      if (flCropHarvest) return

! --- n-p-k 
      if( flCropNut) then
!        initialise and start reading
         filnam = trim(pathcrop)//trim(cropfil(icrop))//'.crp'
         nut = getun2 (10,90,2)
         call rdinit(nut,logf,filnam)

         CALL rdsdou ('LRNR', LRNR)
         CALL rdsdou ('LSNR', LSNR)
         CALL rdsdou ('NLAI', NLAI)
         CALL rdsdou ('NLUE', NLUE)
         CALL rdsdou ('NMAXSO', NMAXSO)
         CALL rdsdou ('NPART', NPART)
         CALL rdsdou ('NFIXF', NFIXF)
         CALL rdsdou ('NSLA', NSLA)
         CALL rdsdou ('RNFLV', RNFLV)
         CALL rdsdou ('RNFRT', RNFRT)
         CALL rdsdou ('RNFST', RNFST)
         CALL rdsdou ('TCNT', TCNT) 
         CALL rdsdou ('DVSNLT', DVSNLT)
         CALL rdsdou ('DVSNT', DVSNT)
         CALL rdsdou ('RDRNS', RDRNS)
         CALL rdsdou ('FNTRT', FNTRT)
         CALL rdsdou ('FRNX', FRNX)
         CALL rdadou ('NMXLV', NMXLV, 30, ILNMXL)
!        Read harvest losses (fractions of leaves, stems, storage organs)
         call rdsdor ('FraHarLosOrm_lv',0.d0,1.0d0,FraHarLosOrm_lv) 
         call rdsdor ('FraHarLosOrm_st',0.d0,1.0d0,FraHarLosOrm_st) 
         call rdsdor ('FraHarLosOrm_so',0.d0,1.0d0,FraHarLosOrm_so) 

! -      close input file
         close(nut)

!        open output files and write header 
         if (icrop.eq.1) then
            call outbalcropOM1(1,pathwork,outfil,project,date,daycrop,  &
     &         t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)
            call outbalcropOM2(1,pathwork,outfil,project,date,daycrop,  &
     &         t,dvs,tsum,storagediff,wlv,wst,wso,wrt,delt,         &
     &         grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)
            call outbalcropN(1,pathwork,outfil,project,date,daycrop,    &
     &         t,dvs,tsum,nuptt,nfixtt,anlvi,ansti,anrti,ansoi,anlv,&
     &         anst,anrt,anso,nlossl,nlossr,nlosss,nbalan,nni)
         endif
      endif

! --- CO2 impact
      filnam = trim(pathcrop)//trim(cropfil(icrop))//'.crp'
      uco2 = getun2 (10,90,2)
      call rdinit(uco2,logf,filnam)
      flco2 = .false.
      if(rdinqr('flco2')) then
        call rdslog ('flco2',flco2)
      endif
      fco2amax = 1.0d0           ! factor to correct AMAX for CO2
      fco2eff = 1.0d0            ! factor to correct EFF for CO2
      fco2tra = 1.0d0            ! factor to correct TRA for CO2
      if(flco2) then
!        Read CO2-crop data
         CALL rdadou ('CO2AMAXTB', CO2AMAXTB, 30, ifnd) 
         CALL rdfdou ('CO2EFFTB', CO2EFFTB, 30, ifnd) 
         CALL rdfdou ('CO2TRATB', CO2TRATB, 30, ifnd) 
      endif
 
!---  Fraction of deceased leaves entering the soil system
      FraDeceasedLvToSoil = 0.0d0
      if(rdinqr('FraDeceasedLvToSoil')) then
         call rdsdor ('FraDeceasedLvToSoil',0.d0,1.0d0,FraDeceasedLvToSoil) 
      endif
      
      close(uco2)

!---  Read CO2-air data from a separate file
      if(flco2) then
         filnam = trim(pathcrop)//'Atmospheric.co2'
         uco2 = getun2 (30,90,2)
         call rdinit(uco2,logf,filnam)
         CALL rdainr ('CO2year', 1000, 3000, CO2year, mayrs, ifnd)
         CALL rdfdor ('CO2ppm',  10.0d0, 1000.0d0, CO2ppm, mayrs, ifnd)
         close(uco2)
      endif

! --- initial values of crop parameters
      if(swcf.ne.3) then
        swinter = 1
      else
        swgc    = 1
        swinter = 3
      endif

! --- maximum rooting depth
      rdm = min(rdmax,rdc)

! --- skip next initialization if crop parameters are read from *.END file
      if (t1900 - tstart .gt. 1.d-3 .or. swinco .ne. 3 .or.             &
     &   (t1900 - cropstart(icrop) .gt. -1.d-3 .and.                    &
     &                t1900 - cropstart(icrop) .lt. 1.d-3)) then
      
! ---   actual rooting depth
        rd = min(rdi,rdm)
!       if(.not.flCropEmergence)  rd = 0.0d0   !!!!!! >>>>> please pay attention
        rdpot = min(rdi,rdm)

        dvs = 0.0d0
        flAnthesis = .false.
        tsum = 0.0d0
        fr = afgen (frtb,30,dvs)
        fl = afgen (fltb,30,dvs)
        fs = afgen (fstb,30,dvs)
        fo = afgen (fotb,30,dvs)
        sla(1) = afgen (slatb,30,dvs)
        lvage(1) = 0.0d0
        ilvold = 1
        slapot(1) = afgen (slatb,30,dvs)
        lvagepot(1) = 0.0d0
        ilvoldpot = 1

! ---   initial state variables of the crop
        wrt = fr*tdwi
        wrtpot = wrt
        tadw = (1.0d0-fr)*tdwi
        tadwpot = tadw
        wst = fs*tadw
        wstpot = wst
        wso = fo*tadw
        wsopot = wso
        wlv = fl*tadw
        wlvpot = wlv
!       KRO-BOO-20160403: intro because comparison with Wofost
        laiem = wlv*sla(1)  
        lv(1) = wlv
        lvpot(1) = wlv
        lasum = laiem     
        lasumpot = laiem     
        laiexp = laiem     
        glaiex = 0.0d0
        laimax = laiem
        lai = lasum+ssa*wst+spa*wso 
        laipot = lai 
        dwrt = 0.0d0
        dwrtpot = 0.0d0
        dwlv = 0.0d0
        dwlvCrop = 0.0d0
        dwlvSoil = 0.0d0
        dwlvpot = 0.0d0
        dwso = 0.0d0
        dwst = 0.0d0
        dwstpot = 0.0d0
        if(flCropNut) then
          WLVt0 = wlv
          WSTt0 = wst
          WSOt0 = wso
          WRTt0 = wrt
        endif
! ---   n-p-k 
        if( flCropNut) then
!******************************************************************
!         initial maximum nutrient concentrations in plant organs 
!         per kg biomass [kg N kg-1 dry biomass] at sowing added IS
!******************************************************************        
          call nutrsow(anlv,anst,anrt,anso)
!******************************************************************
!         initial maximum nutrient concentrations in plant organs 
!         per kg biomass [kg N kg-1 dry biomass] at emergence added IS
!******************************************************************
          call nutremrg(nmxlv,lsnr,lrnr,wlv,wst,wrt,                    &
     &      anlv,anst,anrt,anso,anlvi,ansti,anrti,ansoi,dvs)
        endif

! --- initial summation variables of the crop
        gasst = 0.0d0
        gasstpot = 0.0d0
        mrest = 0.0d0 
        mrestpot = 0.0d0 
        cwdm = 0.0d0
        cwdmpot = 0.0d0
! --- only for vernalisation
        vern = 0.0d0             ! vernalisation state (d)
        flvernalised = .FALSE.   ! crop not vernalised (-)

! --- end skip above initialization if crop parameters are read from *.END file
      endif
      
! ---   set crop height and cropfactor
        if (swcf.ne.3) then
          cf = afgen (cftb,(2*magrs),dvs)
          ch = afgen (chtb,(2*magrs),dvs)
        else
          cf        = afgen (cftb,(2*magrs),lai)
          cfeic     = afgen (cfeictb,(2*magrs),lai)
          ch        = afgen(chtb,(2*magrs),lai)
          siccapact = siccaplai*lai
        endif

! --- initialize matric flux potential and hleaf
      if (swdrought .eq. 2) then                                        
        call MatricFlux(1,h(1),1,dummy)
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,numnod
         twilt(i) = watcon(i,wiltpoint,cofgen,swsophy,numtab,sptab,     &
     &                     ientrytab)
         hroot(i) = h(i)
        enddo
        hleaf = -2000.d0
      endif                                                             

! --- initialize germination routine
      if (initcrp(icrop) .eq. 2) then

        if (t1900 - tstart .gt. 1.d-3 .and. swinco .ne. 3) then
!         Start the crop based on specified sowing date             
          call wofostarablelandgerm(1,0.0d0)
        endif
      endif

! -      n-p-k 
      if( flCropNut) then
        call nutrinit    (nlossl,nlossr,nlosss,                         &
     &                    nuptt,rnlv,rnst,rnrt,rnso,                    &
     &                    rnldlv,rnldst,rnldrt,                         &
     &                    nfixtt,nni,NLOSSLDeceasedLvToSoil)
      endif


      return

      case (2)

! === calculate potential rate and state variables =====================

! --- rates of change of the crop variables ----------------------------

! --- phenological development rate 
      call astro(daynr,lat,rad,dayl,daylp,sinld,cosld,difpp,            &
     &           atmtr,dsinbe)

! --- increase in temperature sum
      dtsum = afgen (dtsmtb,30,tav)

! --- phenological development rate for potential AND actual crops
      if (swsoybean.eq.0) then
! ---   standard crops
        if (dvs.lt.1.0d0) then     
! ---     vegetative phase
          dvred = 1.0d0
          vernfac = 1.0d0
          vernrate = 0.0d0
          if (idsl.ge.1) then
             dvred = max(0.0d0,min(1.0d0,(daylp-dlc)/(dlo-dlc)))
          endif
          if (idsl.eq.2) then
!            vernalisation rate,based on routines from pyWofost (Allard de Wit, 2015)
             if(.not.flvernalised) then
                if(dvs.lt.verndvs) then
                   vernrate = afgen (verntb,30,tav)
                   r = (vern - vernbase) / (vernsat - vernbase)
                   vernfac = interpol(0.0d0,1.0d0,r)
                else
                   flvernalised = .true.
                endif
             endif
          endif
          dvr = vernfac * dvred*dtsum/tsumea
        else
! ---     generative phase
          dvr = dtsum/tsumam
        endif    

      else if (swsoybean.eq.1) then
! ---   soybean
        call mgtemprf(tav,toptdvr,tmindvr,tmaxdvr,rfmgtemp)
        call mgphotoprf(mg,daynr,lat,popt,pcrt,flphenodayl,rfmgphotop)
        if (dvs.lt.1.0d0) then 
! ---     vegetative phase
          if(flrfphotoveg) then
             dvr = dvrmax1 * rfmgphotop * rfmgtemp
          else
             dvr = dvrmax1 * rfmgtemp
          endif
        else
! ---     generative phase
          dvr = dvrmax2 * rfmgphotop * rfmgtemp
        endif
      endif

!     adjust development stage for realistic TSUM
      if (dvs.ge.1.d0 .and. (.not. flAnthesis)) then
        flAnthesis = .true.
        dvs = 1.0d0
      end if
      
! === daily dry matter production 

!***************************************************************     
!     Assimilation correction for CO2 changes in atmosphere 
!     borrowed from Lintul4 added IS
!***************************************************************   
      if(flco2) then
        indexyr = ifindi (CO2year, mayrs, 1, mayrs, iyear)
        if (indexyr.lt.1 .or. indexyr.gt.mayrs) then
          Messag ='Input if CO2year or CO2ppm inconsistent, correct'
          call fatalerr ('wofost',messag)
        endif
        CO2 = CO2ppm(indexyr)
        fco2amax = afgen(CO2AMAXTB,30,CO2)
        fco2eff = afgen(CO2EFFTB,30,CO2)    
        fco2tra = afgen(CO2TRATB,30,CO2)    
      endif
      effc = fco2eff * eff
      amax = fco2amax * afgen (amaxtb,30,dvs) * afgen (tmpftb,30,tavd)
      call totass (dayl,amax,effc,laipot,kdif,rad,difpp,dsinbe,         &
     &             sinld,cosld,dtgapot)
! --- correction for low minimum temperature
      dtgapot = dtgapot * afgen (tmnftb,30,tmnr)
! --- potential assimilation in kg ch2o per ha
      pgasspot = dtgapot * 30.0d0/44.0d0

! --- reduction due to limited attainable maximum (expert option)
      reltr = 1.0d0        
      if (flpotrelmf) reltr = relmf
! --- optional reduction of pgass to gass
      gasspot = pgasspot * reltr

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration
      rmrespot = (rmr*wrtpot+rml*wlvpot+rms*wstpot+rmo*wsopot)*         &
     &            afgen(rfsetb,30,dvs)
      teff = q10**((tav-25.0d0)/10.0d0)
      mrespot = dmin1(gasspot,rmrespot*teff)
      asrcpot = gasspot - mrespot

! --- partitioning factors
      fr = afgen(frtb,30,dvs)
      fl = afgen(fltb,30,dvs)
      fs = afgen(fstb,30,dvs)
      fo = afgen(fotb,30,dvs)
! --- check on partitioning
      call chckprt(dvs,fr,fl,fs,fo)    

! --- dry matter increase
      cvf = 1.0d0/((fl/cvl+fs/cvs+fo/cvo)*(1.0d0-fr)+fr/cvr)
      dmipot = cvf*asrcpot
! --- check on carbon balance
      ccheck = (gasspot-mrespot-(fr+(fl+fs+fo)*(1.0d0-fr))*dmipot/cvf)  &
     &         /max(0.0001d0,gasspot)      
      if (abs(ccheck).gt.0.0001d0) then
        Messag ='The carbon balance is not correct'
        call fatalerr ('wofost',messag)
      endif

! == = growth rate by plant organ

! --- root extension
      rrpot = min (rdm-rdpot,rri)
      if (fr.le.0.0d0.or.pgasspot.lt.1.0d0) rrpot = 0.0d0

! --- growth rate roots and aerial parts
      admipot = (1.0d0-fr)*dmipot
      grrtpot = fr*dmipot
      drrtpot = wrtpot*afgen (rdrrtb,30,dvs)
      gwrtpot = grrtpot - drrtpot

! --- weight of new leaves
      grlvpot = fl*admipot

! --- death of leaves due to water stress or high lai
      laicr = 3.2d0/kdif
      dslvpot = wlvpot*max(0.0d0,min(0.03d0,0.03d0*                     &
     &           (laipot-laicr)/laicr))

! --- death of leaves due to exceeding life span:

! --- first: leaf death due to water stress or high lai is imposed 
! ---        on array until no more leaves have to die or all leaves
! ---        are gone

      restpot = dslvpot*delt
      i1 = ilvoldpot

      do while (restpot.gt.lvpot(max(i1,1)).and.i1.ge.1)
        restpot = restpot - lvpot(i1) 
        i1 = i1-1
      enddo

! --- then: check if some of the remaining leaves are older than span,
! ---       sum their weights

      dalvpot = 0.0d0
      if (lvagepot(max(i1,1)).gt.span .and. restpot.gt.0.0d0            &
     &                   .and.i1.ge.1) then
        dalvpot = lvpot(i1) - restpot
        restpot = 0.0d0
        i1 = i1-1
      endif

      do while (i1.ge.1.and.lvagepot(max(i1,1)).gt.span)
        dalvpot = dalvpot+lvpot(i1)
        i1 = i1-1
      enddo

      dalvpot = dalvpot/delt

! --- finally: calculate total death rate leaves
      drlvpot = dslvpot + dalvpot

! --- physiologic ageing of leaves per time step
      fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! --- specific leaf area valid for current timestep
      slatpot = afgen (slatb,30,dvs)

! --- calculation of specific leaf area in case of exponential growth:
! --- leaf area not to exceed exponential growth curve
      if (laiexp.lt.6.0d0) then
        dteff = max (0.0d0,tav-tbase)
! ---   increase in leaf area during exponential growth
        glaiex = laiexp*rgrlai*dteff
! ---   source-limited increase in leaf area
        glasolpot = grlvpot*slatpot
! ---   actual increase is determined by lowest value
        glapot = min (glaiex,glasolpot)
! ---   slat will be modified in case gla equals glaiex
        if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
      endif  

! --- growth rate stems
      grstpot = fs*admipot
! --- death rate stems
      drstpot = afgen (rdrstb,30,dvs)*wstpot
! --- net growth rate stems
      gwstpot = grstpot - drstpot

! --- growth rate storage organs
      gwsopot = fo*admipot

! ----integrals of the crop --------------------------------------------

! --- leaf death (due to water stress or high lai) is imposed on array 
! --- untill no more leaves have to die or all leaves are gone

      dslvtpot = dslvpot*delt
      i1 = ilvoldpot
      do while (dslvtpot.gt.0.and.i1.ge.1)
        if (dslvtpot.ge.lvpot(i1)) then
          dslvtpot = dslvtpot-lvpot(i1)
          lvpot(i1) = 0.0d0
          i1 = i1-1
        else
          lvpot(i1) = lvpot(i1)-dslvtpot
          dslvtpot = 0.0d0
        endif
      enddo

! --- leaves older than span die
      do while (lvagepot(max(i1,1)).ge.span.and.i1.ge.1)
        lvpot(i1) = 0.0d0
        i1 = i1-1
      enddo

! --- oldest class with leaves
      ilvoldpot = i1

! --- shifting of contents, updating of physiological age
      do i1 = ilvoldpot,1,-1
        lvpot(i1+1) = lvpot(i1)
        slapot(i1+1) = slapot(i1)
        lvagepot(i1+1) = lvagepot(i1)+fysdel*delt
      enddo
      ilvoldpot = ilvoldpot + 1

! --- new leaves in class 1
      lvpot(1) = grlvpot*delt
      slapot(1) = slatpot
      lvagepot(1) = 0.0d0 

! --- calculation of new leaf area and weight
      lasumpot = 0.0d0
      wlvpot = 0.0d0
      do i1 = 1,ilvoldpot
        lasumpot = lasumpot + lvpot(i1)*slapot(i1)
        wlvpot = wlvpot + lvpot(i1)
      enddo

! --- leaf area index in case of exponential growth
      laiexp = laiexp+glaiex*delt

! --- dry weight of living plant organs
      wrtpot = wrtpot + gwrtpot*delt
      wstpot = wstpot + gwstpot*delt
      wsopot = wsopot + gwsopot*delt

! --- total above ground biomass
      tadwpot = wlvpot + wstpot + wsopot
      tadwpot = tadwpot ! for Forcheck

! --- dry weight of dead plant organs (roots,leaves & stems)
      dwrtpot = dwrtpot + drrtpot*delt
      dwlvpot = dwlvpot + drlvpot*delt
      dwstpot = dwstpot + drstpot*delt

! --- dry weight of dead and living plant organs
      twlvpot = wlvpot + dwlvpot
      twstpot = wstpot + dwstpot
      cwdmpot = twlvpot + twstpot + wsopot

! --- total gross assimilation and maintenance respiration
      gasstpot = gasspot + gasstpot
      mrestpot = mrespot + mrestpot

! --- leaf area index
      laipot = lasumpot + ssa*wstpot + spa*wsopot
!     prevent immediate lai reduction at emergence
!     KRO-BOO-20160403: suppressed because deviates from Wofost
!      laipot = max(laipot, laiem)

      ! --- rooting depth
      rdpot = rdpot + rrpot

! --- vernalisation state (d)
      if(idsl.eq.2) then
          vern = vern + vernrate
          if(.not.flvernalised .and. vern.ge.vernsat) then
              flvernalised = .true.
          else
              if(flvernalised .and. vern.lt.vernsat) then
                 write(messag,'(a,i6,a)') ' critical DVS,',             &
     &           ' for vernalised reached, day = ', daycum,             &
     &           ' but vernalization requirements not yet filfilled ',  &
     &           ' forcing vernalization now'
                 call warn ('wofost',messag,logf,swscre)
              endif
          endif
      endif

      return

      case (3)

! === calculate actual rate and state variables =====================
! === with optional calculation of (flCropNut) water AND nutrient stress

!***************************************************************     
!     Assimilation correction for CO2 changes in atmosphere 
!     borrowed from Lintul4 added IS
!***************************************************************   
      if(flco2) then
          indexyr = ifindi (CO2year, mayrs, 1, mayrs, iyear)
          if (indexyr.lt.1 .or. indexyr.gt.mayrs) then
             Messag ='Input if CO2year or CO2ppm inconsistent, correct'
            call fatalerr ('wofost',messag)
          endif
          CO2 = CO2ppm(indexyr)
          fco2amax = afgen(CO2AMAXTB,30,CO2)
          fco2eff = afgen(CO2EFFTB,30,CO2)    
          fco2tra = afgen(CO2TRATB,30,CO2)    
      endif
      effc = fco2eff * eff
      amax = fco2amax * afgen (amaxtb,30,dvs) * afgen (tmpftb,30,tavd)
      ptra = fco2tra*ptra

! --- rates of change of the crop variables ----------------------------
 
! --- gross assimilation
      call totass(dayl,amax,effc,lai,kdif,rad,                          &
     &            difpp,dsinbe,sinld,cosld,dtga)

! --- correction for low minimum temperature
      dtga = dtga * afgen (tmnftb,30,tmnr)
! --- potential assimilation in kg ch2o per ha
      pgass = dtga * 30.0d0/44.0d0

! --- water stress reduction of pgass to gass
      if(abs(ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(0.0d0,min(1.0d0,tra/ptra))
      endif

      if( flCropNut) then
! --      water and nitrogen stress reduction of pgass to gass
          call NUTRIE (reltr,NLUE,WLV,WST,DVS,                          &
     &                   ANLV,ANST,NMXLV,                               &
     &                   NMAXLV,NMAXST,NMAXRT,LRNR,LSNR,                &
     &                   NNI,RNFLV,RNFST,FRNX,FSTR)
          gass = pgass * FSTR
      else
        gass = pgass * reltr
      endif

! --- management factor 
!     (other forms of stress, not accounted for)
      gass = gass * relmf

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration
      rmres = (rmr*wrt+rml*wlv+rms*wst+rmo*wso)*afgen(rfsetb,30,dvs)
      mres = dmin1(gass,rmres*teff)
      asrc = gass-mres

! --- partitioning factors
      fr = afgen(frtb,30,dvs)
      fl = afgen(fltb,30,dvs)
      fs = afgen(fstb,30,dvs)
      fo = afgen(fotb,30,dvs)
! --- check on partitioning
      call chckprt(dvs,fr,fl,fs,fo)    

      if( flCropNut) then
!********************************************************************         
!         partitioning correction as influenced by water and N stress
!         Note: the partioning depends only on the Nitrogen stress,
!         not on the P and K stress. Personal communication Joost Wolf        
!         added IS
!******************************************************************** 
          CALL SUBPAR (reltr,NPART,NNI,FR,FL,FS,FO)
      endif
      
! --- conversion factor 
      CVF = 1.0d0/((FL/CVL+FS/CVS+FO/CVO)*(1.0d0-FR)+FR/CVR)
! --- dry matter increase
      dmi = cvf*asrc
! --- check on carbon balance
      call CHCKCBL(DVS,CVF,DMI,FR,FL,FS,FO,GASS,MRES)

! --- growth rate by plant organ

! --- root extension
      rr = min (rdm-rd,rri)
      if (fr.le.0.0d0.or.pgass.lt.1.0d0) rr = 0.0d0

! --- growth rate roots and aerial parts
      call relgrwt(dmi,fr,fl,fs,fo,grrt,grlv,grst,grso)

! --- death of leaves due to water stress or high lai or nitrogen stress
      call deaths(flcropnut,wlv,kdif,lai,NNI,perdl,rdrns,reltr,dslv)

! --- death of leaves due to exceeding life span:
      call deatha(dslv,delt,ilvold,lv,lvage,span,i1,dalv)

! --- death rate leaves as result of death due to water stress or high lai and 
!                                    death due to exceeding life span
      drlv = dslv+dalv

! --- death rate stems
      drst = wst * afgen (rdrstb,30,dvs)

! --- death rate roots
      drrt = wrt * afgen (rdrrtb,30,dvs)

! --- net growth rate stems, roots, storage organs
      gwst = grst - drst
      gwrt = grrt - drrt
      drso = 0.0d0    ! death rate of storage organs is assumed to be 0
      gwso = grso - drso


! --- specific leaf area valid for current timestep
      if(flCropNut) then
!       nutrient and water stress
        slat = afgen (slatb,30,dvs)*EXP(-NSLA * (1.0d0-NNI))
      else
        slat = afgen (slatb,30,dvs)
      endif
!
!     Do not allow slat higher than slatpot; slatpot can be limited by exponential growth
      slat = min(slat,slatpot)   ! pvw

! --- calculation of specific leaf area in case of exponential growth:
! --- leaf area not to exceed exponential growth curve
! --- FSTR is actual stress: water and nutrient 
      Fstress = reltr
      if (flcropnut) then
         Fstress = FSTR
         if ((DVS .LT. 0.2d0).AND.(LAI .LT. 0.75d0)) then
           Fstress = reltr * EXP(-NLAI* (1.0d0 - NNI))
         endif
      endif
      call GLAI(Fstress,LAIEXP,tav,TBASE,RGRLAI,GRLV,SLAT,GLA)


! ---- UPDATE STATES: integrals of the crop --------------------------------------------

! --- phenological development stage
      dvs = min(dvs+dvr*delt,dvsend)
      tsum = tsum + dtsum*delt

! --- leaf death (due to water stress or high lai) is imposed on array 
! --- untill no more leaves have to die or all leaves are gone
      call lvdth(delt,dslv,span,ilvold,lvage,lv,i1)

! --- oldest class with leaves
      ilvold = i1

! --- shifting of contents, updating of physiological age
      call lvshft(delt,fysdel,ilvold,grlv,slat,lv,lvage,sla)

      ilvold = ilvold+1

! --- calculation of new leaf area and weight
      call lvwgli(ilvold,lv,sla,lasum,wlv)

! --- dry weight of living plant organs
      wrt = wrt+gwrt*delt
      wst = wst+gwst*delt
      wso = wso+gwso*delt

! --- total above ground biomass
      tadw = wlv+wst+wso

! --- dry weight of dead plant organs (roots,leaves & stems)

      dwrt = dwrt + drrt*delt
      dwlv = dwlv + drlv*delt
      dwst = dwst + drst*delt
      dwso = dwso + drso*delt   ! dummy, because drso is assumed to be 0
!     split dwlv
      idwlvCrop = (1.0d0-FraDeceasedLvToSoil) * drlv*delt
      idwlvSoil = FraDeceasedLvToSoil * drlv*delt
      dwlvCrop = dwlvCrop + idwlvCrop
      dwlvSoil = dwlvSoil + idwlvSoil
      
! --- dry weight of dead and living plant organs
!     twrt = wrt+dwrt
      twlv = wlv+dwlv
      twst = wst+dwst
      cwdm = twlv+twst+wso

! --- total gross assimilation and maintenance respiration
      gasst = gass + gasst
      mrest = mres + mrest

! --- leaf area index
      lai = lasum+ssa*wst+spa*wso
! --- determine maximum lai
      laimax = max (lai,laimax)
! --- determine minimum lai to prevent dying straight after 
!       emergence when growth is slowed down due to low temperature
!     KRO-BOO-20160403: suppressed because deviates from Wofost
!      lai = max(lai, laiem)

! --- update state variables for nutrient stress

!     Calling the subroutine for N losses of leaves, roots and stem storage
!     organs (kg N ha-1 d-1)
      if(flCropNut)then 

!        Calling the subroutines for N demand of leaves, roots and stem storage
!        organs (kg N ha-1 d-1)
         CALL NDEMND(WLV,WST,WRT,WSO,NMAXLV,NMAXST,                     &
     &                               NMAXRT,NMAXSO,ANLV,ANST,ANRT,ANSO, &
     &                               TCNT,NDEML,NDEMS,NDEMR,NDEMSO)

!        Total N demand (kg N ha-1)

         NDEMTO = MAX (0.0d0,(NDEML + NDEMS + NDEMR))

!        Nutrient uptake limiting factor (-) at low moisture conditions in the
!        rooted soil layer before anthesis. After anthesis/DVSNLT there is no
!        nutrient uptake from the soil
         NLIMIT = INSW(DVS-DVSNLT,INSW(reltr-0.01d0,0.0d0,1.0d0),0.d0)
         NdemandSoil = (1.d0-NFIXF) * NDEMTO * NLIMIT
         NdemandBioFix =  NFIXF * NDEMTO * NLIMIT

      end if
         
      return

      case (4)

      if(flCropNut) then

!        Total N uptake (kg N ha-1 d-1) from soil and by biological N fixation         
         NUPTR = (MAX(0.d0, MIN(NdemandSoil, NsupplySoil) ))/DELT
         NFIXTR = (MAX(0.d0, NdemandBioFix))/DELT

!        Calling the subroutine to estimate the translocatable nutrients in leaves, stem, roots and
!        storage organs (kg N ha-1)
         CALL NTRLOC(ANLV,ANST,ANRT,WLV,WST,WRT,RNFLV,RNFST,RNFRT,      &
     &                  FNTRT,ATNLV,ATNST,ATNRT,ATN)

!        N supply to the storage organs (kg N ha-1 d-1)      
         NSUPSO = INSW (DVS-DVSNT,0.0d0,ATN/TCNT)

!        Rate of N uptake in grains (kg N ha-1 d-1)
         RNSO =  MIN (NDEMSO,NSUPSO)

!        Calling the subroutine to calculate nutrient translocation from leaves, stem, and roots (kg N ha-1 d-1)
         CALL NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN,RNTLV,RNTST,RNTRT)

!        Calling the subroutine to compute the partitioning of the total
!        nutrient uptake rate (NUPTR) over the leaves, stem and roots (kg N ha-1 d-1)
         CALL RNUSUB(NDEML,NDEMS,NDEMR,NUPTR,                           &
     &                  NFIXTR,NDEMTO, RNULV,RNUST,RNURT)
        
!        Calling routine to calculate nutrient losses due to dying leaves, stems           
!        and roots (kg N ha-1 d-1)    
         CALL RNLD(DRLV,DRRT,DRST,RNFLV,RNFRT,RNFST,                    &
     &                RNLDLV,RNLDRT,RNLDST)


!   - ---Rate of change of N in crop organs   
         RNLV = RNULV - RNTLV - RNLDLV
         RNST = RNUST - RNTST - RNLDST
         RNRT = RNURT - RNTRT - RNLDRT

!-----   Total N  uptake by crop over time (kg N ha-1) from soil and by biological fixation
         NUPTT = NUPTT + NUPTR*DELT
         NFIXTT= NFIXTT+ NFIXTR*DELT

!-----   Actual N amount in various living organs and total living N amount(kg N ha-1)
         ANLV =  max(0.0d0, (ANLV + RNLV*DELT) )
         ANST =  max(0.0d0, (ANST + RNST*DELT) )
         ANRT =  max(0.0d0, (ANRT + RNRT*DELT) )
         ANSO =  ANSO + RNSO*DELT
!!!         NLIVT=  ANLV + ANST + ANRT + ANSO

!-----   N losses from leaves, roots and stems due to senescence and total N loss (kg N ha-1)
         NLOSSL =  NLOSSL + RNLDLV*DELT
         NLOSSR =  NLOSSR + RNLDRT*DELT
         NLOSSS =  NLOSSS + RNLDST*DELT
!!!         NLOSST =  NLOSSL + NLOSSR + NLOSSS

!----    total N  in living and dead roots
!!!         NROOT= ANRT + NLOSSR

!       increment values of dead weight of plant organs,
!       to be used in the soil nutrient submodel
         idwrt = drrt*delt
         idwlv = idwlvSoil
         idwst = 0.0d0
         idwso = 0.0d0
         iNLOSSR =  RNLDRT*DELT
         iNLOSSL =  FraDeceasedLvToSoil * RNLDLV*DELT
         NLOSSLDeceasedLvToSoil =  NLOSSLDeceasedLvToSoil + iNLOSSL
         NLOSSL =  NLOSSL - iNLOSSL
         iNLOSSS =  0.0d0
         iNLOSSO =  0.0d0
         HarLosOrm_rt = 0.0d0; HarLosOrm_lv = 0.0d0; HarLosOrm_st = 0.0d0
         HarLosOrm_dwrt = 0.0d0; HarLosOrm_dwlv = 0.0d0; HarLosOrm_dwst = 0.0d0
         HarLosOrm_so = 0.0d0; HarLosOrm_tot = 0.0d0 
         HarLosNit_rt = 0.0d0; HarLosNit_lv = 0.0d0 
         HarLosNit_st = 0.0d0; HarLosNit_so = 0.0d0 
         HarLosNit_dwrt = 0.0d0; HarLosNit_dwlv = 0.0d0 
         HarLosNit_dwst = 0.0d0; HarLosNit_dwso = 0.0d0
!        during the last day of the crop period: add the weight of living roots 
!        to the dead roots and reset living weight to zero
         if (flHarvestDay .or. (dvs.ge.dvsend) .or.                     &
     &                 abs(t1900-1.0d0-cropend(icrop)).lt.1.0d-3 ) then
            HarLosOrm_rt = wrt
            HarLosOrm_dwlv =  FraHarLosOrm_lv * dwlv
            HarLosOrm_lv   = FraHarLosOrm_lv * wlv + HarLosOrm_dwlv
            HarLosOrm_dwst =  FraHarLosOrm_st * dwst
            HarLosOrm_st   = FraHarLosOrm_st * wst + HarLosOrm_dwst
            HarLosOrm_dwso =  FraHarLosOrm_so * dwso
            HarLosOrm_so   = FraHarLosOrm_so * wso + HarLosOrm_dwso
            HarLosOrm_tot = HarLosOrm_rt + FraHarLosOrm_lv * wlv +      & 
     &             FraHarLosOrm_st * wst + FraHarLosOrm_so * wso
            wrt = wrt - HarLosOrm_rt
            wlv = wlv - FraHarLosOrm_lv * wlv
            wst = wst - FraHarLosOrm_st * wst
            wso = wso - FraHarLosOrm_so * wso
            idwrt = idwrt + HarLosOrm_rt
            idwlv = idwlv + HarLosOrm_lv
            idwst = idwst + HarLosOrm_st
            idwso = idwso + HarLosOrm_so
            HarLosNit_rt = ANRT
            HarLosNit_dwlv = FraHarLosOrm_lv * NLOSSL 
            HarLosNit_lv = FraHarLosOrm_lv * ANLV + HarLosNit_dwlv 
            HarLosNit_dwst = FraHarLosOrm_st * NLOSSS 
            HarLosNit_st = FraHarLosOrm_st * ANST + HarLosNit_dwst 
            HarLosNit_dwso = FraHarLosOrm_so * 0.0d0 
            HarLosNit_so = FraHarLosOrm_so * ANSO + HarLosNit_dwso 
            iNLOSSL = iNLOSSL + HarLosNit_lv
            iNLOSSS = iNLOSSS + HarLosNit_st
            iNLOSSO = iNLOSSO + HarLosNit_so
            iNLOSSR = iNLOSSR + HarLosNit_rt
            ANLV = ANLV - FraHarLosOrm_lv * ANLV
            ANST = ANST - FraHarLosOrm_st * ANST
            ANSO = ANSO - FraHarLosOrm_so * ANSO
            ANRT = ANRT - HarLosNit_rt
        end if

      endif

      if(flCropNut) then

! ----- CHECK and WRITE MASS BALANCE: dry matter of crop

!       output of OM balance1: from air to partitioning (kg/ha DM CH2O)
        call outbalcropOM1(2,pathwork,outfil,project,date,daycrop,      &
     &       t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)

! -     OM balance2: storage difference(kg/ha DM CH2O)
        storagediff = HarLosOrm_rt + HarLosOrm_lv + HarLosOrm_st +      &
     &                HarLosOrm_so 
        storagediff = storagediff + (wlv+wst+wso+wrt) -                 &
     &                          (wlvt0+wstt0+wsot0+wrtt0)
        ombalan = storagediff - ( (grlv+grst+grso+grrt)*delt -          &
     &            (drlv+drst+drso+drrt)*delt ) -                        &
     &            (HarLosOrm_dwlv + HarLosOrm_dwst + HarLosOrm_dwso)
        if (abs(ombalan) .ge. 1.0d0) then
           write(messag,'(a,f8.3)')                                     &
     &      ' Warning Wofost: OM balance2 not 0, OMBAL = ',ombalan
           call warn ('wofost',messag,logf,swscre)
!     &     ' OM balance2 not 0, simulation stopped OMBAL=',ombalan
!           call fatalerr ('wofost',messag)
        endif
!       output of OM balance2
        call outbalcropom2(2,pathwork,outfil,project,date,daycrop,      &
     &         t,dvs,tsum,storagediff,wlv,wst,wso,wrt,delt,             &
     &         grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)

! ----- CHECK and WRITE MASS BALANCE: nitrogen of crop  

        NBALAN =  NUPTT + NFIXTT + (ANLVI+ANSTI+ANRTI+ANSOI)            &
     &      - (ANLV+ANST+ANRT+ANSO) - (NLOSSL+NLOSSR+NLOSSS)            &
     &      - NLOSSLDeceasedLvToSoil                                    &
     &      - (HarLosNit_lv+HarLosNit_st+HarLosNit_so+HarLosNit_rt)     &
     &      +  HarLosNit_dwlv + HarLosNit_dwst + HarLosNit_dwso

!       output of N balance
        call outbalcropN(2,pathwork,outfil,project,date,daycrop,        &
     &         t,dvs,tsum,NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,ANLV,&
     &         ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,nni)
        IF (abs(NBALAN) .GE. 1.0d-03) then
           write(messag,'(1a,i6,a,f8.3)') ' Nitrogen balance not 0,'//  &
     &     ' simulation stopped, day = ', daycum,' NBAL=',NBALAN
!*           call fatalerr ('wofost',messag)
           call warn ('CropGrowth_Wofost',messag,logf,swscre)
        endif
 
        if (flHarvestDay .or. (dvs.ge.dvsend) .or.                      &
     &                 abs(t1900-1.0d0-cropend(icrop)).lt.1.0d-3 ) then
          gwst  = 0.0d0
          gwrt  = 0.0d0
          gwso  = 0.0d0
          grlv  = 0.0d0
          NdemandSoil = 0.0d0
          HarLosOrm_tot = 0.0d0
        endif
      endif

! --- rooting depth
      if(swdmi2rd.eq.0) then
        rd = rd+rr
      elseif(swdmi2rd.eq.1) then
!       rooting depth limitation by relative dry matter increase (dmi/dmipot)
        if(dmi.gt.nihil) then
          if(dmipot .ge. 1.0d-3) then
            rd = rd + rr * dmi/dmipot
          endif
!         maximum rooting depth
          rdm = min(rdmax,rdc)
          rd = min(rd,rdm)
        endif
      endif


! --- crop factor or crop height
      if (swcf.ne.3) then
        cf = afgen (cftb,(2*magrs),dvs)
        ch = afgen (chtb,(2*magrs),dvs)
      else
        cf        = afgen (cftb,72,lai)
        cfeic     = afgen (cfeictb,72,lai)
        ch        = afgen(chtb,72,lai)
        siccapact = siccaplai*lai
      endif


! ----- UPDATE STATES of dry matter organs
      wlvt0 = wlv
      wstt0 = wst
      wsot0 = wso
      wrtt0 = wrt

! --- crop finish conditions based on dvs or lai
!      if ( ((dvs.ge.dvsend) .or.                                         &
!     &     (lai.le.0.002d0.and.dvs.gt.0.5d0)).and. (.not. flCropHarvest))&
!      ckro - 20170316 suppressed to prevent low Ypot
      if (dvs.ge.dvsend) then
        flCropEmergence = .false.
        flCropHarvest = .true.
        flHarvestDay = .true.
        flCropReadFile = .true.
      endif

      case default
         call fatalerr ('Wofost', 'Illegal value for TASK')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine grass(task)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : detailed grass growth routine 
! ----------------------------------------------------------------------
      use variables
      implicit none
 
      integer   i1,task !RB20140317
      integer   idelaypot,idelay,i,swhydrlift, uco2,getun2,ifnd

      real(8)   laicr,lasum,mres,grazlivinglv
      real(8)   admi,afgen,amax,asrc,ccheck,cosld,cvf,rlwtb(22)
      real(8)   dalv,dayl,delt,dmi
      real(8)   drrt,drst,dslv,dslv1,dslv2,dslvt,dteff,dtga
      real(8)   fcheck,fl,fr,fs,drlv
      real(8)   fysdel,gass,gla,glaiex,glasol,grlv,grrt,grst
      real(8)   gwrt,gwst,pgass,rest,rmres
      real(8)   sinld,slat,teff,twlv,twst,glaiexpot
      real(8)   lasumpot,drst1,drst2
      real(8)   dtgapot,drst1pot,drst2pot
      real(8)   pgasspot,gasspot,rmrespot,mrespot,asrcpot,dmipot
      real(8)   admipot,grrtpot,drrtpot,gwrtpot,grlvpot,dslv1pot
      real(8)   dslv2pot,dslvpot,restpot,dalvpot,drlvpot
      real(8)   glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real(8)   dslvtpot,twlvpot,twstpot,tagpspot,tagps
      real(8)   dummy
      real(8)   dmharvest,dmlastharvest,dateharvest(999),dmgrazing
      real(8)   lsdb(100),DaysGrazing(100),UptGrazing(100),lsda(366)
      real(8)   DelayRegrowthTab(2*100),DaysGrazingtab(2*100),UptGrazingtab(2*100)
      real(8)   wrtmax,rdm,LossGrazing(100),LossGrazingtab(2*100)
      real(8)   watcon,uptgraz,tagprest,lossgraz  !,grazingfactor
      integer   daylastharvest,swharvest,daysgrazpot,daysgraz
      logical   fltsum200
      character(len=11) tmp
      character(len=200) filnam,messag
      logical   rdinqr
      real(8)   CO2,effc
      integer   ifindi,indexyr
      real(8)   dmmowtb(20)
      integer   maxdaymow
      logical   flmowingtb
!      integer   iseqgmtmp
      logical   flpotrelmf    ! Flag indicating calculation of limited attainable yield instead of potential yield (due to management factor)
      real(8)   relmf         ! relative Management factor that reduces crop growth (-)
      logical   flGrassGrowth ! flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]
      character(len=11) ::  dateGrassGrowth            ! date of start of GrassGrowth
      
      parameter (delt=1.0d0)
      include 'params.fi'

      save
! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initialization at start of crop =========================================
      call InitializeCrop
      rlwtb = 0.0d0
      dmmowtb = 0.0d0
      dateharvest = 0.0d0
      DelayRegrowthTab = 0.0d0
      DaysGrazingtab = 0.0d0
      UptGrazingtab = 0.0d0
      lsda = 0.0d0

! --- read grass input data
      call readgrass (cropfil(icrop),pathcrop,tdwi,laiem,rgrlai,slatb,  &
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
     &  flCropCalendar,t1900,tstart,tend,swinco,cropstart(icrop),       &
     &  fltsumttd,tsumtime,tsumtemp,tsumdepth,cuptgraz,cuptgrazpot,     &
     &  saltmax,saltslope,salthead)

! --- initial values
      if(swcf.ne.3) then
        swinter = 1
      else
        swgc    = 1
        swinter = 3
      endif
      seqgrazmowpot = seqgrazmow

! --- development stage (not used by Grassland, instead Daynrs are used)
      dvs = -99.99d0

! --- maximum rooting depth
      rdm = min(rdmax,rdc)

! --- skip next initialization if crop parameters are read from *.END file
      if (t1900 - tstart .gt. 1.d-3 .or. swinco .ne. 3 .or.             &
     &   (t1900 - cropstart(icrop) .gt. -1.d-3 .and.                    &
     &                t1900 - cropstart(icrop) .lt. 1.d-3)) then

! ---   actual rooting depth
        rd = min(rdi,rdm)
        rdpot = min(rdi,rdm)

        iseqgm = 1
        iseqgmpot = iseqgm

! ---   initial values of crop parameters
        rid = 1.0d0
        fr = afgen (frtb,30,rid)
        fl = afgen (fltb,30,rid)
        fs = afgen (fstb,30,rid)
        sla(1) = afgen (slatb,30,rid)
        lvage(1) = 0.d0
        ilvold = 1
        idregr = 0
        slapot(1) = afgen (slatb,30,rid)
        lvagepot(1) = 0.d0
        ilvoldpot = 1
        idregrpot = 0

! ---   initial state variables of the crop
        wrt = fr*tdwi
        wrtpot = wrt
        wst = fs*(1.0d0-fr)*tdwi
        wstpot = wst
        wlv = laiem/sla(1)
        wlvpot = wlv
!     KRO-BOO-20160403: intro because comparison with Wofost
        laiem = wlv*sla(1)  ! is not input !
        lv(1) = wlv
        lvpot(1) = lv(1)
        lasum = laiem
        lasumpot = lasum     
        glaiex = 0.0d0
        laiexp = laiem
        laiexppot = laiem
        laimax = laiem
        lai = lasum+ssa*wst
        laipot = lai
        dwrt = 0.d0
        dwrtpot = dwrt
        dwlv = 0.d0
        dwlvpot = dwlv
        dwst = 0.d0
        dwstpot = dwst
        rid = dble(daycrop)
        daymow    = 0
        daymowpot = 0

! ---   initial summation variables of the crop
        tagp = wlv+wst
        tagppot = tagp
        tagpt = 0.0d0
        tagptpot = 0.0d0
        cuptgraz = 0.0d0
        cuptgrazpot = 0.0d0
        tsum200 = 0.0d0
        flGrassGrowth = .true.
        if (fltsumttd) then                               
          call sumttd('initial',flGrassGrowth,dateGrassGrowth)
        endif

! --- end skip above initialization if crop parameters are read from *.END file
      endif

      if (swcf.ne.3) then
        cf = afgen (cftb,(2*magrs),rid)
        ch = afgen (chtb,(2*magrs),rid)
      else
        cf        = afgen (cftb,(2*magrs),lai)
        cfeic     = afgen (cfeictb,(2*magrs),lai)
        ch        = afgen(chtb,(2*magrs),lai)
        siccapact = siccaplai*lai
      endif

! --- initialize matric flux potential and hleaf
      if (swdrought .eq. 2) then                                        
        call MatricFlux(1,h(1),1,dummy)
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,numnod
         twilt(i) = watcon(i,wiltpoint,cofgen,swsophy,numtab,sptab,     &
     &                     ientrytab)
         hroot(i) = h(i)
        enddo
        hleaf = -2000.d0
      endif                                                             

! --- CO2 impact
!     initialise 
      filnam = trim(pathcrop)//trim(cropfil(icrop))//'.crp'
      uco2 = getun2 (10,90,2)
      call rdinit(uco2,logf,filnam)
      flco2 = .false.
      if(rdinqr('flco2')) then
        call rdslog ('flco2',flco2)
      endif
      fco2amax = 1.0d0           ! factor to correct AMAX for CO2
      fco2eff = 1.0d0            ! factor to correct EFF for CO2
      fco2tra = 1.0d0            ! factor to correct TRA for CO2
      if(flco2) then
!        Read CO2 data
         CALL rdadou ('CO2AMAXTB', CO2AMAXTB, 30, ifnd) 
         CALL rdadou ('CO2EFFTB', CO2EFFTB, 30, ifnd) 
         CALL rdadou ('CO2TRATB', CO2TRATB, 30, ifnd) 
      endif
      close(uco2)
!     Read CO2-air data from a separate file
      if(flco2) then
         filnam = trim(pathcrop)//'Atmospheric.co2'
         uco2 = getun2 (30,90,2)
         call rdinit(uco2,logf,filnam)
         CALL rdainr ('CO2year', 1000, 3000, CO2year, mayrs, ifnd)
         CALL rdfdor ('CO2ppm',  10.0d0, 1000.0d0, CO2ppm, mayrs, ifnd)
         close(uco2)
      endif

! --- harvest
!     initialise 
      if (swharvest.eq.2) then
        iharvest = 1
        do while (t1900 .gt. dateharvest(iharvest))
          iharvest = iharvest + 1
        enddo
      endif      
      
      return

      case (2)

! === calculate potential rate and state variables ======================================

! --- rates of change of the grass variables ---------------------------------------------

      rid = dble(daycrop)

! -   grass growth initiated by tsum from 1st day of calendar year
      flGrassGrowth = .true.
      tsum200 = tsum200 + max(0.0d0,tav)
      if (fltsum200.and.tsum200.lt.200.d0) then
        flGrassGrowth = .false.
      endif
! -   grass growth initiated by temperature, time and depth
      if (fltsumttd .and. dateGrassGrowth.eq.'undefined') then                               
!       suppress growth as long as 3 criteria are not met
        call sumttd('dynamic',flGrassGrowth,dateGrassGrowth)
      endif

! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop.ge.idregrpot) then

! ===   daily dry matter production ===

!***************************************************************     
!       Assimilation correction for CO2 changes in atmosphere 
!       borrowed from Lintul4 added IS
!***************************************************************   
        if(flco2) then
            indexyr = ifindi (CO2year, mayrs, 1, mayrs, iyear)
            if (indexyr.lt.1 .or. indexyr.gt.mayrs) then
            Messag ='Input if CO2year or CO2ppm inconsistent, correct'
              call fatalerr ('wofost',messag)
            endif
            CO2 = CO2ppm(indexyr)
            fco2amax = afgen(CO2AMAXTB,30,CO2)
            fco2eff = afgen(CO2EFFTB,30,CO2)    
            fco2tra = afgen(CO2TRATB,30,CO2)    
        endif
        effc = fco2eff * eff
        amax = fco2amax * afgen (amaxtb,30,rid) * afgen (tmpftb,30,tavd)
        call astro(daynr,lat,rad,dayl,daylp,sinld,cosld,difpp,          &
     &             atmtr,dsinbe)
        call totass(dayl,amax,effc,laipot,kdif,rad,                     &
     &              difpp,dsinbe,sinld,cosld,dtgapot)
! ---   correction for low minimum temperature
        dtgapot = dtgapot * afgen (tmnftb,30,tmnr)
! ---   potential assimilation in kg ch2o per ha
        pgasspot = dtgapot * 30.0d0/44.0d0

! ---   reduction due to limited attainable maximum (optional)
        reltr = 1.0d0        
        if (flpotrelmf) reltr = relmf
        gasspot = pgasspot * reltr        

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmrespot=(rmr*wrtpot+rml*wlvpot+rms*wstpot)*afgen(rfsetb,30,rid)
        teff = q10**((tav-25.0d0)/10.0d0)
        mrespot = min (gasspot,rmrespot*teff)
        asrcpot = gasspot-mrespot

! ---   partitioning factors
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)
! ---   check on partitioning
        fcheck = fr+(fl+fs)*(1.0d0-fr) - 1.0d0
        if (abs(fcheck).gt.0.0001d0) then
          write(tmp,'(f6.3)') rid
          tmp = adjustl (tmp)
          Messag ='The sum of partitioning factors for leaves, stems'// &
     &    ' and storage organs is not equal to one at time '            &
     &    //trim(tmp)//'.'
          call fatalerr ('grass_pot',messag)
        endif

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmipot = cvf*asrcpot
! ---   check on carbon balance
        ccheck = (gasspot-mrespot-(fr+(fl+fs)*(1.0d0-fr))*dmipot/cvf)   &
     &         /max(0.0001d0,gasspot)      
        if (abs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr ('grass_pot',messag)
        endif


! ===   growth rate by plant organ ===

! ---   root length (not used, because Rooting depth is (for grassland) 
!       dependent on available root biomass (weight)
!        rrpot = min (rdm-rdpot,rri)
!        if (fr.le.0.or.pgasspot.lt.1.0d0) rrpot = 0.0d0

! ---   growth rate roots and aerial parts
! ---   after reaching a live weight of wrtmax (default 2500 kg), the
! ---   growth of the roots is balanced by the death of root tissue
        grrtpot = fr*dmipot
        if (wrtpot.gt.wrtmax) then
          drrtpot = grrtpot
          drrtpot = max(drrtpot,wrtpot*afgen (rdrrtb,30,rid))
        else
          drrtpot = wrtpot*afgen (rdrrtb,30,rid)
        endif
        gwrtpot = grrtpot-drrtpot

! ---   growth rate leaves

! ---   weight of new leaves
        admipot = (1.0d0-fr)*dmipot
        grlvpot = fl*admipot

! ---   death of leaves due to water stress or high lai
        dslv1pot = 0.0d0
        laicr = 3.2d0/kdif
        dslv2pot=wlvpot*max(0.0d0,                                      &
     &                  min(0.03d0,0.03d0*(laipot-laicr)/laicr))
        dslvpot = max (dslv1pot,dslv2pot) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        restpot = dslvpot*delt
        i1 = ilvoldpot

        do while (restpot.gt.lvpot(max(i1,1)).and.i1.ge.1)
          restpot = restpot-lvpot(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalvpot = 0.0d0
        if (lvagepot(max(i1,1)).gt.span.and.restpot.gt.0.and.           &
     &                          i1.ge.1) then
          dalvpot = lvpot(i1)-restpot
          restpot = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvagepot(max(i1,1)).gt.span)
          dalvpot = dalvpot+lvpot(i1)
          i1 = i1-1
        enddo

        dalvpot = dalvpot/delt

! ---   death rate leaves and growth rate living leaves
        drlvpot   = dslvpot+dalvpot

! ---   leaf area not to exceed exponential growth curve
        slatpot = afgen (slatb,30,rid)
        if (laiexppot.lt.6.0d0) then
          dteff = max (0.0d0,tav-tbase)
          glaiexpot = laiexppot*rgrlai*dteff
! ---   source-limited increase in leaf area
          glasolpot = grlvpot*slatpot
          glapot = min (glaiexpot,glasolpot)
! ---   adjustment of specific leaf area of youngest leaf class
          if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
        endif  

! ---   growth rate stems
        grstpot = fs*admipot
! ---   death of stems due to water stress is zero in case of potential growth
        drst1pot = 0.0d0
! ---   death of stems due to ageing
        drst2pot = afgen (rdrstb,30,rid)*wstpot
        drstpot = (drst1pot+drst2pot)/delt 
        gwstpot = grstpot-drstpot

! ----  integrals of the crop --------------------------------------------

!       after cutting, growth is initialized again and the weight of the sward is stored
        if(seqgrazmowpot(iseqgmpot).eq.2) then ! mowing
          
          flharvestpot = .false.
          
          if (swharvest.eq.1) then ! use dry matter threshold
      
            if (flmowingtb) then   ! use of mowing table
              dmharvest = afgen(dmmowtb,20,rid)
              daymowpot = daymowpot + 1
              if (tagppot .gt. dmharvest .or.                           &
     &                 (daymowpot.gt.maxdaymow.and.iseqgmpot.gt.1)) then
                flharvestpot = .true.
              endif
            else                   ! use of fixed threshold
              if (tagppot .gt. dmharvest .or. (daynr.gt.daylastharvest  &
     &          .and.tagppot.gt.dmlastharvest)) then
                flharvestpot = .true.
              endif
            endif
          
          endif
          
          if (swharvest.eq.2) then  ! use fixed mowing dates
            
            if(t1900.gt.dateharvest(iharvest)) then
              flharvestpot = .true.
            endif
          
          endif

        else if(seqgrazmowpot(iseqgmpot).eq.1) then ! grazing
!       Daily uptake by livestock is assumed to be an extra loss rate 

!         Assumption: no delay in regrowth during and after grazing
          idregrpot = daycrop

          if((swharvest.eq.2 .and. t1900.gt.dateharvest(iharvest)) .or. &
     &       (swharvest.eq.1 .and. tagppot.gt.dmgrazing) .or.           &
     &        flgrazingpot) then

!           Initialise Count nr of days with grazing
            if(.not.flgrazingpot) then
              idaysgrazpot = 0
            endif

!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgraz = lsda(iseqgmpot) *                                 &
     &                         afgen(uptgrazingtab,200,lsda(iseqgmpot))

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgraz = lsda(iseqgmpot) *                                &
     &                         afgen(lossgrazingtab,200,lsda(iseqgmpot))

!           verify if uptake is possible: tagprest should remain after grazing
            if ((tagppot-uptgraz-lossgraz).gt.tagprest) then
              flgrazingpot = .true.
              cuptgrazpot = cuptgrazpot + uptgraz

!           distribute grazing over stems and leaves (living and dead parts)
            wstpot = wstpot - (uptgraz+lossgraz) * wstpot / tagppot
            dwstpot = dwstpot - (uptgraz+lossgraz) * dwstpot / tagppot
            dwlvpot = dwlvpot - (uptgraz+lossgraz) * dwlvpot / tagppot
            grazlivinglv = (uptgraz+lossgraz) * wlvpot / tagppot

!             reduce leave weights
              i1 = ilvoldpot
              do while (grazlivinglv.gt.0.and.i1.ge.1)
                if (grazlivinglv.ge.lvpot(i1)) then
                  grazlivinglv = grazlivinglv-lvpot(i1)
                  lvpot(i1) = 0.0d0
                  i1 = i1-1
                else
                  lvpot(i1) = lvpot(i1)-grazlivinglv
                  grazlivinglv = 0.d0
                endif
              enddo

!             Check number of days with grazing
              daysgrazpot=int(afgen(daysgrazingtab,200,lsda(iseqgmpot)))
              idaysgrazpot = idaysgrazpot + 1
              if(idaysgrazpot .eq. daysgrazpot) then
                flgrazingpot = .false.
                iseqgmpot = iseqgmpot + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (flgrazingpot .or. swharvest.eq.2) then
              flgrazingpot = .false.
              iseqgmpot = iseqgmpot + 1
            endif

          endif
        endif
        
        if (flharvestpot) then
          iseqgmpot = iseqgmpot + 1
          flharvestpot = .false.
          slapot(1) = afgen (slatb,30,rid)
          fl = afgen (fltb,30,rid)
          fs = afgen (fstb,30,rid)
          wlvpot = mowrest / (1.d0 + (fs/fl))
          wstpot = fs/fl*wlvpot
          dwlvpot = 0.0d0
          dwstpot = 0.0d0
          lvagepot(1) = 0.0d0
          ilvoldpot = 1
          lasumpot = wlvpot * slapot(1)
          laiexppot = lasumpot
          lvpot(1) = wlvpot

          gwstpot = 0.0d0
          gwrtpot = 0.0d0
          drlvpot = 0.0d0
          drstpot = 0.0d0
          drrtpot = 0.0d0

          daymowpot = 0
          
          tagpspot =max(0.0d0,(tagppot-(wlvpot+dwlvpot+wstpot+dwstpot)))
          tagptpot = tagptpot + tagpspot

! ---     set regrowth delay
          idelaypot = int(afgen(DelayRegrowthTab,200,tagpspot))
          idregrpot = daycrop + idelaypot

        endif

        if (daycrop.ge.idregrpot) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvtpot = dslvpot*delt
          i1 = ilvoldpot
           do while (dslvtpot.gt.0.and.i1.ge.1)
            if (dslvtpot.ge.lvpot(i1)) then
              dslvtpot = dslvtpot-lvpot(i1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            else
              lvpot(i1) = lvpot(i1)-dslvtpot
              dslvtpot = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (lvagepot(max(i1,1)).ge.span.and.i1.ge.1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          ilvoldpot = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvoldpot,1,-1
            lvpot(i1+1) = lvpot(i1)
            slapot(i1+1) = slapot(i1)
            lvagepot(i1+1) = lvagepot(i1)+fysdel*delt
          enddo
          ilvoldpot = ilvoldpot+1

! ---     new leaves in class 1
          lvpot(1) = grlvpot*delt
          slapot(1) = slatpot
          lvagepot(1) = 0.d0

! ---     calculation of new leaf area and weight
          lasumpot = 0.d0
          wlvpot = 0.d0
          do i1 = 1,ilvoldpot
            lasumpot = lasumpot+lvpot(i1)*slapot(i1)
            wlvpot = wlvpot+lvpot(i1)
          enddo

          laiexppot = laiexppot+glaiexpot*delt

        endif

! ---   dry weight of living plant organs
        wrtpot = wrtpot+gwrtpot*delt
        wstpot = wstpot+gwstpot*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        dwrtpot = dwrtpot+drrtpot*delt
        dwlvpot = dwlvpot+drlvpot*delt
        dwstpot = dwstpot+drstpot*delt

! ---   dry weight of dead and living plant organs
        twlvpot = wlvpot+dwlvpot
        twstpot = wstpot+dwstpot
        tagppot = twlvpot+twstpot

! ---   leaf area index
        laipot = lasumpot+ssa*wstpot
!       prevent immediate lai reduction at emergence
!       KRO-BOO-20160403: suppressed because deviates from Wofost
!       laipot = max(laipot, laiem)

! ---   rooting depth as function of available root weight
        rdpot = afgen (rlwtb,22,wrtpot)
        rdpot = min(rdpot,rdm)
      endif

      return

      case (3)

! === calculate actual rate and state variables ======================================

! --- rates of change of the crop variables ---------------------------------------------

!     correction of potential transpiration in relation to reference crop
!     (default = 1.0, range = 0.8 - 1.2)
      ptra = fco2tra*ptra

! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop.ge.idregr) then

! ===   daily dry matter production ===

!***************************************************************     
!       Assimilation correction for CO2 changes in atmosphere 
!       borrowed from Lintul4 added IS
!***************************************************************   
        if(flco2) then
            indexyr = ifindi (CO2year, mayrs, 1, mayrs, iyear)
            if (indexyr.lt.1 .or. indexyr.gt.mayrs) then
            Messag ='Input if CO2year or CO2ppm inconsistent, correct'
              call fatalerr ('wofost',messag)
            endif
            CO2 = CO2ppm(indexyr)
            fco2amax = afgen(CO2AMAXTB,30,CO2)
            fco2eff = afgen(CO2EFFTB,30,CO2)    
            fco2tra = afgen(CO2TRATB,30,CO2)    
        endif
        effc = fco2eff * eff
        amax = fco2amax * afgen (amaxtb,30,rid) * afgen (tmpftb,30,tavd)
! ---   gross assimilation
        call astro(daynr,lat,rad,dayl,daylp,sinld,cosld,difpp,          &
     &             atmtr,dsinbe)
        call totass(dayl,amax,effc,lai,kdif,rad,                        &
     &              difpp,dsinbe,sinld,cosld,dtga)
! ---   correction for low minimum temperature
        dtga = dtga * afgen (tmnftb,30,tmnr)
! ---   potential assimilation in kg ch2o per ha
        pgass = dtga * 30.0d0/44.0d0

! ---   water stress reduction of pgass to gass and limited attainable maximum
        if(abs(ptra).lt.nihil) then
          reltr = 1.0d0
        else
          reltr = max(0.0d0,min(1.0d0,tra/ptra))
        endif
        gass = pgass * reltr

! ---   relative management factor that reduces crop growth        
        gass = gass * relmf

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmres = (rmr*wrt+rml*wlv+rms*wst)*afgen(rfsetb,30,rid)
        teff = q10**((tav-25.0d0)/10.0d0)
        mres = min (gass,rmres*teff)
        asrc = gass-mres

! ---   partitioning factors (relevant for restart)
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmi = cvf*asrc
! ---   check on carbon balance
        ccheck = (gass-mres-(fr+(fl+fs)*(1.0d0-fr))*dmi/cvf)            &
     &         /max(0.0001d0,gass)      
        if (abs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr ('grass_act',messag)
        endif

! ===   growth rate by plant organ ===

! ---   root length is not used, because for grassland rooting depth depends on available root biomass
!        rr = min (rdm-rd,rri)
!        if (fr.le.0.or.pgass.lt.1.0d0) rr = 0.0d0
!        rr = 0.0d0

! ---   growth rate roots and aerial parts
! ---   after reaching a live weight of wrtmax (default 2500 kg), the
! ---   growth of the roots is balanced by the death of root tissue
        grrt = fr*dmi
        if (wrt.gt.wrtmax) then
!original drrt = wrt - wrtmax
          drrt = grrt
          drrt = max(drrt,wrt*afgen (rdrrtb,30,rid))
!         CO2 loss
!d         co2rootloss = drrt*44.0d0/33.0d0
!original  grrt = 0.0d0
        else
          drrt = wrt*afgen (rdrrtb,30,rid)
        endif
        admi = (1.0d0-fr)*dmi
!original drrt = 0.0d0
        gwrt = grrt-drrt

!       CO2 fixation and loss
!d       co2rootfix = grrt*44.0d0/33.0d0
!d       co2rootloss = co2rootloss + drrt*44.0d0/33.0d0

! ---   growth rate leaves

! ---   weight of new leaves
        grlv = fl*admi

! ---   death of leaves due to water stress or high lai
        dslv1 = wlv*(1.0d0-reltr)*perdl
        laicr = 3.2d0/kdif
        dslv2 = wlv*max(0.0d0,min(0.03d0,0.03d0*(lai-laicr)/laicr))
        dslv = max (dslv1,dslv2) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        rest = dslv*delt
        i1 = ilvold

        do while (rest.gt.lv(max(i1,1)).and.i1.ge.1)
          rest = rest-lv(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalv = 0.0d0
        if (lvage(max(i1,1)).gt.span.and.rest.gt.0.and.i1.ge.1) then
          dalv = lv(i1)-rest
          rest = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvage(max(i1,1)).gt.span)
          dalv = dalv+lv(i1)
          i1 = i1-1
        enddo

        dalv = dalv/delt

! ---   death rate leaves and growth rate living leaves
        drlv   = dslv+dalv

! ---   physiologic ageing of leaves per time step
        slat = afgen (slatb,30,rid)

! ---   leaf area not to exceed exponential growth curve
        if (laiexp.lt.6.0d0) then
          dteff = max (0.0d0,tav-tbase)
          glaiex = laiexp*rgrlai*dteff
! ---     source-limited increase in leaf area
          glasol = grlv*slat
          gla = min (glaiex,glasol)
! ---     adjustment of specific leaf area of youngest leaf class
          if (grlv.gt.0.0d0) slat = gla/grlv
        endif  

! ---   growth rate stems
        grst = fs*admi
! ---   death of stems due to water stress
        drst1 = wst*(1.0d0-reltr)*perdl
! ---   death of stems due to ageing
        drst2 = afgen (rdrstb,30,rid)*wst
        drst = (drst1+drst2)/delt 
        gwst = grst-drst

! ----  integrals of the crop --------------------------------------------

!       after cutting, growth is initialized again and the weight of the sward is stored
        if(seqgrazmow(iseqgm).eq.2) then ! mowing

          flharvest = .false.
          
          if (swharvest.eq.1) then ! use dry matter threshold
      
            if (flmowingtb) then   ! use of mowing table
              dmharvest = afgen(dmmowtb,20,rid)
              daymow = daymow + 1
              if (tagp .gt. dmharvest .or.                              &
     &                 (daymow.gt.maxdaymow.and.iseqgm.gt.1)) then
                flharvest = .true.
              endif
            else                   ! use of fixed threshold
              if (tagp .gt. dmharvest .or. (daynr.gt.daylastharvest     &
     &          .and.tagp.gt.dmlastharvest)) then
                flharvest = .true.
              endif
            endif
          
          endif
          
          if (swharvest.eq.2) then  ! use fixed mowing dates
            
            if(t1900.gt.dateharvest(iharvest)) then
              iharvest = iharvest + 1
              flharvest = .true.
            endif
          
          endif
              

        else if(seqgrazmow(iseqgm).eq.1) then ! grazing
!       Daily uptake by livestock is assumed to be an extra loss rate 

!         Assumption: no delay in regrowth during and after grazing
          idregr = daycrop

          if((swharvest.eq.2 .and. t1900.gt.dateharvest(iharvest)) .or. &
     &       (swharvest.eq.1 .and. tagp.gt.dmgrazing) .or.              &
     &        flgrazing) then

!           Initialise Count nr of days with grazing
            if(.not.flgrazing) then
              idaysgraz = 0
            endif

!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgraz = lsda(iseqgm)*afgen(uptgrazingtab,200,lsda(iseqgm))

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgraz = lsda(iseqgm) *                                   &
     &                           afgen(lossgrazingtab,200,lsda(iseqgm))

!           verify if uptake is possible: tagprest should remain after grazing
            if ((tagp-uptgraz-lossgraz).gt.tagprest) then
              flgrazing = .true.
              cuptgraz = cuptgraz + uptgraz

!             distribute grazing over stems and leaves (living and dead parts)
              wst = wst - (uptgraz+lossgraz) * wst / tagp
              dwst = dwst - (uptgraz+lossgraz) * dwst / tagp
              dwlv = dwlv - (uptgraz+lossgraz) * dwlv / tagp
              grazlivinglv = (uptgraz+lossgraz) * wlv / tagp

!             reduce leave weights
              i1 = ilvold
              do while (grazlivinglv.gt.0.and.i1.ge.1)
                if (grazlivinglv.ge.lv(i1)) then
                  grazlivinglv = grazlivinglv-lv(i1)
                  lv(i1) = 0.0d0
                  i1 = i1-1
                else
                  lv(i1) = lv(i1)-grazlivinglv
                  grazlivinglv = 0.d0
                endif
              enddo

!             Check number of days with grazing
              daysgraz = int(afgen(daysgrazingtab,200,lsda(iseqgm))) 
              idaysgraz = idaysgraz + 1
              if(idaysgraz .eq. daysgraz) then
                flgrazing = .false.
                iseqgm = iseqgm + 1
                iharvest = iharvest + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (flgrazing .or. swharvest.eq.2) then
              flgrazing = .false.
              iseqgm = iseqgm + 1
              iharvest = iharvest + 1
            endif

          endif
        endif

        if (flharvest) then
          iseqgm = iseqgm + 1
          flharvest = .false.
          sla(1) = afgen (slatb,30,rid)
          fl = afgen (fltb,30,rid)
          fs = afgen (fstb,30,rid)
          wlv = mowrest / (1.d0 + (fs/fl))
          wst = fs/fl*wlv
          dwlv = 0.0d0
          dwst = 0.0d0
          lvage(1) = 0.0d0
          ilvold = 1
          lasum = wlv * sla(1)
          laiexp = lasum
          lv(1) = wlv

          gwst = 0.0d0
          gwrt = 0.0d0
          drlv = 0.0d0
          drst = 0.0d0
          drrt = 0.0d0

          tagps = max (0.0d0,(tagp-(wlv+dwlv+wst+dwst)))
          tagpt = tagpt + tagps

          daymow = 0
          
! ---     regrowth delay
          idelay = int(afgen(DelayRegrowthTab,200,tagps))
          idregr = daycrop + idelay

        endif

        if (daycrop.ge.idregr) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvt = dslv*delt
          i1 = ilvold
          do while (dslvt.gt.0.and.i1.ge.1)
            if (dslvt.ge.lv(i1)) then
              dslvt = dslvt-lv(i1)
              lv(i1) = 0.0d0
              i1 = i1-1
            else
              lv(i1) = lv(i1)-dslvt
              dslvt = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (lvage(max(i1,1)).ge.span.and.i1.ge.1)
              lv(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          ilvold = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvold,1,-1
            lv(i1+1) = lv(i1)
            sla(i1+1) = sla(i1)
            lvage(i1+1) = lvage(i1)+fysdel*delt
          enddo
          ilvold = ilvold+1

! ---     new leaves in class 1
          lv(1) = grlv*delt
          sla(1) = slat
          lvage(1) = 0.d0 

! ---     calculation of new leaf area and weight
          lasum = 0.d0
          wlv = 0.d0
          do i1 = 1,ilvold
            lasum = lasum+lv(i1)*sla(i1)
            wlv = wlv+lv(i1)
          enddo

          laiexp = laiexp+glaiex*delt

        endif

! ---   dry weight of living plant organs
        wrt = wrt+gwrt*delt
        wst = wst+gwst*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        dwrt = dwrt+drrt*delt
        dwlv = dwlv+drlv*delt
        dwst = dwst+drst*delt

! ---   dry weight of dead and living plant organs
        twlv = wlv+dwlv
        twst = wst+dwst
        tagp = twlv+twst

! ---   leaf area index
        lai = lasum+ssa*wst
        laimax = max (lai,laimax)

! ---   rooting depth as function of root weight
        rd = afgen (rlwtb,22,wrt)
        rd= min(rd,rdm)

      endif

      case default
         call fatalerr ('Grass', 'Illegal value for TASK')
      end select

      return
      end
      
!** PROGRAM:
!**    WOFOSTARABLELANDGERM
!** 
!** DESCRIPTION:
!**    Subroutine WOFOSTARABLELANDGERM simulates the germination of arable 
!**               land crop
!** 
!************************************************************************
!c
!c
      subroutine wofostarablelandgerm(task,prz1)
!c
      use        Variables
!c
      implicit   none
!c
!c     arguments
      integer(4) task     ! 0=init, 1= simulation of 1 day
      real(8)    prz1     ! Soil water pressure head in root zone at depth z1 (L)
      real(8)    pFz1
!c
!c     Locals
      real(8)    tsumemesub
!c
      if (task .eq. 1) then !-------start germination simulation--start task 1
!c
        tsumgerm = 0.0d0
        dvs      = 0.0d0
      else !------------------------------update temperature sum--start task 2
!c
        pFz1 = DLOG10(MAX(1.0d0,-prz1))
        if (prz1 .lt. hdrygerm) then
!c
!c         Dry situation

          tsumemesub = agerm * pFz1 + cgerm
        elseif (prz1 .ge. hdrygerm .and. prz1 .le. hwetgerm) then       
!c
!c         Optimal situation
          tsumemesub = tsumemeopt
        else
!c
!c         Wet situation, avoid DLOG10-error due to positive prz1
          tsumemesub = -agerm * pFz1 + bgerm
        endif
!c
!c       Update of temperature sum, for the time step of 1 day

        if (tav .gt. TBASEM)then
           if( tav .lt. TEFFMX) then
              if(tsumemesub.lt.0.1d0) then
                  tsumgerm = tsumgerm + (tav-TBASEM)
              else
                  tsumgerm = tsumgerm +(tsumemeopt/tsumemesub)*(tav-TBASEM)
              endif
           else
              if(tsumemesub.lt.0.1d0) then
                  tsumgerm = tsumgerm + (TEFFMX-TBASEM)
              else
                  tsumgerm = tsumgerm +(tsumemeopt/tsumemesub)*             &
     &                                             (TEFFMX-TBASEM)
              endif
           endif
        endif

      endif
!c
      return
      end

! ----------------------------------------------------------------------
      subroutine totass (dayl,amax,eff,lai,kdif,avrad,difpp,            &
     &                   dsinbe,sinld,cosld,dtga)

!*  Purpose: This routine calculates the daily total gross CO2
!*           assimilation by performing a Gaussian integration over
!*           time. At three different times of the day, irradiance is
!*           computed and used to calculate the instantaneous canopy
!*           assimilation, whereafter integration takes place. More
!*           information on this routine is given by Spitters et al.
!*           (1988).

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  DAYL    R8  Astronomical daylength (base = 0 degrees)     h      O
!*  AMAX    R8  Assimilation rate at light saturation      kg CO2/   I
!*                                                        ha leaf/h   
!*  EFF     R8  Initial light use efficiency              kg CO2/J/  I
!*                                                        ha/h m2 s   
!*  LAI     R8  Leaf area index                             ha/ha    I
!*  KDIF    R8  Extinction coefficient for diffuse light             I
!*  AVRAD   R8  Daily shortwave radiation                  J m-2 d-1 I
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of
!*              light                                      J m-2 s-1 I
!*  DSINBE  R8  Daily total of effective solar height         s      I
!*  SINLD   R8  Seasonal offset of sine of solar height       -      I
!*  COSLD   R8  Amplitude of sine of solar height             -      I
!*  DTGA    R8  Daily total gross assimilation           kg CO2/ha/d O

!*  FATAL ERROR CHECKS: none
!*  SUBROUTINES and FUNCTIONS called : ASSIM
!*  FILE usage : none

!*  Authors: Daniel van Kraalingen 
!*  Date   : April 1991

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

      implicit none

!*     formal parameters
      real(8) dayl,amax,eff,lai,kdif,avrad,difpp,dsinbe,sinld,cosld,dtga

!*     local parameters
      integer i1
      real(8) hour,pi,sinb,par,pardif,pardir,fgros
      real(8) xgauss(3),wgauss(3)

      parameter (pi=3.1415926d0)
      save

!**
!*     gauss points and weights are stored in an array
      data xgauss /0.1127017d0, 0.5000000d0, 0.8872983d0/
      data wgauss /0.2777778d0, 0.4444444d0, 0.2777778d0/

!*     calculation of assimilation is done only when it will not be zero
!*     (AMAX >0, LAI >0)
      dtga  = 0.0d0
      if (amax.gt.0.0d0.and.lai.gt.0.0d0) then
         do 10 i1=1,3
            hour   = 12.0d0+0.5d0*dayl*xgauss(i1)
            sinb   = max(0.0d0,                                         &
     &                   sinld+cosld*cos(2.0d0*pi*(hour+12.0d0)/24.0d0))
            par    = 0.5d0*avrad*sinb*(1.0d0+0.4d0*sinb)/dsinbe
            pardif = min(par,sinb*difpp)
            pardir = par-pardif
            call assim (amax,eff,lai,kdif,sinb,pardir,pardif,fgros)
            dtga = dtga+fgros*wgauss(i1)
10       continue
         dtga = dtga*dayl
      end if

      return
      end

      subroutine assim (amax,eff,lai,kdif,sinb,pardir,pardif,fgros)

!*     Chapter 13 in documentation WOFOST Version 4.1 (1988)

!*     This routine calculates the gross CO2 assimilation rate of
!*     the whole crop, FGROS, by performing a Gaussian integration
!*     over depth in the crop canopy. At three different depths in
!*     the canopy, i.e. for different values of LAI, the
!*     assimilation rate is computed for given fluxes of photosynthe-
!*     tically active radiation, whereafter integration over depth
!*     takes place. More information on this routine is given by
!*     Spitters et al. (1988). The input variables SINB, PARDIR
!*     and PARDIF are calculated in routine TOTASS.

!*     Subroutines and functions called: none.
!*     Called by routine TOTASS.

!*     Author: D.W.G. van Kraalingen, 1986

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  AMAX    R8  Maximum CO2 assimilation rate              kg/ha/hr  I
!*  EFF     R8  Light use efficiency of a leaf         kg CO2 / J adsorbed  I
!*  LAI     R8  Leaf area index                               -      I
!*  KDIF    R8  Extinction coefficient for diffuse visible light -   I
!*  SINB    R8  ...........nog invullen ..............               -      I
!*  PARDIR  R8  ...........nog invullen ..............                 -      I
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of     
!*              light                                      J m-2 s-1 I
!*  PARDIF  R8  ...........nog invullen ..............             -      I
!*  FGROS   R8  ...........nog invullen ..............           s      O
!**    
!*13.1 declarations
      implicit none

!*     formal parameters
      real(8) amax,eff,lai,kdif,sinb,pardir,pardif,fgros

!*     local parameters
      integer i
      real(8) scv,refh,refs,kdirbl,kdirt,laic,visdf,vist,visd,visshd
      real(8) fgrsh,vispp,fgrsun,fslla,fgl
      real(8) xgauss(3),wgauss(3)

      save

!*     initialize GAUSS array and scattering coefficient
      data xgauss /0.1127017d0, 0.5000000d0, 0.8872983d0/
      data wgauss /0.2777778d0, 0.4444444d0, 0.2777778d0/
      data scv /0.2d0/

!*13.2 extinction coefficients KDIF,KDIRBL,KDIRT
      refh   = (1.0d0-sqrt(1.0d0-scv))/(1.0d0+sqrt(1.0d0-scv))
      refs   = refh*2.0d0/(1.0d0+1.6d0*sinb)
      kdirbl = (0.5d0/sinb)*kdif/(0.8d0*sqrt(1.0d0-scv))
      kdirt  = kdirbl*sqrt(1.0d0-scv)

!*13.3 three-point Gaussian integration over LAI
      fgros  = 0.0d0
      do 10 i=1,3
         laic   = lai*xgauss(i)
!*        absorbed diffuse radiation (VISDF),light from direct
!*        origine (VIST) and direct light(VISD)
         visdf  = (1.0d0-refs)*pardif*kdif  *exp (-kdif  *laic)
         vist   = (1.0d0-refs)*pardir*kdirt *exp (-kdirt *laic)
         visd   = (1.0d0-scv) *pardir*kdirbl*exp (-kdirbl*laic)
!*        absorbed flux in W/m2 for shaded leaves and assimilation
         visshd = visdf+vist-visd
         fgrsh  = amax*(1.0d0-exp(-visshd*eff/max(2.0d0,amax)))
!*        direct light absorbed by leaves perpendicular on direct
!*        beam and assimilation of sunlit leaf area
         vispp  = (1.0d0-scv)*pardir/sinb
         if (vispp.le.0.0d0) then
            fgrsun = fgrsh
         else
            fgrsun = amax*(1.0d0-(amax-fgrsh)                           &
     &          *(1.0d0-exp (-vispp*eff/max(2.0d0,amax)))/ (eff*vispp))
         end if
!*        fraction of sunlit leaf area (FSLLA) and local
!*        assimilation rate (FGL)
         fslla  = exp (-kdirbl*laic)
         fgl    = fslla*fgrsun+(1.0d0-fslla)*fgrsh
!*        integration
         fgros  = fgros+fgl*wgauss(i)
10    continue

      fgros  = fgros*lai
      return
      end

      subroutine astro (iday,lat,avrad,                                 &
     &                  dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe)

!*  Purpose: This subroutine calculates astronomic daylength,
!*           diurnal radiation characteristics such as the atmospheric
!*           transmission, diffuse radiation etc.. This routine has
!*           been modified so that it uses arrays to hold some input
!*           output variables for faster processing 

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  IDAY    I4  Day number (Jan 1st = 1)                      -      I
!*  LAT     R8  Latitude of the site                       degrees   I
!*  AVRAD   R8  Daily shortwave radiation                  J m-2 d-1 I
!*  DAYL    R8  Astronomical daylength (base = 0 degrees)     h      O
!*  DAYLP   R8  Astronomical daylength (base =-4 degrees)     h      O
!*  SINLD   R8  Seasonal offset of sine of solar height       -      O
!*  COSLD   R8  Amplitude of sine of solar height             -      O
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of     
!*              light                                      J m-2 s-1 O
!*  ATMTR   R8  Daily atmospheric transmission                -      0
!*  DSINBE  R8  Daily total of effective solar height         s      O

!*  FATAL ERROR CHECKS: none
!*  SUBROUTINES and FUNCTIONS called : none
!*  FILE usage : none

!*  Authors: Daniel van Kraalingen
!*  Date   : April 1991

!*  Modification: Include checks for 0<=daylength<=24 hour
!*                Remove caching of results
!*  Author      : Allard de Wit
!*  Date        : January 2011

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

      implicit none
!*     formal parameters
      integer iday
      real(8) lat,avrad,dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe

!*     local parameters
      real(8) pi,angle,rad
      real(8) dec,sc,aob,aob_corr,angot,dsinb,frdif

      parameter (pi=3.1415926d0, angle=-4.0d0, rad=0.0174533d0)

!*     Error check on latitude
      if (abs(lat).gt.90.d0) call fatalerr                              &
     &   ('astro','lat > 90 or lat < -90')

!*     Declination and solar constant for this day
      dec = -asin(sin(23.45d0*rad)*cos(2.d0*pi*dble(iday+10)/365.0d0))
      sc  = 1370.d0*(1.d0+0.033d0*cos(2.d0*pi*dble(iday)/365.d0))

!*     calculation of daylength from intermediate variables
!*     SINLD, COSLD and AOB
      sinld = sin(rad*lat)*sin(dec)
      cosld = cos(rad*lat)*cos(dec)
      aob = sinld/cosld

!*     For very high latitudes and days in summer and winter a limit is  
!*     inserted to avoid math errors when daylength reaches 24 hours in 
!*     summer or 0 hours in winter.

!*     Calculate solution for base=0 degrees
      if (abs(aob).le.1.0d0) then
         dayl  = 12.0d0*(1.d0+2.d0*asin(aob)/pi)
!*        integrals of sine of solar height
         dsinb  = 3600.d0*                                              &
     &            (dayl*sinld+24.d0*cosld*sqrt(1.d0-aob**2)/pi)
         dsinbe = 3600.d0*                                              &
     &            (dayl*(sinld+0.4d0*(sinld**2+cosld**2*0.5d0))+  &
     &     12.d0*cosld*(2.d0+3.d0*0.4d0*sinld)*sqrt(1.d0-aob**2)/pi)
      else
         if (aob.gt.1.0d0)  dayl = 24.0d0
         if (aob.lt.-1.0d0) dayl =  0.0d0
!*        integrals of sine of solar height      
         dsinb  = 3600.d0*(dayl*sinld)
         dsinbe = 3600.d0*                                              &
     &            (dayl*(sinld+0.4d0*(sinld**2+cosld**2*0.5d0)))
      endif

!*     Calculate solution for base=-4 (ANGLE) degrees
      aob_corr = (-sin(angle*rad)+sinld)/cosld
      if (abs(aob_corr).le.1.0d0) then 
         daylp = 12.0d0*(1.d0+2.d0*asin(aob_corr)/pi)
      else
         if (aob_corr.gt.1.0d0)  daylp = 24.0d0
         if (aob_corr.lt.-1.0d0) daylp =  0.0d0
      endif

!*     extraterrestrial radiation and atmospheric transmission
      angot  = sc*dsinb
!*     Check for DAYL=0 as in that case the angot radiation is 0 as well
      if (dayl.gt.0.0d0) then
          atmtr = avrad/angot
      else
          atmtr = 0.0d0
      endif

!*     estimate fraction diffuse irradiation
      if (atmtr.gt.0.75d0) frdif = 0.23d0
      if (atmtr.le.0.75d0.and.atmtr.gt.0.35d0)                          &
     &  frdif = 1.33d0-1.46d0*atmtr
      if (atmtr.le.0.35d0.and.atmtr.gt.0.07d0)                          &
     &  frdif = 1.d0-2.3d0*(atmtr-0.07d0)**2
      if (atmtr.le.0.07d0) frdif = 1.d0

      difpp = frdif*atmtr*0.5d0*sc

      RETURN
      END

! ----------------------------------------------------------------------
      subroutine outbalcropN(task,pathwork,outfil,project,date,daycrop, &
     &             t,dvs,tsum,NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,ANLV,&
     &               ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,NNI) 
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output N balance files 
! ----------------------------------------------------------------------
      implicit none

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum   !,laipot,lai,cf,rdpot,rd,ch,crt0,crt1
!      real(8) cwdmpot,cwdm,wsopot,wso,wstpot,wst,wlvpot,wlv,wrtpot,wrt
      real(8) NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,ANLV
      real(8) ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,NNI
! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   getun,nba
      character(len=160) filnam,filtext

      save    nba

      comma = ',' 
    
      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.nba'
      nba = getun (20,90)
      call fopens(nba,filnam,'new','del')
      filtext = 'output of N-balance of detailed crop growth model'
      call writehead (nba,1,filnam,filtext,project)

      write (nba,100)
 100    format ('*',/,                                                  &
     & '*             day     day      -    grC   kg/ha   kg/ha   kg/', &
     & 'ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   ', &
     & 'kg/ha   kg/ha   kg/ha  kg/ha',/,                                &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,  NUPTT, NFIXTT,  ANLVI'&
     & ,',  ANSTI,  ANRTI,  ANSOI,   ANLV,   ANST,   ANRT',             &
     & ',   ANSO, NLOSSL, NLOSSR, NLOSSS, NBALAN, NNI')
      return

      case (2)
! --- write dynamic data ----------------------------------------------------
      write (nba,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,NUPTT,comma,NFIXTT,comma,ANLVI,comma,ANSTI,comma,     &
     & ANRTI,comma,ANSOI,comma,ANLV,comma,ANST,comma,ANRT,comma,        &
     & ANSO,comma,NLOSSL,comma,NLOSSR,comma,NLOSSS,comma,NBALAN,        &
     & comma,NNI
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 15(a1,f7.2) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (nba)

      case default
         call fatalerr ('OutbalCropN', 'Illegal value for TASK')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine outbalcropOM1(task,pathwork,outfil,project,date,       &
     &      daycrop,t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output OM balance files 
! ----------------------------------------------------------------------
      implicit none
      include 'params.fi'

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck
! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   getun,om1
      character(len=160) filnam,filtext
      real(8)   OMroot,OMleaves,OMstems,OMstorage  !,Cccheck

      save    om1

      comma = ',' 

      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.om1'
      om1 = getun (20,90)
      call fopens(om1,filnam,'new','del')
      filtext = 'output of OM1-balance (kg/ha DM increase per time'//   &
     & 'step) of detailed crop growth model'
      call writehead (om1,1,filnam,filtext,project)

      write (om1,100)
 100    format ('*',/,                                                  &
     & '*             day     day      -   grCd   kg/ha   kg/ha',       &
     & '   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha       -   kg/ha',/,    &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,   gass,   mres,OMroots'&
     & ,',OMleaves,OMstems,OMstorage, dmi,    cvf,OMcheck')
      return

      case (2)

! --- write dynamic data ----------------------------------------------------
      if(cvf.lt.nihil) then
        OMroot    = 0.0d0 
        OMleaves  = 0.0d0
        OMstems   = 0.0d0
        OMstorage = 0.0d0
      else
        OMroot    = fr*dmi/cvf
        OMleaves  = fl*(1.0d0-fr)*dmi/cvf
        OMstems   = fs*(1.0d0-fr)*dmi/cvf
        OMstorage = fo*(1.0d0-fr)*dmi/cvf
      endif  

      write (om1,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,gass,comma,mres,comma,OMroot,comma,OMleaves,comma,    &
     & OMstems,comma,OMstorage,comma,dmi,comma,cvf,comma,ccheck
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 9(a1,f7.2) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (om1)

      case default
         call fatalerr ('OutbalCropOM1', 'Illegal value for TASK')
      end select

      return
      end 

      subroutine outbalcropom2(task,pathwork,outfil,project,date,       &
     &         daycrop,t,dvs,tsum,storagediff,wlv,wst,wso,wrt,      &
     &         delt,grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output OM balance files 
! ----------------------------------------------------------------------
      implicit none

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum,storagediff,wlv,wst,wso,wrt,delt
      real(8) grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan

! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   getun,om2
      character(len=160) filnam,filtext

      save    om2

      comma = ',' 

      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.om2'
      om2 = getun (20,90)
      call fopens(om2,filnam,'new','del')
      filtext = 'output of OM2-balance (kg/ha DM, cumulative and '//    &
     & 'increments) of detailed crop growth model'
      call writehead (om2,1,filnam,filtext,project)

      write (om2,100)
 100    format ('*',/,                                                  &
     & '*             day     day      - degday   kg/ha   kg/ha ',      &
     & '  kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha ',      &
     & '  kg/ha kg/ha/d kg/ha/d kg/ha/d kg/ha/d',/,                     &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,storagediff,wlv,',      &
     & '    wst,    wso,    wrt,   grlv,   grst,   grso,   grrt,',      &
     & 'ombalan,   drlv,   drst,   drso,   drrt')

      return

      case (2)

! --- write dynamic data ----------------------------------------------------
      write (om2,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,storagediff,comma,wlv,comma,wst,comma,wso,comma,      &
     & wrt,comma,grlv*delt,comma,grst*delt,comma,grso*delt,comma,       &
     & grrt*delt,comma,drlv*delt,comma,drst*delt,comma,drso*delt,       &
     & comma,drrt*delt,comma,ombalan
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 14(a1,f7.1) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (om2)

      case default
         call fatalerr ('OutbalCropOM2', 'Illegal value for TASK')
      end select

      return
      end 


      SUBROUTINE CHCKPRT(DVS,FR,FL,FS,FO)       
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Checks the partitioning factors, and interrupt in case of error
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    R8  FR        Fraction of total dry matter partitioned to the roots (-)
!       I    R8  FL        Fraction of total dry matter partitioned to the leaves (-)
!       I    R8  FS        Fraction of total dry matter partitioned to the stems (-)
!       I    R8  FO        Fraction of total dry matter partitioned to the storage organs (-)
! ----------------------------------------------------------------------
      implicit none
! --- global
      real(8)   DVS,FR,FL,FS,FO
! --- local
      character(len=300) messag
      real(8)   FCHECK

!*     check on partitioning
      fcheck = fr+(fl+fs+fo)*(1.0d0-fr) - 1.0d0
      if (abs (fcheck).gt.0.0001d0) then
!        write (messag,'(a,f5.2,/,3(a,g12.5),/,2(a,g12.5))')             &
        write (messag,'(a,f5.2,3(a,g12.5),2(a,g12.5))')             &
     &      ' error in partitioning functions, dvs= ',dvs,             &
     &      ' fcheck = ',fcheck,' fr = ',fr,' fl = ',fl,                &
     &      ' fs = ',fs,' fo = ',fo
        call fatalerr ('wofost',messag)
      end if
      
      return
      end
 
      SUBROUTINE CHCKCBL(DVS,CVF,DMI,FR,FL,FS,FO,GASS,MRES)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Checks the carbon balance, and interrupt in case of error
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    R8  DVS       development stage (-)
!       I    R8  CVF        (-)
!       I    R8  DMI       dry matter increase (-)
!       I    R8  FR        Fraction of total dry matter partitioned to the roots (-)
!       I    R8  FL        Fraction of total dry matter partitioned to the leaves (-)
!       I    R8  FS        Fraction of total dry matter partitioned to the stems (-)
!       I    R8  FO        Fraction of total dry matter partitioned to the storage organs (-)
!       I    R8  GASS      Gross assimilation (-)
!       I    R8  MRES      maintenance respiration (-)
! ----------------------------------------------------------------------
      implicit none
! --- global
      real(8)   DVS,CVF,DMI,FR,FL,FS,FO,GASS,MRES
! --- local
      character(len=300) messag
      real(8)   CCHECK

!*     check on C-balance
      CCHECK = (GASS-MRES-(FR+(FL+FS+FO)*(1.0d0-FR))*DMI/CVF)           &
     &       /MAX (0.0001d0,GASS)
      IF (ABS (CCHECK).GT.0.0001d0) THEN
        WRITE (messag,'(A,I3,/,3(A,G12.5),/,A,4G12.5,/,2(A,G12.5))')    &
     &     ' Carbon flows nog balanced on day ',DVS,                    &
     &     ' CCHECK = ',CCHECK,' GASS = ',GASS,' MRES = ',MRES,         &
     &     ' FR,L,S,O = ',FR,FL,FS,FO,' DMI = ',DMI,' DVF = ',CVF
        call fatalerr ('wofost',messag)
      END IF    
      RETURN
      END


      SUBROUTINE RELGRWT(DMI,FR,FL,FS,FO,GRRT,GRLV,GRST,GRSO)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : To calculate relative growth rate of roots, stems leaves
!                          and storage organs  
! Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DSLV     Death rate leaves                              kg ha-1 d-1 DM
!   I    R8  DMI      total dry matter increase                      - 
!   I    R8  FR       Fraction of total dry matter partitioned to the roots (-)
!   I    R8  FL       Fraction of total dry matter partitioned to the leaves (-)
!   I    R8  FS       Fraction of total dry matter partitioned to the stems (-)
!   I    R8  FO       Fraction of total dry matter partitioned to the storage organs (-)
!   O    R8  GRRT     Growth rate roots                              kg ha-1 d-1 DM
!   O    R8  GRLV     Growth rate leaves                             kg ha-1 d-1 DM
!   O    R8  GRST     Growth rate stems                              kg ha-1 d-1 DM
!   O    R8  GRSO     Growth rate storage organs                     kg ha-1 d-1 DM
! ===== ==== =======  =============================================  ==============
!   O    R8  ADMI      above ground dry matter increase (-)
      implicit none
! --- global
      real(8)   DMI,FR,FL,FS,FO,GRRT,GRLV,GRST,GRSO
! --- local
      real(8)   ADMI
!     save
            
      ADMI = (1.0d0-FR)*DMI
      GRRT = FR*DMI
      GRLV = FL*ADMI
      GRST = FS*ADMI
      GRSO = FO*ADMI
      
      RETURN
      END

      subroutine deaths(flcropnut,wlv,kdif,lai,NNI,perdl,rdrns,         &
     &                  reltr,dslv)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Compute the relative death rate leaves due
!                          to stress (kg DM ha-1 d-1)   
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    L   flCropNut Flag indicating simulation of nutrient stress -
!   I    R8  WLV       Dry weight of living leaves                   kg ha-1 DM
!   O    R8  reltr
!   I    R8  KDIF      Extinction coeff. for diffuse visible light   -
!   I    R8  LAI       Leaf Area Index
!   I    R8  NNI  
!   I    R8  PERDL     Max.rel.death rate of leaves due to water strs -
!   I    R8  RDRNS     Max.rel.death rate of leaves due to nitrogen strs -
!   O    R8  DSLV      Death rate leaves                             kg ha-1 d-1 DM
! ===== ==== =======  =============================================  ==============
!   -    R8  DSLV1     Death rate leaves due to water stress         kg ha-1 d-1 DM
!   -    R8  DSLV1     Death rate leaves due to self-shading         kg ha-1 d-1 DM
!   -    R8  LAICR     
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      logical   flCropNut
      real(8)   DSLV,KDIF,LAI,NNI,PERDL,RDRNS,reltr,WLV
! --- local
      real(8)   DSLV1,DSLV2,LAICR
!      save     
      
!     death rate of leaves due to water stress
      DSLV1 = WLV*(1.d0-reltr)*PERDL
      
!     death rate of leaves due high LAI
      LAICR = 3.2d0/KDIF
      DSLV2 = WLV*max(0.0d0, min(0.03d0, (0.03d0*(LAI-LAICR)/LAICR)))
      DSLV  = MAX (DSLV1, DSLV2)
      
!     death rate increase due to nutrient shortage
      IF(flCropNut .AND. NNI.LT.1.0d0) THEN
         DSLV = DSLV + WLV*RDRNS * (1.0d0-NNI)
      END IF 
      
      RETURN
      END
      
      subroutine deatha(dslv,delt,ilvold,lv,lvage,span,i1,dalv)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : To compute the relative death rate leaves due              * 
!                          to ageing (kg DM ha-1 d-1)   

!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT      
!   I    R8  REST   
!   I    I   i1vold      
!   I    R8  SPAN  
!   I    R8  LVAGE
!   I/O  I   11      
!   O    R8  DALV      Death rate leaves due to ageing (kg/ha/d DM)
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   i1,ilvold
      real(8)   DSLV,DELT,SPAN,DALV,lv(366),lvage(366)
! --- local
      real(8)   REST
            
! --- first: leaf death due to water stress or high lai is imposed on array
! ---        until no more leaves have to die or all leaves are gone

      rest = dslv*delt
      i1 = ilvold

      do while (rest.gt.lv(max(i1,1)).and.i1.ge.1)
        rest = rest-lv(i1) 
        i1 = i1-1
      enddo

! --- then: check if some of the remaining leaves are older than span,
! ---       sum their weights

      dalv = 0.0d0
      if (lvage(max(i1,1)).ge.span .and. rest.gt.0.0d0 .and.i1.ge.1)then
        dalv = lv(i1)-rest
        rest = 0.0d0
        i1 = i1-1
      endif

      do while (i1.ge.1.and.lvage(max(i1,1)).gt.span)
        dalv = dalv+lv(i1)
        i1 = i1-1
      enddo

      dalv = dalv/delt
      
      RETURN
      END      

      SUBROUTINE GLAI(Fstress,LAIEXP,TEMP,TBASE,RGRLAI,GRLV,SLAT,GLA)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : exponential, sink limited leave increase
!                          and adjustment of specific leaf area of youngest 
!                          leaf class. Adjust for water and nitrogen stress
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  Fstress
!   I    R8  LAIEXP
!   O    R8  GLA       
!   I/O  R8  SLAT     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      real(8)   Fstress,LAIEXP,TEMP,TBASE,RGRLAI,GRLV,SLAT,GLA
! --- local
      real(8)   DTEFF,GLAIEX,GLASOL
!      save     

      IF (LAIEXP.LT.6.0d0) THEN
         DTEFF  = MAX (0.d0,TEMP-TBASE)
         GLAIEX = Fstress * LAIEXP*RGRLAI*DTEFF
!*        source-limited increase in leaf area
         GLASOL = GRLV*SLAT
!*        sink-limited increase in leaf area
         GLA    = MIN (GLAIEX, GLASOL)
!*        adjustment of specific leaf area of youngest leaf class
         IF (GRLV.GT.0.d0) SLAT = GLA/GRLV
      END IF
      
      RETURN
      END
      SUBROUTINE LVDTH(DELT,DSLV,SPAN,ILVOLD,LVAGE,LV,I1)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Impose leave death on LV array 
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT
!   I    R8  DSLV       
!   I    I   SPAN
!   I    I   ILVOLD
!   I    R8  LVAGE    
!   I/O  R8  LV     
!   O    I   I1
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold,i1
      real(8)   delt,dslv,lv(366),lvage(366),span
! --- local
      real(8)   dslvt

      save     

! --- remaining leaves
      dslvt = dslv*delt
      i1 = ilvold
      do while ((dslvt.gt.0.0d0) .and. (i1.ge.1))
        if (dslvt.ge.lv(i1)) then
          dslvt = dslvt-lv(i1)
          lv(i1) = 0.0d0
          i1 = i1-1
        else
          lv(i1) = lv(i1)-dslvt
          dslvt = 0.0d0
        endif
      enddo

! --- leaves older than span die
      do while (lvage(max(i1,1)).ge.span.and.i1.ge.1)
        lv(i1) = 0.0d0
        i1 = i1-1
      enddo

      RETURN
      END

      SUBROUTINE LVSHFT(DELT,FYSDEL,ILVOLD,grlv,slat,LV,LVAGE,SLA)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Shift contents of LV, LVAGE, SLA tables with one day  
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT
!   I    R8  FYSDEL       
!   I    I   ILVOLD
!   I/O  R8  LV     
!   I/O  R8  LVAGE     
!   I/O  R8  SLA     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold
      real(8)   delt,fysdel,grlv,slat,lv(366),lvage(366),sla(366)
! --- local
      integer   i1

!      save     

!     Shift contents of LV, LVAGE, SLA tables with one day
      do i1 = ilvold,1,-1
        lv(i1+1) = lv(i1)
        sla(i1+1) = sla(i1)
        lvage(i1+1) = lvage(i1)+fysdel*delt
      enddo

!     new leaves in class 1
      lv(1) = grlv*delt
      sla(1) = slat
      lvage(1) = 0.0d0 
     
      return
      end
      SUBROUTINE LVWGLI(ILVOLD,LV,SLA,LASUM,WLV)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Calculate new wlv and lai 
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    I   ILVOLD
!   I    R8  LV     
!   I    R8  SLA     
!   O    R8  LASUM     
!   O    R8  WLV     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold
      real(8)   lv(366),sla(366),lasum,wlv
! --- local
      integer    i1
!      save     

! --- calculation of new leaf area and weight
      lasum = 0.0d0
      wlv = 0.0d0
      do i1 = 1,ilvold
        lasum = lasum+lv(i1)*sla(i1)
        wlv = wlv+lv(i1)
      enddo
  
      return
      end

! --- for soybean (swsoybean=1): temperature and photoperiodicity
      subroutine mgtemprf(tav,toptdvr,tmindvr,tmaxdvr,rfmgtemp)
! ----------------------------------------------------------------------
!     Last modified      : Sept 2015
!       based on routines needed from pyWofost for soybean (Allard de Wit, 2015)
!
!     Purpose            : temperature reduction factor for soybean (short day)
!       approach and parameters based on Setiyono et al. doi 10.1016/j.fcr.2006.07.011
!       http://digitalcommons.unl.edu/agronomyfacpub/112
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  tav        
!   I    R8  toptdvr       
!   I    R8  tmindvr
!   I    R8  tmaxdvr     
!   O    R8  rfmgtemp    
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      real(8)   tav, toptdvr, tmindvr, tmaxdvr, rfmgtemp
! --- local
      real(8)   alpha, p1, p2, p3, p4

      alpha = log(2.0d0)/(log((tmaxdvr-tmindvr)/(toptdvr-tmindvr)))
      if(tav.lt.tmindvr .or. tav.gt.tmaxdvr) then
        rfmgtemp = 0.0d0
      else
        p1 = 2.0d0 * (tav - tmindvr)**alpha
        p2 = (toptdvr - tmindvr)**alpha
        p3 = (tav - tmindvr)**(2.0d0*alpha)
        p4 = (toptdvr - tmindvr)**(2.0d0*alpha)
        rfmgtemp = (p1 * p2 - p3) / p4

      endif

      return
      end

      subroutine mgphotoprf(mg,iday,lat,popt,pcrt,flphenodayl,rfmgphotop)
! ----------------------------------------------------------------------
!     Last modified      : Sept 2015
!       based on routines needed from pyWofost for soybean (Allard de Wit, 2015)
!
!     Purpose            : Photoperiod reduction factor for soybean (short day)
!       approach and parameters based on Setiyono et al. doi 10.1016/j.fcr.2006.07.011
!       http://digitalcommons.unl.edu/agronomyfacpub/112
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  mg       maturity group 
!   I    I   iday     daynr      
!   I    R8  lat      lattitude
!   I    L   flphenodayl Flag to allow input of POPT and PCRT or using 
!                        empirical relation from Setiyono et al
!   I    R8  popt     optimal daylength for phenological developm.   hr
!   I    R8  pcrt     critical daylength for phenological developm.  hr
!   O    R8  rfmgphotop  reduction factor for photoperiodicity       -  
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   iday
      logical   flphenodayl
      real(8)   lat, popt, pcrt
      real(8)   mg, rfmgphotop
! --- local
      real(8)   alpha,m,p0,p1,p2
      real(8)   dec,daylp,pi,rad,sinld,cosld,aob
      parameter (pi=3.1415926d0, rad=0.0174533d0)

! astronomic daylength according to solar elevation angle of -0.833 day
!*     Declination and solar constant for this day
      dec = -asin(sin(23.45d0*rad)*cos(2.d0*pi*dble(iday+10)/365.0d0))
      sinld = sin(rad*lat)*sin(dec)
      cosld = cos(rad*lat)*cos(dec)
      aob = sinld/cosld
      if (abs(aob).le.1.0d0) then 
         daylp = 12.0d0*(1.d0+2.d0*asin(aob)/pi)
      else
         if (aob.gt.1.0d0)  daylp = 24.0d0
         if (aob.lt.-1.0d0) daylp =  0.0d0
      endif

! First determine Popt and Pcrt based on maturity group rating
      m = 3.0d0
      if(flphenodayl) then
        continue            ! popt and pcrt  are input
      else
        popt = 12.759d0 - 0.388d0*mg - 0.058d0*mg**2
        pcrt = 27.275d0 - 0.493d0*mg - 0.066d0*mg**2
      endif
      alpha = log(2.0d0)/log(((pcrt - popt)/m) + 1.0d0)
      p0 = (pcrt - popt)/m

      if(daylp.lt.popt) then
        rfmgphotop = 1.0d0
      else
        if (daylp.gt.pcrt) then
          rfmgphotop = 0.0d0
        else
          p1 = (daylp - popt)/m + 1.0d0
          p2 = (pcrt - daylp)/(pcrt - popt)
          rfmgphotop = (p1*(p2**p0))**alpha
        endif
      endif

!      write(99,*)  iday, daylp, rfmgphotop
      
      return
      end

      subroutine sumttd(task,flGrassGrowth,dateGrassGrowth)
! ----------------------------------------------------------------------
!     Last modified      : Jan 2016
!     Author             : Joop Kroes
!
!     Purpose            : Suppress grass growth as long as 3 criteria are not met: 
!                          temperature, time and depth
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    C   condition case: 'initial' or 'dynamic' (-)
!       I    R8  tsumtemp  temperature limit to initiate grass growth  [0.0..20.0 grC, R]
!       I    I   tsumtime  time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
!       I    R8  tsumdepth depth at which temp above tsumtemp for grass growth [0.0..100.0 cm below soil surface, R]
!       I    R8  z         depth of a node (L)
!       I    R8  tsoil     Array with soil temperatures (�C) for each compartment
!       O    L   flGrassGrowth flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]
! ----------------------------------------------------------------------
      use Variables
      implicit none

! --- global
      character(len=*), intent(in) :: task
      logical, intent(out)         :: flGrassGrowth      ! flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]

! --- local
      integer    :: tsumtimecum      ! cumulative, from 1-jan, time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
      logical    :: fltsumtemp       ! flag to indicate if temperature criteria is met 
      logical    :: fltsimprev       ! flag to indicate if temperature criterion is met during simulated previous day
      logical    :: fltsimcount      ! flag to indicate nr of contiuous simulated days that temperature criteria is met
      integer    :: cmpcrit, node
      character(len=11), intent(out) ::  dateGrassGrowth            ! date of start of GrassGrowth
!     output
      integer    :: uo, getun   !, idum, ios
      character(len=160)  :: filnam, filtext
      character(len=1)  :: comma
      logical    :: flexist, flopened
      
      save

      comma = ',' 

      
      select case(task)

      case('initial')
          tsumtimecum = 0
          fltsumtemp = .false.
          fltsimprev = .false.
          ! find depth of compartment with critical soil temperature
          cmpcrit = 1
          do node = 1, numnod
             if (z(node) .le. (-1.0d0*tsumdepth)) then
                cmpcrit = node
                exit
             endif
          enddo
          ! no growth as long as 3 criteria are not met
          flGrassGrowth = .false.
          dateGrassGrowth = 'undefined'
          
          ! === open output file and write headers
          filnam = trim(pathwork)//trim(outfil)//'.ttd'
          flopened = .false.
          if(uo.gt.0) then
              !  Inquiry by Unit
              inquire (uo, opened=flopened, exist=flexist)
          endif
          if (.not.flexist .and. .not.flopened) then
              uo = getun (20,90)
              call fopens(uo,filnam,'new','del')
              filtext = 'output of subr sumttd'
              call writehead (uo,1,filnam,filtext,project)
              write (uo,100)
 100          format (' Date,z(cmpcrit),tsoil(cmpcrit),',               &
     &           'fltsumtemp,fltsimprev,fltsimcount')
          endif

      case('dynamic')
          ! temperature and depth criterium
          if(tsoil(cmpcrit).ge.tsumtemp) then
              fltsumtemp = .true.
          else
              fltsumtemp = .false.
          endif
          ! timing criterium: set flag for continuous days that exceed critical temperature
          if(fltsumtemp .and. fltsimprev) then
              tsumtimecum = tsumtimecum + 1
          else
              tsumtimecum = 0
              fltsimprev = .false.
          endif
          !  timing criterium: count nr of days exceeding critical temperature
          if(tsumtimecum.ge.tsumtime) then
              fltsimcount = .true.
          else
              fltsimcount = .false.
          endif
          ! no growth as long as 3 criteria are not met
          flGrassGrowth = .false.
          if(fltsumtemp .and. fltsimprev .and. fltsimcount) then
              call dtdpst ('year-month-day',t1900,dateGrassGrowth)
              flGrassGrowth = .true.
          endif

          ! timing criteria for next timestep
          if(fltsumtemp) then
              fltsimprev = .true.
          endif
          
          ! === write output 
          write (uo,200) Date,comma,z(cmpcrit),comma,tsoil(cmpcrit),    &
     &           comma,fltsumtemp,comma,fltsimprev,comma,fltsimcount
 200      format (a11,2(a1,f7.2),3(a1,i3))

      end select

      return
      end
