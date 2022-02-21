! File VersionID:
!   $Id: swapoutput.f90 330 2017-06-12 10:22:55Z kroes006 $
! ----------------------------------------------------------------------
      subroutine swapoutput(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : open and write general swap output files 
! ----------------------------------------------------------------------

      use Variables
      implicit none

      integer task

      select case (task)
      case (1)

! === open output files ===============================

! --  bal file
      call outbal (task,bal,logf,swscre,swdra,numnod,nrlevs,swsolu,     &
     &    ioutdat,cevap,cgird,cgrai,cqbot,tstart,cqrot,crunoff,crunon,  &
     &    cQMpOutDrRap,cqdra,cqdrain,dz,z,samini,sampro,samcra,sqprec,  &
     &    sqirrig,sqbot,dectot,rottot,sqrap,sqdra,pond,volact,volini,   &
     &    t1900,outdat,outfil,pathwork,project,                         &
     &    caintc,csubl,PondIni,WaSrDm1,WaSrDm2,WaSrDm1Ini,WaSrDm2Ini,   &
     &    swsnow,cgsnow,csnrai,snowinco,ssnow)

! --  blc file
      if (swblc .eq. 1)                                                 &
     &  call outblc(task,blc,outfil,pathwork,cgrai,cnrai,cgird,         &
     &     cnird,cqrot,cevap,volact,volini,z,dz,swmacro,nrlevs,swirfix, &
     &     schedule,numnod,snowinco,ssnow,cgsnow,cmelt,caintc,csnrai,   &
     &     cqprai,ioutdat,t1900,outdat,tstart,project,pond,pondini,     &
     &     cqdrainin,cqdrainout,cinund,crunoff,cqtdo,cqtup,cqbotdo,     &
     &     cqbotup,csubl,crunon,                                        &
     &     (CQMpInIntSatDm1+CQMpInIntSatDm2+CQMpInMtxSatDm1+            &
     &     CQMpInMtxSatDm2),(CQMpOutMtxSatDm1+CQMpOutMtxSatDm2+         &
     &     CQMpOutMtxUnsDm1+CQMpOutMtxUnsDm2),CQMpInTopPreDm1,          &
     &     CQMpInTopPreDm2,CQMpInTopLatDm1,CQMpInTopLatDm2)

      return

      case (2)

! === write actual data ===============================

      if (flbaloutput) then

! --    bal file
        call outbal (task,bal,logf,swscre,swdra,numnod,nrlevs,swsolu,   &
     &    ioutdat,cevap,cgird,cgrai,cqbot,tstart,cqrot,crunoff,crunon,  &
     &    cQMpOutDrRap,cqdra,cqdrain,dz,z,samini,sampro,samcra,sqprec,  &
     &    sqirrig,sqbot,dectot,rottot,sqrap,sqdra,pond,volact,volini,   &
     &    t1900,outdat,outfil,pathwork,project,                         &
     &    caintc,csubl,PondIni,WaSrDm1,WaSrDm2,WaSrDm1Ini,WaSrDm2Ini,   &
     &    swsnow,cgsnow,csnrai,snowinco,ssnow)

! --  blc file
        if (swblc .eq. 1)                                               &
     &  call outblc(task,blc,outfil,pathwork,cgrai,cnrai,cgird,         &
     &     cnird,cqrot,cevap,volact,volini,z,dz,swmacro,nrlevs,swirfix, &
     &     schedule,numnod,snowinco,ssnow,cgsnow,cmelt,caintc,csnrai,   &
     &     cqprai,ioutdat,t1900,outdat,tstart,project,pond,pondini,     &
     &     cqdrainin,cqdrainout,cinund,crunoff,cqtdo,cqtup,cqbotdo,     &
     &     cqbotup,csubl,crunon,                                        &
     &     (CQMpInIntSatDm1+CQMpInIntSatDm2+CQMpInMtxSatDm1+            &
     &     CQMpInMtxSatDm2),(CQMpOutMtxSatDm1+CQMpOutMtxSatDm2+         &
     &     CQMpOutMtxUnsDm1+CQMpOutMtxUnsDm2),CQMpInTopPreDm1,          &
     &     CQMpInTopPreDm2,CQMpInTopLatDm1,CQMpInTopLatDm2)

      endif

      return

      case (3)

! === close output files ================================

      close (bal)
      if (swblc .eq. 1) close (blc)

! --- final message log file
      write(logf,'(/,a)') ' Swap simulation okay!'
      close (logf) 

      case default
         call fatalerr ('SWAPoutput', 'Illegal value for Task')
      end select
      
      return
      end 

! ----------------------------------------------------------------------
      subroutine soilwateroutput(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : open and write soil water output files 
! ----------------------------------------------------------------------

      use Variables
      implicit none

      integer task

      select case (task)

      case (1)
! === open output files and write headers ===============================

! --  wba file
      call outwba (1,wba,daynr,daycum,swscre,cevap,cgird,               &
     &    cgrai,csnrai,cnird,cnrai,cpeva,cptra,cqbot,cqdra,cqrot,crunon,&
     &    crunoff,gwl,cQMpOutDrRap,pond,t1900,date,volact,volini,       &
     &    wbalance,outfil,pathwork,project,flprintshort,floutput,       &
     &    swsnow,cqprai,ssnow,snowinco,PondIni,flheader)

! --  inc file
      call outinc (1,inc,daynr,daycum,igrai,isnrai,igsnow,              &
     &   igird,iintc,irunon,iruno,iptra,iqrot,ipeva,ievap,isubl,iqdra,  &
     &   iQMpOutDrRap,iqbot,t1900,date,outfil,pathwork,project,flheader,&
     &   gwl,volact,volini,pond,PondIni,ssnow,snowinco,flprintshort)

! --  str file
      call outstr (1,str,daynr,daycum,ies0,iet0,iew0,ipeva,             &
     &  iptra,iqrot,iqredwet,iqreddry,iqredsol,iqredfrs,t1900,          &
     &  date,outfil,pathwork,project,flheader,flprintshort)

! --  vap file
      if (swvap.eq.1) call outvap (1,vap,daynr,numnod,daycum,z,         &
     &    cml,t1900,theta,h,k,tsoil,q,outfil,pathwork,project,          &
     &    swheader,isqtop,isqbot,qdraincomp,qrot,cmsy,dz,date,          &
     &    flprintshort)

! --  rot file, only when drought stress according to De Jong van Lier
     if (swdrought.eq.2) then
          call outrot(1,rot,daynr,noddrz,daycum,z,hxylem,t1900,         &
     &    theta,hm1,q,outfil,pathwork,project,swheader,qrot,            &
     &    dz,date,flprintshort,hroot,inq,inqrot,mroot,mflux,            &
     &    rootrho,rootphi,hleaf)
      endif
          
! --  capillary rise output file
      if(swcapriseoutput) call capriseoutput(task)

! --  extensive formatted output file for solute studies
      if (swafo.ge.1) then
        call outafo (task,afo,outfil,pathwork,numnod,outper,period,     &
     &  pond,gwl,theta,h,inq,inqrot,nrlevs,inqdra,ievap,ipeva,iptra,    &
     &  iruno,numlay,botcom,dz,thetas,cofgen,kdif,kdir,swafo,igrai,     &
     &  igird,inrai,inird,iintc,gc,lai,tav,tsoil,daycum,                &
     &  rd,cf,wbalance,project,tstart,tend,swdiscrvert,numnodnew,dznew, &
     &  SwMacro,Ssnow,igSnow,isnrai,iSubl,irunon, FrArMtrx,IPondBeg,    &
     &  ISsnowBeg,IThetaBeg,VlMpStDm1,VlMpStDm2,DiPoCp,WalevDm1,VlMpDm1,&
     &  WaSrDm1,VlMpDm2,WaSrDm2,IQInTopPreDm1,IQInTopLatDm1,            &
     &  IQExcMtxDm1Cp,IQOutDrRapCp,IAvFrMpWlWtDm1,IQInTopPreDm2,        &
     &  IQInTopLatDm2,IQExcMtxDm2Cp,IAvFrMpWlWtDm2,IWaSrDm1Beg,         &
     &  IWaSrDm2Beg,CritDevMasBal,tcum,nod1lay,swsophy,numtab,sptab,    &
     &  ientrytab,                                                      &
     &  out_etr,out_hum,out_rad,out_tmn,out_tmx,out_wet,out_win)
      endif

! --  generate file with soil physical parameters
      call outsoilphys ()

! --  extensive unformatted output file for solute studies
      if (swaun.ge.1) then
        call outaun (task,aun,outfil,pathwork,numnod,outper,period,     &
     &  pond,gwl,theta,h,inq,inqrot,nrlevs,inqdra,ievap,ipeva,iptra,    &
     &  iruno,numlay,botcom,dz,thetas,cofgen,kdif,kdir,swaun,igrai,     &
     &  igird,inrai,inird,iintc,gc,lai,tav,tsoil,daycum,                &
     &  rd,cf,wbalance,project,tstart,tend,swdiscrvert,numnodnew,dznew, &
     &  SwAfo,SwMacro,Ssnow,igSnow,isnrai,iSubl,irunon,FrArMtrx,        &
     &  IPondBeg,ISsnowBeg,IThetaBeg,VlMpStDm1,VlMpStDm2,DiPoCp,        &
     &  WalevDm1,VlMpDm1,WaSrDm1,VlMpDm2,WaSrDm2,IQInTopPreDm1,         &
     &  IQInTopLatDm1,IQExcMtxDm1Cp,IQOutDrRapCp,IAvFrMpWlWtDm1,        &
     &  IQInTopPreDm2,IQInTopLatDm2,IQExcMtxDm2Cp,IAvFrMpWlWtDm2,       &
     &  IWaSrDm1Beg,IWaSrDm2Beg,CritDevMasBal,tcum,nod1lay,             &
     &  swsophy,numtab,sptab,ientrytab)
      endif

! --  output of recharge/storage for modflow
      if(swoutputmodflow.eq.1) call OutputModflow(1)

      return
      
      case (2)
! === write actual data ===============================

! --  wba file
      call outwba (2,wba,daynr,daycum,swscre,cevap,cgird,               &
     &    cgrai,csnrai,cnird,cnrai,cpeva,cptra,cqbot,cqdra,cqrot,crunon,&
     &    crunoff,gwl,cQMpOutDrRap,pond,t1900,date,volact,volini,       &
     &    wbalance,outfil,pathwork,project,flprintshort,floutput,       &
     &    swsnow,cqprai,ssnow,snowinco,PondIni,flheader)

! --  inc file
      call outinc (2,inc,daynr,daycum,igrai,isnrai,igsnow,              &
     &   igird,iintc,irunon,iruno,iptra,iqrot,ipeva,ievap,isubl,iqdra,  &
     &   iQMpOutDrRap,iqbot,t1900,date,outfil,pathwork,project,flheader,&
     &   gwl,volact,volini,pond,PondIni,ssnow,snowinco,flprintshort)

! --  str file
      if (flCropCalendar) then
        if (.not. flCropHarvest) then
          call outstr (2,str,daynr,daycum,ies0,iet0,iew0,ipeva,         &
     &    iptra,iqrot,iqredwet,iqreddry,iqredsol,iqredfrs,t1900,        &
     &    date,outfil,pathwork,project,flheader,flprintshort)
        endif
      endif

! --  vap file
      if (swvap.eq.1) call outvap (2,vap,daynr,numnod,daycum,z,         &
     &    cml,t1900,theta,h,k,tsoil,q,outfil,pathwork,project,          &
     &    swheader,isqtop,isqbot,qdraincomp,qrot,cmsy,dz,date,          &
     &    flprintshort)

! --  rot file
      if (swdrought.eq.2 .and. ptra .gt. 1.0d-10)                       & 
     &    call outrot (2,rot,daynr,noddrz,daycum,z,                     &
     &    hxylem,t1900,theta,hm1,q,outfil,pathwork,project,swheader,    &
     &    qrot,dz,date,flprintshort,hroot,inq,inqrot,mroot,mflux,       &
     &    rootrho,rootphi,hleaf)

! --  capillary rise output file
      if(swcapriseoutput) call capriseoutput(task)

! --  extensive formatted output file for solute studies
      if (swafo.ge.1) then
        call outafo (task,afo,outfil,pathwork,numnod,outper,period,     &
     &  pond,gwl,theta,h,inq,inqrot,nrlevs,inqdra,ievap,ipeva,iptra,    &
     &  iruno,numlay,botcom,dz,thetas,cofgen,kdif,kdir,swafo,igrai,     &
     &  igird,inrai,inird,iintc,gc,lai,tav,tsoil,daycum,                &
     &  rd,cf,wbalance,project,tstart,tend,swdiscrvert,numnodnew,dznew, &
     &  SwMacro,Ssnow,igSnow,isnrai,iSubl,irunon, FrArMtrx,IPondBeg,    &
     &  ISsnowBeg,IThetaBeg,VlMpStDm1,VlMpStDm2,DiPoCp,WalevDm1,VlMpDm1,&
     &  WaSrDm1,VlMpDm2,WaSrDm2,IQInTopPreDm1,IQInTopLatDm1,            &
     &  IQExcMtxDm1Cp,IQOutDrRapCp,IAvFrMpWlWtDm1,IQInTopPreDm2,        &
     &  IQInTopLatDm2,IQExcMtxDm2Cp,IAvFrMpWlWtDm2,IWaSrDm1Beg,         &
     &  IWaSrDm2Beg,CritDevMasBal,tcum,nod1lay,swsophy,numtab,sptab,    &
     &  ientrytab,                                                      &
     &  out_etr,out_hum,out_rad,out_tmn,out_tmx,out_wet,out_win)
      endif

! --  extensive unformatted output file for solute studies
      if (swaun.ge.1) then
        call outaun (task,aun,outfil,pathwork,numnod,outper,period,     &
     &  pond,gwl,theta,h,inq,inqrot,nrlevs,inqdra,ievap,ipeva,iptra,    &
     &  iruno,numlay,botcom,dz,thetas,cofgen,kdif,kdir,swaun,igrai,     &
     &  igird,inrai,inird,iintc,gc,lai,tav,tsoil,daycum,                &
     &  rd,cf,wbalance,project,tstart,tend,swdiscrvert,numnodnew,dznew, &
     &  SwAfo,SwMacro,Ssnow,igSnow,isnrai,iSubl,irunon,FrArMtrx,        &
     &  IPondBeg,ISsnowBeg,IThetaBeg,VlMpStDm1,VlMpStDm2,DiPoCp,        &
     &  WalevDm1,VlMpDm1,WaSrDm1,VlMpDm2,WaSrDm2,IQInTopPreDm1,         &
     &  IQInTopLatDm1,IQExcMtxDm1Cp,IQOutDrRapCp,IAvFrMpWlWtDm1,        &
     &  IQInTopPreDm2,IQInTopLatDm2,IQExcMtxDm2Cp,IAvFrMpWlWtDm2,       &
     &  IWaSrDm1Beg,IWaSrDm2Beg,CritDevMasBal,tcum,nod1lay,             &
     &  swsophy,numtab,sptab,ientrytab)
      endif

! --  output of recharge/storage for modflow
      if(swoutputmodflow.eq.1) call OutputModflow(2)

      return

      case (3)
! === write final values and close output files ===========================

! --- write final pressure heads, solute concentrations and soil temperatures
      call outend (                                                     &
     &          numnod,h,flSolute,flAgeTracer,cml,z,fltemperature,tsoil,&
     &          ssnow,pond,dt,icrop,croptype,cropname,flSurfaceWater,   &
     &          wls,swredu,ldwet,spev,outfil,pathwork,project,          &
     &  rd,rdpot,dvs,flanthesis,tsum,ilvold,ilvoldpot,wrt,wrtpot,tadw,  &
     &  tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,lai,laipot,dwrt,&
     &  dwrtpot,dwlv,dwlvpot,dwst,dwstpot,gasst,gasstpot,mrest,mrestpot,&
     &  cwdm,cwdmpot,sla,slapot,lvage,lvagepot,lv,                      &
     &  lvpot,daycrop,nofd,atmin7,tsumgerm,rid,flgrazingpot,idregr,     &
     &  idregrpot,laiexppot,laimax,tagp,tagppot,idaysgraz,idaysgrazpot, &
     &  daymow,daymowpot,tagpt,tagptpot,tsum200,flCropHarvest,iharvest, &
     &  iseqgm,iseqgmpot,flgrazing,flCropEmergence,                     &
     &  flCropCalendar,slw,cuptgraz,cuptgrazpot)

! --- close output files
      close (wba)
      close (inc)
      close (str)
      if (swvap.eq.1) close (vap)
      if (swafo.ge.1) close (afo)
      if (swaun.ge.1) close (aun)
      if(swoutputmodflow.eq.1) call OutputModflow(3)
      if(swcapriseoutput) call capriseoutput(3)

      case default
         call fatalerr ('SoilWaterOutput', 'Illegal value for Task')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine outwba (task,wba,daynr,daycum,swscre,cevap,cgird,      &
     &    cgrai,csnrai,cnird,cnrai,cpeva,cptra,cqbot,cqdra,cqrot,crunon,&
     &    crunoff,gwl,cQMpOutDrRap,pond,t1900,date,volact,volini,       &
     &    wbalance,outfil,pathwork,project,flprintshort,floutput,       &
     &    swsnow,cqprai,ssnow,snowinco,PondIni,flheader)
! ----------------------------------------------------------------------
!     date               : July 2002
!     purpose            : write water balance data to outnam.wba file
! ---------------------------------------------------------------------
      implicit none

! -   global
      integer   wba,daynr,daycum,swscre,swsnow,task

      real(8)   cevap,cgird,cgrai,csnrai,cnird,cnrai,cpeva,cptra,cqbot
      real(8)   cqdra,cqrot,crunoff,gwl,cQMpOutDrRap,pond
      real(8)   volact,volini,wbalance,PondIni
      real(8)   cqprai,ssnow,snowinco,crunon,t1900
      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      character(len=11) date
      logical   flheader,flprintshort,floutput

! -   local
      integer   getun
      real(8)   dstor
      character(len=80) filnam,filtext
      character(len=10) gwlout
      character(len=19) datexti
      character(len=1)  comma

! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)
      
! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.wba'
      wba = getun (20,90)
      call fopens(wba,filnam,'new','del')
      filtext = 'cumulative water balance components (cm)'
      call writehead (wba,1,filnam,filtext,project)
      if (swscre .eq. 1) then
! ---   header of screen
        filtext = 'quick view water balance components'
        call writehead (5,1,'screen',filtext,project)
      endif

! --- write header in wba file
      if (flprintshort) then
        write (wba,10)
      else
        write (wba,12)
      endif
 10     format ('*',/,                                                  &
     &  '       Date,    Time,Day,  Dcum,Rain_g,Rain_n,Irr_g,Irr_n,',   &
     &  '  RunOn, RunOff,  Tpot,  Tact,  Epot,  Eact,  Drain,    Bot,', &
     &  '  DStor,      Gwl ,  Pond, Wbal,      Date2')
 12     format ('*',/,                                                  &
     &  '       Date,Day,  Dcum,Rain_g,Rain_n,Irr_g,Irr_n,  RunOn,',    &
     &  ' RunOff,  Tpot,  Tact,  Epot,  Eact,   Drain,     Bot,  ',     &
     &  ' DStor,       Gwl,  Pond, Wbal,      Date2')


! --- write header screen output
      if (swscre.eq.1) then
        write (*,20)
 20     format(/,t20,'cumulative water balance components (cm)')
        write (*,22)
 22     format (/,'        date   rain  irrig runoff transp evapor',    &
     &          '  drain   qbot    gwl   wbal'/,                        &
     &          '              gross  gross        actual actual',      &
     &          '    net    net           cum')
      endif
      return

      case (2)

! === write output water balance components ===============================

! --- write header in case of new balance period
      if (flheader) then
        if (flprintshort) then
          write (wba,10)
        else
          write (wba,12)
        endif
      endif

! --- determine date and date-time
      call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)

! --- write output record wba file
      gwlout = "          "
      if (gwl.lt.998.0d0)  write(gwlout,'(f9.1)') gwl
      if (swsnow.eq.0) then
        dstor = (volact + pond) - (volini + PondIni)
        if (flprintshort) then
          write (wba,25) datexti,comma,daynr,comma,daycum,comma,        &
     &    cgrai+csnrai,comma,cnrai,comma,cgird,comma,cnird,comma,crunon,&
     &    comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,comma,cevap,&
     &    comma,(cqdra+cQMpOutDrRap),comma,cqbot,comma,dstor,comma,     &
     &    gwlout,comma,pond,comma,wbalance,comma,date
        else
          write (wba,30) date,comma,daynr,comma,daycum,comma,           &
     &    cgrai+csnrai,comma,cqprai,comma,cgird,comma,cnird,comma,      &
     &    crunon,comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,     &
     &    comma,cevap,comma,(cqdra+cQMpOutDrRap),comma,cqbot,comma,     &
     &    dstor,comma,gwlout,comma,pond,comma,wbalance,comma,date
        endif
      else 
        dstor = (volact + pond + ssnow) - (volini + PondIni + snowinco)
        if (flprintshort) then
          write (wba,25) datexti,comma,daynr,comma,daycum,comma,        &
     &    cgrai+csnrai,comma,cnrai,comma,cgird,comma,cnird,comma,crunon,&
     &    comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,comma,cevap,&
     &    comma,(cqdra+cQMpOutDrRap),comma,cqbot,comma,dstor,comma,     &
     &    gwlout,comma,pond,comma,wbalance,comma,date
        else
          write (wba,30) date,comma,daynr,comma,daycum,comma,           &
     &    cgrai+csnrai,comma,cqprai,comma,cgird,comma,cnird,comma,      &
     &    crunon,comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,     &
     &    comma,cevap,comma,(cqdra+cQMpOutDrRap),comma,cqbot,comma,     &
     &    dstor,comma,gwlout,comma,pond,comma,wbalance,comma,date
        endif
      endif

 25   format (a19,a1,i3,a1,i6,2(a1,f7.2),2(a1,f5.1),a1,f7.2,a1,         &
     &     f7.2,a1,f7.2,3(a1,f8.2),3(a1,f7.2),2a,a1,f6.2,a1,f5.2,a1,a11)

 30   format (a11,a1,i3,a1,i6,2(a1,f7.2),2(a1,f5.1),a1,f7.2,a1,f7.2,    &
     &      4(a1,f6.2),3(a1,f8.3),2a,a1,f6.2,a1,f5.2,a1,a11)

! --- write output record screen

      if (swscre.eq.1 .and. floutput) then
        if (flheader) write (*,22)
        write (unit=*, fmt=40)                            &
     &                date,cgrai,cgird,crunoff,cqrot,cevap,             &
     &               (cqdra+cQMpOutDrRap),cqbot,gwl,wbalance
 40     format(1x,a11,f8.1,6f7.2,f7.1,f7.2)
      endif

      case default
         call fatalerr ('outwba', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outinc (task,inc,daynr,daycum,igrai,isnrai,igsnow,     &
     &   igird,iintc,irunon,iruno,iptra,iqrot,ipeva,ievap,isubl,iqdra,  &
     &   iQMpOutDrRap,iqbot,t1900,date,outfil,pathwork,project,flheader,&
     &   gwl,volact,volini,pond,PondIni,ssnow,snowinco,flprintshort)
! ----------------------------------------------------------------------
!     date               : july 2002
!     purpose            : write water balance increments to outnam.inc file
! ---------------------------------------------------------------------
      implicit none

! --- global
      integer   inc,daynr,daycum,task

      real(8)   igrai,igird,iintc,iruno,iptra,iqrot,ipeva,ievap,isubl
      real(8)   iqdra
      real(8)   isnrai,igsnow,irunon,iQMpOutDrRap,iqbot,t1900,gwl
      real(8)   volact,volini,pond,PondIni,ssnow,snowinco

      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      character(len=11) date

      logical   flheader,flprintshort

! --- local
      integer   getun
      real(8)  baldev,dstor
      character(len=80) filnam,filtext
      character(len=1)  comma
      character(len=10) gwlout
      character(len=19) datexti
      
      real(8), save ::   VolOld,PondOld,SnowOld
      
      
! ----------------------------------------------------------------------
      comma = ','
      
      select case (task)
      case (1)

      t1900 = t1900

! --- open output file once
      filnam = trim(pathwork)//trim(outfil)//'.inc'
      inc = getun (20,90)
      call fopens(inc,filnam,'new','del')
      filtext = 'water balance increments (cm/day)'
      call writehead (inc,1,filnam,filtext,project)

! --- write header of inc file
      if (flprintshort) then
        write (inc,10)
      else
        write (inc,12)
      endif
  10  format('*',/,                                                     &
     & '      Date,    Time,Day,  Dcum,      Rain,     Snow,',          &
     & '      Irrig,    Interc,     Runon,    Runoff,      Tpot,',      &
     & '      Tact,      Epot,      Eact,  Drainage,   QBottom,',       &
     & '       Gwl,  dstorage,    baldev')
  12  format('*',/,                                                     &
     & '       Date,Day,  Dcum,      Rain,     Snow,',                  &
     & '      Irrig,    Interc,     Runon,    Runoff,      Tpot,',      &
     & '      Tact,      Epot,      Eact,  Drainage,   QBottom,',       &
     & '       Gwl,  dstorage,    baldev')

      
      VolOld  = volini
      PondOld = PondIni
      SnowOld = snowinco
     
      return

      case (2)

! --- write header in case of new balance period
      if (flheader) then
        if (flprintshort) then
          write (inc,10)
        else
          write (inc,12)
        endif
      endif

! --- determine date and date-time
      call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)

! --- write output record .inc file
      gwlout = "          "
      if (gwl.lt.998.0d0)  write(gwlout,'(f9.1)') gwl
      dstor = (volact + pond + ssnow) - (VolOld + PondOld + SnowOld)
      baldev = (igrai+isnrai+igsnow+igird+irunon) - dstor -             &
     & (iintc+iruno+iqrot+ievap+isubl+iQMpOutDrRap+iqdra+(-1.0d0*iqbot))
      if (flprintshort) then
        write (inc,20) datexti,comma,daynr,comma,daycum,comma,          &
     &    igrai+isnrai,comma,igsnow,comma,igird,comma,iintc,comma,      &
     &    irunon,comma,iruno,comma,iptra,comma,iqrot,comma,ipeva,comma, &
     &    ievap,comma,(iQMpOutDrRap+iqdra),comma,iqbot,                 &
     &    comma,gwlout,comma,dstor,comma,baldev            !comma,storage
      else
        write (inc,22) date,comma,daynr,comma,daycum,comma,             &
     &    igrai+isnrai,comma,igsnow,comma,igird,comma,iintc,comma,      &
     &    irunon,comma,iruno,comma,iptra,comma,iqrot,comma,ipeva,comma, &
     &    ievap,comma,(iQMpOutDrRap+iqdra),comma,iqbot,                 &
     &    comma,gwlout,comma,dstor,comma,baldev           !comma,storage
      endif
 20   format (a19,a1,i3,a1,i6,12(a1,f10.5),2a,2(a1,f10.5))     !,(a1,e12.5)
 22   format (a11,a1,i3,a1,i6,12(a1,f10.5),2a,2(a1,f10.5))     !,(a1,e12.5)

      VolOld = volact
      PondOld = pond
      SnowOld = ssnow
      
      case default
         call fatalerr ('outinc', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outstr (task,str,daynr,daycum,ies0,iet0,iew0,ipeva,    &
     &  iptra,iqrot,iqredwet,iqreddry,iqredwsol,iqredfrs,t1900,         &
     &  date,outfil,pathwork,project,flheader,flprintshort)
! ----------------------------------------------------------------------
!     date               : March 2008
!     purpose            : write ETpot and stress factors to outfil.str file
! ---------------------------------------------------------------------
      implicit none

! --- global
      integer   str,daynr,daycum,task
      real(8)   ipeva,iptra,iqrot
      real(8)   iqredwet,iqreddry,iqredwsol,iqredfrs,t1900
      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      character(len=11) date
      logical   flheader,flprintshort
      real(8)   ies0  !  potential evaporation rate from a wet bare soil [cm/d]
      real(8)   iet0  !  potential transpiration rate from a dry crop [cm/d]
      real(8)   iew0  !  potential transpiration rate from a wet crop [cm/d]
!

! --- local
      integer   getun
      character(len=80)  filnam
      character(len=132) filtext
      character(len=1)   comma
      character(len=19)  datexti
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

      t1900 = t1900

! --- open output file once
      filnam = trim(pathwork)//trim(outfil)//'.str'
      str = getun (20,90)
      call fopens(str,filnam,'new','del')
      filtext = 'ES0, ET0, EW0, Epot, Tpot, Tact and 4 Tact-stress '//  &
     & 'values for wetness, drought, salinity and frost (cm/period)'
      call writehead (str,1,filnam,filtext,project)

! --- write header of str file
      if (flprintshort) then
        write (str,10)
      else
        write (str,12)
      endif
  10  format('*',/,                                                     &
     & '       Date,    Time,Day,  Dcum,ESoilWet,Tcropdry,Tcropwet,',   &
     & '    Epot,    Tpot,    Tact, Tredwet, Treddry, Tredsol, Tredfrs')
  12  format('*',/,                                                     &
     & '       Date,Day,  Dcum,ESoilWet,Tcropdry,Tcropwet,',            &
     & '    Epot,    Tpot,    Tact, Tredwet, Treddry, Tredsol, Tredfrs')
      return

      case (2)

! --- write header in case of new balance period
      if (flheader) then
        if (flprintshort) then
          write (str,10)
        else
          write (str,12)
        endif
      endif

! --- determine date and date-time
      call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)

! --- write output record .str file
      if (flprintshort) then
        write (str,20) datexti,comma,daynr,comma,daycum,                &
     &                 comma,ies0,comma,iet0,comma,iew0,comma,ipeva,    &
     &                 comma,iptra,comma,iqrot,comma,iqredwet,          &
     &                 comma,iqreddry,comma,iqredwsol,comma,iqredfrs
      else
        write (str,22) date,comma,daynr,comma,daycum,                   &
     &                 comma,ies0,comma,iet0,comma,iew0,comma,ipeva,    &
     &                 comma,iptra,comma,iqrot,comma,iqredwet,          &
     &                 comma,iqreddry,comma,iqredwsol,comma,iqredfrs
      endif
 20   format (a19,a1,i3,a1,i6,10(a1,f8.4))
 22   format (a11,a1,i3,a1,i6,10(a1,f8.4))

      case default
         call fatalerr ('outstr', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outvap (task,vap,daynr,numnod,daycum,z,cml,            &
     &    t1900,theta,h,k,tsoil,q,outfil,pathwork,project,swheader,     &
     &    isqtop,isqbot,qdraincomp,qrot,cmsy,dz,date,flprintshort)
! ----------------------------------------------------------------------
!     date               : December 2007
!     purpose            : write output of soil profile date
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer   vap,daynr,numnod,daycum,task,swheader
      real(8)   z(macp),cml(macp),t1900,isqtop,isqbot,dz(macp)
      real(8)   theta(macp),h(macp),k(macp+1),tsoil(macp),q(macp+1)
      real(8)   qdraincomp(macp),qrot(macp),cmsy(macp)
      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      character(len=11) date
      character(len=19) datexti
      character(len=11) inidate
      logical   flprintshort
! --- local
      integer   node,getun
      real(8)   sflux
      character(len=80) filnam,filtext
      character(len=1)  comma

! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)
! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.vap'
      vap = getun (20,90)
      call fopens(vap,filnam,'new','del')
      filtext = 'soil profile data'
      call writehead (vap,1,filnam,filtext,project)
      write (vap,100)

! --- write header in vap file
      if (flprintshort) then
        write (vap,200)
      else
        write (vap,210)
      endif

! --- write initial profile data to vap file

      if (flprintshort) then
! ---   determine date and date-time
        call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)
        do node = 1,numnod
          if (node.eq.1) then
             sflux = isqtop
              else
             sflux = 0.5d0 * (cml(node) + cml(node-1)) * q(node)
           end if
           write (vap,300) datexti,comma,z(node),comma,theta(node),     &
     &       comma,h(node),comma,k(node),comma,qdraincomp(node),comma,  &
     &       qrot(node),comma,q(node),comma,tsoil(node),comma,cml(node),&
     &      comma,cmsy(node),comma,sflux,comma,(z(node)+0.5d0*dz(node)),&
     &       comma,(z(node)-0.5d0*dz(node)),comma,daynr,comma,daycum
        end do
        sflux = isqbot
        write (vap,400) datexti,comma,(z(numnod)-0.5d0*dz(numnod)),     &
     &    comma,comma,comma,comma,comma,                                &
     &    comma,q(numnod+1),comma,comma,comma,comma,sflux,comma,        &
     &    (z(numnod)-0.5d0*dz(numnod)),comma,comma,daynr,comma,daycum

      else
! ---   determine initial date 
        call dtdpst ('year-month-day',t1900-0.1d0,inidate)
        do node = 1,numnod
          if (node.eq.1) then
             sflux = isqtop
              else
             sflux = 0.5d0 * (cml(node) + cml(node-1)) * q(node)
           end if
           write (vap,310) inidate,comma,z(node),comma,theta(node),     &
     &       comma,h(node),comma,k(node),comma,qdraincomp(node),comma,  &
     &       qrot(node),comma,q(node),comma,tsoil(node),comma,cml(node),&
     &      comma,cmsy(node),comma,sflux,comma,(z(node)+0.5d0*dz(node)),&
     &       comma,(z(node)-0.5d0*dz(node)),comma,daynr,comma,daycum
        end do
        sflux = isqbot
        write (vap,410) inidate,comma,(z(numnod)-0.5d0*dz(numnod)),     &
     &    comma,comma,comma,comma,comma,                                &
     &    comma,q(numnod+1),comma,comma,comma,comma,sflux,comma,        &
     &    (z(numnod)-0.5d0*dz(numnod)),comma,comma,daynr,comma,daycum

      endif

      return

      case (2)

! === write actual profile data ===========================================


! --- write header in vap file
      if (swheader .eq. 1) then
        if (flprintshort) then
          write (vap,200)
        else
          write (vap,210)
        endif
      endif

      if (flprintshort) then
! ---   determine date and date-time
        call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)
        do node = 1,numnod
          if (node.eq.1) then
             sflux = isqtop
              else
             sflux = 0.5d0 * (cml(node) + cml(node-1)) * q(node)
           end if
           write (vap,300) datexti,comma,z(node),comma,theta(node),     &
     &       comma,h(node),comma,k(node),comma,qdraincomp(node),comma,  &
     &       qrot(node),comma,q(node),comma,tsoil(node),comma,cml(node),&
     &      comma,cmsy(node),comma,sflux,comma,(z(node)+0.5d0*dz(node)),&
     &       comma,(z(node)-0.5d0*dz(node)),comma,daynr,comma,daycum
        end do
        sflux = isqbot
        write (vap,400) datexti,comma,(z(numnod)-0.5d0*dz(numnod)),     &
     &    comma,comma,comma,comma,comma,                                &
     &    comma,q(numnod+1),comma,comma,comma,comma,sflux,comma,        &
     &    (z(numnod)-0.5d0*dz(numnod)),comma,comma,daynr,comma,daycum

      else
        do node = 1,numnod
          if (node.eq.1) then
             sflux = isqtop
              else
             sflux = 0.5d0 * (cml(node) + cml(node-1)) * q(node)
           end if
           write (vap,310) date,comma,z(node),comma,theta(node),        &
     &       comma,h(node),comma,k(node),comma,qdraincomp(node),comma,  &
     &       qrot(node),comma,q(node),comma,tsoil(node),comma,cml(node),&
     &      comma,cmsy(node),comma,sflux,comma,(z(node)+0.5d0*dz(node)),&
     &       comma,(z(node)-0.5d0*dz(node)),comma,daynr,comma,daycum
        end do
        sflux = isqbot
        write (vap,410) date,comma,(z(numnod)-0.5d0*dz(numnod)),        &
     &    comma,comma,comma,comma,comma,                                &
     &    comma,q(numnod+1),comma,comma,comma,comma,sflux,comma,        &
     &    (z(numnod)-0.5d0*dz(numnod)),comma,comma,daynr,comma,daycum

      endif

 100  format(                                                           &
     & '* Explanation:   instantaneous fluxes of drainage, ',           &
     & 'root extraction, water and solute ',/                           &
     & '*                fluxes of water and solute apply to ',         &
     & 'top of compartment',/                                           &
     & '*                solute1 = concentration in soil water',/       &
     & '*                solute2 = total concentration ',               &
     & '(dissolved + adsorbed)')

 200  format(/,                                                         &
     & '                         cm,  cm3/cm3,         cm,       cm/d,',&
     & '       cm/d,       cm/d,       cm/d,     ºC,     mg/cm3,     ', &
     &'mg/cm3,   mg/cm2/d,     cm,     cm,  nr,   nr',/                 &
     & '      date,    time,  depth, wcontent,      phead,    hconduc,',&
     & '   drainage,    rootext,  waterflux,   temp,    solute1,    ',  &
     & 'solute2, soluteflux,    top, bottom, day, dcum')

 210  format(/,                                                         &
     & '                 cm,  cm3/cm3,         cm,       cm/d,',        &
     & '       cm/d,       cm/d,       cm/d,     ºC,     mg/cm3,     ', &
     &'mg/cm3,   mg/cm2/d,     cm,     cm,  nr,   nr',/                 &
     & '       date,  depth, wcontent,      phead,    hconduc,',        &
     & '   drainage,    rootext,  waterflux,   temp,    solute1,    ',  &
     & 'solute2, soluteflux,    top, bottom, day, dcum')

 300  format(a19,a1,f7.1,a1,f9.3,1p,5(a1,e11.3),a1,0p,f7.2,1p,          &
     &       3(a1,e11.3),a1,0p,2(f7.1,a1),i4,a1,i5)

 310  format(a11,a1,f7.1,a1,f9.3,1p,5(a1,e11.3),a1,0p,f7.2,1p,          &
     &       3(a1,e11.3),a1,0p,2(f7.1,a1),i4,a1,i5)

 400  format(a19,a1,f7.1,a1,t38,a1,t50,a1,t62,a1,t74,a1,t86,a1,1p,e11.3,&
     &       a1,t106,a1,t118,a1,t130,a1,e11.3,a1,0p,f7.1,a1,t158,a1,i4, &
     &       a1,i5)

 410  format(a11,a1,f7.1,a1,9x,4(a1,11x),a1,1p,e11.3,a1,7x,2(a1,11x),a1,&
     &       e11.3,a1,0p,f7.1,a1,7x,a1,i4,a1,i5)

      case default
         call fatalerr ('outvap', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outrot (task,rot,daynr,noddrz,daycum,z,hxylem,         &
     &    t1900,theta,hm1,q,outfil,pathwork,project,swheader,qrot,      &
     &    dz,date,flprintshort,hroot,inq,inqrot,mroot,mflux,            &
     &    rootrho,rootphi,hleaf)

! ----------------------------------------------------------------------
!     date               : February 2012
!     purpose            : write output of microscopic root water uptake
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer   rot,daynr,noddrz,daycum,task,swheader
      real(8)   z(macp),t1900,dz(macp),hxylem,hroot(macp)
      real(8)   theta(macp),hm1(macp),q(macp+1),inq(macp+1)
      real(8)   qrot(macp),inqrot(macp),mroot(macp),mflux(macp)
      real(8)   rootrho(macp),rootphi(macp),hleaf
      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      character(len=11) date
      character(len=19) datexti
!      character(len=11) inidate
      logical   flprintshort

! --- local
      integer   node,getun
      character(len=80) filnam,filtext
      character(len=1)  comma

! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)


! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.rot'
      rot = getun (20,90)
      call fopens(rot,filnam,'new','del')
      filtext = 'microscopic root water uptake'
      call writehead (rot,1,filnam,filtext,project)
      write (rot,100)

! --- write header in rot file
      if (swheader .eq. 0) then
        if (flprintshort) then
          write (rot,200)
        else
          write (rot,210)
        endif
      endif

      return

      case (2)

! === write actual profile data ===========================================

! --- write header in rot file
      if (swheader .eq. 1) then
        if (flprintshort) then
          write (rot,200)
        else
          write (rot,210)
        endif
      endif

      if (flprintshort) then
! ---   determine date and date-time
        call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)
        do node = 1,noddrz
           write (rot,300) datexti,comma,z(node),comma,hleaf,comma,     &
     &       hxylem,comma,                                              &
     &       hroot(node),comma,hm1(node),comma,inqrot(node),comma,      &
     &       qrot(node),comma,inq(node),comma,q(node),comma,mroot(node),&
     &       comma,mflux(node),comma,rootrho(node),comma,rootphi(node), &
     &       comma,theta(node),comma,(z(node)+0.5d0*dz(node)),comma,    &
     &       (z(node)-0.5d0*dz(node)),comma,daynr,comma,daycum
        end do

      else
        do node = 1,noddrz
           write (rot,310) date,comma,z(node),comma,hleaf,comma,        &
     &       hxylem,comma,                                              &
     &       hroot(node),comma,hm1(node),comma,inqrot(node),comma,      &
     &       qrot(node),comma,inq(node),comma,q(node),comma,mroot(node),&
     &       comma,mflux(node),comma,rootrho(node),comma,rootphi(node), &
     &       comma,theta(node),comma,(z(node)+0.5d0*dz(node)),comma,    &
     &       (z(node)-0.5d0*dz(node)),comma,daynr,comma,daycum
        end do
 
      endif

 100  format(                                                           &
     & '* Explanation:   fluxes of soil water (qsoilw and iqsoilw)',    &
     & ' apply to top of compartment;',/                                &
     & '*                both instantaneous (qroot and qsoilw) and',    &
     & ' incremental fluxes (qroot and qsoilw) are listed')

 200  format(/,                                                         &
     & '                         cm,         cm,         cm,',          &
     & '         cm,',                                                  &
     & '         cm,         cm,       cm/d,         cm,       cm/d,',  &
     & '      cm2/d,      cm2/d,       /cm2,        d/m,  cm3/cm3,',    &
     & '     cm,     cm,  nr,   nr',/                                   &
     & '       date,  depth,      hleaf,     hxylem,      hroot,',      &
     & '      hsoil,     iqroot,      qroot,    iqsoilw,     qsoilw,',  &
     & '      Mroot,      Msoil,    RootRho,    RootPhi, wcontent,',    &
     & '    top, bottom, day, dcum')

 210  format(/,                                                         &
     & '                 cm,         cm,         cm,         cm,',      &
     & '         cm,         cm,       cm/d,         cm,       cm/d,',  &
     & '      cm2/d,      cm2/d,       /cm2,        d/m,  cm3/cm3,',    &
     & '     cm,     cm,  nr,   nr',/                                   &
     & '       date,  depth,      hleaf,     hxylem,      hroot,',      &
     & '      hsoil,     iqroot,      qroot,    iqsoilw,     qsoilw,',  &
     & '      Mroot,      Msoil,    RootRho,    RootPhi, wcontent,',    &
     & '    top, bottom, day, dcum')

 300  format(a19,a1,f7.1,4(a1,f11.0),8(a1,e11.3),a1,f9.3,a1,2(f7.1,a1), &
     &       i4,a1,i5)

 310  format(a11,a1,f7.1,4(a1,f11.0),8(a1,e11.3),a1,f9.3,a1,2(f7.1,a1), &
     &       i4,a1,i5)

      case default
         call fatalerr ('outrot', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outbal (task,bal,logf,swscre,swdra,numnod,nrlevs,      &
     & swsolu,ioutdat,cevap,cgird,cgrai,cqbot,tstart,cqrot,crunoff,     &
     & crunon,cQMpOutDrRap,cqdra,cqdrain,dz,z,samini,sampro,samcra,     &
     & sqprec,sqirrig,sqbot,dectot,rottot,sqrap,sqdra,pond,volact,      &
     & volini,t1900,outdat,outfil,pathwork,project,                     &
     & caintc,csubl,PondIni,WaSrDm1,WaSrDm2,WaSrDm1Ini,WaSrDm2Ini,      &
     & swsnow,cgsnow,csnrai,snowinco,ssnow)
! ----------------------------------------------------------------------
!     date               : July 2002
!     purpose            : write overview balances to bal file
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer   bal,swdra,numnod,nrlevs,swsolu,ioutdat
      integer   logf,swscre,task,swsnow
      real(8)   cevap,cgird,cgrai,cqbot,caintc,tstart,crunon
      real(8)   cqrot,crunoff,cQMpOutDrRap,cqdra,cqdrain(5),dz(macp)
      real(8)   z(macp),samini,sampro,samcra,sqprec,sqirrig,sqbot,dectot
      real(8)   rottot,sqrap,sqdra,pond,volact,volini,t1900
      real(8)   outdat(maout),csubl,PondIni
      real(8)   WaSrDm1,WaSrDm2,WaSrDm1Ini,WaSrDm2Ini
      real(8)   precip,cgsnow,csnrai,snowinco,ssnow
      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project

! --- local
      integer   i,getun
      character(len=80) filnam,filtext
      character(len=11) datbegin,datend
! ----------------------------------------------------------------------

      select case (task)
      case (1)
      logf = logf ! for Forcheck
      swscre = swscre
! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.bal'
      bal = getun (20,90)
      call fopens(bal,filnam,'new','del')
      filtext='overview of actual water and solute balance components'
      call writehead (bal,1,filnam,filtext,project)

      return

      case (2)

! --- begin date of balance period
      if (ioutdat .eq. 2) then
        call dtdpst ('year-month-day',tstart+0.1d0,datbegin)
      else
        call dtdpst ('year-month-day',                                &
     &                                outdat(ioutdat-2)+1.1d0,datbegin)
      endif

! --- end date of balance period
      call dtdpst ('year-month-day',t1900-0.9d0,datend)

! --- write output record
      write (bal,20) datbegin,datend
      write (bal,22) (-z(numnod) + 0.5*dz(numnod))

      if (swsolu .eq. 1) then
          write (bal,24) (volact+pond+WaSrDm1+WaSrDm2+ssnow),           &
     &                                            (sampro+samcra),      &
     &                  (volini+PondIni+WaSrDm1Ini+WaSrDm2Ini+snowinco),&
     &                                             samini,              &
     &                  (volact+pond+WaSrDm1+WaSrDm2+ssnow-             &
     &                   volini-PondIni-WaSrDm1Ini-WaSrDm2Ini-snowinco),&
     &                                            (sampro+samcra-samini)
      else
          write (bal,25) (volact+pond+WaSrDm1+WaSrDm2+ssnow),           &
     &                  (volini+PondIni+WaSrDm1Ini+WaSrDm2Ini+snowinco),&
     &                  (volact+pond+WaSrDm1+WaSrDm2+ssnow-             &
     &                   volini-PondIni-WaSrDm1Ini-WaSrDm2Ini-snowinco)
      endif

      precip = cgrai
      if (swsnow.eq.1) then
         precip = precip + cgsnow + csnrai
      endif

      if (swsnow.eq.1) then
         write (bal,26) precip,caintc,crunon,crunoff,cgird,cqrot,cqbot, &
     &               (cevap+csubl),cQMpOutDrRap
      else
         write (bal,27) precip,caintc,crunon,crunoff,cgird,cqrot,cqbot, &
     &               (cevap+csubl),cQMpOutDrRap
      endif 
 
      if (swdra .ne. 0) then
        do i = 1,nrlevs
          write (bal,28) i,cqdrain(i)
        end do
      endif
      write(bal,30) (precip+cgird+cqbot+crunon),                        &
     &   (caintc+crunoff+cqrot+cevap+csubl+cQMpOutDrRap+cqdra)

      if (swsolu .eq. 1) then
        write (bal,34) sqprec,dectot,sqirrig,rottot,sqbot,sqrap,sqdra
        write (bal,36) (sqprec+sqirrig+sqbot),                          &
     &    (dectot+rottot+sqrap+sqdra)
      endif

 20   format(/'Period',t20,':',t23,a11,' until  ',a11)
 22   format('Depth soil profile',t20,':',f8.2,' cm')
 24   format(/,T13,'Water storage',t36,'Solute storage',/,              &
     &   'Final',t9,':',t15,f8.2,' cm',t31,e12.4,' mg/cm2',/,           &
     &   'Initial',t9,':',t15,f8.2,' cm',t31,e12.4,' mg/cm2',/,         &
     &    t13,13('='),t33,17('='),/,                                    &
     &   'Change',t15,f8.2,' cm',t31,e12.4,' mg/cm2')
 25   format(/,T13,'Water storage',/,'Final',t9,':',t15,f8.2,' cm',/,   &
     &   'Initial',t9,':',t15,f8.2,' cm',/,                             &
     &    t13,13('='),/,'Change',t15,f8.2,' cm')
 26   format(//,'Water balance components (cm)',//,                     &
     &   'In',t30,'Out',/,25('='),t30,28('='),/,                        &
     &   'Rain',t16,':',f9.2,t30,'Interception',t48,':',f9.2,/,         &
     &   'Runon',t16,':',f9.2,t30,'Runoff',t48,':',f9.2,/,              &
     &   'Irrigation',t16,':',f9.2,t30,'Transpiration',t48,':',f9.2,/,  &
     &   'Bottom flux',t16,':',f9.2,t30,'Soil evaporation',t48,':',f9.2,&
     &   /,t30,'Crack flux',t48,':',f9.2)
 27   format(//,'Water balance components (cm)',//,                     &
     &   'In',t30,'Out',/,25('='),t30,28('='),/,                        &
     &   'Rain + snow',t16,':',f9.2,t30,'Interception',t48,':',f9.2,/,  &
     &   'Runon',t16,':',f9.2,t30,'Runoff',t48,':',f9.2,/,              &
     &   'Irrigation',t16,':',f9.2,t30,'Transpiration',t48,':',f9.2,/,  &
     &   'Bottom flux',t16,':',f9.2,t30,'Soil evaporation',t48,':',f9.2,&
     &   /,t30,'Crack flux',t48,':',f9.2)
 28   format(t30,'Drainage level',i2,t48,':',f9.2)
 30   format(25('='),t30,28('='),/,                                     &
     &   'Sum',t16,':',f9.2,t30,'Sum',t48,':',f9.2)
 34   format(//,'Solute balance components (mg/cm2)',//,                &
     &   'In',t30,'Out',/,25('='),t30,28('='),/,                        &
     &   'Rain',t13,':',e12.4,t30,'Decomposition',t45,':',e12.4,/,      &
     &   'Irrigation',t13,':',e12.4,t30,'Root uptake',t45,':',e12.4,/,  &
     &   'Bottom flux',t13,':',e12.4,t30,'Cracks',t45,':',e12.4,/,      &
     &   t30,'Drainage',t45,':',e12.4)
 36   format(25('='),t30,28('='),/,                                     &
     &   'Sum',t13,':',e12.4,t30,'Sum',t45,':',e12.4,/)

      case default
         call fatalerr ('outbal', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine OutCropFixed(task,date,t,daycrop,dvs,tsum,lai,cf,rd,   &
     &   crp,ch)
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : open and write fixed crop output files 
! ----------------------------------------------------------------------
      implicit none

! --- global variables ------------------
      integer task,daycrop,crp
      real(8) t,dvs,tsum,lai,cf,rd,ch
      character(len=11) date
! --- local
      character(len=1) comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! --- write header of new crop ----------------------------------------------
      write (crp,100)
 100  format ('*',/,                                                    &
     & '*             day     day      -    grC       -      -      cm',&
     & '       -       cm    cm ',                                      &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '    kg/ha    kg/ha    kg/ha    kg/ha',                          &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha',/ &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM, LAIpot,    LAI,',      &
     & ' Height,CrpFac,RootdPot, Rootd,    PWLV,     WLV,    ',         &
     & 'PWST,     WST,    PWRT,     WRT,   CPWDM,    CWDM,   CPWSO,',   &
     & '    CWSO,PGRASSDM, GRASSDM,  PMOWDM,   MOWDM, PGRAZDM,  GRAZDM',&
     & ', DWLVCROP, DWLVSOIL,     DWST,     DWRT,     DWSO,HarLosOrm')

      return

      case (2)

! --- write actual data ------------------------------------------------------

! --- write output record
      write (crp,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,"       ",comma,lai,comma,ch,comma,cf,                &
     & comma,"       ",comma,nint(rd),                                  &
     & comma,comma,comma,comma,comma,comma,comma,comma,comma,comma,     &
     & comma,comma,comma,comma,comma,comma,                             &
     & comma,comma,comma,comma,comma,comma
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0,a1,a,2(a1,f7.2),(a1,f6.2),&
     &      (a1,a),(a1,i7),16(a1,'        ') ,                          &
     &      6(a1,'         ') )

      case default
         call fatalerr ('outcropfixed', 'Illegal value for Task')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine outend (                                               &
     &          numnod,h,flSolute,flAgeTracer,cml,z,fltemperature,tsoil,&
     &          ssnow,pond,dt,icrop,croptype,cropname,flSurfaceWater,   &
     &          wls,swredu,ldwet,spev,outfil,pathwork,project,          &
     &  rd,rdpot,dvs,flanthesis,tsum,ilvold,ilvoldpot,wrt,wrtpot,tadw,  &
     &  tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,lai,laipot,dwrt,&
     &  dwrtpot,dwlv,dwlvpot,dwst,dwstpot,gasst,gasstpot,mrest,mrestpot,&
     &  cwdm,cwdmpot,sla,slapot,lvage,lvagepot,lv,                      &
     &  lvpot,daycrop,nofd,atmin7,tsumgerm,rid,flgrazingpot,idregr,     &
     &  idregrpot,laiexppot,laimax,tagp,tagppot,idaysgraz,idaysgrazpot, &
     &  daymow,daymowpot,tagpt,tagptpot,tsum200,flCropHarvest,iharvest, &
     &  iseqgm,iseqgmpot,flgrazing,flCropEmergence,                     &
     &  flCropCalendar,slw,cuptgraz,cuptgrazpot)
! ----------------------------------------------------------------------
!     Update             : September 2015
!     date               : July 2002
!     purpose            : write final result to .end file
! ---------------------------------------------------------------------
      implicit none
      include  'arrays.fi'

! --- global
      integer   numnod,count,ilvold,ilvoldpot,daycrop,nofd,idaysgrazpot
      integer   idregr,idregrpot,iseqgm,iseqgmpot,iharvest,idaysgraz
      integer   daymow,daymowpot
      real(8)   h(macp),cml(macp),z(macp),tsoil(macp), Ssnow, Pond, wls
      integer   icrop                ! Current crop number
      integer   croptype(macrop)     ! Switch for type of crop model: 
      integer   swredu               ! Switch for reduction of soil evaporation
      real(8)   ldwet                ! Length of dry period (L) used in Black's model
      real(8)   spev                 ! Cumulative potential evaporation (L) used in Boesten/Stroosnijder
      real(8)   dt                   ! Final Time step (T)
      real(8) rd,rdpot,dvs,tsum,wrt,wrtpot,tadw,tsumgerm,rid,slw
      real(8) tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,lai,laipot
      real(8) dwrtpot,dwlv,dwlvpot,dwst,dwstpot,gasst,gasstpot,mrest
      real(8) cwdm,cwdmpot,dwrt,mrestpot
      real(8) laiexppot,laimax,tagp,tagppot,tagpt,tagptpot,tsum200
      real(8) sla(366),slapot(366),lvage(366),lvagepot(366),lv(366)
      real(8) lvpot(366),atmin7(7),cuptgraz,cuptgrazpot
      logical  flsurfacewater        ! flag to indicate simulation of separate surface water system
      logical   flsolute             ! Flag indicating simulation of solute transport
      logical   flAgeTracer          ! Flag indicating simulation of Ageing
      logical   flanthesis           ! Flag indicating anthesis stage of a crop
      logical   fltemperature        ! Flag indicating simulation of soil heat flow
!      logical   flcropoutput         ! Flag indicating that crop output can and should be written
      logical   flCropHarvest,flgrazing,flCropEmergence
      logical   flCropCalendar,flgrazingpot !,flCultivate
      character(len=40) cropname(macrop)  ! Array with crop names
      character(len=16) outfil
      character(len=*)  pathwork,project

! --- local
      integer   i,fin,getun,swanthesis,swcropharvest,swgrazing
      integer   swCropEmergence,swgrazingpot
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','
! ----------------------------------------------------------------------

! --- open output file for final profile data and write heading
      filnam = trim(pathwork)//trim(outfil)//'.end'
      fin = getun (20,90)
      call fopens(fin,filnam,'new','del')
      filtext = 'final state variables'
      call writehead (fin,1,filnam,filtext,project)


! --- write snow- and ponding-layer
      write (fin,101)  ssnow
 101  format (/,'*  Snow layer (Ssnow in cm)',/,' Ssnow = ',e12.5)
      write (fin,102)  slw
 102  format (/,'*  Liquid water in snow layer (cm)',/,' Slw = ',e12.5)
      write (fin,103)  pond
 103  format (/,'*  Ponding layer (Pond in cm)',/,' Pond = ',e12.5)

! --- write soil water pressure heads
      write (fin,113)
 113  format (/,'*  Soil water pressure heads  (z in cm, h in cm)')
      write (fin,114) 
 114  format (t6,'  z_h,h')
      do i = 1, numnod
        write (fin,110) z(i),comma,h(i)
      end do
 110  format (f10.1,a1,1p,e12.5) 

! --- write solute concentrations
      if(flSolute .or. flAgeTracer) then
        write (fin,125)
 125    format (/,                                                      &
     &          '*  Solute concentrations (z in cm, Cml in mg/cm3)')
        write (fin,120)
 120    format (t6,'z_Cml,Cml')
        do i = 1, numnod
          write (fin,110) z(i),comma,cml(i)
        end do
      endif

! --- write soil temperatures
      if(fltemperature) then
        write (fin,135)
 135    format (/,'*  Soil temperatures  (z in cm, Tsoil in C)')
        write (fin,130)
 130    format (t4,'z_Tsoil,Tsoil')
        do i = 1, numnod
          write (fin,110) z(i),comma,tsoil(i)
        end do
      endif

! --- write surface water level
      if(flSurfaceWater) then
        write (fin,159)
 159    format (/,'*  Surface water')
        write (fin,151) wls
 151    format ('*  Surface water level (cm)',/,' wls = ',e12.5)
      endif

! --- write soil evaporation reservoirs
      if(swredu.eq.1) then
        write (fin,162)
 162    format (/,'*  Soil evaporation reservoir (Black)')
        write (fin,163) ldwet
 163    format ('*  Time after significant rainfall (d)',/,             &
     &          ' ldwet = ',e12.5)
      endif
      if(swredu.eq.2) then
        write (fin,164)
 164    format (/,                                                      &
     &          '*  Soil evaporation reservoir (Boesten/Stroosnijder)')
        write (fin,165) spev
 165    format ('*  Rainfall excess (cm)',/,                            &
     &          ' spev = ',e12.5)
      endif

! --- write length of final timestep (d)
      write (fin,176)
 176  format (/,'*  Timing parameters')
      write (fin,177)  dt
 177  format ('*  Length of final time step (d)',/,' dt = ',e12.5)

! --- write weather data
      write (fin,180)
 180  format (/,                                                        &
     &          '*  Weather data'/,                                     &
     &          '*  Minimum temperatures of last week')
      write (fin,182)  (atmin7(i), i=1,7)
 182  format (' atmin7 = ',/,7e14.5)

      if (flCropCalendar) then
         
! ---   write crop growth data

        if (croptype(icrop) .eq. 1) then
          write (fin,184) trim(cropname(icrop))
 184      format (/,'*  Crop growth data of last crop: ',a)
           ! ---     write fixed crop data
          write (fin,186)  nofd
 186      format (' nofd = ',i4)
          write (fin,187)  daycrop
 187      format (' daycrop = ',i4)
          write (fin,188)  dvs
 188      format (' dvs = ',e15.8)
          write (fin,190)  tsum
 190      format (' tsum = ',e15.8)
        else
! ---     write crop data of wofost or simulated grass
          
          if (.not. flCropHarvest) then          
            write (fin,184) trim(cropname(icrop))
            write (fin,202)  rd
 202        format (' rd = ',e15.8)
            write (fin,204)  rdpot
 204        format (' rdpot = ',e15.8)
            write (fin,206)  dvs
 206        format (' dvs = ',e15.8)
            write (fin,207)  daycrop
 207        format (' daycrop = ',i4)
            if (flanthesis) then
              swanthesis = 1
            else
              swanthesis = 0
            endif 
            write (fin,208)  swanthesis
 208        format (' swanthesis = ',i4)
            write (fin,210)  tsum
 210        format (' tsum = ',e15.8)
            write (fin,212)  ilvold
 212        format (' ilvold = ',i4)
            write (fin,214)  ilvoldpot
 214        format (' ilvoldpot = ',i4)
            write (fin,216)  wrt
 216        format (' wrt = ',e15.8)
            write (fin,218)  wrtpot
 218        format (' wrtpot = ',e15.8)
            write (fin,220)  tadw
 220        format (' tadw = ',e15.8)
            write (fin,222)  tadwpot
 222        format (' tadwpot = ',e15.8)
            write (fin,224)  wst
 224        format (' wst = ',e15.8)
            write (fin,226)  wstpot
 226        format (' wstpot = ',e15.8)
            write (fin,228)  wso
 228        format (' wso = ',e15.8)
            write (fin,230)  wsopot
 230        format (' wsopot = ',e15.8)
            write (fin,232)  wlv
 232        format (' wlv = ',e15.8)
            write (fin,234)  wlvpot
 234        format (' wlvpot = ',e15.8)
            write (fin,236)  laiexp
 236        format (' laiexp = ',e15.8)
            write (fin,238)  lai
 238        format (' lai = ',e15.8)
            write (fin,240)  laipot
 240        format (' laipot = ',e15.8)
            write (fin,242)  dwrt
 242        format (' dwrt = ',e15.8)
            write (fin,244)  dwrtpot
 244        format (' dwrtpot = ',e15.8)
            write (fin,246)  dwlv
 246        format (' dwlv = ',e15.8)
            write (fin,248)  dwlvpot
 248        format (' dwlvpot = ',e15.8)
            write (fin,250)  dwst
 250        format (' dwst = ',e15.8)
            write (fin,252)  dwstpot
 252        format (' dwstpot = ',e15.8)
            write (fin,254)  gasst
 254        format (' gasst = ',e15.8)
            write (fin,256)  gasstpot
 256        format (' gasstpot = ',e15.8)
            write (fin,258)  mrest
 258        format (' mrest = ',e15.8)
            write (fin,260)  mrestpot
 260        format (' mrestpot = ',e15.8)
            write (fin,270)  cwdm
 270        format (' cwdm = ',e15.8)
            write (fin,272)  cwdmpot
 272        format (' cwdmpot = ',e15.8)
            write (fin,274)  tsumgerm
 274        format (' tsumgerm = ',e15.8)
            write (fin,276)  nofd
 276        format (' nofd = ',i4)
            write (fin,277)  iharvest
 277        format (' iharvest = ',i4)
            write (fin,278)  rid
 278        format (' rid = ',e15.8)
            write (fin,280)  idregr
 280        format (' idregr = ',i4)
            write (fin,282)  idregrpot
 282        format (' idregrpot = ',i4)
            write (fin,284)  laiexppot
 284        format (' laiexppot = ',e15.8)
            write (fin,286)  laimax
 286        format (' laimax = ',e15.8)
            write (fin,287)  tagp
 287        format (' tagp = ',e15.8)
            write (fin,288)  tagppot
 288        format (' tagppot = ',e15.8)
            write (fin,289)  tagpt
 289        format (' tagpt = ',e15.8)
            write (fin,290)  tagptpot
 290        format (' tagptpot = ',e15.8)
            write (fin,291)  cuptgraz
 291        format (' cuptgraz = ',e15.8)
            write (fin,292)  cuptgrazpot
 292        format (' cuptgrazpot = ',e15.8)
            write (fin,293)  tsum200
 293        format (' tsum200 = ',e15.8)
            if (flCropHarvest) then
              swcropharvest = 1
            else
              swcropharvest = 0
            endif 
            write (fin,294)  daymow
 294        format (' daymow = ',i4)
            write (fin,295)  daymowpot
 295        format (' daymowpot = ',i4)
            write (fin,298)  swcropharvest
 298        format (' swcropharvest = ',i4)
            write (fin,302)  iseqgm
 302        format (' iseqgm = ',i4)
            write (fin,304)  iseqgmpot
 304        format (' iseqgmpot = ',i4)
            if (flgrazing) then
              swgrazing = 1
            else
              swgrazing = 0
            endif 
            write (fin,306)  swgrazing
 306        format (' swgrazing = ',i4)
            if (flgrazingpot) then
              swgrazingpot = 1
            else
              swgrazingpot = 0
            endif 
            write (fin,307)  swgrazingpot
 307        format (' swgrazingpot = ',i4)
            if (flCropEmergence) then
              swCropEmergence = 1
            else
              swCropEmergence = 0
            endif 
            write (fin,308)  swCropEmergence
 308        format (' swCropEmergence = ',i4)
            write (fin,310)  idaysgraz
 310        format (' idaysgraz = ',i4)
            write (fin,312)  idaysgrazpot
 312        format (' idaysgrazpot = ',i4)
            
            write (fin,314)
 314       format (/' day            sla         slapot          lvage',&
     &         '       lvagepot             lv          lvpot')
            count = max(ilvold,ilvoldpot)
            do i = 1,count
              write (fin,316)  i,sla(i),slapot(i),lvage(i),lvagepot(i), &
     &        lv(i),lvpot(i)
 316          format (i4,6e15.8)
            enddo
          else

            write (fin,184) trim(cropname(icrop))
          
            if (flCropEmergence) then
              swCropEmergence = 1
            else
              swCropEmergence = 0
            endif 
            write (fin,308) swCropEmergence
            if (flCropHarvest) then
              swcropharvest = 1
            else
              swcropharvest = 0
            endif 
            write (fin,298)  swcropharvest
          
          endif
        endif
      endif

! --- end write crop growth data
      close (fin)

      return
      end

! ----------------------------------------------------------------------
      subroutine OutWofost(task,date,daycrop,crp,t,dvs,tsum,laipot,lai, &
     &     rdpot,rd,ch,cf,cwdmpot,cwdm,wsopot,wso,                      &
     &     wlvpot,wlv,wstpot,wst,wrtpot,wrt,                            &
     &     dwlvCrop,dwlvSoil,dwst,dwrt,dwso,HarLosOrm_tot)
! ----------------------------------------------------------------------
!     UpDate             : Aug 2014
!     Date               : Oct 2004   
!     Purpose            : Write detailed crop growth output files 
! ----------------------------------------------------------------------
      implicit none

! --- global variables ------------------
      integer task,daycrop,crp
      real(8) t,dvs,tsum,laipot,lai,cf,rdpot,rd,ch
      real(8) cwdmpot,cwdm,wsopot,wso,wstpot,wst,wlvpot,wlv,wrtpot,wrt
      real(8) dwlvCrop,dwlvSoil,dwst,dwrt,dwso,HarLosOrm_tot
      character(len=11) date
! --- local variables ------------------
      character(len=1) comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! --- write header of new crop ----------------------------------------------
      write (crp,100)
 100  format ('*',/,                                                    &
     & '*             day     day      -    grC       -      -      cm',&
     & '       -       cm    cm ',                                      &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '    kg/ha    kg/ha    kg/ha    kg/ha',                          &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha',/ &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM, LAIpot,    LAI,',      &
     & ' Height,CrpFac,RootdPot, Rootd,    PWLV,     WLV,    ',         &
     & 'PWST,     WST,    PWRT,     WRT,   CPWDM,    CWDM,   CPWSO,',   &
     & '    CWSO,PGRASSDM, GRASSDM,  PMOWDM,   MOWDM, PGRAZDM,  GRAZDM',&
     & ', DWLVCROP, DWLVSOIL,     DWST,     DWRT,     DWSO,HarLosOrm')

      return

      case (2)

! --- write actual data ------------------------------------------------------

! --- write output record
      write (crp,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,laipot,comma,lai,comma,ch,comma,cf,comma,nint(rdpot), &
     & comma,nint(rd),                                                  &
     & comma,nint(wlvpot),comma,nint(wlv),comma,nint(wstpot),           &
     & comma,nint(wst),comma,nint(wrtpot),comma,nint(wrt),comma,        &
     & nint(cwdmpot),comma,nint(cwdm),comma,nint(wsopot),comma,nint(wso)&
     & ,comma,comma,comma,comma,comma,comma,                            &
     &  comma,dwlvCrop,comma,dwlvSoil,comma,dwst,comma,dwrt,comma,dwso, &
     &  comma,HarLosOrm_tot
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0,3(a1,f7.2),(a1,f6.2),     &
     &        2(a1,i7),  10(a1,i8), 6(a1,'        ') ,                  &
     &        6(a1,f9.2) )

      case default
         call fatalerr ('outWofost', 'Illegal value for Task')
      end select

      return
      end
! ----------------------------------------------------------------------
      subroutine OutGrass(task,date,daycrop,crp,t,dvs,tsum200,laipot,   &
     &         lai,rdpot,rd,ch,cf,tagppot,tagp,tagptpot,tagpt,          &
     &         wlvpot,wlv,wstpot,wst,wrtpot,wrt,cuptgraz,cuptgrazpot)
! ----------------------------------------------------------------------
!     UpDate             : Aug 2014
!     Date               : Oct 2004   
!     Purpose            : Write detailed grass simulation output files 
! ----------------------------------------------------------------------
      implicit none

! --- global variables ------------------
      integer task,daycrop,crp
      real(8) t,dvs,laipot,lai,cf,rdpot,rd,ch
      real(8) tagppot,tagp,tagptpot,tagpt,tsum200
      real(8) wlvpot,wlv,wstpot,wst,wrtpot,wrt,cuptgraz,cuptgrazpot
      character(len=11) date
! --- local variables ------------------
      character(len=1) comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! --- write header of new crop ----------------------------------------------

      write (crp,100)
 100  format ('*',/,                                                    &
     & '*             day     day      -    grC       -      -      cm',&
     & '       -       cm    cm ',                                      &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '    kg/ha    kg/ha    kg/ha    kg/ha',                          &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha',/ &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM, LAIpot,    LAI,',      &
     & ' Height,CrpFac,RootdPot, Rootd,    PWLV,     WLV,    ',         &
     & 'PWST,     WST,    PWRT,     WRT,   CPWDM,    CWDM,   CPWSO,',   &
     & '    CWSO,PGRASSDM, GRASSDM,  PMOWDM,   MOWDM, PGRAZDM,  GRAZDM',&
     & ', DWLVCROP, DWLVSOIL,     DWST,     DWRT,     DWSO,HarLosOrm')

      return

      case (2)

! --- write actual data ------------------------------------------------------

! --- write output record
      write (crp,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum200,comma,laipot,comma,lai,comma,ch,comma,cf,                &
     & comma,nint(rdpot),comma,nint(rd),                                &
     & comma,nint(wlvpot),comma,nint(wlv),comma,nint(wstpot),           &
     & comma,nint(wst),comma,nint(wrtpot),comma,nint(wrt),              &
     & comma,comma,comma,comma,comma,nint(tagppot),                     &
     & comma,nint(tagp),comma,nint(tagptpot),comma,nint(tagpt),         &
     & comma,nint(cuptgrazpot),comma,nint(cuptgraz),                    &
     & comma,comma,comma,comma,comma,comma
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0,3(a1,f7.2),(a1,f6.2),     &
     &  2(a1,i7), 6(a1,i8), 4(a1,'        '), 6(a1,i8),                 &
     &  6(a1,'         '))

      case default
         call fatalerr ('outGrass', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine SoluteOutput(task) 
! ----------------------------------------------------------------------
!     Date               : November 2004   
!     Purpose            : open and write solute output files 
! ----------------------------------------------------------------------

      use Variables
      implicit none

      integer task

      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  sba file
      call outsba (1,sba,daynr,daycum,sampro,sqbot,project,             &
     &             sqdra,solbal,dectot,rottot,sqprec,                   &
     &             date,sqirrig,outfil,pathwork,flheader)

      return

      case (2)

! === write actual data ===============================

! --  sba file
      call outsba (2,sba,daynr,daycum,sampro,sqbot,project,             &
     &             sqdra,solbal,dectot,rottot,sqprec,                   &
     &             date,sqirrig,outfil,pathwork,flheader)

      return

      case (3)

! === close output files ===========================

! --- close sba file
      close (sba)

      case default
         call fatalerr ('SoluteOutput', 'Illegal value for Task')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine outsba (task,sba,daynr,daycum,sampro,sqbot,            &
     &    project,sqdra,solbal,dectot,rottot,sqprec,                    &
     &    date,sqirrig,outfil,pathwork,flheader)
! ----------------------------------------------------------------------
!     date               : July 2002
!     purpose            : output of salt balance
! ---------------------------------------------------------------------
      implicit none

! --- global
      integer   sba,daynr,daycum,task
      real(8)   sampro,sqbot,sqdra
      real(8)   solbal,dectot,rottot,sqprec,sqirrig
      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      character(len=11) date
      logical   flheader
! --- local variables ------------------
      integer   getun
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','


      select case (task)
      case (1)

! --- open output file -------------------------------------------------
      filnam = trim(pathwork)//trim(outfil)//'.sba'
      sba = getun (20,90)
      call fopens(sba,filnam,'new','del')
      filtext = 'cumulative solute balance components (mg/cm2)'
      call writehead (sba,1,filnam,filtext,project)

! --- write header in sba file
      write (sba,10)
10    format ('*',/,                                                    &
     & '       Date, Day,  Dcum,      Flux top,   Root uptake, ',       &
     & 'Decomposition,      Drainage,   Flux bottom,       Storage,',   &
     & '       Balance,        Date2')

      return

      case (2)

! === write actual data =================================================

! --- write header in case of new balance period
      if (flheader) write (sba,10)

! --- write output solute balance components ----------------------------

      write(sba,15) date,comma,daynr,comma,daycum,comma,(sqprec+sqirrig)&
     & ,comma,rottot,comma,dectot,comma,sqdra,comma,sqbot,comma,        &
     & sampro,comma,solbal,comma,date

 15   format(a11,a1,i4,a1,i6,6(a1,e14.5),a1,e14.2,a1,a13)

      case default
         call fatalerr ('outsba', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine AgeTracerOutput(task) 
! ----------------------------------------------------------------------
!     Date               : October 2010
!     Purpose            : open and write Groundwater Ageing output files 
! ----------------------------------------------------------------------
      use Variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer agep,agee,ageq

      save   agee,agep,ageq

      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  age files
      call outage (1,agep,agee,ageq,daynr,daycum,date,outper,           &
     &          project,nrlevs,outfil,pathwork,numnod,z,cml,            &
     &          AgeGwl1m,icAgeBot,icAgeDra,icAgeRot,icAgeSur,inqdra)

      return

      case (2)

! === write actual data ===============================

! --  age files
      call outage (2,agep,agee,ageq,daynr,daycum,date,outper,           &
     &          project,nrlevs,outfil,pathwork,numnod,z,cml,            &
     &          AgeGwl1m,icAgeBot,icAgeDra,icAgeRot,icAgeSur,inqdra)

      return

      case (3)

! === close output files ===========================

! --- close sba file
      close (agep)
      close (agee)

      case default
         call fatalerr ('AgeTracerOutput', 'Illegal value for Task')
      end select

      return
      end 
! ----------------------------------------------------------------------
      subroutine outage(task,agep,agee,ageq,daynr,daycum,date,outper,   &
     &          project,nrlevs,outfil,pathwork,numnod,z,cml,            &
     &          AgeGwl1m,icAgeBot,icAgeDra,icAgeRot,icAgeSur,inqdra)
! ----------------------------------------------------------------------
!     date               : October 2010
!     purpose            : output of groundwater age
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer   agep,agee,ageq,daynr,daycum,task,numnod,nrlevs
      real(8)   outper             ! Length of actual output interval (T)
      real(8)   z(macp)            ! Depth of a node (L)
      real(8)   cml(macp)          ! Age (d) of groundwater
      real(8)   AgeGwl1m           ! Age (d) of groundwater in upper 1 meter of saturated zone
      real(8)   icAgeBot           ! Incremental (over output interval) age (d) of groundwater leaving bottom comp.
      real(8)   icAgeDra(madr)     ! Incremental (over output interval) age (d) of groundwater leaving bottom comp.
      real(8)   icAgeRot           ! Incremental (over output interval) age (d) of groundwater leaving by root uptake
      real(8)   icAgeSur           ! Incremental (over output interval) age (d) of groundwater leaving by surface runoff
      real(8)   inqdra(madr,macp)  ! Array with intermediate amounts of lateral drainage for each level and compartment (L)
      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      character(len=11) date
! --- local variables ------------------
      integer   getun,reclngth,node,level
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      real(8)   iqdrainout(madr)   ! Cumulative (over 1 output timestep) drainage flux (L) for each drainage level
! ----------------------------------------------------------------------
      comma = ','


      select case (task)
      case (1)

! --- open output files -------------------------------------------------
!     age of groundwater as profile
      filnam = trim(pathwork)//trim(outfil)//'.ageProfile.csv'
      agep = getun (20,90)
      reclngth = 50 + 12*numnod 
      open(unit=agep,file=filnam,status='unknown',recl=reclngth)
      filtext = 'Groundwater age profiles (all age-values in days)'
      call writehead (agep,1,filnam,filtext,project)
!     age of groundwater in effluents: drains, transpiration, leaching, runoff
      filnam = trim(pathwork)//trim(outfil)//'.ageEffluent.csv'
      agee = getun (20,90)
      call fopens(agee,filnam,'new','del')
      filtext = 'Groundwater age effluent (all age-values in days)'
      call writehead (agee,1,filnam,filtext,project)
!     effluent drain water fluxes
      filnam = trim(pathwork)//trim(outfil)//'.ageEffluentqDrain.csv'
      ageq = getun (20,90)
      call fopens(ageq,filnam,'new','del')
      filtext = 'Drain water effluent (mm/day)'
      call writehead (ageq,1,filnam,filtext,project)

! --- write headers of files
      write (agep,9) (z(node),node=1,numnod)
      if (numnod.le.9) then
         write (agep,10) (node,node=1,numnod)
      else
         write (agep,11) (node,node=1,9), (node,node=10,numnod)
      endif
  9   format('*',t8,' NodeDepth (cm) =,,',18(',',f7.2) ,997(',',f8.2) )
 10   format('*',/, t8,'Date,Day,Daycum',   9(',Node',i3.3) )
 11   format('*',/, t8,'Date,Day,Daycum',   9(',Node',i3.3),            &
     &                                   1015(',Node',i4.4) )
 
      write (agee,'(3a)') ' Date,daynr,daycum,AgeGwl1m,AgeBottom,',     &
     &'AgeRootUpt,AgeRunoff,AgeDrainSys1,AgeDrainSys2,AgeDrainSys3,',   &
     &'AgeDrainSys4,AgeDrainSys5'

      write (ageq,'(3a)') ' Date,daynr,daycum,',                        &
     &'qDrainSys1,qDrainSys2,qDrainSys3,qDrainSys4,qDrainSys5'

      return

      case (2)

! === write actual data =================================================

!     age of groundwater as profile
      write(agep,15) date,comma,daynr,comma,daycum,                     &
     &               (comma,cml(node),node=1,numnod)
 15   format(a11,a1,i4,a1,i6,1p,1024(a1,e10.3))

!     age of groundwater in effluents: drains, transpiration, leaching, runoff
!     and age (d) of groundwater in upper 1 meter of saturated zone
      write(agee,16) date,comma,daynr,comma,daycum,comma,AgeGwl1m,comma,&
     &               icAgeBot/outper,comma,icAgeRot/outper,comma,       &
     &               icAgeSur/outper,                                   &
     &               (comma,icAgeDra(level)/outper,level=1,nrlevs)
 16   format(a11,a1,i4,a1,i6,1p,9(a1,e10.3))

!     qdrain discharge-effluent (without infiltration!)
      do level = 1,nrlevs
        iqdrainout(level) = 0.0d0
        do node = 1,numnod
          if (inqdra(level,node).gt.0.0d0) then
           iqdrainout(level) = iqdrainout(level) + inqdra(level,node)
          endif      
        enddo
      enddo
      write(ageq,16) date,comma,daynr,comma,daycum,                     &
     &               (comma,iqdrainout(level),level=1,nrlevs)

      case default
         call fatalerr ('outage', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine TemperatureOutput(task) 
! ----------------------------------------------------------------------
!     Date               : November 2004   
!     Purpose            : open and write soil temperature output files 
! ----------------------------------------------------------------------

! --- global variables ------------------
      use Variables

      implicit none
      integer task


      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  heat params file
      if (swcalt.eq.2)  call outheapar()

! --  tem file
      if (swtem .eq. 1) call outtem (task,numnod,date,daynr,tem,daycum, &
     &    tav,tebot, tsoil,tetop,outfil,pathwork,flheader,project)

      return

      case (2)

! === write actual data ===============================

! --  tem file
      if (swtem .eq. 1) call outtem (task,numnod,date,daynr,tem,daycum, &
     &    tav,tebot, tsoil,tetop,outfil,pathwork,flheader,project)

      return

      case (3)

! === close output files ===========================

! --- close tem file
      if (swtem .eq. 1) close (tem)

      case default
         call fatalerr ('TemperatureOutput', 'Illegal value for Task')
      end select

      return
      end 



! ----------------------------------------------------------------------
      subroutine outheapar ()
! ----------------------------------------------------------------------
!     date               : February 2005
!     purpose            : Output of soil heat conductivity and capacity
! ---------------------------------------------------------------------
! --- global
      use variables

! --- local variables ------------------
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      integer   getun,hea,lay,node,j
      real(8)   heacap(macp),heacnd(macp),thetadum(numnod)
! ---------------------------------------------------------------------
      comma = ',' 

! === open output file =================================================
      filnam = trim(pathwork)//'heatparam.csv'
      hea = getun (20,90)
      call fopens(hea,filnam,'new','del')
      filtext = 'soil heat conductivity and capacity'
      call writehead (hea,1,filnam,filtext,project)

! --- write header
      write(hea,'(a)') 'layer, theta,heacap(J/cm3/K),heacnd(J/cm/K/d)'

! --- write thermal properties: heat capacity and thermal conductivity 
      do node = 1,numnod
        thetadum(node) = thetas(node)
      enddo
      do lay = 1,numlay
!         find first Node of the Layer
          Node= 1
          do while(Layer(Node).ne.lay)
            Node= Node + 1
          enddo
          do j = 1,21
            thetadum(node) = thetar(node) + dfloat(j-1) *               &
     &                    (thetas(node)-thetar(node)) / 20.0d0
            call devries(node,thetadum,thetas,heacap,heacnd,fquartz,    &
     &        fclay,forg)
            write(hea,22) lay, comma, thetadum(node), comma,            &
     &                 heacap(node), comma, heacnd(node)
22          format(i4,a1,f8.5,2(a1,e14.5))
          enddo
      enddo

! --- close heat params file
      close (hea)

      return
      end

! ----------------------------------------------------------------------
      subroutine outsoilphys ()
! ----------------------------------------------------------------------
!     date               : March 2008
!     purpose            : Output of soil physical parameters
!        cofgen(1,node) = ores
!        cofgen(2,node) = osat
!        cofgen(3,node) = ksatfit
!        cofgen(4,node) = alfa
!        cofgen(5,node) = lexp
!        cofgen(6,node) = npar
!        cofgen(7,node) = 1.d0 - (1.d0 / npar)
!        cofgen(8,node) = dummy
!        cofgen(9,node) = h_enpr
!        cofgen(10,node)= ksatexm
! ---------------------------------------------------------------------
! --- global
      use variables
      implicit none 

! --- local variables ------------------
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      integer   getun,soi,lay,swfrostx,hh,numhead,node
      real(8)   watcon,moiscap,hconduc,relsat, dhx
      real(8)   dimocax,kx,thetax, FrArMtrx1,rfcpx
      parameter (numhead = 341)
      real(8)   hx(numhead)
! ---------------------------------------------------------------------
      comma = ',' 
      FrArMtrx1 = 1.d0
      rfcpx = 1.d0
      swfrostx = 0

! === open output file =================================================
      filnam = trim(pathwork)//'SoilPhysParam.csv'
      soi = getun (20,90)
      call fopens(soi,filnam,'new','del')
      filtext = 'soil physical parameters '
      call writehead (soi,1,filnam,filtext,project)

! --- write header
      write(soi,'(a)')                                                  &
     &    'layer,head(cm),theta(cm3.cm3),C(cm-1),RelSat(-),k(cm.d-1)'

! --- determine range of heads for which soil physical params are calculated
      if(swsophy.eq.0)then
         hx(1) = 0.0d0
         hx(2) = -1.0d-7
         hx(numhead) = -1.0d07
         do hh = 3,numhead-1
           hx(hh) = hx(hh-1)*1.1d0
         enddo
      end if

! --- soil physical params: calculate and write 
      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         if(swsophy.eq.1)then
            hx(1) = 0.0d0
            hx(2) = -1.0d-1
            hx(numhead) = sptab(1,node,1)
            dhx = 10.d0**((log10(-sptab(1,node,1))+1.d0)/(numhead-2))
            do hh = 3,numhead
              hx(hh) = hx(hh-1)* dhx
            enddo
         end if

         do hh = 1, numhead
            thetax = watcon(node,hx(hh),cofgen,swsophy,numtab,sptab,    &
     &               ientrytab)
            dimocax = moiscap(node,hx(hh),cofgen,dt,swsophy,numtab,     &
     &               sptab,ientrytab)
            kx = hconduc (node,hx(hh),thetax,cofgen,swfrostx,rfcpx,     &
     &                    swsophy,numtab,sptab,ientrytab)
!           in case of static macropores FrArMtrx < 1
            if(swmacro.eq.1)  kx = FrArMtrx1 * kx
!           write output
            relsat = (thetax-cofgen(1,Node)) /                          &
     &                 (cofgen(2,Node)-cofgen(1,Node))
            write (soi,22) lay, comma, hx(hh), comma, thetax, comma,    &
     &                     dimocax, comma, relsat, comma, kx
 22         format(i10,5(a,1pe15.7))
          enddo
      enddo

! --- close soil physical params file
      close (soi)

      return
      end

! ----------------------------------------------------------------------
      subroutine outtem (task,numnod,date,daynr,tem,daycum,tav,tebot,   &
     &    tsoil,tetop,outfil,pathwork,flheader,project)

! ----------------------------------------------------------------------
!     date               : November 2004
!     purpose            : Output of soil temperatures           
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer   numnod,daynr,tem,daycum,task
      real(8)   tav,tebot,tsoil(macp),tetop
      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      character(len=11) date
      logical   flheader
! --- local variables ------------------
      integer      getun,i, reclngth
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.tem'
      tem = getun (20,90)
!      reclngth = 36 + 7*numnod 
      reclngth = 50 + 7*numnod 
      open(unit=tem,file=filnam,status='unknown',recl=reclngth)
      filtext = 'soil temperature profiles (oC)'
      call writehead (tem,1,filnam,filtext,project)

! --- write header
      if (numnod.le.9) then
         write (tem,10) (i,i=1,numnod)
      else
         write (tem,11) (i,i=1,9), (i,i=10,numnod)
      endif

 10   format('*',/,                                                     &
     & t8,'Date,Day,Daycum,   Tav, Tetop',9(',    T',i1),               &
     &', TeBot')

 11   format('*',/,                                                     &
     & t8,'Date,Day,Daycum,   Tav, Tetop',9(',    T',i1),               &
     &1024(',   T',i2),', TeBot')

      write (tem,'(a11,a1,i3,a1,i6,1024(a1,f6.1:))') '    Initial'      &
     &      ,comma,daynr,comma,daycum,comma,tav,comma,tetop,            &
     &      (comma,tsoil(i),i=1,numnod),comma,tebot

      return

      case (2)

! === write actual soil temperature data ================================

! --- write header in case of new balance period
      if (flheader) write (tem,10) 

! --- write soil temperature profile
!     PWB: idem
!      write (tem,'(a11,a1,i3,a1,i6,<numnod+3>(a1,f6.1:))') date
!     &      ,comma,daynr,comma,daycum,comma,tav,comma,tetop,
!     &      (comma,tsoil(i),i=1,numnod),comma,tebot
      write (tem,'(a11,a1,i3,a1,i6,1024(a1,f6.1:))') date               &
     &      ,comma,daynr,comma,daycum,comma,tav,comma,tetop,            &
     &      (comma,tsoil(i),i=1,numnod),comma,tebot

      case default
         call fatalerr ('outtem', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine IrrigationOutput(task) 
! ----------------------------------------------------------------------
!     Date               : November 2004   
!     Purpose            : open and write irrigation output files 
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- global variables ------------------
      integer   task
! --- local variables ------------------
      integer   getun
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! === open output file and write header ===============================

      filnam = trim(pathwork)//trim(outfil)//'.irg'
      irg = getun (20,90)
      call fopens(irg,filnam,'new','del')
      filtext = 'irrigation data'
      call writehead (irg,1,filnam,filtext,project)

      write (irg,100)
 100  format ('*',/,                                                    &
     &  '*      Crop name        Date   Day   Day  Irrigation',         &
     &  '     Solute',/,                                                &
     &  '*                             Crop   Cum          cm',         &
     &  '     mg/cm3',/,                                                &
     &  '*<=============><==========><====><====><==========>           &
     &<=========>',/,                                                   &
     &  '       Crop name,       Date,  Day,  Day, Irrigation,',        &
     &  '    Solute')

      flIrg1Start = .false.

      return

      case (2)

! === write actual data ===============================

! --- write header in case of new balance period
      if (flheadirg) then
        write (irg,100)
       flheadirg = .false.
      endif

! --- write output record irrigation
      if (irrigevent.eq.1 .and. flCropCalendar                          &
     &                                   .and. .not. flCropHarvest) then
        write (irg,'(a16,a1,a11,a1,i4,a1,i6,a1,2x,f8.1,a1,e11.3)')      &
     &     cropname(icrop),comma,date,comma,daycrop,                    &
     &     comma,daycum,comma,gird,comma,cirr
      elseif (irrigevent.eq.1) then
        write (irg,'(t1,a,t18,a11,a1,i4,a1,i6,a1,2x,f8.1,a1,e11.3)')    &
     &     'bare soil,',date,comma,daycrop,                             &
     &     comma,daycum,comma,gird,comma,cirr

      elseif (irrigevent.eq.2) then
        write (irg,'(a16,a1,a11,a1,i4,a1,i6,a1,2x,f8.1,a1,e11.3)')      &
     &     cropname(icrop),comma,date,comma,daycrop,                    &
     &     comma,daycum,comma,gird,comma,cirr
      endif

      flIrrigationOutput = .false.

      return

      case (3)

! === close output files ===========================

! --- close irg file
      close (irg)

      case default
         call fatalerr ('IrrigationOutput', 'Illegal value for Task')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine SnowOutput(task) 
! ----------------------------------------------------------------------
!     Date               : December 2004   
!     Purpose            : open and write snow pack data
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer   getun
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.snw'
      snw = getun (20,90)
      call fopens(snw,filnam,'new','del')
      filtext = 'snow pack output data (cm/period)'
      call writehead (snw,1,filnam,filtext,project)

! --- write header
      write (snw,10)
 10   format ('*',/,                                                    &
     &   '    date,      dcum,  rainfall,  snowfall,snowstorage, ',     &
     &   'meltflux,sublimation')

      return

      case (2)

! === write actual soil temperature data ================================

! --- write header in case of new balance period
      if (flheader) write (snw,10) 

! --- write actual data
      write (snw,20) date,comma,daycum,comma,snrai,comma,gsnow,         &
     &               comma,ssnow,comma,melt,comma,subl
20    format (a11,a1,i6,1x,5(a1,f10.4))

      return

      case (3)

! === close output file ===========================

! --- close snw file
      close (snw)

      case default
         call fatalerr ('SnowOutput', 'Illegal value for Task')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine outblc(task,blc,outfil,pathwork,cgrai,cnrai,cgird,     &
     &     cnird,cqrot,cevap,volact,volini,z,dz,swmacro,nrlevs,swirfix, &
     &     schedule,numnod,snowinco,ssnow,cgsno,cmelt,caintc,csnrai,    &
     &     cqprai,ioutdat,t1900,outdat,tstart,project,pond,pondini,     &
     &     cqdrainin,cqdrainout,cinund,crunoff,cqtdo,cqtup,cqbotdo,     &
     &     cqbotup,csubl,crunon,CQMpExfMtx,                             &
     &     CQMpInfMtx,CQMpInTopPreDm1,CQMpInTopPreDm2,                  &
     &     CQMpInTopLatDm1,CQMpInTopLatDm2)
! ----------------------------------------------------------------------
!     date               : December 2004
!     purpose            : Write detailed overview of water balance
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

!     global
      integer   blc,numnod,ioutdat,swmacro,nrlevs,swirfix,schedule,task

      real(8)   cevap,cgird,cgrai,cnird,cnrai
      real(8)   cqrot,dz(macp),z(macp), pond, pondini
      real(8)   volact,volini,tstart,t1900,outdat(maout)
      real(8)   snowinco,ssnow,cgsno,cmelt,caintc,csnrai          
      real(8)   cqprai,cinund,crunoff,csubl,crunon
      real(8)   cqdrainout(nrlevs),cqdrainin(nrlevs)
      real(8)   cqtdo,cqtup,cqbotdo,cqbotup  
      real(8)   CQMpInTopPreDm1, CQMpInTopPreDm2
      real(8)   CQMpInTopLatDm1, CQMpInTopLatDm2

      character(len=16) outfil
      character(len=*)  pathwork
      character(len=80) project
      
!     local
      integer   getun,level

      real(8)   CQMpInfMtx, CQMpExfMtx, CQMpInTop,soilout
      real(8)   plantin,plantout,snowin,snowout,pondin,pondout,soilin

      character(len=300) filnam
      character(len=11)  datbegin,datend
      character(len=80)  filtext
! ----------------------------------------------------------------------

      select case (task)
      case (1)
! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.blc'
      blc = getun (50,90)
      call fopens(blc,filnam,'new','del')
      filtext='detailed overview of water balance components (cm)'
      call writehead (blc,1,filnam,filtext,project)

      return

      case (2)

! --- begin date of balance period
      if (ioutdat .eq. 2) then
        call dtdpst ('year-month-day',tstart+0.1d0,datbegin)
      else
        call dtdpst ('year-month-day',                                &
     &                                outdat(ioutdat-2)+1.1d0,datbegin)
      endif

! --- end date of balance period
      call dtdpst ('year-month-day',t1900-0.9d0,datend)

! --- initialize macropore variables
      if (swmacro.ne.1) then
         CQMpInfMtx = 0.0d0
         CQMpExfMtx = 0.0d0
         CQMpInTop = 0.0d0
      else if (swmacro.eq.1) then    
         CQMpInTop = CQMpInTopPreDm1 + CQMpInTopPreDm2 +                &
     &               CQMpInTopLatDm1 + CQMpInTopLatDm2
      endif

! --- write output record
      write (blc,20) datbegin,datend
      write (blc,22) (-z(numnod) + 0.5*dz(numnod))
      write (blc,40)      
      write (blc,41) snowinco,pondini,volini,ssnow,pond,volact 
      write (blc,42) cgrai
      write (blc,44) csnrai,cqprai,cnrai
      if (swirfix.eq.1 .or. schedule.eq.1) then
        write (blc,46) cgird
        write (blc,67) cnird,cnird
      endif
      write (blc,45) caintc
      write (blc,43) cgsno
      write (blc,47) cmelt,cmelt
      write (blc,68) csubl
      write (blc,48) cqrot
      write (blc,491) cevap
      write (blc,492) crunon,crunoff
      write (blc,50) cinund
      write (blc,52) cqtdo,cqtdo
      write (blc,53) cqtup,cqtup
      if (swmacro.eq.1)  write(blc,54) CQMpInTop, CQMpInfMtx, CQMpExfMtx

      if (nrlevs .ge. 1) then 
        write (blc,55)
        do level=1,nrlevs
          write (blc,56) level,cqdrainin(level),level, cqdrainout(level)
        enddo
      endif
      write (blc,61) cqbotup,cqbotdo

! --- sum
      plantin = cgrai+cgird
      snowin = snowinco+cgsno+csnrai
      pondin = pondini+cqprai+cnird+cmelt+cinund+cqtup+crunon
      soilin = volini+cqtdo+cqbotup+CQMpInfMtx
      do level = 1,nrlevs
        soilin = soilin + cqdrainin(level)
      enddo
      plantout = cnrai+caintc+cnird
      snowout = cmelt+ssnow+csubl
      pondout = pond+crunoff+cqtdo+cevap+CQMpInTop 
      soilout = volact+cqrot+cqtup+cqbotdo+CQMpExfMtx
      do level = 1,nrlevs
        soilout = soilout + cqdrainout(level)
      enddo

      write (blc,63) plantin,snowin,pondin,soilin,plantout,snowout,     &
     &       pondout,soilout
      write (blc,64) (ssnow-snowinco),(pond-pondini),(volact-volini)
      write (blc,65) (plantout-plantin),(snowout-snowin),               &
     &     (pondout-pondin),(soilout-soilin)

 20   format(/'Period',t20,':',t23,a11,' until  ',a11)
 22   format('Depth soil profile',t20,':',f8.2,' cm')
 40   FORMAT(49('='),'+',49('='),/,'INPUT',t50,'|',t52,'OUTPUT',/,t17,  &
     &'   PLANT    SNOW    POND    SOIL',t50,'|',t67,                   &
     &'   PLANT    SNOW    POND    SOIL',/,49('='),'+',49('='))
 41   FORMAT('Initially Present',t25,3f8.2,t50,'|',t52,                 &
     &       'Finally present',t75,3f8.2)
 42   FORMAT('Gross Rainfall',t17,f8.2,t50,'|')
 43   FORMAT('Snowfall',t25,f8.2,t50,'|')
 44   FORMAT('Nett Rainfall',t25,2f8.2,t50,'|',t52,'Nett Rainfall',     &
     &        t67,f8.2)
 45   FORMAT(t50,'|',t52,'Interception',t67,f8.2)
 46   FORMAT('Gross Irrigation',t17,f8.2,t50,'|')
 67   FORMAT('Nett Irrigation',t33,f8.2,t50,'|',t52,'Nett Irrigation',  &
     &        t67,f8.2)
 47   FORMAT('Snowmelt',t33,f8.2,t50,'|',t52,'Snowmelt',t75,f8.2)
 68   FORMAT(t50,'|',t52,'Sublimation',t75,f8.2)
 48   FORMAT(t50,'|',t52,'Plant Transpiration',t91,f8.2) 
 491  FORMAT(t50,'|',t52,'Soil Evaporation',t83,f8.2)
 492  FORMAT('Runon',t33,f8.2,t50,'|',t52,'Runoff',t83,f8.2)
 50   FORMAT('Inundation',t33,f8.2,t50,'|')
 52   FORMAT('Infiltr. Soil Surf.',t41,f8.2,t50,'|',t52,                &
     &       'Infiltr. Soil Surf.',t83,f8.2)
 53   FORMAT('Exfiltr. Soil Surf.',t33,f8.2,t50,'|',t52,                &
     &       'Exfiltr. Soil Surf.',t91,f8.2)
 54   FORMAT(t50,'|',t52,'Inflow top macroprs',t83,f8.2,/,              &
     &       'Infiltr. macropores',t41,f8.2,t50,'|',t52,                &
     &       'Exfiltr. macropores',t91,f8.2)
 55   FORMAT('Infiltr. subsurf.',t50,'|',t52,'Drainage')
 56   FORMAT('- system',i2,t41,f8.2,t50,'|',t52,'- system',i2,t91,f8.2)
 61   FORMAT('Upward seepage',t41,f8.2,t50,'|',t52,'Downward seepage',  &
     &        t91,f8.2,/,49('='),'+',49('='))
 63   FORMAT('Sum',t17,4f8.2,t50,'|',t52,'Sum',t67,4f8.2,/,49('='),'+', &
     &        49('='))
 64   FORMAT('Storage Change',t25,3f8.2)
 65   FORMAT('Balance Deviation',t18,f7.2,3f8.2,/,99('='),/)
  
      case default
         call fatalerr ('outbLC', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outafo (task,afo,outfil,pathwork,numnod,outper,period, &
     &  pond,gwl,theta,h,inq,inqrot,nrlevs,inqdra,ievap,ipeva,iptra,    &
     &  iruno,numlay,botcom,dz,thetas,cofgen,kdif,kdir,swafo,igrai,     &
     &  igird,inrai,inird,iintc,gc,lai,tav,tsoil,daycum,rd,cf,wbalance, &
     &  project,tstart,tend,swdiscrvert,numnodnew,dznew,SwMacro,Ssnow,  &
     &  igSnow,isnrai,iSubl,irunon, FrArMtrx,IPondBeg,ISsnowBeg,        &
     &  IThetaBeg,VlMpStDm1,VlMpStDm2,DiPoCp,WalevDm1,VlMpDm1,WaSrDm1,  &
     &  VlMpDm2,WaSrDm2,IQInTopPreDm1,IQInTopLatDm1,InQExcMtxDm1Cp,     &
     &  InQOutDrRapCp,IAvFrMpWlWtDm1,IQInTopPreDm2,IQInTopLatDm2,       &
     &  InQExcMtxDm2Cp,IAvFrMpWlWtDm2,IWaSrDm1Beg,IWaSrDm2Beg,          &
     &  CritDevMasBal,tcum,nod1lay,swsophy,numtab,sptab,ientrytab,      &
     &  out_etr,out_hum,out_rad,out_tmn,out_tmx,out_wet,out_win)
! ----------------------------------------------------------------------
!     Date               : 29-jan-2003
!     Purpose            : ANIMO/PEARL output: formatted hydrological data
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! -   global
      integer   swdiscrvert,numnodnew,task
      integer   ientrytab(macp,0:matabentries)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment

      real(8)   dznew(macp)
      integer   afo,botcom(maho),datea(6),daycum,numlay,numnod,nrlevs
      integer   period,swafo,nod1lay(maho)
      real(8)   gwl,inqdra(Madr,macp),pond,rd,cf,tstart,tend,timjan1
      real(8)   ievap,ipeva,iptra,iruno,inq(macp+1),inqrot(macp)
      real(8)   igrai,igird,inrai,inird,isnrai,gc,lai,tav,irunon,iintc
      real(8)   theta(macp),h(macp),thetas(macp),watcon,dz(macp)
      real(8)   cofgen(12,macp),kdif,kdir
      real(8)   dval(maho),tsoil(macp),wbalance,tcum,outper
      character(len=*)  outfil,pathwork
      character(len=80) project

      integer   SwMacro
      real(8)   FrArMtrx(MaCp),IPondBeg, Ssnow,ISsnowBeg,igSnow,iSubl
      real(8)   VlMpStDm1(macp),VlMpStDm2(macp),DiPoCp(macp)
      real(8)   WalevDm1,VlMpDm1,WaSrDm1,VlMpDm2,WaSrDm2
      real(8)   IQInTopPreDm1,IQInTopLatDm1,IQInTopPreDm2,IQInTopLatDm2
      real(8)   InQOutDrRapCp(macp)
      real(8)   InQExcMtxDm1Cp(macp),IAvFrMpWlWtDm1(macp)
      real(8)   InQExcMtxDm2Cp(macp),IAvFrMpWlWtDm2(macp)
      real(8)   IThetaBeg(MaCp),CritDevMasBal
      real(8)   IWaSrDm1Beg, IWaSrDm2Beg
      integer   swsophy,numtab(macp)
      real(8)   sptab(5,macp,matab)
! -   CFO file for PEARL        
      real(4)   out_etr,out_hum,out_rad,out_tmn,out_tmx,out_wet,out_win

! -   local
      integer   getun,ftype,node,bruny,eruny,brund,erund,lay,level
      integer   swop,botcomNew(maho)
      real(4)   fsec
      real(8)   hNew(macp),inqNew(macp+1)
      real(8)   inqdraNew(Madr,macp),thetaNew(macp),inqrotNew(macp)
      real(8)   TsoilNew(0:macp),tsoili(0:macp)
      real(8)   DiPoCpNew(macp), IAvFrMpWlWtDm1New(macp)
      real(8)   IAvFrMpWlWtDm2New(macp), InQExcMtxDm1CpNew(macp)
      real(8)   InQExcMtxDm2CpNew(macp), InQOutDrRapCpNew(macp)
      real(8)   IThetaBegNew(MaCp)
      real(8)   VlMpStDm1New(macp), VlMpStDm2New(macp)
      character(len=80)  filtext
      character(len=300) filnam
      character(len=4)   afoext
      logical   FlOpenFileDev

!     save local variables
      save      brund, FlOpenFileDev, swop


! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initial output ====================================================

! -   verify reduced grid for output
      if(SwDiscrvert.eq.1) then
         call checkDiscrVert(numnod,dz,numnodNew,dzNew)
      endif

! -   set switch Swop
      Swop = 1
! -   BFO and CFO file for PEARL        
      if (Swafo.ge.2 .and. SwMacro.eq.1) Swop = 2

! -   convert vertical discretization: initial values
      do node = 1,numnod
         tsoili(node) = tsoil(node)
      enddo
      call convertdiscrvert(SwDiscrvert,1,Swop,nrlevs,numlay,           &
     &         botcom,numnod,dz,h,theta,inq,inqrot,inqdra,              &
     &         IThetaBeg,Tsoili,cofgen,                                 &
     &         numnodNew,dzNew,botcomNew,hNew,thetaNew,inqNew,          &
     &         inqrotNew,inqdraNew,IThetaBegNew,TsoilNew,               &
     &         DiPoCp,FrArMtrx,IAvFrMpWlWtDm1,IAvFrMpWlWtDm2,           &
     &         InQExcMtxDm1Cp,InQExcMtxDm2Cp,InQOutDrRapCp,VlMpStDm1,   &
     &         VlMpStDm2,DiPoCpNew,IAvFrMpWlWtDm1New,IAvFrMpWlWtDm2New, &
     &         InQExcMtxDm1CpNew,InQExcMtxDm2CpNew,InQOutDrRapCpNew,    &
     &         VlMpStDm1New,VlMpStDm2New,swsophy,numtab,sptab,ientrytab)
! -   open file
      if (swafo.eq.1) then
         afoext = '.afo'
      elseif (swafo.eq.2) then
         afoext = '.bfo'
      elseif (swafo.eq.3) then
         afoext = '.cfo'
      endif
      filnam = trim(pathwork)//trim(outfil)//afoext
      afo = getun (20,90)
      call fopens(afo,filnam,'new','del')


! ---   write initial part
      if (swafo.ge.2) then
         ftype = 1
         if (swafo.eq.2) then
            filtext = 'formatted hydrological data; BFO version'
         else if (swafo.eq.3) then
            filtext = 'formatted hydrological data; CFO version'
         endif
         call writehead(afo,ftype,filnam,filtext,project)
      endif

      call dtdpar (tstart,datea,fsec)
      bruny = datea(1)
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      brund = nint (tstart - timjan1 + 1.0d0)
      call dtdpar (tend,datea,fsec)
      eruny = datea(1)
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      erund = nint (tend - timjan1 + 1.0d0)
      if (swafo.eq.1) then
         write (afo,4000) bruny,eruny,                                  &
     &                 real(brund-1),real(erund),real(period)
      else if (swafo.ge.2) then
         write (afo,4001) swop
         write (afo,4000) bruny,eruny,real(brund-1),real(erund)
      endif

      write (afo,4010) numnodNew,numlay,nrlevs
      write (afo,4010) (botcomNew(lay),lay=1,numlay)
!  -    ThetaS should be known for new soil layers; this works but is not very nice
      write (afo,4020) (real(thetas(botcom(lay))), lay=1,numlay)

      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         dval(lay) = watcon(node,-100.d0,cofgen,swsophy,numtab,sptab,   &
     &                      ientrytab)
      enddo
      write (afo,4020) (real(dval(lay)),lay=1,numlay) 

      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         dval(lay) = watcon(node,-15849.0d0,cofgen,swsophy,numtab,sptab,&
     &               ientrytab)
      enddo
      write (afo,4020) (real(dval(lay)),lay=1,numlay) 

      write (afo,4020) (0.01*real(dzNew(node)),node=1,numnodNew)

      if (swafo.ge.2 .and. swop.eq.2) then
         write (afo,4020)                                               &
     &         (0.01*real(VlMpStDm1New(node)),node=1,numnodNew)
         write (afo,4020)                                               &
     &         (0.01*real(VlMpStDm2New(node)),node=1,numnodNew)
         write (afo,4020) (0.01*real(DiPoCpNew(node)),node=1,numnodNew)
      endif

      write (afo,4020) (real(thetaNew(node)),node=1,numnodNew)
      write (afo,4020) -0.01*real(gwl),0.01*real(pond)

      if (swafo.ge.2) then
         write (afo,4020) 0.01*real(Ssnow)
         write (afo,4020) (real(tsoilNew(node)),node = 1,numnodNew) 
      endif

      if (swafo.ge.2 .and. swop.eq.2) then
         write (afo,4020) 0.01*real(WaLevDm1), 0.01*real(VlMpDm1),      &
     &       0.01*real(WaSrDm1), 0.01*real(VlMpDm2), 0.01*real(WaSrDm2)
      endif

4000  format (2(1x,i8),3(1x,f8.0),1x,i8) 
4001  format (1x,i8)
4010  format (5(1x,i8))
4020  format (8(1x,f10.6))

      FlOpenFileDev = .false.

      return
      
      case (2)

! === output during simulation ===============================================

! -     convert vertical discretization: dynamic part
        do node = 1,numnod
          tsoili(node) = tsoil(node)
        enddo
        call ConvertDiscrVert(SwDiscrvert,2,Swop,nrlevs,numlay,         &
     &         botcom,numnod,dz,h,theta,inq,inqrot,inqdra,              &
     &         IThetaBeg,Tsoili,cofgen,                                 &
     &         numnodNew,dzNew,botcomNew,hNew,thetaNew,inqNew,          &
     &         inqrotNew,inqdraNew,IThetaBegNew,TsoilNew,               &
     &         DiPoCp,FrArMtrx,IAvFrMpWlWtDm1,IAvFrMpWlWtDm2,           &
     &         InQExcMtxDm1Cp,InQExcMtxDm2Cp,InQOutDrRapCp,VlMpStDm1,   &
     &         VlMpStDm2,DiPoCpNew,IAvFrMpWlWtDm1New,IAvFrMpWlWtDm2New, &
     &         InQExcMtxDm1CpNew,InQExcMtxDm2CpNew,InQOutDrRapCpNew,    &
     &         VlMpStDm1New,VlMpStDm2New,swsophy,numtab,sptab,ientrytab)

        if (swafo.eq.1) then
          write (afo,30) real(tcum),                                    &
     &    0.01*real((igrai+igird)/outper),                              &
     &    0.01*real(iintc/outper),                                      &
     &    0.01*real(ievap/outper),                                      &
     &    0.0,                                                          &
     &    0.01*real(ipeva/outper),0.01*real(iptra/outper),              &
     &    0.01*real(iruno/outper),                                      &
     &    -0.01*real(gwl),0.01*real(pond)
        elseif (swafo.ge.2) then
!          write (afo,*) real(t), real(outper),
          write (afo,33) (brund*1.0d0-1.0d0+tcum), real(outper),        &
     &        0.01*real((igrai+isnrai)/outper),                         &
     &        0.01*real(igsnow/outper),                                 &
     &        0.01*real(igird/outper),                                  &
     &        0.01*real( (iintc-(igird-inird))/outper ),                &
     &        0.01*real( (igird-inird)/outper ),                        &
     &        0.01*real(isubl/outper),                                  &
     &        0.01*real(ievap/outper),                                  &
     &        0.0,                                                      &
     &        0.01*real(ipeva/outper),                                  &
     &        0.01*real(iptra/outper),                                  &
     &        0.01*real(irunon/outper),                                 &
     &        0.01*real(iruno/outper),                                  &
     &        -0.01*real(gwl), 0.01*real(pond),                         &
     &        0.01*real(ssnow), 0.01*real(wbalance)
        endif
        write (afo,40) (real(hNew(node))           ,node=1,numnodNew)
        write (afo,50) (real(thetaNew(node))       ,node=1,numnodNew)
        write (afo,50) (0.01*real(inqrotNew(node)/outper),node=1,       &
     &   numnodNew)
        write (afo,50) (-0.01*real(inqNew(node)/outper),                &
     &   node=1,numnodNew+1)
        if (nrlevs.gt.0) then
          do level=1,nrlevs
            write (afo,50) (0.01*real(inqdraNew(level,node)/outper),    &
     &        node=1,numnodNew)
          enddo
        end if
        if (swafo.ge.2) then
!ckro-20141112   exponential relation between soil cover and lai
          gc = min( (1.0d0-exp(-1.0d0*kdif*kdir*lai)), 1.0d0)
          write (afo,'(4(1x,f10.4))') real(gc),real(lai),               &
     &     abs(real(0.01d0*rd)),real(cf)
! -   CFO file for PEARL        
         if (swafo.eq.2) then
            write (afo,'(1x,f7.2)') real(tav)
         else if (swafo.eq.3) then
! -   if swafo = 3, CFO file with extra output concerning meteo:
!     RAD(kJ.m-2), Tmin(grC), Tmax(grC), HUM(kPA), WIND(m.s-1), ETref (m.d-1), WET (d)
!     and missing values voor WET: -1.0 if SWRAIN =/= 2      
            write (afo,'(3x,f8.1,f9.2,2(2x,f9.2),f9.2,3x,2(3x,f8.5))')  &
     &     out_rad, out_tmn, out_tmx, out_hum, out_win, out_etr, out_wet
         endif
          write (afo,4020) (real(tsoilNew(node)),node = 1,numnodNew) 
        endif

        if (swafo.ge.2 .and. swop.eq.2) then
          write (AFO,4020) 0.01*real(WaLevDm1), 0.01*real(VlMpDm1),     &
     &             0.01*real(WaSrDm1), 0.01*real(IQInTopPreDm1/outper), &
     &             0.01*real(IQInTopLatDm1/outper)
          write (AFO,4020) (0.01*real(InQExcMtxDm1CpNew(node)/outper),  &
     &                      node=1,numnodNew)
          write (AFO,4020) (0.01*real(InQOutDrRapCpNew(node)/outper),   &
     &                      node=1,numnodNew)
          write (AFO,4020) (real(IAvFrMpWlWtDm1New(node)/outper),       &
     &                      node=1,numnodNew)
          write (AFO,4020) 0.01*real(VlMpDm2), 0.01*real(WaSrDm2),      &
     &                     0.01*real(IQInTopPreDm2/outper),             &
     &                     0.01*real(IQInTopLatDm2/outper)
          write (AFO,4020) (0.01*real(InQExcMtxDm2CpNew(node)/outper),  &
     &                      node=1,numnodNew)
          write (AFO,4020) (real(IAvFrMpWlWtDm2New(node)/outper),       &
     &                      node=1,numnodNew)
        endif

        if (swafo.ge.2 .and. CritDevMasBal.gt.1.d-30) then
! --- Checking of mass balances of sub systems per period OutPer
           call checkmassbal (DayCum,nrlevs,NumNodNew,                  &
     &              SwMacro,outfil,pathwork,FlOpenFileDev,DZNew,        &
     &              CritDevMasBal,ievap,igird,igrai,                    &
     &              igSnow,inird,inrai,inqdraNew,InQExcMtxDm1CpNew,     &
     &              InQExcMtxDm2CpNew,inqNew,InQOutDrRapCpNew,inqrotNew,&
     &              IPondBeg,IQInTopPreDm1,IQInTopLatDm1,IQInTopPreDm2, &
     &              IQInTopLatDm2,ISsnowBeg,IThetaBegNew,               &
     &              iruno,irunon,isnrai,iSubl,pond,Ssnow,thetaNew,      &
     &              IWaSrDm1Beg,IWaSrDm2Beg,WaSrDm1,WaSrDm2)
        end if

30    format (1x,f6.0,7(1x,f9.6),3(1x,f7.4))
33    format (1(1x,f12.6),12(1x,f9.6),4(1x,f9.4))
40    format (8(1x,e10.3e2))
50    format (10(1x,f9.6))
!30    format (1x,f6.0,7(1x,f9.6),3(1x,f7.4))   ! 7
!33    format (1(1x,f12.6),12(1x,f9.6),4(1x,f9.4))
!40    format (8(1x,e10.3e2))  ! 9
!50    format (10(1x,f9.6))    ! 8

      case default
         call fatalerr ('outafo', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outaun (task,aun,outfil,pathwork,numnod,outper,period, &
     &  pond,gwl,theta,h,inq,inqrot,nrlevs,inqdra,ievap,ipeva,iptra,    &
     &  iruno,numlay,botcom,dz,thetas,cofgen,kdif,kdir,swaun,igrai,     &
     &  igird,inrai,inird,iintc,gc,lai,tav,tsoil,daycum,rd,cf,wbalance, &
     &  project,tstart,tend,swdiscrvert,numnodnew,dznew,SwAfo,SwMacro,  &
     &  Ssnow,igSnow,isnrai,iSubl,irunon,FrArMtrx,IPondBeg,             &
     &  ISsnowBeg,IThetaBeg,VlMpStDm1,VlMpStDm2,DiPoCp,WalevDm1,VlMpDm1,&
     &  WaSrDm1,VlMpDm2,WaSrDm2,IQInTopPreDm1,IQInTopLatDm1,            &
     &  InQExcMtxDm1Cp,InQOutDrRapCp,IAvFrMpWlWtDm1,IQInTopPreDm2,      &
     &  IQInTopLatDm2,InQExcMtxDm2Cp,IAvFrMpWlWtDm2,IWaSrDm1Beg,        &
     &  IWaSrDm2Beg,CritDevMasBal,tcum,nod1lay,swsophy,numtab,sptab,    &
     &  ientrytab)
! ----------------------------------------------------------------------
!     Date               : 29-jan-2003
!     Purpose            : ANIMO/PEARL output: unformatted hydrological data
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! -   global
      integer   swdiscrvert,numnodnew,task
      integer   ientrytab(macp,0:matabentries)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment

      real(8)   dznew(macp)
      integer   aun,botcom(maho),datea(6),daycum,numlay,numnod,nrlevs
      integer   period,swaun,nod1lay(maho)
      real(8)   gwl,inqdra(Madr,macp),pond,rd,cf,tstart,tend,timjan1
      real(8)   ievap,ipeva,iptra,iruno,inq(macp+1),inqrot(macp)
      real(8)   igrai,igird,inrai,inird,isnrai,gc,lai,tav,irunon,iintc
      real(8)   theta(macp),h(macp),thetas(macp),watcon,dz(macp)
      real(8)   cofgen(12,macp),kdif,kdir,tcum,outper
      real(8)   dval(maho),tsoil(macp),wbalance,CritDevMasBal
      character(len=*)  outfil,pathwork
      character(len=80) project

      integer   SwAfo, SwMacro
      real(8)   FrArMtrx(MaCp),IPondBeg, Ssnow,ISsnowBeg,igSnow,iSubl
      real(8)   VlMpStDm1(macp),VlMpStDm2(macp),DiPoCp(macp)
      real(8)   WalevDm1,VlMpDm1,WaSrDm1,VlMpDm2,WaSrDm2
      real(8)   IQInTopPreDm1,IQInTopLatDm1,IQInTopPreDm2,IQInTopLatDm2
      real(8)   InQOutDrRapCp(macp), IThetaBeg(MaCp)
      real(8)   InQExcMtxDm1Cp(macp),IAvFrMpWlWtDm1(macp)
      real(8)   InQExcMtxDm2Cp(macp),IAvFrMpWlWtDm2(macp)
      real(8)   IWaSrDm1Beg, IWaSrDm2Beg

      integer   swsophy,numtab(macp)
      real(8)   sptab(5,macp,matab)

! -   local
      integer   swop
      integer   getun,ftype,node,bruny,eruny,brund,erund,lay,level
      integer   botcomNew(maho)
      real(4)   fsec
      real(8)   hNew(macp),inqNew(macp+1)
      real(8)   inqdraNew(Madr,macp),thetaNew(macp),inqrotNew(macp)
      real(8)   TsoilNew(0:macp),tsoili(0:macp)
      real(8)   DiPoCpNew(macp), IAvFrMpWlWtDm1New(macp)
      real(8)   IAvFrMpWlWtDm2New(macp), InQExcMtxDm1CpNew(macp)
      real(8)   InQExcMtxDm2CpNew(macp), InQOutDrRapCpNew(macp)
      real(8)   IThetaBegNew(MaCp),VlMpStDm1New(macp),VlMpStDm2New(macp)
      character(len=80)  filtext
      character(len=300) filnam
      character(len=4)   aunext

      logical   FlOpenFileDev

!     save local variables
      save      brund, FlOpenFileDev, swop

! ----------------------------------------------------------------------
      
      select case (task)
      case (1)

! === output initial condition ==========================================

! -   verify reduced grid for output
      if(SwDiscrvert.eq.1) then
         call checkDiscrVert(numnod,dz,numnodNew,dzNew)
      endif

! -   set switch Swop
      Swop = 1
      if (Swaun.eq.2 .and. SwMacro.eq.1) Swop = 2

! -     convert vertical discretization: initial values
      do node = 1,numnod
        tsoili(node) = tsoil(node)
      enddo
      call convertdiscrvert(SwDiscrvert,1,Swop,nrlevs,numlay,           &
     &         botcom,numnod,dz,h,theta,inq,inqrot,inqdra,              &
     &         IThetaBeg,Tsoili,cofgen,                                 &
     &         numnodNew,dzNew,botcomNew,hNew,thetaNew,inqNew,          &
     &         inqrotNew,inqdraNew,IThetaBegNew,TsoilNew,               &
     &         DiPoCp,FrArMtrx,IAvFrMpWlWtDm1,IAvFrMpWlWtDm2,           &
     &         InQExcMtxDm1Cp,InQExcMtxDm2Cp,InQOutDrRapCp,VlMpStDm1,   &
     &         VlMpStDm2,DiPoCpNew,IAvFrMpWlWtDm1New,IAvFrMpWlWtDm2New, &
     &         InQExcMtxDm1CpNew,InQExcMtxDm2CpNew,InQOutDrRapCpNew,    &
     &         VlMpStDm1New,VlMpStDm2New,swsophy,numtab,sptab,ientrytab)

! -   open file
      if (swaun.eq.1) then
        aunext = '.aun'
      elseif (swaun.eq.2) then
        aunext = '.bun'
      endif
      filnam = trim(pathwork)//trim(outfil)//aunext
      aun = getun (20,90)
      open(unit=aun,file=filnam,status='unknown',form='unformatted')

! --- write initial part
      if (swaun.eq.2) then
          ftype = 2
          filtext = 'unformatted hydrological data'
          call writehead(aun,ftype,filnam,filtext,project)
      endif

      call dtdpar (tstart,datea,fsec)
      bruny = datea(1)
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      brund = nint (tstart - timjan1 + 1.0d0)
      call dtdpar (tend,datea,fsec)
      eruny = datea(1)
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      erund = nint (tend - timjan1 + 1.0d0)
      if (swaun.eq.1) then
        write (aun) bruny,eruny,                                        &
     &                   real(brund-1),real(erund),real(period)
      else if (swaun.eq.2) then
        write (aun) swop
        write (aun) bruny,eruny,real(brund-1),real(erund)
      endif

      write (aun) numnodNew,numlay,nrlevs
      write (aun) (botcomNew(lay),lay=1,numlay)
!  -  ThetaS should be known for new soil layers; this works but is not very nice
      write (aun) (real(thetas(botcom(lay))), lay=1,numlay)

      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         dval(lay) = watcon(node,-100.d0,cofgen,swsophy,numtab,sptab,   &
     &               ientrytab)
      enddo
      write (aun) (real(dval(lay)),lay=1,numlay) 

      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         dval(lay) = watcon(node,-15849.0d0,cofgen,swsophy,numtab,sptab,&
     &               ientrytab)
      enddo
      write (aun) (real(dval(lay)),lay=1,numlay) 

      write (aun) (0.01*real(dzNew(node)),node=1,numnodNew)

      if (swaun.eq.2 .and. swop.eq.2) then
        write (aun) (0.01*real(VlMpStDm1New(node)),node=1,numnodNew)
        write (aun) (0.01*real(VlMpStDm2New(node)),node=1,numnodNew)
        write (aun) (0.01*real(DiPoCpNew(node)),node=1,numnodNew)
      endif

      write (aun) (real(thetaNew(node)),node=1,numnodNew)
      write (aun) -0.01*real(gwl),0.01*real(pond)

      if (swaun.eq.2) then
        write (aun) 0.01*real(Ssnow)
        write (aun) (real(tsoilNew(node)),node = 1,numnodNew) 
      endif

      if (swaun.eq.2 .and. swop.eq.2) then
        write (aun) 0.01*real(WaLevDm1), 0.01*real(VlMpDm1),            &
     &        0.01*real(WaSrDm1), 0.01*real(VlMpDm2), 0.01*real(WaSrDm2)
      endif

      FlOpenFileDev = .false.

      return

      case (2)

! === output during simulation ===============================================

! -   convert vertical discretization: dynamic part
      do node = 1,numnod
        tsoili(node) = tsoil(node)
      enddo
      call convertdiscrvert(SwDiscrvert,2,Swop,nrlevs,numlay,           &
     &         botcom,numnod,dz,h,theta,inq,inqrot,inqdra,              &
     &         IThetaBeg,Tsoili,cofgen,                                 &
     &         numnodNew,dzNew,botcomNew,hNew,thetaNew,inqNew,          &
     &         inqrotNew,inqdraNew,IThetaBegNew,TsoilNew,               &
     &         DiPoCp,FrArMtrx,IAvFrMpWlWtDm1,IAvFrMpWlWtDm2,           &
     &         InQExcMtxDm1Cp,InQExcMtxDm2Cp,InQOutDrRapCp,VlMpStDm1,   &
     &         VlMpStDm2,DiPoCpNew,IAvFrMpWlWtDm1New,IAvFrMpWlWtDm2New, &
     &         InQExcMtxDm1CpNew,InQExcMtxDm2CpNew,InQOutDrRapCpNew,    &
     &         VlMpStDm1New,VlMpStDm2New,swsophy,numtab,sptab,ientrytab)

      if (swaun.eq.1) then
        write (aun) real(tcum),                                         &
     &    0.01*real((igrai+igird)/outper),                              &
     &    0.01*real(iintc/outper),                                      &
     &    0.01*real(ievap/outper),                                      &
     &    0.0,                                                          &
     &    0.01*real(ipeva/outper),                                      &
     &    0.01*real(iptra/outper),                                      &
     &    0.01*real(iruno/outper),                                      &
     &    -0.01*real(gwl),0.01*real(pond)
      elseif (swaun.eq.2) then
        write (aun) real(brund*1.0d0-1.0d0+tcum), real(outper),         &
     &        0.01*real((igrai+isnrai)/outper),                         &
     &        0.01*real(igsnow/outper),                                 &
     &        0.01*real(igird/outper),                                  &
     &        0.01*real( (iintc-(igird-inird))/outper ),                &
     &        0.01*real( (igird-inird)/outper ),                        &
     &        0.01*real(isubl/outper),                                  &
     &        0.01*real(ievap/outper),                                  &
     &        0.0,                                                      &
     &        0.01*real(ipeva/outper),                                  &
     &        0.01*real(iptra/outper),                                  &
     &        0.01*real(irunon/outper),                                 &
     &        0.01*real(iruno/outper),                                  &
     &        -0.01*real(gwl), 0.01*real(pond),                         &
     &        0.01*real(ssnow), 0.01*real(wbalance)
      endif
      write (aun) (real(hNew(node))                 ,node=1,numnodNew)
      write (aun) (real(thetaNew(node))             ,node=1,numnodNew)
      write (aun) (0.01*real(inqrotNew(node)/outper),node=1,numnodNew)
      write (aun) (-0.01*real(inqNew(node)/outper),node=1,numnodNew+1)
      if (nrlevs.gt.0) then
        do level=1,nrlevs
          write (aun) (0.01*real(inqdraNew(level,node)/outper),         &
     &        node=1,numnodNew)
        enddo
      end if
      if (swaun.eq.2) then
!ckro-20141112   exponential relation between soil cover and lai
        gc = min( (1.0d0-exp(-1.0d0*kdif*kdir*lai)), 1.0d0)
        write (aun) real(gc),real(lai),abs(real(0.01d0*rd)),real(cf)
        write (aun) real(tav)
        write (aun) (real(tsoilNew(node)),node = 1,numnodNew) 
      endif
      if (swaun.eq.2 .and. swop.eq.2) then
        write (AUN) 0.01*real(WaLevDm1), 0.01*real(VlMpDm1),            &
     &             0.01*real(WaSrDm1), 0.01*real(IQInTopPreDm1/outper), &
     &             0.01*real(IQInTopLatDm1/outper)
        write (AUN) (0.01*real(InQExcMtxDm1CpNew(node)/outper),         &
     &                      node=1,numnodNew)
        write (AUN) (0.01*real(InQOutDrRapCpNew(node)/outper),          &
     &                      node=1,numnodNew)
        write (AUN) (real(IAvFrMpWlWtDm1New(node)/outper),              &
     &                      node=1,numnodNew)
        write (AUN) 0.01*real(VlMpDm2), 0.01*real(WaSrDm2),             &
     &                     0.01*real(IQInTopPreDm2/outper),             &
     &                     0.01*real(IQInTopLatDm2/outper)
        write (AUN) (0.01*real(InQExcMtxDm2CpNew(node)/outper),         &
     &                      node=1,numnodNew)   
        write (AUN) (real(IAvFrMpWlWtDm2New(node)/outper),              &
     &                      node=1,numnodNew)
      end if

      if ((swaun.ge.1.or.swafo.ge.1) .and. CritDevMasBal.gt.1.d-30) then
 
! --- Checking of mass balances of sub systems per period OutPer
        call checkmassbal (DayCum,nrlevs,NumNodNew,                     &
     &              SwMacro,outfil,pathwork,FlOpenFileDev,DZNew,        &
     &              CritDevMasBal,ievap,igird,igrai,                    &
     &              igSnow,inird,inrai,inqdraNew,InQExcMtxDm1CpNew,     &
     &              InQExcMtxDm2CpNew,inqNew,InQOutDrRapCpNew,inqrotNew,&
     &              IPondBeg,IQInTopPreDm1,IQInTopLatDm1,IQInTopPreDm2, &
     &              IQInTopLatDm2,ISsnowBeg,IThetaBegNew,               &
     &              iruno,irunon,isnrai,iSubl,pond,Ssnow,thetaNew,      &
     &              IWaSrDm1Beg,IWaSrDm2Beg,WaSrDm1,WaSrDm2)
 
      endif

      case default
         call fatalerr ('outaun', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine WriteSwapOk(Project)
!-----------------------------------------------------------------------
!      Date    : December 2004
!      Purpose : Create file Swap.ok at end of simulation
!-----------------------------------------------------------------------
      implicit none
      
      integer   cexf,getun
      character(len=80) project

!     create file Swap.ok to let environment programs verify termination
      cexf = getun (20,90)
      call fopens(cexf,'Swap.ok','new','del')
      call writehead (cexf,1,'Swap.ok',                                 &
     &  'this header only: simulation succesfully terminated',project)
      close (cexf)

      return
      end


! ----------------------------------------------------------------------
      subroutine SurfaceWaterOutput(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : open and write surface water output files 
! ----------------------------------------------------------------------

      use variables
      implicit none

      integer task

      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  drf file
      if (swdrf.eq.1) call outdrf (task,outfil,drf,pathwork,daynr,date, &
     &  nrpri,nrlevs,cqdrain,cqdrd,crunoff,cQMpOutDrRap,flheader)

      if (swswb.eq.1) call outswb(task,outfil,vtair,pathwork,daynr,     &
     & daycum,hbweir,overfl,gwl,pond,wlstar,wls,swstini,swst,cqdrd,     &
     & crunoff,cQMpOutDrRap,cwsupp,swb,numadj,hwlman,cwout,swsec,swman, &
     & nmper,impend,imper,project,logf,swscre,date,t1900,t,outper,iyear)

      return

      case (2)

! === write actual data ===============================

! --  drf file
      if (swdrf.eq.1) call outdrf (task,outfil,drf,pathwork,daynr,date, &
     &  nrpri,nrlevs,cqdrain,cqdrd,crunoff,cQMpOutDrRap,flheader)

      if (swswb.eq.1) call outswb(task,outfil,vtair,pathwork,daynr,     &
     & daycum,hbweir,overfl,gwl,pond,wlstar,wls,swstini,swst,cqdrd,     &
     & crunoff,cQMpOutDrRap,cwsupp,swb,numadj,hwlman,cwout,swsec,swman, &
     & nmper,impend,imper,project,logf,swscre,date,t1900,t,outper,iyear)

      return

      case (3)

! === close output files ===========================

      close (drf)
      close (swb)

      case default
         call fatalerr ('SurfaceWaterOutput', 'Illegal value for Task')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine outdrf (task,outfil,drf,pathwork,daynr,date,           &
     &  nrpri,nrlevs,cqdrain,cqdrd,crunoff,cQMpOutDrRap,flheader)
! ----------------------------------------------------------------------
!     Date               : 10/6/99
!     Purpose            : write drainage fluxes, surface runoff, rapid
!                          drainage to  OUTNAM.DRF file
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer   drf,daynr,nrpri,nrlevs,task
      real(8)   cqdrain(Madr),cqdrd,crunoff,cQMpOutDrRap
      character(len=16) outfil
      character(len=80) pathwork
      character(len=11) date
      logical   flheader

! --- local
      integer   level,getun
      real(8)   c1qdrain(Madr),c1qdrd,c1runoff,c1qdrar
      character(len=80) filnam
      character(len=1)  comma

      save
      comma = ','
! ----------------------------------------------------------------------

      select case (task)
      case (1)
! === open output file and write header ================================

! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.drf'
      drf = getun (20,90)
      call fopens(drf,filnam,'new','del')

! --- write header
      if (nrpri .eq. 1) then
        write (DRF,100) nrlevs
      else
        write (DRF,110) nrlevs
      endif
      write (DRF,120)

! --- Output format
 100    format (                                                        &
     &'* CUMULATIVE DRAINAGE FLUXES, SURFACE RUNOFF AND RAPID',         &
     &' DRAINAGE:',/,'*',/,                                             &
     &'* Total number of levels for drainage fluxes :',I3,/,'*',/       &
     &'* First level (primary system) NOT included in sw-reservoir',/,  &
     &'*',/)

 110    format (                                                        &
     &'* CUMULATIVE DRAINAGE FLUXES, SURFACE RUNOFF AND RAPID',         &
     &' DRAINAGE:',/,'*',/,                                             &
     &'* Total number of levels for drainage fluxes :',I3,/,'*',/       &
     &'* First level (primary system) IS included in sw-reservoir',/,   &
     &'*')

 120    format (                                                        &
     &'* Meaning of symbols:',/,                                        &
     &'* - CQDinc(1-5): Drainage fluxes all level (1-5)',/,             &
     &'* - CQDRDinc   : Total of drainage fluxes into surface',         &
     &            ' water reservoir (SWSRF>=1)',/,                      &
     &'* - CRUNOFFinc : Surface runoff (increment can be <0.0)',/,      &
     &'* - CQDRARinc  : Rapid drainage (increment always > = 0.0)',/,   &
     &'* - CQDcum(1-5): Drainage fluxes all level (1-5)',/,             &
     &'* - CQDRDcum   : Total of drainage fluxes into surface',         &
     &            ' water reservoir (SWSRF>=1)',/,                      &
     &'* - CRUNOFFcum : Surface runoff (increment can be <0.0)',/,      &
     &'* - CQDRARcum  : Rapid drainage (increment always > = 0.0)',/,   &
     &'* - cum        : cumulative value',/,                            &
     &'* - inc        : incremental value',/,'*',/,                     &
     &'     Date, Daynr,CQDinc1,CQDinc2,CQDinc3,CQDinc4,CQDinc5,',      &
     &'CQDRDinc,CRUNOFFinc,CQDRARinc,',                                 &
     &'CQDcum1,CQDcum2,CQDcum3,CQDcum4,CQDcum5,',                       &
     &'CQDRDcum,CRUNOFFcum,CQDRARcum')
 
      do level=1,nrlevs
        c1qdrain(level)=0.0d0
      enddo
      c1qdrd = 0.0d0
      c1runoff = 0.0d0
      c1qdrar = 0.0d0

      return

      case (2)

! === output of data ===========================================================

! --- write header and reset cumulative fluxes in case of new balance period
      if (flheader) then
        if (nrpri .eq. 1) then
          write (DRF,100) nrlevs
        else
          write (DRF,110) nrlevs
        endif
        write (DRF,120)

        do level=1,nrlevs
          c1qdrain(level)=0.0d0
        enddo
        c1qdrd = 0.0d0
        c1runoff = 0.0d0
        c1qdrar = 0.0d0
      endif

! --- write output record
      write (DRF,20) date,comma,daynr,                                  &
     &  comma,(cqdrain(1)-c1qdrain(1)),                                 &
     &  comma,(cqdrain(2)-c1qdrain(2)),                                 &
     &  comma,(cqdrain(3)-c1qdrain(3)),                                 &
     &  comma,(cqdrain(4)-c1qdrain(4)),                                 &
     &  comma,(cqdrain(5)-c1qdrain(5)),                                 &
     &  comma,(cqdrd-c1qdrd),comma,(crunoff-c1runoff),                  &
     &  comma,(cQMpOutDrRap-c1qdrar),                                   &
     &  comma,cqdrain(1),comma,cqdrain(2),                              &
     &  comma,cqdrain(3),comma,cqdrain(4),comma,cqdrain(5),             &
     &  comma,cqdrd,comma,crunoff,comma,cQMpOutDrRap

 20   format (A11,a1,I4,8(a1,f8.2),8(a1,f8.1))

! --- store cumulative values
      do level=1,nrlevs
          c1qdrain(level)=cqdrain(level)
      enddo
      c1qdrd = cqdrd
      c1runoff = crunoff
      c1qdrar = cQMpOutDrRap

      case default
         call fatalerr ('outdrf', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outswb(task,outfil,vtair,pathwork,daynr,daycum,hbweir, &
     & overfl,gwl,pond,wlstar,wls,swstini,swst,cqdrd,crunoff,           &
     & cQMpOutDrRap,cwsupp,swb,numadj,hwlman,cwout,swsec,swman,nmper,   &
     & impend,imper,project,logf,swscre,date,t1900,t,outper,iyear)
! ----------------------------------------------------------------------
!     Date               : 21/08/99
!     Purpose            : write surface water balance data to 
!                          OUTNAM.SWB file, surface water management
!                          data to OUTNAM.MAN. The files overlap
! ---------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- Global
      integer   swb,man,daynr,daycum,swsec,swman(mamp),nmper
      integer   imper,numadj,logf,swscre,task,iyear
      real(8)   gwl,pond,wlstar,swstini,wls,swst,cqdrd,crunoff
      real(8)   impend(mamp),cQMpOutDrRap,t,outper
      real(8)   cwsupp,cwout,hwlman,vtair,hbweir(nmper),t1900
      character(len=16) outfil
      character(len=80) pathwork,project
      character(len=11) date
      logical   overfl

! --- Local
      integer   getun,nrOfDays
      real(8)   gwlev,cqdrf1,c1wsupp,c1wout,delbal,small,zero
      character(len=1)   spc,comma
      character(len=200) messag
      character(len=132) filnam
      character(len=80)  filtext
      character(len=19) datetime
      logical   dtleap

      save

      data      small     /0.0001d0/
      data      zero      /0.0d0/
      comma = ','
! ----------------------------------------------------------------------

      select case (task)
      case (1)
! === open output files and write headers ==============================

! --- open output file once
      filnam = trim(pathwork)//trim(outfil)//'.swb'
      swb = getun (20,90)
      call fopens(swb,filnam,'new','del')
      filtext = 'Surface water balance increments (cm/period)'
      call writehead (swb,1,filnam,filtext,project)

! --- write header
      write (swb,10)
      write (swb,11)
 10   format ('*',/,                                                    &
     &'* Surface water system:',/,                                      &
     &'* - GWL           : groundwater level + ponding (<0. = bel',     &
     &'ow s.s) (cm)',/,                                                 &
     &'* - WLST          : target level,  aut. weir (<0. = b.s.s.',     &
     &') (cm)',/,                                                       &
     &'* - WLST          : weir crest,   fixed weir (<0. = b.s.s.',     &
     &')  (cm)',/,                                                      &
     &'* - WLS           : surface water level (<0. = b.s.s.)  (c',     &
     &'m)',/,                                                           &
     &'* - swst           : storage in sw-reservoir, per unit area',    &
     &' of* subcatchment (cm)',/,                                       &
     &'* - DRORR         : drainage fluxes (>0. = into sw) + runo',     &
     &'ff +* rapid drainage (cm)',/,                                    &
     &'* - QSUPP         : external supply to sw-reservoir (cm)',/,     &
     &'* - QOUT          : outflow from sw-reservoir (cm)',/,           &
     &'* - DRORRcum      : drainage fluxes (>0. = into sw) + runo',     &
     &'ff +* rapid drainage (cumulative, cm)')
 11   format (                                                          &
     &'* - QSUPPcum      : external supply to sw-reservoir (cumul',     &
     &'ative, cm)',/,                                                   &
     &'* - QOUTcum       : outflow from sw-reservoir (cumulative,',     &
     &' cm)',/,'*')

! ---   determine which management period the model is in:
      if (swsec .eq. 2) then
        imper = 0
 100    imper = imper + 1
!       Error handling
        if (imper .gt. nmper) then
          messag = 'error sw-management(IMPER), more than defined'
          call fatalerr ('Outswb',messag)
        endif
        if (t1900-1.d0-0.1d-10 .gt. impend(imper)) goto 100
      endif

! --- add ponding to groundwater level, if gwl at soil surface:
      if (abs(gwl) .lt. 1.0d-7) then
          gwlev = pond
      else
          gwlev = gwl
      endif

! --- write record for initial state
      write (SWB,20) comma,gwlev,comma,wlstar,comma,wls,comma,swst,     &
     & comma,(0.),comma,(0.),comma,(0.),comma,(0.),comma,(0.),comma,(0.)
 20   format ('* initial  ,   - ,    - ',                               &
     & 3(a1,f7.1),a1,f6.1,6(a1,f7.2),/,'*',/,                           &
     &'      Date,Daynr,DayCum,    GWL,   WLST,    WLS,   swst,  ',     &
     &'DRORR,  QSUPP,  QOUT,DRORRcum,QSUPPcum,QOUTcum')

! --- compose filename for management data
      if (SWSEC .eq. 2) then

! --    open output file once
        filnam = trim(pathwork)//trim(outfil)//'.man'
        man = getun (20,90)
        call fopens(man,filnam,'new','del')
        filtext = 'Surface water management (cm/period)'
        call writehead (man,1,filnam,filtext,project)

! --    write header
        write (MAN,30)
 30     format ('*',/,                                                  &
     &'* Surface water management:',/,                                  &
     &'* - SWMAN         : type of weir (f = fixed; a = automatic)',/,  &
     &'* - GWL           : groundwater level + ponding (<0. = belo',    &
     &'w s.s) (cm)',/,                                                  &
     &'* - HWLMAN        : pressure head used for target level (cm)',/, &
     &'* - VTAIR         : total air volume in soil profile (cm)',/,    &
     &'* - WLST          : target level of autom. weir (<0. = b.s.',    &
     &'s.)  (cm)',/,                                                    &
     &'*              or : crest  level of fixed  weir (<0. = b.s.',    &
     &'s.)  (cm)',/,                                                    &
     &'* - WLS           : surface water level (<0. = b.s.s.) (cm)',/,  &
     &'* - QOUT          : surface water outflow (<0. = supply) (cm)',/,&
     &'* - NUMADJ        : number of adjustments of target level (',    &
     &'cm)',/,                                                          &
     &'* - OVERFL        : flag for overflow of aut. weir (o)',/,       &
     &'* - CREST         : Crest level (<0. = b.s.s.) (cm)',/,          &
     &'*',/,                                                            &
     &'     Date,DayNr,DayCum,SWMan,    GWL,HWLMAN, VTAIR,   WLST,',    &
     &'    WLS,   QOUT, NUMADJ,OVERFL,  CREST')
      endif

      return

      case (2)

! === actual output data =================================================

! --- reset cumulative fluxes each year
      if (abs(t-outper).lt.small) then
        cqdrf1 = zero
        c1wsupp = zero
        c1wout = zero
      endif

! --- water balance error
      delbal = (swst + cwout) - (swstini + cqdrd+crunoff + cwsupp +     &
     &                          cQMpOutDrRap)
      if (delbal .gt. 0.05d0) then
        call dtdpst                                                     &
     &        ('year-month-day,hour:minute:seconds',t1900,datetime)
        write(messag,1001) delbal,datetime
 1001   format('  Error in cumulative water balance of surfwater:',     &
     &         f5.2,'   at : ',a19)
        call warn ('Outswb',messag,logf,swscre)
      endif

! --- add ponding to groundwater level, if gwl at soil surface:
      if (abs(gwl) .lt. 1.0d-7) then
          gwlev = pond
      else
          gwlev = gwl
      endif

! --- write output record OUTNAM.SWB
      write (SWB,40) date,comma,daynr,comma,daycum,comma,gwlev,comma,   &
     &  wlstar,comma,wls,comma,swst,comma,                              &
     &  (cqdrd+crunoff+cQMpOutDrRap-cqdrf1),comma,(cwsupp-c1wsupp),     &
     &  comma,(cwout-c1wout),comma,(cqdrd+crunoff+cQMpOutDrRap),        &
     &  comma,cwsupp,comma,cwout
 40   format (a11,a1,I4,a1,I6,3(a1,f7.1),a1,f6.1,6(a1,f7.2))

! --- write output record OUTNAM.MAN
      if (swsec .eq. 2) then
        if (overfl) then
          spc = 'o'
        else
          spc = '-'
        endif
        if (swman(imper) .eq. 1) then
          write(man,50) date,comma,daynr,comma,daycum,comma,gwlev,comma,&
     &      wlstar,comma,wls,comma,((cwout-c1wout)-(cwsupp-c1wsupp)),   &
     &      comma,numadj,comma,spc,comma,hbweir(imper)
        else
          write(man,60) date,comma,daynr,comma,daycum,comma,gwlev,comma,&
     &      hwlman,comma,vtair,comma,wlstar,comma,wls,comma,            &
     &      ((cwout-c1wout)-(cwsupp-c1wsupp)),comma,numadj,comma,spc,   &
     &      comma,hbweir(imper)
        endif
 50     format (a11,a1,i4,a1,I6,',    f',a1,f7.1,2(',     -'),          &
     &          3(a1,f7.1),a1,i7,a1,5x,a1,a1,f7.1)

 60     format (a11,a1,i4,a1,I6,',    a',a1,f7.1,2(a1,f7.1),            &
     &          3(a1,f7.1),a1,i7,a1,5x,a1,a1,f7.1)
      endif

! --- store cumulative values
      cqdrf1 = (cqdrd+crunoff+cQMpOutDrRap)
      c1wsupp = cwsupp
      c1wout = cwout

! --- saves storage as initial values for next year
      nrOfDays = 365
      if (dtleap(iyear)) nrOfDays = nrOfDays + 1
      if (daynr.eq.nrOfDays) swstini = swst

      case default
         call fatalerr ('outswb', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine warn (modul,messag,logf,swscre)
      implicit none
! ----------------------------------------------------------------------
! --- global
      integer       logf,swscre
      character(len=*) modul,messag
 
! --- local
      character(len=400)  messages
! ----------------------------------------------------------------------
 
      messages = 'Warning from module '//modul//' : '//trim(messag)

      write (logf,'(a)') trim(messages)
      if (swscre .gt. 0) write (*,'(2x,a)') trim(messages)

      return
      end

      subroutine writehead(outf,ftype,filnam,filtext,project)
!-----------------------------------------------------------------------
!      Date    : January 2006
!      Purpose : Writes a header to output files of SWAP

!-----------------------------------------------------------------------
      implicit none

!     global
      integer       outf,ftype
      character(len=*) project,filnam,filtext

!     local
      integer       date_time(6)
      real(8)       dpactualtime
      real(4)       dum
      character(len=80)  model_id,dtstring,String
      character(len=132) Version
!-----------------------------------------------------------------------

! --- Version nr of model, to appear in all output files
!     together with revision nr
      include 'description.fi'

      model_id = 'Swap '//trim(Version)

! --- get actual time
      call dtnow (date_time)
      dum = 0.0
      call dtardp (date_time, dum, dpactualtime)
      call dtdpst ('year-month-day hour:min:sec', dpactualtime,       &
     &             dtstring)

      if (outf.eq.5) then
! ---   write to screen
        write (*,16)  trim(project)
        write (*,17)  trim(filtext)
        write (*,19)  trim(model_id)
        write (*,20)  trim(dtstring)
      else if (ftype.eq.1) then
!       write formatted output file
        write (outf,16)  trim(project)
        write (outf,17)  trim(filtext)
        write (outf,18)  trim(filnam)
        write (outf,19)  trim(model_id)
        write (outf,20)  trim(dtstring)
      else
!       write unformatted output file (header has fixed length of 80 characters)
        write (String,'(a80)')  '* Project:       '//trim(project)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* File content:  '//trim(filtext)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* File name:     '//trim(filnam)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* Model version: '//trim(model_id)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* Generated at:  '//trim(dtstring)
        write (outf)  adjustl(String)
      endif

 16   format('* Project:       ',a)
 17   format('* File content:  ',a)
 18   format('* File name:     ',a)
 19   format('* Model version: ',a)
 20   format('* Generated at:  ',a)

      return
      end

! ----------------------------------------------------------------------
      subroutine CloseTempFil
! ----------------------------------------------------------------------
!     date               : Aug 2004             
!     purpose            : delete temporary files of TTUTIL
! ----------------------------------------------------------------------
      implicit  none

      integer cexf,getun
      logical fileopen

! --- delete temporary files
      cexf = getun (10,90)
      call rddtmp(cexf)
      inquire(unit=20,opened=fileopen)
      if(fileopen)close(20, STATUS = 'DELETE')

      return
      end


! ----------------------------------------------------------------------
      subroutine checkDiscrVert(numnod,dz,numnodNew,dzNew)
!     date               : 20081105
!     purpose            : verify reduced vertical discretizationface,
! global   formal parameters  : (i = input, o = output)
!     numnod       ! Number of nodes or compartments....................... i
!     numnodnew    ! Number of desired nodes for soil water quality models..i
!     dz(macp)     ! Compartment thickness (L) ............................ i
!     dznew(macp)  ! Desired dz for soil water quality models (L) ......... i
! local
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'
! global
      integer   numnod,numnodNew
      real(8)   dz(macp),dznew(macp)
! local
      integer   in,io,iotmp
      real(8)   cumdzN(macp),cumdzOld,cumdzNew
      character(len=200) messag
! local parameters
      real(8)   small
      data      small     /0.0001d0/
! ----------------------------------------------------------------------


! --  boundary of new compartments must equal boundaries of old compartmts
!     or: new size is sum of old sizes
      iotmp = 0
      do in = 1,numnodNew
        cumdzN(in) = 0.0d0
        do io = 1,numnod
          if(io.gt.iotmp .and. cumdzN(in).lt.dzNew(in)) then 
             cumdzN(in) = cumdzN(in) + dz(io)
             iotmp = io
          endif
        enddo
      enddo
      do in = 1,numnodNew
        if (abs(dzNew(in)-cumdzN(in)).gt.small) then
          write(messag,'(a,i5)')                                        &
     &    'New discretization has error in size of new compartment ',in
          call fatalerr ('checkDiscrVert',messag)
        endif
      enddo

! --  cumulative thickness
      cumdzOld = 0.0d0
      do io = 1,numnod
         cumdzOld = cumdzOld + dz(io)
      enddo
      cumdzNew = 0.0d0
      do in = 1,numnodNew
         cumdzNew = cumdzNew + dznew(in)
      enddo
      if (abs(cumdzNew-cumdzOld).gt.small) then
        write(messag,'(a,f10.5)')                                       &
     &    'New discretization is wrong, total length differs >',small
        call fatalerr ('checkDiscrVert',messag)
      endif


      return
      end

! ----------------------------------------------------------------------
      subroutine OutputModflow(task) 
! ----------------------------------------------------------------------
!     Date               : July 2009
!     Purpose            : open and write data to explore storage/recharge
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer   getun,sto,nod1m, nod
      real(8)   vsat,vt0,gwlt0,vt1,qre,qv1m,gwlt1,stocoav,stocot1
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      integer   i,swBotbtmp
      real(8)   vair(2),dgwl(2),gwltmp,pondtmp,qbottmp, vt2
      real(8)   thetatmp(macp),htmp(macp),xd(6)


      save      sto,vsat,vt0,gwlt0, stocot1

! ----------------------------------------------------------------------
      comma = ','


      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.sto'
      sto = getun (20,90)
      call fopens(sto,filnam,'new','del')
      filtext = 'Storage and recharge output data '//                   &
     &'(q in cm/period; gwl in cm; sto in cm)'
      call writehead (sto,1,filnam,filtext,project)

! --- write header
      write (sto,10)
 10   format ('*',/,                                                    &
     &        'date,cumtime,qvss,qtr,qevs,qrun,qdr,qbot,qre,qv1m,'//    &
     &        'gwlt0,gwlt1,vt0,vt1,stocoav,stocot1')

      vsat = 0.0d0
      do nod = 1,numnod
         vsat = vsat + cofgen(2,nod) * dz(nod)
      end do

      return

      case (2)

! === write actual data ================================

! --- write header in case of new balance period
      if (flheader) write (sto,10) 
! --- write actual data
      vt1 = 0.0d0
      do nod = 1,numnod
         if(nod.gt.1)then
           if(z(nod)-0.5d0*dz(nod).le.-100.d0                           &
     &                   .and. z(nod-1)-0.5d0*dz(nod-1).gt.-100.d0)then
              nod1m = nod
           end if
         end if
         vt1 = vt1 + theta(nod) * dz(nod)
      end do
      vt1 = vsat - vt1
      qre = iqdra - iqbot
      qv1m = q(nod1m)
      gwlt1 = gwl
      if(dabs(gwlt1 - gwlt0) .gt.1.0d-6)then
         stocoav = - (vt1 - vt0) / (gwlt1 - gwlt0)
      else
         stocoav = 1.0d+6
      end if   


      vt0 = vt1
      gwlt0 = gwlt1


! --- initilization
      xd(1) = q(1)
      xd(2) = iqrot
      xd(3) = ievap
      xd(4) = iruno
      xd(5) = iqdra
      xd(6) = iqbot
      
      swBotbtmp = swbotb
      qbottmp = qbot
      gwltmp = gwl
      pondtmp = pond
      do nod =1,numnod
         thetatmp(nod) = theta(nod)
         htmp(nod)     = h(nod)
      end do


      swbotb =1
      dgwl(1) = 1.0d0
      dgwl(2) = -1.0d0
! ---   determine SoilWater bottom boundary conditions

      do i=1,2

        gwlinp = gwltmp + dgwl(i)

!        call BoundBottom

        fldtreduce = .true.
        do while(fldtreduce)
           fldtreduce = .false.

! ---   calculate drainage fluxes
           if (fldrain)                           call Drainage
           if (.not.fldecdt .and. flSurfaceWater) call SurfaceWater(2)
           if (SwFrost.eq.1)                      call FrozenBounds
      
! ---   calculate SoilWater 
           if (.not.fldecdt) then 
! --- save state variables of time = t
              call SoilWaterStateVar(1)
     
! --- calculate new soil water state variables
              call headcalc

! ---   calculate surface water balace
              if (flSurfaceWater) call SurfaceWater(3)
           end if

! ---   update time variables and switches/flags
           if(fldecdt .or. (flMacroPore .and. FlDecMpRat))then
              call SoilWaterStateVar(2)
              call TimeControl(3)
              fldtreduce = .true.
           end if

        end do

        vt2 = 0.0d0
        do nod = 1,numnod
           vt2 = vt2 + theta(nod) * dz(nod)
        end do
        vair(i) = vsat - vt2

      end do


      stocot1 = - (vair(1) - vair(2))  / (dgwl(1) - dgwl(2))


      write (sto,20) date,comma,t1900,(comma,xd(i),i=1,6),comma,        &
     &               qre,comma,qv1m,comma,gwlt0,comma,gwlt1,comma,      &
     &               vt0,comma,vt1,comma,stocoav,comma,stocot1
20    format (a11,a1,f12.4,14(a1,f10.4))


      swbotb = swBotbtmp
      qbot = qbottmp
      gwl = gwltmp
      pond = pondtmp
      do nod =1,numnod
         theta(nod) = thetatmp(nod)
         h(nod)     = htmp(nod)
      end do


      return

      case (3)

! === close output file ===========================

! --- close snw file
      close (sto)

      case default
         call fatalerr ('OutputModFlow', 'Illegal value for Task')
      end select

      return
      end


! ----------------------------------------------------------------------
      subroutine capriseoutput(task) 
! ----------------------------------------------------------------------
!     Date               : September 2016
!     Purpose            : open and write data to explore capillary rise to/from rootzone
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer   getun,sto
      character(len=300) filnam
      character(len=200) filtext
      character(len=1)   comma

      save      sto

! ----------------------------------------------------------------------
      comma = ','


      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.crz'
      sto = getun (20,90)
      call fopens(sto,filnam,'new','del')
      filtext = 'Capillary rise q to/from rootzone when roots '//      &
     &' are present (rd>0) (q in cm/period, pos=upward, neg=downward)'
      call writehead (sto,1,filnam,filtext,project)

! --- write header
      write (sto,10)
 10   format ('*',/,                                                    &
     &        'date,tcum,noddrz,rd,inq(noddrz)')
      return

      case (2)

! === write actual data ================================
      if(rd.gt.0.0d0) then
        write (sto,20) date,comma,tcum,comma,noddrz,comma,rd,comma,&
    &                  inq(noddrz+1)
 20     format (a11,a1,f12.4,a1,i4,a1,f10.4,a1,f20.12)
      endif

      return

      case (3)

! === close output file ===========================

! --- close snw file
      close (sto)

      case default
         call fatalerr ('capriseoutput', 'Illegal value for Task')
      end select

      return
      end
