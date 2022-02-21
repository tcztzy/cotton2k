! File VersionID:
!   $Id: management_soil.f90 329 2017-06-12 08:09:14Z kroes006 $
! ----------------------------------------------------------------------
      subroutine SoilManagement(task) 
! ----------------------------------------------------------------------
!     Date               : March 2015
!     Purpose            : save and reset soil water state variables
! ----------------------------------------------------------------------

! --- global variables
      use Variables
      use Wofost_Soil_Declarations
      use Wofost_Soil_Interface

      implicit none

! --- local variables
      character(len=300) filnam, filnamce
!      character(len=100) MatName(masme)
      integer   task, sme, smm, snp, getun2, nsmm, nsme, oup
      Integer :: i,j, idum, getun
      real(8) :: dum, t4, t3, help
      real(8)   MatAmount(masme) ! Array with amount of applied material (kg/ha)
      real(8)   smedate(masme)   ! Array with soilmanagement event dates (-)

!0.3  Local variables
      integer :: le, fn
      real(8):: ProdRate0, ProducPot, ProducAct, TCSF
      real(8):: cNH4_t_Ndemand_rate_limited
      real(8):: cNH4_av_Ndemand_rate_limited
      real(8):: cNH4_t_Nsupply_rate_limited
      real(8):: cNH4_av_Nsupply_rate_limited
      real(8):: cNO3_t_Ndemand_rate_limited
      real(8):: cNO3_av_Ndemand_rate_limited
      real(8):: cNO3_t_Nsupply_rate_limited
      real(8):: cNO3_av_Nsupply_rate_limited
      real(8):: Nsupply_Ndemand_rate_limited
      real(8):: Nsupply_Nsupply_rate_limited
      real(8):: dum1, dum2, dum3, dum4 !, dum5, dum6, dum7
      character(len=1) comma
      character(len=160) filtext, line

      integer :: stat      ! ,imat
      real(8) :: xAmend, xAppAge, xOrgMatFrac, xOrgNFrac, xNH4NFrac,    &
     &           xNO3NFrac

      real(8) :: FactNuptJuvenil
      logical :: FlNuptJuvenil

! optional input for Experts
      integer :: swexpertN      
      real(8) :: AppAgeArableRt, AppAgeArableLv
      real(8) :: AppAgeArableSt, AppAgeArableSo
      real(8) :: AppAgeGrassRt, AppAgeGrassLv, AppAgeGrassSt
      logical :: rdinqr
      
      save

      select case (task)
      case (1)

! --- open output file once
      filnam = trim(pathwork)//trim(outfil)//'_nut.csv'
      nut = getun (20,90)
      call fopens(nut,filnam,'new','del')
      filtext = 'nutrient balance increments (kg/ha)'
      call writehead (nut,1,filnam,filtext,project)
! --- write header of nut file
      write (nut,'(a,/,a)')'*','Date,Day,Dcum,'//                       &
     &    'FOM_old,FOM_end,FOM_dif,FOM_add,FOM_cres,FOM2Bio,FOM2Hum,'// &
     &    'FOM_dis,Bio_old,Bio_end,Bio_dif,Bio2Bio,Bio2Hum,Bio_dis,'//  &
     &    'Hum_old,Hum_end,Hum_dif,Hum_add,Hum_cres,Hum2Bio,Hum2Hum,'// &
     &    'Hum_dis,cDissi,NFOM_old,NFOM_end,NFOM_dif,NFOM_add,'//       &
     &    'NFOM_cres,NFOM2Bio,NFOM2Hum,NFOM_min,'//                     &
     &    'NBio_old,NBio_end,NBio_dif,NBio2Bio,NBio2Hum,NBio_min,'//    &
     &    'NHum_old,NHum_end,NHum_dif,NHum_add,NHum_cres,NHum2Bio,'//   &
     &    'NHum2Hum,NHum_min,Nminer,'//                                 &
     &    'NH4_old,NH4_end,NH4_dif,NH4N_amend,NH4N_cres,'//             &
     &    'NH4_intop,NH4_inlat,NH4_inbot,'//                            &
     &    'NH4_upt,NH4_out,NH3N_volat,NH4_nitrif,NO3_old,NO3_end,'//    &
     &    'NO3_dif,NO3N_amend,NO3N_cres,NO3_intop,NO3_inlat,'//         &
     &    'NO3_inbot,NO3_upt,NO3_out,NO3_denitr,Ndemand,Nsupply,'//     &
     &    'WVol_old,WVol_end,WVol_dif,WFl_inTop,WFl_inLat,WFl_inBot,'// &
     &    'WFl_Eva,WFl_Tra,WFl_Out,'//                                  &
     &    'idwrt,idwlv,idwst,iNLOSSL,iNLOSSR,iNLOSSS,idwso,iNLOSSO,'//  &
     &    'WFPS,red_T,red_W,red_W_Nit,red_W_Den,red_Resp'

      filnamce = trim(pathwork)//trim(outfil)//'_crop_ext.csv'
      cropext = getun (20,90)
      call fopens(cropext,filnamce,'new','del')
      write(cropext,'(2a)')'    nishmi     nishma     niromi     nirom',&
     &                   'a     poshmi     poshma     poromi     poroma'
      write(cropext,'(a)')' to be replaced'
      write(cropext,'(2a)')'year-mo-da,Ntotuptake,Ptotuptake,DMcressur',&
     &        'f,N-cressurf,P-cressurf,DMcresbott,N-cresbott,P-cresbott'
      Ntotuptake = 0.0d0; Ptotuptake = 0.0d0; DMcressur  = 0.0d0
      Ncressurf  = 0.0d0; Pcressurf  = 0.0d0; DMcresbott = 0.0d0
      Ncresbott  = 0.0d0; Pcresbott  = 0.0d0; nishmi = 1.0d0
      nishma = 0.0d0; niromi = 1.0d0; niroma = 0.0d0; poshmi = 1.0d0
      poshma = 0.0d0; poromi = 1.0d0; poroma = 0.0d0

      call Wofost_SoilParameters

! -   Soil management events:  

! -   open file with soil management material definitions
      filnam = trim(pathwork)//trim(project)//'.smm'
      smm = getun2 (10,90,2)
      call rdinit(smm,logf,filnam)
      call rdainr ('MatNum',1,maxmat,MatNum,maxmat,nsmm)
      call rdfcha ('MatName',MatName,maxmat,nsmm)
      call rdfdor ('AppAge',0.0d0,5000.0d0,AppAge,maxmat,nsmm)
      call rdfdor ('OrgMatFrac',0.0d0,1.0d0,OrgMatFrac,maxmat,nsmm)
      call rdfdor ('OrgNFrac',0.0d0,1.0d0,OrgNFrac,maxmat,nsmm)
      call rdfdor ('NH4NFrac',0.0d0,1.0d0,NH4NFrac,maxmat,nsmm)
      call rdfdor ('NO3NFrac',0.0d0,1.0d0,NO3NFrac,maxmat,nsmm)
      close(smm)

!     open file with soil management events
      filnam = trim(pathwork)//trim(project)//'.sme'
      sme = getun2 (10,90,2)
      call rdinit(sme,logf,filnam)
      call rdatim ('smedate',smedate,masme,nsme)
      call rdfinr ('MatNum',1,masme,MatNum,masme,nsme)
      call rdfdor ('Dosagekgha',0.0d0,500000.0d0,MatAmount,masme,nsme)
      call rdfdor ('VolatFraction',0.0d0,500000.0d0,VolaFrac,masme,nsme)
      close(sme)

!     Sorting of dosages according to time of application

      Do i=1,nsme-1
         Do j=i+1,nsme
            If(smedate(i) .Gt. smedate(j))Then
               Dum          = smedate(j)
               smedate(j)   = smedate(i)
               smedate(i)   = Dum
               idum         = MatNum(j)
               MatNum(i)     = MatNum(i)
               MatNum(j)     = idum
               Dum          = MatAmount(j)
               MatAmount(j) = MatAmount(i) 
               MatAmount(i) = Dum
               Dum          = VolaFrac(j)
               VolaFrac(j)  = VolaFrac(i)
               VolaFrac(i)  = Dum
            End If
         End Do
      End Do

!     Converting kg/ha -> kg/m2

      Do i=1,nsme
         Amend(i) = 1.0d-4 * MatAmount(i)
      End Do

!     Grouping of dosages per date/time

      j = 1
      NuAmend(j) = 1
      TimeAmend(j) = smedate(1)
      iAmend(1,1) = 1
      do isme=2,nsme
         if(smedate(isme)-smedate(isme-1).lt.0.d-3)then
            NuAmend(j) = NuAmend(j) + 1
         else
            j = j + 1
            NuAmend(j) = 1
            TimeAmend(j) = smedate(isme)
         end if
         iAmend(j,NuAmend(j)) = isme
      end do
      

! -   open file with soil nutrient parameters
      filnam = trim(pathwork)//trim(project)//'.snp'
      snp = getun2 (10,90,2)
      call rdinit(snp,logf,filnam)
!     initial states
      call rdsdor ('FOM1_t',0.0d0,1000.0d0,FOM_t(1))
      call rdsdor ('FOM2_t',0.0d0,1000.0d0,FOM_t(2))
      call rdsdor ('FOM3_t',0.0d0,1000.0d0,FOM_t(3))
      call rdsdor ('FOM4_t',0.0d0,1000.0d0,FOM_t(4))
      call rdsdor ('FOM5_t',0.0d0,1000.0d0,FOM_t(5))
      call rdsdor ('FOM6_t',0.0d0,1000.0d0,FOM_t(6))
      call rdsdor ('FOM7_t',0.0d0,1000.0d0,FOM_t(7))
      call rdsdor ('FOM8_t',0.0d0,1000.0d0,FOM_t(8))
      call rdsdor ('Bio_t',0.0d0,1000.0d0,Bio_t)
      call rdsdor ('Hum_t',0.0d0,1000.0d0,Hum_t)
      call rdsdor ('cNH4_t',0.0d0,1000.0d0,cNH4_t)
      call rdsdor ('cNO3_t',0.0d0,1000.0d0,cNO3_t)
!     boundary concentrations
      call rdsdor ('cNH4N_top',0.0d0,1000.0d0,cNH4N_top)
      call rdsdor ('cNH4N_lat',0.0d0,1000.0d0,cNH4N_lat)
      call rdsdor ('cNH4N_seep',0.0d0,1000.0d0,cNH4N_seep)
      call rdsdor ('cNO3N_top',0.0d0,1000.0d0,cNO3N_top)
      call rdsdor ('cNO3N_lat',0.0d0,1000.0d0,cNO3N_lat)
      call rdsdor ('cNO3N_seep',0.0d0,1000.0d0,cNO3N_lat)
!     Coefficients and rate constants
      call rdsdor ('Temp_ref',0.0d0,1000.0d0,Temp_ref)
      call rdsdor ('SorpCoef',0.0d0,1000.0d0,SorpCoef)
      call rdsdor ('RateConNitrif_ref',0.0d0,1000.0d0,RateConNitrif_ref)
      call rdsdor ('RateConDenitr_ref',0.0d0,1000.0d0,RateConDenitr_ref)
!     Response function parameters
      call rdsdor ('WFPSCrit',0.0d0,1000.0d0,WFPSCrit)
      call rdsdor ('WFPScrit2',0.0d0,1000.0d0,WFPScrit2)
      call rdsdor ('CdissiHalf',0.0d0,1000.0d0,CdissiHalf)
!     Soil supply uptake parameters
      call rdsdor ('TCSF_N',0.0d0,1000.0d0,TCSF_N)
      LaiCritNupt = 0.0d0
      if(rdinqr('LaiCritNupt')) then
        call rdsdor ('LaiCritNupt',0.0d0,10.0d0,LaiCritNupt)
      endif
!     Effective depth of soil layer
      call rdsdor ('dz_WSN',0.0d0,1000.0d0,dz_WSN)
!     for experts: apparent age of crop specific residues roots, leaves, shoots, storageOrgans 
      swexpertN = 0
      AppAgeArableRt = 1.27d0
      AppAgeArableLv = 0.99d0
      AppAgeArableSt = 1.27d0
      AppAgeArableSo = 1.27d0
      AppAgeGrassRt = 1.20d0
      AppAgeGrassLv = 0.92d0
      AppAgeGrassSt = 0.99d0
      RateconHum_exp = -1.0d0
      if(rdinqr('swexpertN')) then
        call rdsinr ('swexpertN',0,1,swexpertN)
        if (swexpertN.eq.1) then
          call rdsdor ('AppAgeArableRt',0.0d0,100.0d0,AppAgeArableRt)
          call rdsdor ('AppAgeArableLv',0.0d0,100.0d0,AppAgeArableLv)
          call rdsdor ('AppAgeArableSt',0.0d0,100.0d0,AppAgeArableSt)
          call rdsdor ('AppAgeArableSo',0.0d0,100.0d0,AppAgeArableSo)
          call rdsdor ('AppAgeGrassRt',0.0d0,100.0d0,AppAgeGrassRt)
          call rdsdor ('AppAgeGrassLv',0.0d0,100.0d0,AppAgeGrassLv)
          call rdsdor ('AppAgeGrassSt',0.0d0,100.0d0,AppAgeGrassSt)
          call rdsdor ('RateconHum_exp',0.0d0,100.0d0,RateconHum_exp)
          RateconHum_exp = RateconHum_exp / 365.2425d0
          call rdstim('tstartHumexp',tstartHumexp)
        endif
      endif
      close(snp)

! --- vertical discretization of soil profile

      dum1=0.0; dum2=0.0; idum = 0
      do i=1,numnod
         if(dum1 + 1.0d-2 * dz(i) .lt. dz_WSN)then
            dum1 = dum1 + 1.0d-2 * dz(i)
            dum2 = dum2 + cofgen(2,i) * 1.0d-2 * dz(i)
         end if
      end do
      WFrac_sat = dum2 / dum1


!     initialize counter for amendments
      isme = 1
      t_WSNold = tstart
      Hum_add    = 0.0d0
      NHum_add   = 0.0d0
      FOM_add    = 0.0d0
      NFOM_add   = 0.0d0
      NH4N_amend = 0.0d0
      NO3N_amend = 0.0d0
      Hum_cres   = 0.0d0
      NHum_cres  = 0.0d0
      FOM_cres   = 0.0d0
      NFOM_cres  = 0.0d0
      NH4N_cres  = 0.0d0
      NO3N_cres  = 0.0d0
      NH4N_volat = 0.0d0
      FOM2Bio    = 0.0d0
      FOM2Hum    = 0.0d0
      Bio2Bio    = 0.0d0
      Bio2Hum    = 0.0d0
      NFOM2Bio   = 0.0d0
      NFOM2Hum   = 0.0d0
      NBio2Bio   = 0.0d0
      Bio2Hum    = 0.0d0

      FOM_old    = 0.0d0
      NFOM_old   = 0.0d0
      do fn = 1,nf
         FOM_old    = FOM_old + FOM_t(fn) * dz_WSN                      !kg/m2
         NFOM_old   = NFOM_old + NFracFOM(fn) * FOM_t(fn) * dz_WSN      !kg/m2
      end do
      Bio_old       = Bio_t  * dz_WSN                                   !kg/m2
      NBio_old      = NFracBio * Bio_t * dz_WSN                         !kg/m2
      Hum_old       = Hum_t * dz_WSN                                    !kg/m2
      NHum_old      = NFracHum * Hum_t * dz_WSN                         !kg/m2
      dum1=0.0; dum2=0.0; idum = 0
      do i=1,numnod
         if(dum1 + 1.0d-2 * dz(i) .lt. dz_WSN)then
            dum1 = dum1 + 1.0d-2 * dz(i)
            dum2 = dum2 + theta(i) * 1.0d-2 * dz(i)
         end if
      end do
      WFrac_t      = dum2 / dum1
      WFrac_t0     = WFrac_t
      NH4_old = (WFrac_t0 + DryBD * SorpCoef ) * cNH4_t * dz_WSN        !kg/m2
      NO3_old = WFrac_t0 * cNO3_t * dz_WSN                              !kg/m2

      FOM_end = FOM_old                                                 !kg/m2
      Bio_end = Bio_old                                                 !kg/m2
      Hum_end = Hum_old                                                 !kg/m2
      NFOM_end = NFOM_old                                               !kg/m2
      NBio_end = NBio_old                                               !kg/m2
      NHum_end = NHum_old                                               !kg/m2
      NH4_end = NH4_old                                                 !kg/m2
      NO3_end = NO3_old                                                 !kg/m2

      return

      case (2)
!     update varaibles
      do fn=1,nf
         FOM_t0(fn) = FOM_t(fn)
      end do
      Bio_t0  = Bio_t
      Hum_t0  = Hum_t
      cNH4_t0 = cNH4_t
      cNO3_t0 = cNO3_t
      WFrac_t0 = WFrac_t
      t1900Soil = t1900

      return

      case (3)
! --- timed soil management events
      if (abs(TimeAmend(isme) +1.0d0 - t1900) .lt. 1.d-3) then

         call Wofost_SoilAmendents
         isme = isme + 1

      endif

      return

      case (4)
! --- mineralisation (production of mineral nitrogen)

      dum1=0.0; dum2=0.0; dum3=0.0; dum4=0.0
      idum = 0
      do i=1,numnod
         if(dum1 + 1.0d-2 * dz(i) .lt. dz_WSN)then
            idum = idum + 1
            dum1 = dum1 + 1.0d-2 * dz(i)
            dum2 = dum2 + tsoil(i) * 1.0d-2 * dz(i)
!            dum3 = dum3 + thetm1(i) * 1.0d-2 * dz(i)
            dum4 = dum4 + theta(i) * 1.0d-2 * dz(i)
         end if
      end do
      Temp          = dum2 / dum1 
!      WFrac_t0      = dum3 / dum1
      WFrac_t       = dum4 / dum1
      dt_WSN        = t1900 - t_WSNold
      t_WSNold      = t1900

      call Wofost_SoilRateConstants(1)

      call Wofost_SoilOrgMatN

! --- transformation and transport processes of soluble nitrogen

!     results from potential crop growth
      Ndemand =    1.0d-4 * NdemandSoil

!     water balance items of soil layer for which nutriemnts are simulated
      dum1=0.0;dum2=0.0;dum3=0.0;dum4=0.0
      idum = 0
      
      do i=1,numnod
         dum2 = dum2 + inqrot(i)/outper 
         if(dum1 + 1.0d-2 * dz(i) .lt. dz_WSN)then
            dum1 = dum1 + 1.0d-2 * dz(i)
            do le=1,5
               dum3 = dum3 - min(0.0d0, (inqdra(le,i)/outper) )
            end do
         end if
      end do
      help = 1.0d-2 * (igrai+isnrai+igsnow+igird-iintc+irunon-iruno) /  &
     &                 outper
      SoilEvap = 1.0d-2 * ievap / outper
      if( (help + SoilEvap) .lt. 0.0d0 ) help = - SoilEvap
      Wflux_inTop   = help + SoilEvap 
      Wflux_inLat   = 1.0d-2 * dum3 
      Wflux_transp  = 1.0d-2 * dum2 
      dum4          = help + 1.0d-2 * (  dum3 - dum2 ) +                &
     &                (WFrac_t0 - WFrac_t) * dz_WSN/dt_WSN 
      Wflux_out     = max(0.0d0,dum4)
      Wflux_inBot   = -min(0.0d0,dum4)

      call Wofost_SoilRateConstants(2)

!     ammonium



! Mineralisation, but Nminer is negative in case of immobilisation
      ProdRate0 = Nminer / dt_WSN


! Crop uptake
      FactNuptJuvenil = 0.0d0
      FlNuptJuvenil = .false.
      if(flCropcalendar .and.dvs.lt. 1.0d0 .and. LaiCritNupt.gt.1.0d-02 &
     &                                     .and. lai.lt.LaiCritNupt)then
         FactNuptJuvenil = ( LaiCritNupt - lai ) / LaiCritNupt
         FlNuptJuvenil = .true.
      end if
      FactNuptJuvenil = max(FactNuptJuvenil,0.0d0) 

!     N demand rate limiting
      ProducPot = ProdRate0 - Ndemand / dz_WSN                          ! kg/m3
      TCSF       = 0.0d0
      call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,         &
     &                 Wflux_out, Wflux_transp,Wflux_inBot, Wflux_inTop,&
     &                 Wflux_inLat, TCSF, RateConNitrif, ProducPot,     &
     &                 ProducAct, DryBD, SorpCoef, cNH4N_seep,cNH4N_top,&
     &                 cNH4N_lat, cNH4_t0, cNH4_t, cNH4_av)
      cNH4_t_Ndemand_rate_limited  = cNH4_t
      cNH4_av_Ndemand_rate_limited = cNH4_av
      Nsupply_Ndemand_rate_limited = ProdRate0 - ProducAct

!     Combined N demand rate and NH4 supply rate limiting
      ProducPot = ProdRate0 - FactNuptJuvenil * Ndemand / dz_WSN 
      TCSF      = TCSF_N*(0.5*(WFrac_t+WFrac_t0) + DryBD*SorpCoef)/     &
     &            (0.5*(WFrac_t+WFrac_t0)) * (1.0d0 - FactNuptJuvenil)
!        limitation for rooting depth
!      TCSF      = TCSF * max(1.0d0, (0.01* rd/dz_WSN))
      call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,         &
     &                 Wflux_out, Wflux_transp,Wflux_inBot, Wflux_inTop,&
     &                 Wflux_inLat, TCSF, RateConNitrif, ProducPot,     &
     &                 ProducAct, DryBD, SorpCoef, cNH4N_seep,cNH4N_top,&
     &                 cNH4N_lat, cNH4_t0, cNH4_t, cNH4_av)
      cNH4_t_Nsupply_rate_limited  = cNH4_t
      cNH4_av_Nsupply_rate_limited = cNH4_av
      Nsupply_Nsupply_rate_limited = TCSF*Wflux_transp * cNH4_av /      &
     &     dz_WSN + FactNuptJuvenil * Ndemand / dz_WSN                  ! kg/m3

      if(cNH4_t_Ndemand_rate_limited.lt.cNH4_t_Nsupply_rate_limited)then
         cNH4_t      = cNH4_t_Nsupply_rate_limited
         cNH4_av     = cNH4_av_Nsupply_rate_limited
         NsupplyNH4N = Nsupply_Nsupply_rate_limited
      else
         cNH4_t      = cNH4_t_Ndemand_rate_limited
         cNH4_av     = cNH4_av_Ndemand_rate_limited
         NsupplyNH4N = Nsupply_Ndemand_rate_limited
      end if

!     nitrification
      ProdRate0 = 0.5d0 * (WFrac_t + WFrac_t0) * RateConNitrif * cNH4_av

!     N demand rate limiting
      ProducPot = ProdRate0 - ( Ndemand /dz_WSN - NsupplyNH4N )
      TCSF       = 0.0
      call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,         &
     &                 Wflux_out, Wflux_transp,Wflux_inBot, Wflux_inTop,&
     &                 Wflux_inLat, TCSF, RateConDenitr, ProducPot,     &
     &                 ProducAct,0.0d0, 0.0d0, cNO3N_seep, cNO3N_top,   &
     &                 CNO3N_lat, cNO3_t0, cNO3_t, cNO3_av)
      cNO3_t_Ndemand_rate_limited  = cNO3_t
      cNO3_av_Ndemand_rate_limited = cNO3_av
      Nsupply_Ndemand_rate_limited = ProdRate0 -ProducAct 

!     Combined N demand rate and NO3 supply rate limiting
      ProducPot = ProdRate0 - FactNuptJuvenil * ( Ndemand /dz_WSN -     &
     &                                                     NsupplyNH4N )
      TCSF      = TCSF_N * (1.0d0 - FactNuptJuvenil)
!     limitation for rooting depth
!      TCSF      = TCSF * max(1.0d0, (rd/dz_WSN))
      call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,         &
     &                 Wflux_out, Wflux_transp,Wflux_inBot, Wflux_inTop,&
     &                 Wflux_inLat, TCSF, RateConDenitr, ProducPot,     &
     &                 ProducAct,0.0d0, 0.0d0, cNO3N_seep, cNO3N_top,   &
     &                 CNO3N_lat, cNO3_t0, cNO3_t, cNO3_av)
      cNO3_t_Nsupply_rate_limited  = cNO3_t
      cNO3_av_Nsupply_rate_limited = cNO3_av
      Nsupply_Nsupply_rate_limited = TCSF *Wflux_transp *cNO3_av /      &
     &     dz_WSN + FactNuptJuvenil * ( Ndemand /dz_WSN - NsupplyNH4N ) ! kg/m3 

      if(cNO3_t_Ndemand_rate_limited.lt.cNO3_t_Nsupply_rate_limited)then
         cNO3_t      = cNO3_t_Nsupply_rate_limited
         cNO3_av     = cNO3_av_Nsupply_rate_limited
         NsupplyNO3N = Nsupply_Nsupply_rate_limited
      else
         cNO3_t      = cNO3_t_Ndemand_rate_limited
         cNO3_av     = cNO3_av_Ndemand_rate_limited
         NsupplyNO3N = Nsupply_Ndemand_rate_limited
      end if

      Nsupply = ( NsupplyNH4N + NsupplyNO3N ) * dz_WSN                  ! kg/m2

      Call Wofost_SoilBalanceCheck

      NsupplySoil = 1.0d+04 * Nsupply                                         ! kg/ha


      return

      case (5)

!     No MatName                      AppAge   OrgMatFrac OrgNFrac  NH4NFrac  NO3NFrac
!     11 'Green leaves'               0.92     1.0        t.b.d.    0.0       0.0
!     12 'Overground crop residues'   0.99     1.0        t.b.d.    0.0       0.0
!     13 'Root and stubble residues'  1.57     1.0        t.b.d.    0.0       0.0
!     14 'Grass shoots'               0.92     1.0        t.b.d.    0.0       t.b.d.
!     15 'Grass roots'                1.20     1.0        t.b.d.    0.0       t.b.d.
!     16 'Tree leaves'                2.25     1.0        t.b.d.    0.0       0.0
!     17 'Spruce needles'             3.34     1.0        t.b.d.    0.0       0.0

      if( idwrt .le.1.0d-8 .and. idwst .le.1.0d-8 .and.                 &
     &    idwlv .le.1.0d-8 .and. idwso .le.1.0d-8)return

!      if(CropType(icrop).eq.2)Then ! Arable crop  ! not yet active because grassland is not implemented

         xOrgMatFrac = 1.0d0
         xNH4NFrac   = 0.0d0
         xNO3NFrac   = 0.0d0
         if( idwrt .gt.1.0d-8 )Then
!            imat     = 13
            xAmend   = 1.0d-4 * idwrt
            xAppAge  = AppAgeArableRt
            xOrgNFrac= iNLOSSR / idwrt 
            call Wofost_Soil_CropResidues(xAmend,xAppAge,               &
     &                    xOrgMatFrac,xOrgNFrac,xNH4NFrac,xNO3NFrac)
         end if
         if( idwlv .gt.1.0d-8 )Then
!            imat     = 12
            xAmend   = 1.0d-4 * idwlv
            xAppAge  = AppAgeArableLv
            xOrgNFrac= iNLOSSL / idwlv
            call Wofost_Soil_CropResidues(xAmend,xAppAge,               &
     &                    xOrgMatFrac,xOrgNFrac,xNH4NFrac,xNO3NFrac)
         end if
         if( idwst .gt.1.0d-8 )Then
!            imat     = 13
            xAmend   = 1.0d-4 * idwst
            xAppAge  = AppAgeArableSt
            xOrgNFrac= iNLOSSS / idwst
            call Wofost_Soil_CropResidues(xAmend,xAppAge,               &
     &                    xOrgMatFrac,xOrgNFrac,xNH4NFrac,xNO3NFrac)
         end if
         if( idwso .gt.1.0d-8 )Then
!            imat     = 12
            xAmend   = 1.0d-4 * idwso
            xAppAge  = AppAgeArableSo
            xOrgNFrac= iNLOSSO / idwso
            call Wofost_Soil_CropResidues(xAmend,xAppAge,               &
     &                    xOrgMatFrac,xOrgNFrac,xNH4NFrac,xNO3NFrac)
         end if
         
! not yet active because grassland is not implemented

!      Else If (CropType(icrop).eq.3)Then ! Grassland Detailed
!
!         xOrgMatFrac = 1.0d0
!         xNH4NFrac   = 0.0d0
!         if( idwrt .gt.1.0d-8 )Then
!!            imat     = 15
!            xAmend   = 1.0d-4 * idwrt
!            xAppAge  = AppAgeGrassRt
!            xOrgNFrac= iNLOSSR / idwrt 
!            xNO3NFrac= 0.0d0
!            call Wofost_Soil_CropResidues(xAmend,xAppAge,               &
!     &                    xOrgMatFrac,xOrgNFrac,xNH4NFrac,xNO3NFrac)
!         end if
!         if( idwlv .gt.1.0d-8 )Then
!!            imat     = 14
!            xAmend   = 1.0d-4 * idwlv
!            xAppAge  = AppAgeGrassLv
!            xOrgNFrac= iNLOSSL / idwlv
!            xNO3NFrac= 0.0
!            call Wofost_Soil_CropResidues(xAmend,xAppAge,               &
!     &                    xOrgMatFrac,xOrgNFrac,xNH4NFrac,xNO3NFrac)
!         end if
!         if( idwst .gt.1.0d-8 )Then
!!            imat     = 12
!            xAmend   = 1.0d-4 * idwst
!            xAppAge  = AppAgeGrassSt
!            xOrgNFrac= iNLOSSS / idwst
!            xNO3NFrac= 0.0d0
!            call Wofost_Soil_CropResidues(xAmend,xAppAge,               &
!     &                    xOrgMatFrac,xOrgNFrac,xNH4NFrac,xNO3NFrac)
!         end if
!
!      End If

      idwrt = 0.0d0
      idwlv = 0.0d0
      idwst = 0.0d0
      idwso = 0.0d0
      
      iNLOSSL = 0.0d0
      iNLOSSR = 0.0d0
      iNLOSSS = 0.0d0
      iNLOSSO = 0.0d0

      return

      case (6)
! --- output and wrting

      comma = ','
      t4    = 1.0d+4
      t3    = 1.0d+3
      pnratio = 0.1d0
      
      Ntotuptake = Ntotuptake + t4 * ( NH4_upt + NO3_upt )
      Ptotuptake = Ptotuptake + pnratio * t4 * ( NH4_upt + NO3_upt )
      DMcressur  = DMcressur  + idwlv_1 + idwst_1 + idwso_1
      Ncressurf  = Ncressurf  + iNLOSSL_1 + iNLOSSS_1 + iNLOSSO_1
      Pcressurf  = Pcressurf  + pnratio *(iNLOSSL_1+iNLOSSS_1+iNLOSSO_1)
      DMcresbott = DMcresbott + idwrt_1
      Ncresbott  = Ncresbott  + iNLOSSR_1
      Pcresbott  = Pcresbott  + pnratio * iNLOSSR_1

      write (nut,'(a11,a1,i3,a1,i6,99(a1,1pe11.4:))')                   &
     &    date,comma,daynr,comma,daycum,comma,                          &
     &    t4*FOM_old,comma,t4*FOM_end,comma,t4*(FOM_old-FOM_end),comma, &
     &    t4*FOM_add,comma,t4*FOM_cres,comma,                           &
     &    t4*FOM2Bio,comma,t4*FOM2Hum,comma,t4*FOM_dis,comma,           &
     &    t4*Bio_old,comma,t4*Bio_end,comma,t4*(Bio_old-Bio_end),comma, &
     &    t4*Bio2Bio,comma,t4*Bio2Hum,comma,t4*Bio_dis,comma,           &
     &    t4*Hum_old,comma,t4*Hum_end,comma,t4*(Hum_old-Hum_end),comma, &
     &    t4*Hum_add,comma,t4*Hum_cres,comma,                           &
     &    t4*Hum2Bio,comma,t4*Hum2Hum,comma,t4*Hum_dis,comma,           &
     &    t4*cDissi,comma, t4*NFOM_old,comma,t4*NFOM_end,comma,         &
     &    t4*(NFOM_old-NFOM_end),comma,t4*NFOM_add,comma,               &
     &    t4*NFOM_cres,comma,t4*NFOM2Bio,comma,t4*NFOM2Hum,comma,       &
     &    t4*NFOM_min,comma,t4*NBio_old,comma,t4*NBio_end,comma,        &
     &    t4*(NBio_old-NBio_end),comma,                                 &
     &    t4*NBio2Bio,comma,t4*NBio2Hum,comma,t4*NBio_min,comma,        &
     &    t4*NHum_old,comma,t4*NHum_end,comma,                          &
     &    t4*(NHum_old-NHum_end),comma,t4*NHum_add,comma,               &
     &    t4*NHum_cres,comma,t4*NHum2Bio,comma,t4*NHum2Hum,comma,       &
     &    t4*NHum_min,comma, t4*NH4_miner,comma,                        &
     &    t4*NH4_old,comma,t4*NH4_end,comma,t4*(NH4_old-NH4_end),comma, &
!     &    t4*cNH4_t0*(WFrac_t0 + DryBD * SorpCoef)*dz_WSN, comma,       &
!     &    t4*cNH4_t *(WFrac_t + DryBD * SorpCoef)*dz_WSN, comma,        &
!     &    t4*( cNH4_t0*(WFrac_t0 + DryBD * SorpCoef) -                  &
!     &         cNH4_t *(WFrac_t  + DryBD * SorpCoef) )*dz_WSN, comma,   &
     &    t4*NH4N_amend,comma,t4*NH4N_cres,comma,t4*NH4_intop,comma,    &
     &    t4*NH4_inlat,comma,                                           &
     &    t4*NH4_inbot,comma,t4*NH4_upt,comma,t4*NH4_out,comma,         &
     &    t4*NH4N_volat,comma, t4*NH4_nitrif,comma,                     &
     &    t4*NO3_old,comma,t4*NO3_end,comma,t4*(NO3_old-NO3_end),comma, &
!     &    t4 * cNO3_t0 * WFrac_t0 * dz_WSN,  comma,                     &
!     &    t4 * cNO3_t * WFrac_t * dz_WSN, comma,                        &
!     &    t4 * (cNO3_t0 * WFrac_t0 - cNO3_t * WFrac_t) * dz_WSN, comma, &
     &    t4*NO3N_amend,comma,t4*NO3N_cres,comma,t4*NO3_intop,          &
     &    comma,t4*NO3_inlat,comma,t4*NO3_inbot,comma,t4*NO3_upt,comma, &
     &    t4*NO3_out,comma,t4*NO3_denitr, comma,t4*Ndemand, comma,      &
     &    t4*Nsupply, comma,                                            &
     &    t3*WFrac_t0*dz_WSN, comma, t3*WFrac_t*dz_WSN, comma,          &
     &    t3*(WFrac_t0 - WFrac_t)*dz_WSN, comma,                        &
     &    t3*Wflux_inTop * dt_WSN ,comma,                               &
     &    t3*Wflux_inLat * dt_WSN,comma,                                &
     &    t3*Wflux_inBot * dt_WSN, comma,t3*SoilEvap * dt_WSN, comma,   &
     &    t3*Wflux_transp * dt_WSN, comma,t3*Wflux_out * dt_WSN, comma, &
     &    idwrt_1, comma, idwlv_1, comma, idwst_1, comma,               &
     &    iNLOSSL_1, comma, iNLOSSR_1, comma, iNLOSSS_1, comma,         &
     &    idwso_1, comma,  iNLOSSO_1, comma,                            &
     &    WFPS, comma, red_T, comma, red_W, comma, red_W_Nit, comma,    &
     &    red_W_Den, comma, red_Resp

         FOM_add    = 0.0d0
         FOM_cres   = 0.0d0
         FOM2Bio    = 0.0d0
         FOM2Hum    = 0.0d0
         FOM_dis    = 0.0d0 

         Bio2Bio    = 0.0d0
         Bio2Hum    = 0.0d0
         Bio_dis    = 0.0d0

         Hum_add    = 0.0d0
         Hum_cres   = 0.0d0
         Hum2Bio    = 0.0d0
         Hum2Hum    = 0.0d0
         Hum_dis    = 0.0d0
 
         NFOM_add    = 0.0d0
         NFOM_cres   = 0.0d0
         NFOM2Bio    = 0.0d0
         NFOM2Hum    = 0.0d0
         NFOM_min    = 0.0d0 

         NBio2Bio    = 0.0d0
         NBio2Hum    = 0.0d0
         NBio_min    = 0.0d0

         NHum_add    = 0.0d0
         NHum_cres   = 0.0d0
         NHum2Bio    = 0.0d0
         NHum2Hum    = 0.0d0
         NHum_min    = 0.0d0
 
         NH4N_amend = 0.0d0
         NO3N_amend = 0.0d0
         NH4N_cres  = 0.0d0
         NO3N_cres  = 0.0d0
         NH4N_volat = 0.0d0
         idwrt_1 = idwrt
         idwlv_1 = idwlv
         idwst_1 = idwst
         idwso_1 = idwso
         iNLOSSL_1 = iNLOSSL
         iNLOSSR_1 = iNLOSSR
         iNLOSSS_1 = iNLOSSS
         iNLOSSO_1 = iNLOSSO

! outout to CROP_EXT file

      if ( floutput ) then

         write(cropext,'(a10,8(a1,1pe11.4))')                           &
     &         date,comma,Ntotuptake,comma,Ptotuptake,comma,DMcressur,  &
     &         comma,Ncressurf,comma,Pcressurf,comma,DMcresbott,        &
     &         comma,Ncresbott,comma,Pcresbott

         if(DMcressur.ge.1.0d-06)then
            nishmi = min(nishmi,Ncressurf/DMcressur)
            nishma = max(nishma,Ncressurf/DMcressur)
            poshmi = min(poshmi,Pcressurf/DMcressur)
            poshma = max(poshma,Pcressurf/DMcressur)
         end if
         if(DMcresbott.ge.1.0d-6)then
            niromi = min(niromi,Ncresbott/DMcresbott)
            niroma = max(niroma,Ncresbott/DMcresbott)
            poromi = min(poromi,Pcresbott/DMcresbott)
            poroma = max(poroma,Pcresbott/DMcresbott)
         end if
         Ntotuptake = 0.0d0
         Ptotuptake = 0.0d0
         DMcressur  = 0.0d0
         Ncressurf  = 0.0d0
         Pcressurf  = 0.0d0
         DMcresbott = 0.0d0
         Ncresbott  = 0.0d0
         Pcresbott  = 0.0d0
      end if

      return

      case (7)

      close(cropext)

      open(cropext,file=filnamce,status = 'old')
      iline = 0
      do while (.not. IS_IOSTAT_END(stat) )
         iline = iline + 1
         read(cropext,'(a)',IOSTAT=stat) cstring(iline)
      end do
      close(cropext)
      stat = 0

      open(cropext,file=filnamce,status = 'unknown')
      write(cropext,'(a)') cstring(1)
      write(cropext,'(8e11.4)')nishmi,nishma,niromi,niroma,             &
     &                        poshmi,poshma,poromi,poroma
      do i = 3,iline
         write(cropext,'(a)') cstring(i)
      end do
      close(cropext)

! -   open file with soil nutrient parameters
      filnam = trim(pathwork)//trim(project)//'.snp'
      snp = getun2 (10,90,2)
      open(unit=snp,file=filnam,status='old')

! -   open file to write soil nutrient parameters
      filnam = trim(pathwork)//trim(project)//'_nut.end'
      oup = getun2 (10,90,2)
      open(unit=oup,file=filnam,status='unknown')
!      do while (.not. eof(snp) )
      do while (.not. IS_IOSTAT_END(stat) )
         read(snp,'(a)',IOSTAT=stat)line
         if(line(1:6).eq.'FOM1_t')then
            write(oup,'("FOM1_t =",f12.6)')FOM_t(1)
         else if(line(1:6).eq.'FOM2_t')then
            write(oup,'("FOM2_t =",f12.6)')FOM_t(2)
         else if(line(1:6).eq.'FOM3_t')then
            write(oup,'("FOM3_t =",f12.6)')FOM_t(3)
         else if(line(1:6).eq.'FOM4_t')then
            write(oup,'("FOM4_t =",f12.6)')FOM_t(4)
         else if(line(1:6).eq.'FOM5_t')then
            write(oup,'("FOM5_t =",f12.6)')FOM_t(5)
         else if(line(1:6).eq.'FOM6_t')then
            write(oup,'("FOM6_t =",f12.6)')FOM_t(6)
         else if(line(1:6).eq.'FOM7_t')then
            write(oup,'("FOM7_t =",f12.6)')FOM_t(7)
         else if(line(1:6).eq.'FOM8_t')then
            write(oup,'("FOM8_t =",f12.6)')FOM_t(8)
         else if(line(1:5).eq.'Bio_t')then
            write(oup,'("Bio_t =",f12.6)')BIO_t
         else if(line(1:5).eq.'Hum_t')then
            write(oup,'("Hum_t =",f12.6)')HUM_t
         else if(line(1:6).eq.'cNH4_t')then
            write(oup,'("cNH4_t =",f12.6)')cNH4_t
         else if(line(1:6).eq.'cNO3_t')then
            write(oup,'("cNO3_t =",f12.6)')cNO3_t
         else
            write(oup,'(a)')line
         end if
      end do
      close(snp)
      close(oup)

      case default
         call fatalerr ('SoilManagement', 'Illegal value for TASK')
      end select

      return
      end

