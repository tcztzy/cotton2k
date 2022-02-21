! File VersionID:
!   $Id: wofost_soil_balancecheck.f90 298 2016-07-25 20:07:31Z kroes006 $
! ----------------------------------------------------------------------
      Subroutine Wofost_SoilBalanceCheck

!0    Declarations
      use Wofost_Soil_Declarations
!0.3  intermediate local variables
      integer :: fn

      FOM_old = FOM_end                                                 ! kg/m2
      Bio_old = Bio_end                                                 ! kg/m2
      Hum_old = Hum_end                                                 ! kg/m2
      NFOM_old = NFOM_end                                               ! kg/m2
      NBio_old = NBio_end                                               ! kg/m2
      NHum_old = NHum_end                                               ! kg/m2
      NH4_old = NH4_end                                                 ! kg/m2
      NO3_old = NO3_end


      FOM_end    = 0.0d0                                                ! kg/m2
      NFOM_end   = 0.0d0                                                ! kg/m2
      do fn = 1,nf
         FOM_end    = FOM_end + FOM_t(fn) * dz_WSN
         NFOM_end   = NFOM_end + NFracFOM(fn) * FOM_t(fn) * dz_WSN
      end do
      Bio_end       = Bio_t  * dz_WSN                                   ! kg/m2
      NBio_end      = NFracBio * Bio_t * dz_WSN                         ! kg/m2
      Hum_end       = Hum_t * dz_WSN                                    ! kg/m2
      NHum_end      = NFracHum * Hum_t * dz_WSN                         ! kg/m2

      NH4_end       = cNH4_t * (WFrac_t + DryBD * SorpCoef) * dz_WSN    ! kg/m2
      NH4_miner     = Nminer * dz_WSN                                   ! kg/m2
      NH4_intop     = cNH4N_top  * Wflux_inTop * dt_WSN                 ! kg/m2
      NH4_inlat     = cNH4N_lat  * Wflux_inLat * dt_WSN                 ! kg/m2
      NH4_inbot     = cNH4N_seep * Wflux_inBot * dt_WSN                 ! kg/m2
      NH4_upt       = NsupplyNH4N * dz_WSN                              ! kg/m2
      NH4_out       = cNH4_av * Wflux_out * dt_WSN                      ! kg/m2
      NH4_nitrif    = 0.5 * (WFrac_t+WFrac_t0) * RateConNitrif *        &  ! kg/m2
     &                cNH4_av * dt_WSN * dz_WSN

      NO3_end       = cNO3_t * WFrac_t * dz_WSN                         ! kg/m2
      NO3_intop     = cNO3N_top  * Wflux_inTop * dt_WSN                 ! kg/m2
      NO3_inlat     = cNO3N_lat  * Wflux_inLat * dt_WSN                 ! kg/m2
      NO3_inbot     = cNO3N_seep * Wflux_inBot * dt_WSN                 ! kg/m2
      NO3_upt       = NsupplyNO3N * dz_WSN                              ! kg/m2
      NO3_out       = cNO3_av * Wflux_out * dt_WSN                      ! kg/m2
      NO3_denitr    = 0.5 * (WFrac_t+WFrac_t0) * RateConDenitr *        &  ! kg/m2
     &                cNO3_av  * dt_WSN * dz_WSN

!     Organic matter and organic N balance
      Inp_orgm = FOM_add  +  Hum_add + FOM_cres  +  Hum_cres            ! kg/m2
      Inp_orgN = NFOM_add +  NHum_add + NFOM_cres  +  NHum_cres         ! kg/m2
      Out_orgm = FOM_dis  +  Bio_dis  +  Hum_dis 
      Out_orgN = NFOM_min +  NBio_min +  NHum_min 
      dST_orgm = (FOM_end - FOM_old) + (Bio_end - Bio_old) +            &
     &           (Hum_end - Hum_old)
      dST_orgN = (NFOM_end - NFOM_old) + (NBio_end - NBio_old) +        &
     &           (NHum_end - NHum_old)
      dif_orgm =  dST_orgm - Out_orgm + Inp_orgm   
      dif_orgN =  dST_orgN - Out_orgN + Inp_orgN

!    Ammonium N balance
      Inp_NH4N = NH4N_amend + NH4N_cres + NH4_miner + NH4_intop +       &
     &           NH4_inlat + NH4_inbot
      Out_NH4N = NH4N_volat + NH4_upt + NH4_out + NH4_nitrif
      dST_NH4N = NH4_end - NH4_old
      dif_NH4N = dST_NH4N - Out_NH4N + Inp_NH4N 

!    Nitrate N balance
      Inp_NO3N = NO3N_amend + NO3N_cres + NO3_intop + NO3_inlat +       &
     &           NO3_inbot + NH4_nitrif 
      Out_NO3N = NO3_upt + NO3_out + NO3_denitr 
      dST_NO3N = NO3_end - NO3_old
      dif_NO3N = dST_NO3N - Out_NO3N + Inp_NO3N 
   
      return
      End Subroutine Wofost_SoilBalanceCheck
