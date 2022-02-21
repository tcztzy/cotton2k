! File VersionID:
!   $Id: wofost_soil_cropresidues.f90 323 2017-03-03 15:17:25Z kroes006 $
! ----------------------------------------------------------------------
      Subroutine Wofost_Soil_CropResidues(xAmend,xAppAge,               &
     &                    xOrgMatFrac,xOrgNFrac,xNH4NFrac,xNO3NFrac)

!0    Declarations
      use Wofost_Soil_Declarations
!      integer :: imat
      real(8) :: xAmend, xAppAge, xOrgMatFrac, xOrgNFrac, xNH4NFrac,    &
     &           xNO3NFrac
!0.3  intermediate local variables
      integer :: fn  !, im, matno
      real(8) :: Am_NH4, Am_NO3, Am_OM, Age, fDPM, fRPM, Asfa,          &
     &           OrgNFr, fOrgN1, fOrgN2, FHum, fAsfa1, fAsfa2
      real(8) :: cf, ff(8), help
!     Materials

      Am_NH4 = xNH4NFrac * xAmend                                       ! kg/m2
      Am_NO3 = xNO3NFrac * xAmend                                       ! kg/m2

      cf = WFrac_t0 + DryBD * SorpCoef
      cNH4_t0 = (cf * dz_WSN * cNH4_t0 + Am_NH4) /(cf* dz_WSN )
      cNO3_t0 = (WFrac_t0 * dz_WSN*cNO3_t0+Am_NO3)/(WFrac_t0*dz_WSN)

      Am_OM = xOrgMatFrac * xAmend
      if( Am_OM .ge.1.0d-12 )then
          Age    = xAppAge
          fDPM   = exp(-0.59d0*(Age-0.67d0))
          fHum   = min(1.0d0,max(0.0d0,0.137d0*(Age-2.5d0)))
          fHum   = min(fHum,(xOrgNFrac-NFracFOMmin)/NFracHum)
          fRPM   = min(1.0d0,max(0.0d0,1.0d0 - fDPM - fHum))
          Asfa   = 0.25d0 /(1.0d0+exp(-2.7d0*(Age-2.0d0))) + 0.03d0
          fAsfa1 = (AsfaMax - Asfa) / (AsfaMax-AsfaMin) 
          fAsfa2 = 1.0d0 - fAsfa1
             
          OrgNFr = (xOrgNFrac - fHum * NFracHum) / (1.0d0 - fHum)
          fOrgN1 = (OrgNFr- NFracFOMmin)/(NFracFOMmax - NFracFOMmin)
          fOrgN2 = 1.0d0 - fOrgN1
      
          ff(1) = fDPM * fAsfa1 * fOrgN1
          ff(2) = fDPM * fAsfa1 * fOrgN2
          ff(3) = fDPM * fAsfa2 * fOrgN1
          ff(4) = fDPM * fAsfa2 * fOrgN2
          ff(5) = fRPM * fAsfa1 * fOrgN1
          ff(6) = fRPM * fAsfa1 * fOrgN2
          ff(7) = fRPM * fAsfa2 * fOrgN1
          ff(8) = fRPM * fAsfa2 * fOrgN2

          help      = fHum * Am_OM / dz_WSN
          Hum_t0    = Hum_t0   + help                                   ! kg/m3
          Hum_cres  = Hum_cres + fHum * Am_OM                           ! kg/m2
          NHum_cres = NHum_cres + NFracHum * fHum * Am_OM               ! kg/m2
          do fn = 1,nf
             help        = ff(fn) * Am_OM / dz_WSN
             FOM_t0(fn)  = FOM_t0(fn) + help                            ! kg/m2
             FOM_cres    = FOM_cres + ff(fn) * Am_OM                    ! kg/m2
             NFOM_cres   = NFOM_cres + NFracFOM(fn) * ff(fn) * Am_OM    ! kg/m2
          end do
      end if
      NH4N_cres = NH4N_cres + xNH4NFrac * xAmend                        ! kg/m2
      NO3N_cres = NO3N_cres + xNO3NFrac * xAmend                        ! kg/m2

      return
      End Subroutine Wofost_Soil_CropResidues

