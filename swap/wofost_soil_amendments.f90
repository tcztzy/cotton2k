! File VersionID:
!   $Id: wofost_soil_amendments.f90 298 2016-07-25 20:07:31Z kroes006 $
! ----------------------------------------------------------------------
      Subroutine Wofost_SoilAmendents

!0    Declarations
      use Wofost_Soil_Declarations
!0.3  intermediate local variables
      integer :: im, matno, fn, ii
      real(8) :: Am_NH4, Am_NO3, Am_OM, Age, fDPM, fRPM, Asfa,          &
     &           OrgNFr, fOrgN1, fOrgN2, FHum, fAsfa1, fAsfa2
      real(8) :: cf, ff(8), help
!     Materials

      do ii = 1,NuAmend(isme)
         im = iamend(isme,ii)
             matno = MatNum(im)
         Am_NH4 = NH4NFrac(matno) * (1.0d0 - VolaFrac(im)) * Amend(im)    ! kg/m2
         Am_NO3 = NO3NFrac(matno) * Amend(im)                           ! kg/m2

         cf = WFrac_t0 + DryBD * SorpCoef
         cNH4_t0 = (cf * dz_WSN * cNH4_t0 + Am_NH4) /(cf* dz_WSN )
         cNO3_t0 = (WFrac_t0 * dz_WSN*cNO3_t0+Am_NO3)/(WFrac_t0*dz_WSN)

         Am_OM = OrgMatFrac(matno) * Amend(im)
         if( Am_OM .ge.1.0d-6 )then
            Age    = AppAge(matno)
            fDPM   = exp(-0.59d0*(Age-0.67d0))
            fHum   = min(1.0d0,max(0.0d0,0.137d0*(Age-2.5d0)))
            fHum   = min(fHum,(OrgNFrac(matno)-NFracFOMmin)/NFracHum)
            fRPM   = min(1.0d0,max(0.0d0,1.0d0 - fDPM - fHum))
            Asfa   = 0.25d0 /(1.0d0+exp(-2.7d0*(Age-2.0d0))) + 0.03d0
            fAsfa1 = (AsfaMax - Asfa) / (AsfaMax-AsfaMin)
            fAsfa2 = 1.0d0 - fAsfa1
             
            OrgNFr = (OrgNFrac(matno) - fHum * NFracHum) / (1.0d0-fHum)
            fOrgN1 = (OrgNFr - NFracFOMmin)/(NFracFOMmax - NFracFOMmin)
            fOrgN2 = 1.0d0 - fOrgN1
      
            ff(1) = fDPM * fAsfa1 * fOrgN1
            ff(2) = fDPM * fAsfa1 * fOrgN2
            ff(3) = fDPM * fAsfa2 * fOrgN1
            ff(4) = fDPM * fAsfa2 * fOrgN2
            ff(5) = fRPM * fAsfa1 * fOrgN1
            ff(6) = fRPM * fAsfa1 * fOrgN2
            ff(7) = fRPM * fAsfa2 * fOrgN1
            ff(8) = fRPM * fAsfa2 * fOrgN2

            help     = fHum * Am_OM / dz_WSN
            Hum_t0   = Hum_t0   + help                                  ! kg/m3
            Hum_add  = Hum_add + fHum * Am_OM                           ! kg/m2
            NHum_add = NHum_add + NFracHum * fHum * Am_OM               ! kg/m2
            do fn = 1,nf
               help       = ff(fn) * Am_OM / dz_WSN
               FOM_t0(fn) = FOM_t0(fn) + help                           ! kg/m3
               FOM_add    = FOM_add + ff(fn) * Am_OM                    ! kg/m2
               NFOM_add   = NFOM_add + NFracFOM(fn) * ff(fn) * Am_OM    ! kg/m2
            end do
         end if
         NH4N_amend = NH4N_amend + NH4NFrac(matno) * Amend(im)          ! kg/m2
         NO3N_amend = NO3N_amend + NO3NFrac(matno) * Amend(im)          ! kg/m2
         NH4N_volat = NH4N_volat +  NH4NFrac(matno) * VolaFrac(im)      &
     &                           * Amend(im)                            ! kg/m2
      end do

      return
      End Subroutine Wofost_SoilAmendents

