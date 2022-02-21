! File VersionID:
!   $Id: wofost_soil_rateconstants.f90 323 2017-03-03 15:17:25Z kroes006 $
! ----------------------------------------------------------------------
      Subroutine Wofost_SoilRateConstants(task)

!0    Declarations
      use Wofost_Soil_Declarations

!0.3  intermediate local variables
      Integer :: fn, task
      real(8)   :: Fred_T, Fred_W

      select case (Task)
      case (1)
 
!1    Assign reference values
      do fn = 1,nf
         RateconFOM(fn) = RateconFOM_ref(fn)
      end do
      RateconBio    = RateconBio_ref 
      if(RateconHum_exp .gt. 1.0d-6 .and. t1900Soil.ge.tstartHumexp)then
         RateconHum = RateconHum_exp
      else
         RateconHum = RateconHum_ref
      end if
      RateConNitrif = RateConNitrif_ref
      RateConDenitr = RateConDenitr_ref

!2    Correction for soil temperature
      red_T = Fred_T(Temp,Temp_ref)
      do fn = 1,nf
         RateconFOM(fn) = RateconFOM(fn) * red_T
      end do
      RateconBio = RateconBio * red_T
      RateconHum = RateconHum * red_T


!3    Correction for moisture conditions
      WFPS = 0.5d0 * (WFrac_t + WFrac_t0)/ WFrac_sat 
      red_W = Fred_W(WFPS, WFPScrit)
      do fn = 1,nf
         RateconFOM(fn) = RateconFOM(fn) * red_W
      end do

      RateconBio = RateconBio * red_W
      RateconHum = RateconHum * red_W
      RateConNitrif = RateConNitrif_ref
      RateConDenitr = RateConDenitr_ref
      RateConNitrif = RateConNitrif * red_T
      RateConDenitr = RateConDenitr * red_T
      red_W_Nit     = 0.9d0/(1.0d0+exp(-15.0d0*(WFPS-0.45d0)))          &
     &              + 0.1d0 - 1.0d0 /(1.0d0+exp(-50.0d0*(WFPS-0.95d0)))
      RateConNitrif = RateConNitrif * red_W_Nit
      red_W_Den     = (max(WFPS-WFPScrit2,0.0d0)/(1.0d0-WFPScrit2))**2
      RateConDenitr = RateConDenitr * red_W_Den

      return

      case (2)
 
!     account for respiration activity
      red_Resp      = Cdissi / (CdissiHalf + Cdissi)
      RateConDenitr = RateConDenitr * red_Resp

      case default
         call fatalerr ('Wofost_SoilRateConstants', 'Illegal value for TASK')
      end select

      return
      End Subroutine Wofost_SoilRateConstants

      Function Fred_T(Temp,Temp_ref)
      implicit none
!     response to temperature according to Rijtema et al, 1999
      real(8):: Fred_T, Temp, Temp_ref, r1, r2
      r1 = 1.0d0/(1.0d0+exp(-0.26d0*(Temp-17.0d0))) -                   &
     &     1.0d0/(1.0d0+exp(-0.77d0*(Temp-41.9d0)))
      r2 = 1.0d0/(1.0d0+exp(-0.26d0*(Temp_ref-17.0d0))) -               &
     &     1.0d0/(1.0d0+exp(-0.77d0*(Temp_ref-41.9d0)))
      Fred_T = r1/r2
      return
      End Function Fred_T

      Function Fred_W(WFPS, WFPScrit)
      implicit none
!     response to moisture conditions according to Bril et al (SONICG),
!     adjusted for near saturation conditions according to ANIMO (Groenendijk et al, 2005)
      real(8):: Fred_W, WFPS, WFPScrit, Denumer, A2, A1, A0, Recfanaer
      data Recfanaer/ 0.01d0 /

      If(Wfps .Le. Wfpscrit)Then
         Fred_W  = 6.d0*Wfps**2/(1.d0+9.d0*Wfps**4)
      Else
!     Reduce OM-transformation at (near) saturation circumstances
         Denumer = 81.d0*Wfpscrit**10-162.d0*Wfpscrit**9+81.d0          &
     &               *Wfpscrit**8+18.d0*Wfpscrit**6-36.d0*Wfpscrit**5   &
     &               +18.d0*Wfpscrit**4+Wfpscrit**2-2.d0*Wfpscrit+1.d0
         A2 = ((18.d0*Wfpscrit**4+81.d0*Wfpscrit**8+1.d0)*Recfanaer-    &
     &           162.d0*Wfpscrit**6+108.d0*Wfpscrit**5+6.d0*            &
     &           Wfpscrit**2-12.d0*Wfpscrit)/denumer
         A1 = -2.d0*Wfpscrit*(81.d0*Wfpscrit**8*Recfanaer-              &
     &           108.d0*Wfpscrit**6+54.d0*Wfpscrit**4                   &
     &           +18.d0*Wfpscrit**4*Recfanaer+Recfanaer-6.d0)/denumer
         A0 = Wfpscrit**2*(81.d0*Wfpscrit**8*Recfanaer-                 &
     &           216.d0*Wfpscrit**5+162.d0*Wfpscrit**4+                 &
     &           18.d0*Wfpscrit**4*Recfanaer+Recfanaer-6.d0)/denumer

         Fred_W = A2*Wfps**2+A1*Wfps+A0
      End If

      return
      End Function Fred_W

