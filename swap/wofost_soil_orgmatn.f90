! File VersionID:
!   $Id: wofost_soil_orgmatn.f90 298 2016-07-25 20:07:31Z kroes006 $
! ----------------------------------------------------------------------
      Subroutine Wofost_SoilOrgMatN

!0    Declarations
      use Wofost_Soil_Declarations

!0.3  intermediate local variables
      Integer :: fn
      real(8) :: p1, p2, p3, p4, p5, p111, p112, p121, p122, p211, p212,&
     &           p221, p222, Eval1, Eval2, help, dt,                    &
     &           Bio_av, Hum_av, RHS_bio, RHS_hum

      dt = dt_WSN

!1    Residue of Fresh Organic Matter fractions
      do fn = 1,nf
         FOM_t(fn) = FOM_t0(fn) * exp(-RateconFOM(fn) * dt) 
      end do

!2    Mutual dependent Biomass & Humus pools
!2.1  intermediate variables
      p1 = (1.0d0 - AsfaBio) * RateconBio
      p2 = AsfaBio * RateconHum
      p3 = AsfaHum * RateconBio
      p4 = (1.0d0 - AsfaHum) * RateconHum
!2.2  Eigenvalues of the mutual interacting Biomass - Humus system
      p5 = sqrt((p1-p4)**2+4*p2*p3)
      Eval1 = -(p1+p4+p5)/2.d0
      Eval2 = -(p1+p4-p5)/2.d0
!2.3  Coefficients of the Fundamental Matrix, third subscript refers to resp. Eival1 and Eval2
      p111 = -(p4+Eval1)/p5
      p112 =  (p4+Eval2)/p5
      p121 = -p2/p5
      p122 =  p2/p5
      p211 = -p3/p5
      p212 =  p3/p5
      p221 = -(p1+Eval1)/p5
      p222 =  (p1+Eval2)/p5
!2.4  Residue of Biomass and Humus
      Bio_t = (p111 * exp(Eval1*dt) + p112 * exp(Eval2*dt)) * Bio_t0 +  &
     &        (p121 * exp(Eval1*dt) + p122 * exp(Eval2*dt)) * Hum_t0
      Hum_t = (p211 * exp(Eval1*dt) + p212 * exp(Eval2*dt)) * Bio_t0 +  &
     &        (p221 * exp(Eval1*dt) + p222 * exp(Eval2*dt)) * Hum_t0
      do fn = 1,nf
         Bio_t = Bio_t + AsfaFOM_Bio(fn)* FOM_t0(fn) * RateconFOM(fn) * &
     &           ( p111 / (Eval1 + RateconFOM(fn) ) *                   &
     &           ( exp(Eval1*dt) - exp(-RateconFOM(fn)*dt) ) +          &
     &             p112 / (Eval2 + RateconFOM(fn) ) *                   &
     &           ( exp(Eval2*dt) - exp(-RateconFOM(fn)*dt) ) )
         Bio_t = Bio_t + AsfaFOM_Hum(fn)* FOM_t0(fn) * RateconFOM(fn) * &
     &           ( p121 / (Eval1 + RateconFOM(fn) ) *                   &
     &           ( exp(Eval1*dt) - exp(-RateconFOM(fn)*dt) ) +          &
     &             p122 / (Eval2 + RateconFOM(fn) ) *                   &
     &           ( exp(Eval2*dt) - exp(-RateconFOM(fn)*dt) ) )
         Hum_t = Hum_t + AsfaFOM_Bio(fn)* FOM_t0(fn) * RateconFOM(fn) * &
     &           ( p211 / (Eval1 + RateconFOM(fn) ) *                   &
     &           ( exp(Eval1*dt) - exp(-RateconFOM(fn)*dt) ) +          &
     &             p212 / (Eval2 + RateconFOM(fn) ) *                   &
     &           ( exp(Eval2*dt) - exp(-RateconFOM(fn)*dt) ) )
         Hum_t = Hum_t + AsfaFOM_Hum(fn)* FOM_t0(fn) * RateconFOM(fn) * &
     &           ( p221 / (Eval1 + RateconFOM(fn) ) *                   &
     &           ( exp(Eval1*dt) - exp(-RateconFOM(fn)*dt) ) +          &
     &             p222 / (Eval2 + RateconFOM(fn) ) *                   &
     &           ( exp(Eval2*dt) - exp(-RateconFOM(fn)*dt) ) )
      end do

!     Balance variables in kg/m2
      do fn = 1, nf

!        Gross dissociation of FOM
         help =  (FOM_t0(fn)-FOM_t(fn)) * dz_WSN

!        Incorporation of dissociated FOM into BIO
         FOM2Bio = FOM2Bio + AsfaFOM_Bio(fn) * help
         NFOM2Bio = NFOM2Bio + NFracBio * AsfaFOM_Bio(fn) * help

!        Incorporation of dissociated FOM into Hum
         FOM2Hum = FOM2Hum + AsfaFOM_Hum(fn) * help
         NFOM2Hum = NFOM2Hum + NFracHum * AsfaFOM_Hum(fn) * help

!        Nett dissociation of FOM
         FOM_dis = FOM_dis+(1.0d0-AsfaFOM_Bio(fn)-AsfaFOM_Hum(fn))*help

!        Nett N-mineralisation of FOM
         NFOM_min = NFOM_min + (NFracFOM(fn) - AsfaFOM_Bio(fn)*NFracBio &
     &                       -  AsfaFOM_Hum(fn)*NFracHum) * help 
      end do

!        calculate time averaged Bio and Hum values
      RHS_bio = (Bio_t - Bio_t0) / dt - FOM2Bio / dt / dz_WSN
      RHS_hum = (Hum_t - Hum_t0) / dt - FOM2Hum / dt / dz_WSN
      Bio_av = (p4*RHS_bio+p2*RHS_hum)/(p2*p3-p1*p4)
      Hum_av = (p3*RHS_bio+p1*RHS_hum)/(p2*p3-p1*p4)

!        Conversion of Bio into Bio
      help = dt * dz_WSN
      Bio2Bio  = Bio2Bio  + help * RateconBio * AsfaBio * Bio_av 
      NBio2Bio = NBio2Bio + NFracBio *help *RateconBio *AsfaBio *Bio_av

!        Conversion of Bio into Hum
      Bio2Hum  = Bio2Hum  + help * RateconBio * AsfaHum * Bio_av
      NBio2Hum = NBio2Hum + NFracHum *help *RateconBio *AsfaHum *Bio_av

!        Nett dissimilation of Bio
      Bio_dis  = Bio_dis+(1.0d0-AsfaBio-AsfaHum)*help*RateconBio*Bio_av
      NBio_min = NBio_min + (NFracBio - NFracBio*AsfaBio -              &
     &           NFracHum*AsfaHum) *help *RateconBio *Bio_av

!        Conversion of Hum into Bio
      Hum2Bio  = Hum2Bio  + help * RateconHum * AsfaBio * Hum_av
      NHum2Bio = NHum2Bio + NFracBio *help *RateconHum *AsfaBio *Hum_av

!        Conversion of Hum into Hum
      Hum2Hum  = Hum2Hum  + help * RateconHum * AsfaHum * Hum_av
      NHum2Hum = NHum2Hum + NFracHum *help *RateconHum *AsfaHum *Hum_av

!        Nett dissimilation of Hum
      Hum_dis  = Hum_dis+(1.0d0-AsfaBio-AsfaHum)*help*RateconHum*Hum_av
      NHum_min = NHum_min + (NFracHum - NFracBio*AsfaBio -              &
     &           NFracHum*AsfaHum) * help * RateconHum * Hum_av

!     Dissimilation of carbon
      Cdissi = 0.0d0
      do fn = 1,nf
         Cdissi = Cdissi + CFracFOM(fn) * (FOM_t0(fn) - FOM_t(fn)) *    &
     &            dz_WSN * (1.0d0-AsfaFOM_Bio(fn)-AsfaFOM_Hum(fn))
      end do
      Cdissi = Cdissi + CFracBio * dt * dz_WSN *                        &
     &                 (1.0d0-AsfaBio-AsfaHum) * RateconBio*Bio_av      &
     &                + CFracHum * dt * dz_WSN *                        &
     &                 (1.0d0-AsfaBio-AsfaHum) * RateconHum*Hum_av
 
!     Mineralisation of Nitrogen
      Nminer = 0.0d0
      do fn = 1,nf
         help   =  (FOM_t0(fn)-FOM_t(fn)) 
         Nminer = Nminer + (NFracFOM(fn) - AsfaFOM_Bio(fn)*NFracBio     &
     &                       -  AsfaFOM_Bio(fn)*NFracBio) * help 
      end do
      help = dt 
      Nminer = Nminer + (NFracBio - NFracBio*AsfaBio - NFracHum*AsfaHum)&
     &                * help * RateconBio * Bio_av
      Nminer = Nminer + (NFracHum - NFracBio*AsfaBio - NFracHum*AsfaHum)&
     &                * help * RateconHum * Hum_av

      return
      End Subroutine Wofost_SoilOrgMatN

