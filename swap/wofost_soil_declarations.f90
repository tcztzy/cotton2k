! File VersionID:
!   $Id: wofost_soil_declarations.f90 329 2017-06-12 08:09:14Z kroes006 $
! ----------------------------------------------------------------------
      module  Wofost_Soil_Declarations

      Integer, Parameter :: maxfn = 8
      Integer, Parameter :: maxmat = 20

! --- maximum number of soil management events, should be equal to 
!      PARAMETER  (MASME = 1000) in Arrays.Fi
      Integer, Parameter :: maxamn = 1000

!     Variables resulting from SWAP-simulation of soil water flow
      real(8) :: dz_WSN, WFrac_t, WFrac_t0, Wflux_out, Wflux_transp,    &
     &           Wflux_inBot, Wflux_inTop, Wflux_inLat, SoilEvap,       &
     &           Temp, t_WSNold

!     Parameters used in the response function for temperature and moisture conditions
      real(8) :: Temp_ref, WFrac_sat

!     Parameters used in the Organic matter (nitrogen) model
      Integer :: nf
      real(8) :: RateconFOM_ref(1:maxfn), RateconBio_ref, RateconHum_ref
      real(8) :: RateconHum_exp, tstartHumexp, t1900Soil
      real(8) :: RateconFOM(1:maxfn),     RateconBio,     RateconHum
      real(8) :: AsfaFOM_Bio(1:maxfn), AsfaFOM_Hum(1:maxfn), AsfaBio,   &
     &           AsfaHum
      real(8) :: CFracFOM(1:maxfn),       CFracBio,       CFracHum
      real(8) :: NFracFOM(1:maxfn),       NFracBio,       NFracHum
      real(8) :: NFracFOMmin, NFracFOMmax, AsfaMin, AsfaMax



!     Rate and State variables of the Organic matter model
      real(8) :: FOM_t0(1:maxfn), Bio_t0, Hum_t0
      real(8) :: FOM_t(1:maxfn),  Bio_t,  Hum_t, Cdissi
 
!     Parameters used by the mineral nitrogen model
      real(8) :: RateConNitrif_ref, RateConDenitr_ref
      real(8) :: RateConNitrif,     RateConDenitr,  TCSF_N



!     Rate and State variables of the mineral nitrogen model
      real(8) :: Nminer, Ratecon

!     Ammonium and nitrate solute concentration parameters and variables

      real(8) :: DryBD, SorpCoef, cNH4_t0, cNH4_t, cNH4_av, cNO3_t0,    &
     &           cNO3_t, cNO3_av

!     Concentrations at system boundaries
      real(8) :: cNH4N_top, cNH4N_lat, cNH4N_seep
      real(8) :: cNO3N_top, cNO3N_lat, cNO3N_seep
      real(8) :: Cseep, Ctop, Clat

!     Time control
      real(8) :: dt_WSN

!     Materials involved in soil amendments
      integer      :: iamend(maxamn,maxamn), nuamend(maxamn), namend
      integer      :: MatNum(maxamn)
      integer      :: iAmendTime

      character(len=25) :: MatName(maxmat)
      real(8)      :: Amend(maxamn)
      real(8)      :: AppAge(maxmat), OrgMatFrac(maxmat),               &
     &              OrgNFrac(maxmat), NH4NFrac(maxmat), NO3NFrac(maxmat)
      real(8)      :: VolaFrac(maxamn), NH4N_volat
      real(8)      ::  NH4N_amend, NO3N_amend
      real(8)      ::  NH4N_cres, NO3N_cres
      real(8)      ::  TimeAmend(maxamn)
      integer      :: isme

!     Crop residues amendment
      real(8):: iNLOSSL_1, iNLOSSR_1, iNLOSSS_1, iNLOSSO_1
      real(8):: idwrt_1,   idwlv_1,   idwst_1,   idwso_1


!     Response function parameters
      real(8)::     WFPSCrit, WFPScrit2, CdissiHalf
      real(8)::     WFPS, red_T, red_W, red_W_Nit, red_W_Den, red_Resp


!     Balance check variables
      real(8)::  FOM_add, NFOM_add, Hum_add, NHum_add
      real(8)::  FOM_cres, NFOM_cres, Hum_cres, NHum_cres
      real(8)::  FOM_end, NFOM_end, Hum_end, NHum_end, Bio_end, NBio_end
      real(8)::  NH4_end, NH4_intop, NH4_inlat, NH4_inbot, NH4_upt,     &
     &           NH4_out, NH4_nitrif, NH4_miner
      real(8)::  NO3_end, NO3_intop, NO3_inlat, NO3_inbot, NO3_upt,     &
     &           NO3_out, NO3_denitr
      real(8)::  Inp_orgm, Inp_orgN, Out_orgm, Out_orgN, dST_orgm,      &
     &           dST_orgN, dif_orgm, dif_orgN
      real(8)::  Inp_NH4N, Out_NH4N, dST_NH4N, dif_NH4N
      real(8)::  Inp_NO3N, Out_NO3N, dST_NO3N, dif_NO3N

!        Gross and nett dissimilation / mineralisation
      real(8)::  FOM2Bio, NFOM2Bio, FOM2Hum, NFOM2Hum, FOM_dis, NFOM_min
      real(8)::  Bio2Bio, NBio2Bio, Bio2Hum, NBio2Hum, Bio_dis, NBio_min
      real(8)::  Hum2Bio, NHum2Bio, Hum2Hum, NHum2Hum, Hum_dis, NHum_min

      real(8)::  NsupplyNH4N, NsupplyNO3N
      real(8)::  FOM_old, Bio_old, Hum_old, NFOM_old, NBio_old, NHum_old
      real(8)::  NH4_old, NO3_old


!        unit number for printing nutrient results

!     variables for ANIMO crop_ext file
      Integer :: nut, cropext, iline
      real(8) :: Ntotuptake, Ptotuptake
      real(8) :: DMcressur, Ncressurf, Pcressurf
      real(8) :: DMcresbott,Ncresbott, Pcresbott
      real(8) :: nishmi,nishma,niromi,niroma
      real(8) :: poshmi,poshma,poromi,poroma
      real(8) :: pnratio
      character(len=120) :: cstring(10000)

      save

      end module  Wofost_Soil_Declarations
