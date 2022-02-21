! File VersionID:
!   $Id: wofost_soil_parameters.f90 323 2017-03-03 15:17:25Z kroes006 $
! ----------------------------------------------------------------------
      Subroutine Wofost_SoilParameters
!0    Declarations
      use variables
      use Wofost_Soil_Declarations
!0.3  intermediate local variables
      Integer :: fn


      nf = 8
      RateconFOM_ref(1) = 3.d0 / 365.2425d0
      RateconFOM_ref(2) = RateconFOM_ref(1)
      RateconFOM_ref(3) = RateconFOM_ref(1)
      RateconFOM_ref(4) = RateconFOM_ref(1)
      RateconFOM_ref(5) = 0.3d0 / 365.2425d0
      RateconFOM_ref(6) = RateconFOM_ref(5)
      RateconFOM_ref(7) = RateconFOM_ref(5)
      RateconFOM_ref(8) = RateconFOM_ref(5)
      RateconBio_ref    = 0.66d0 / 365.2425d0
      RateconHum_ref    = 0.02d0 / 365.2425d0

      AsfaMin      = 0.04d0
      AsfaMax      = 0.28d0
      AsfaFOM_Bio(1) = AsfaMin*0.46d0; AsfaFOM_Hum(1) = AsfaMin*0.54d0
      AsfaFOM_Bio(2) = AsfaMin*0.46d0; AsfaFOM_Hum(2) = AsfaMin*0.54d0
      AsfaFOM_Bio(3) = AsfaMax*0.46d0; AsfaFOM_Hum(3) = AsfaMax*0.54d0
      AsfaFOM_Bio(4) = AsfaMax*0.46d0; AsfaFOM_Hum(4) = AsfaMax*0.54d0
      AsfaFOM_Bio(5) = AsfaMin*0.46d0; AsfaFOM_Hum(5) = AsfaMin*0.54d0
      AsfaFOM_Bio(6) = AsfaMin*0.46d0; AsfaFOM_Hum(6) = AsfaMin*0.54d0
      AsfaFOM_Bio(7) = AsfaMax*0.46d0; AsfaFOM_Hum(7) = AsfaMax*0.54d0
      AsfaFOM_Bio(8) = AsfaMax*0.46d0; AsfaFOM_Hum(8) = AsfaMax*0.54d0
      AsfaBio = 0.46d0 * 0.2d0
      AsfaHum = 0.54d0 * 0.2d0

      do fn =1,nf
         CFracFOM(fn) = 0.45d0
      end do
      CFracBio = 0.58d0
      CFracHum = 0.58d0

      NFracFOMmin  = 0.002d0
      NFracFOMmax  = 0.06d0
      NFracFOM(1) = NFracFOMmax
      NFracFOM(2) = NFracFOMmin
      NFracFOM(3) = NFracFOMmax
      NFracFOM(4) = NFracFOMmin
      NFracFOM(5) = NFracFOMmax
      NFracFOM(6) = NFracFOMmin
      NFracFOM(7) = NFracFOMmax
      NFracFOM(8) = NFracFOMmin

      NFracBio = 0.58d0 / 10.d0
      NFracHum = 0.58d0 / 14.d0

      MatName(1)  = 'Cattle manure'            
      MatName(2)  = 'Pig manure'                
      MatName(3)  = 'Poultry manure'           
      MatName(4)  = 'Cattle slurry'            
      MatName(5)  = 'Pig slurry'               
      MatName(6)  = 'Poultry slurry'           
      MatName(7)  = 'Compost'                  
      MatName(8)  = 'Champost'                 
      MatName(9)  = 'Urea'                     
      MatName(10) = 'Mineral N fertilizer'     
      MatName(11) = 'Green leaves'             
      MatName(12) = 'Overground crop residues' 
      MatName(13) = 'Root and stubble residues'
      MatName(14) = 'Grass shoots'             
      MatName(15) = 'Grass roots'              
      MatName(16) = 'Tree leaves'              
      MatName(17) = 'Spruce needles'           

      AppAge(1)  =  3.16d0
      AppAge(2)  =  1.36d0
      AppAge(3)  =  1.36d0
      AppAge(4)  =  3.16d0
      AppAge(5)  =  1.36d0
      AppAge(6)  =  1.36d0
      AppAge(7)  =  1.96d0
      AppAge(8)  =  1.36d0
      AppAge(9)  =  0.00d0
      AppAge(10) =  0.00d0
      AppAge(11) =  0.92d0
      AppAge(12) =  0.99d0
      AppAge(13) =  1.57d0
      AppAge(14) =  0.92d0
      AppAge(15) =  1.20d0
      AppAge(16) =  2.25d0
      AppAge(17) =  3.34d0

      OrgMatFrac(1)  = 0.150d0
      OrgMatFrac(2)  = 0.161d0
      OrgMatFrac(3)  = 0.376d0
      OrgMatFrac(4)  = 0.064d0
      OrgMatFrac(5)  = 0.060d0
      OrgMatFrac(6)  = 0.093d0
      OrgMatFrac(7)  = 0.190d0
      OrgMatFrac(9)  = 0.460d0
      OrgMatFrac(10) = 0.000d0
      OrgMatFrac(11) = 1.000d0
      OrgMatFrac(12) = 1.000d0
      OrgMatFrac(13) = 1.000d0
      OrgMatFrac(14) = 1.000d0
      OrgMatFrac(15) = 1.000d0
      OrgMatFrac(16) = 1.000d0
      OrgMatFrac(17) = 1.000d0

      OrgNFrac(1)  = 0.035d0
      OrgNFrac(2)  = 0.037d0
      OrgNFrac(3)  = 0.058d0
      OrgNFrac(4)  = 0.034d0
      OrgNFrac(5)  = 0.050d0
      OrgNFrac(6)  = 0.047d0
      OrgNFrac(7)  = 0.041d0
      OrgNFrac(8)  = 0.025d0
      OrgNFrac(9)  = 0.000d0
      OrgNFrac(10) = 0.000d0
      OrgNFrac(11) = 0.000d0
      OrgNFrac(12) = 0.000d0
      OrgNFrac(13) = 0.000d0
      OrgNFrac(14) = 0.000d0
      OrgNFrac(15) = 0.000d0
      OrgNFrac(16) = 0.000d0
      OrgNFrac(17) = 0.000d0

      NH4NFrac(1)  = 0.0012d0
      NH4NFrac(2)  = 0.0015d0
      NH4NFrac(3)  = 0.0024d0
      NH4NFrac(4)  = 0.0022d0
      NH4NFrac(5)  = 0.0042d0
      NH4NFrac(6)  = 0.0058d0
      NH4NFrac(7)  = 0.0008d0
      NH4NFrac(8)  = 0.0003d0
      NH4NFrac(9)  = 0.4600d0
      NH4NFrac(10) = 0.5000d0
      NH4NFrac(11) = 0.0000d0
      NH4NFrac(12) = 0.0000d0
      NH4NFrac(13) = 0.0000d0
      NH4NFrac(14) = 0.0000d0
      NH4NFrac(15) = 0.0000d0
      NH4NFrac(16) = 0.0000d0
      NH4NFrac(17) = 0.0000d0

      NO3NFrac(1)  = 0.0000d0
      NO3NFrac(2)  = 0.0000d0
      NO3NFrac(3)  = 0.0000d0
      NO3NFrac(4)  = 0.0000d0
      NO3NFrac(5)  = 0.0000d0
      NO3NFrac(6)  = 0.0000d0
      NO3NFrac(7)  = 0.0000d0
      NO3NFrac(8)  = 0.0000d0
      NO3NFrac(9)  = 0.0000d0
      NO3NFrac(10) = 0.5000d0
      NO3NFrac(11) = 0.0000d0
      NO3NFrac(12) = 0.0000d0
      NO3NFrac(13) = 0.0000d0
      NO3NFrac(14) = 0.0000d0
      NO3NFrac(15) = 0.0000d0
      NO3NFrac(16) = 0.0000d0
      NO3NFrac(17) = 0.0000d0

      iAmendTime = 1

      Temp_ref      = 9.5d0
      DryBD         = BDENS(1)

      return
      End Subroutine Wofost_SoilParameters
