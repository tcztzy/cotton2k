!   File VersionID:
!   $Id: wofostnut.f90 323 2017-03-03 15:17:25Z kroes006 $
! ----------------------------------------------------------------------
!*
!*     Copyright 1988, 2011 Alterra, Wageningen-UR
!*
!*     Licensed under the EUPL, Version 1.1 or as soon they
!*     will be approved by the European Commission - subsequent
!*     versions of the EUPL (the "Licence");
!*     You may not use this work except in compliance with the
!*     Licence.
!*     You may obtain a copy of the Licence at:
!*
!*     http://www.osor.eu/eupl
!*
!*     Unless required by applicable law or agreed to in
!*     writing, software distributed under the Licence is
!*     distributed on an "AS IS" basis,
!*     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
!*     express or implied.
!*     See the Licence for the specific language governing
!*     permissions and limitations under the Licence.
!*
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
! ----------------------------------------------------------------------
      SUBROUTINE NUTRIE (TRANRF,NLUE,WLV,WST,DVS,                       &
     &                   ANLV,ANST,NMXLV,                               &
     &                   NMAXLV,NMAXST,NMAXRT,LRNR,LSNR,                &
     &                   NNI,RNFLV,RNFST,FRNX,FSTR)


!* 3.1 declarations
      implicit none
!     global
!      integer  ILNMXL
      real(8)  TRANRF,NMXLV(30),NLUE                  ! ,TRA,TRAMX
      real(8)  WLV,WST,DVS,ANLV,ANST                  ! ,WRT,WSO,ANRT
      real(8)  NMAXLV,NMAXST,NMAXRT,LRNR,LSNR
      real(8)  NNI,RNFLV,RNFST,FRNX,FSTR              ! ,RNFRT,FNTRT
!     local
      real(8)  afgen
      real(8)  TBGMR,NOPTMR,NUPGMR,NRMR,NFGMR
      real(8)  TINY,NREF
      SAVE        

!     vegetative living above-ground biomass (kg DM ha-1)
      TBGMR = WLV + WST     
        
!     Maximum N concentration in the leaves, from which the N conc. in the
!     stems and roots are derived, as a function of development stage (kg N kg-1 DM)  
      NMAXLV = afgen (NMXLV, 30, DVS)
      
!---- Maximum N concentrations in stems and roots (kg N kg-1 DM)
      NMAXST = LSNR * NMAXLV
      NMAXRT = LRNR * NMAXLV

      CALL NOPTM(WLV,WST,FRNX,NMAXLV,NMAXST,NOPTMR)
           
!     Total N in vegetative living above-ground biomass (kg N ha-1)
      NUPGMR = ANLV + ANST
      
!     N concentration in total vegetative living per kg above-ground biomass  (kg N kg-1 DM) 
      NFGMR  = 1.0d-5
      if(TBGMR.gt.1.0d-5)  NFGMR = NUPGMR/TBGMR

!     Residual N concentration in total vegetative living above-ground biomass  (kg N kg-1 DM) 
      NRMR = 1.0d-5
      if(TBGMR.gt.1.0d-5)  NRMR = (WLV*RNFLV+WST*RNFST)/TBGMR

!     nutrient stress factor           
      TINY=0.001d0
      NNI = TINY
      if((NOPTMR-NRMR).gt.1.0d-5)  then
        NNI = max(TINY, min(1.0d0, ((NFGMR-NRMR)/(NOPTMR-NRMR))))
      endif

!     Nutrient reduction factor     
      NREF= max(0.0d0, min(1.0d0, (1.d0-NLUE*(1.0001d0-NNI)**2)))
      
      IF(TRANRF .LE. NREF) THEN
!     Water stress is more severe than nutrient stress 
         FSTR = TRANRF
      ELSE
!     Nutrient stress is more severe than water stress 
         FSTR = NREF
      ENDIF   
      
      RETURN
      END
    
!*----------------------------------------------------------------------*
!*  SUBROUTINE SUBPAR                                                   *
!*  Purpose: Modification of dry matter partitioning to leaves, stems,  *
!*        roots and storage organs in dependence of water and N stress  *
!*----------------------------------------------------------------------*
      SUBROUTINE SUBPAR (TRANRF,NPART,NNI,FR,FL,FS,FO)

!     declarations
      implicit none
!     global
      real(8)  NPART,TRANRF,NNI,FR,FL,FS,FO
!     local
      real(8)  FRTMOD,FLOLD,FLVMOD

      SAVE        
     
      IF(TRANRF .LT. NNI) THEN

!*     Water stress is more severe as compared to nutrient stress and
!*     partitioning will follow the original assumptions of LINTUL2

        FRTMOD = MAX( 1.0d0, 1.0d0/(TRANRF+0.5d0))
        FR    = MIN(0.6d0, FR * FRTMOD)
        FL    = FL 
        FS    = FS
        FO    = FO 

      ELSE

!*    Nutrient stress is more severe as compared to water stress resulting in
!*    less partitioning to leaves and more to stems

        FLOLD = FL
        FLVMOD = EXP(-NPART* (1.0d0-NNI))
        FL    = FLOLD * FLVMOD
        FS    = FS + FLOLD - FL 
        FR    = FR
        FO    = FO 

      ENDIF
      
      RETURN
      END


!* ---------------------------------------------------------------------*
!*  SUBROUTINE NOPTM                                                    *
!*  Purpose: To compute the maximum N crop concentration            *
!*  organs (kg N kg-1 DM)                                               *
!* ---------------------------------------------------------------------*

      SUBROUTINE NOPTM(WLV,WST,FRNX,NMAXLV,NMAXST,NOPTMR)

!     declarations
      implicit none
!     global
      real(8)  WLV,WST,FRNX,NMAXLV,NMAXST,NOPTMR
!     local
      real(8)  NOPTL,NOPTS,TBGMR
      SAVE        
      
!    * Total vegetative living above-ground biomass (kg DM ha-1)     
      TBGMR = WLV + WST 
      
!     Optimal N amount in vegetative above-ground living biomass and its N concentration
      NOPTL  = FRNX * NMAXLV * WLV
      NOPTS  = FRNX * NMAXST * WST

      NOPTMR  = 1.0d-5
      if(TBGMR.gt.1.0d-5)  NOPTMR = (NOPTL + NOPTS)/TBGMR
      
      RETURN
      END

!* -------------------------------------------------------------------------*
!*  SUBROUTINE NDEMND                                                       *
!*  Purpose: To compute the nutrient demand of crop organs (kg N ha-1)  *
!* -------------------------------------------------------------------------*

      SUBROUTINE NDEMND(WLV,WST,WRT,WSO,                                &
     &               NMAXLV,NMAXST,NMAXRT,NMAXSO,                       &
     &               ANLV,ANST,ANRT,ANSO,TCNT,NDEML,NDEMS,              &
     &               NDEMR,NDEMSO)
!     declarations
      implicit none
!     global
      real(8)  WLV,WST,WRT,WSO,NMAXLV,NMAXST,NMAXRT,NMAXSO
      real(8)  ANLV,ANST,ANRT,ANSO,TCNT,NDEML,NDEMS,NDEMR,NDEMSO
      SAVE        

      NDEML  =  MAX (NMAXLV*WLV  - ANLV, 0.0d0)
      NDEMS  =  MAX (NMAXST*WST  - ANST, 0.0d0)
      NDEMR  =  MAX (NMAXRT*WRT  - ANRT, 0.0d0)
      NDEMSO =  MAX (NMAXSO*WSO  - ANSO, 0.0d0)/TCNT

      RETURN
      END      
      
      
!* ------------------------------------------------------------------------------------------------*
!*  SUBROUTINE NTRLOC                                                                              *
!*  Purpose: To compute the amount of nutrients in the organs that can be translocated(kg N ha-1)  *
!* ------------------------------------------------------------------------------------------------*

      SUBROUTINE NTRLOC(ANLV,ANST,ANRT,WLV,WST,WRT,RNFLV,RNFST,RNFRT,   &
     &                  FNTRT,ATNLV,ATNST,ATNRT,ATN)
!     declarations
      implicit none
!     global
      real(8)  ANLV,ANST,ANRT,WLV,WST,WRT,RNFLV,RNFST,RNFRT
      real(8)  FNTRT,ATNLV,ATNST,ATNRT,ATN
      SAVE        

      ATNLV = MAX (0.d0 , ANLV-WLV*RNFLV)
      ATNST = MAX (0.d0 , ANST-WST*RNFST)
      ATNRT = MAX((ATNLV + ATNST) * FNTRT, ANRT-WRT*RNFRT)
      ATN   = ATNLV +  ATNST + ATNRT
            
      RETURN
      END
      
!* ------------------------------------------------------------------------------------------------*
!*  SUBROUTINE NTRANS                                                                              *
!*  Purpose: To compute the amount of nutrients translocated from different organs (kg N ha-1 d-1) *
!* ------------------------------------------------------------------------------------------------*

      SUBROUTINE NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN,RNTLV,RNTST,RNTRT)
!     declarations
      implicit none
!     global
      real(8)  RNSO,ATNLV,ATNST,ATNRT,ATN,RNTLV,RNTST,RNTRT
!     local
      SAVE        

      RNTLV = 0.d0
      if(ATN.gt.1.0d-5)  RNTLV= RNSO* ATNLV/ ATN
      RNTST = 0.d0
      if(ATN.gt.1.0d-5)  RNTST= RNSO* ATNST/ ATN
      RNTRT = 0.d0
      if(ATN.gt.1.0d-5)  RNTRT= RNSO* ATNRT/ ATN
           
      RETURN
      END
      

!* -------------------------------------------------------------------------*
!*  SUBROUTINE RNUSUB                                                       *
!*  Purpose: To compute the partitioning of the total N uptake rate     *
!*           (NUPTR) over leaves, stem, and roots (kg N ha-1 d-1)       *
!* -------------------------------------------------------------------------*

      SUBROUTINE RNUSUB(NDEML,NDEMS,NDEMR,NUPTR,                        &
     &                  NFIXTR,NDEMTO,RNULV,RNUST,RNURT)
!     declarations
      implicit none
!     global
      real(8)  NDEML,NDEMS,NDEMR,NUPTR,NFIXTR,NDEMTO,RNULV,RNUST,RNURT
!     local

      SAVE        
      
      RNULV = 0.d0
      RNUST = 0.d0
      RNURT = 0.d0
      if(NDEMTO.gt.1.0d-5) then
        RNULV = (NDEML / NDEMTO)* (NUPTR+NFIXTR)
        RNUST = (NDEMS / NDEMTO)* (NUPTR+NFIXTR)
        RNURT = (NDEMR / NDEMTO)* (NUPTR+NFIXTR)
      endif

      RETURN
      END      
      
      
!* ------------------------------------------------------------------------------*
!*  SUBROUTINE RNLD                                                              *
!*  Purpose: To compute the nutrient losses due to dying leaves, stems           *
!*               and roots (kg N ha-1 d-1)                                   *
!* ------------------------------------------------------------------------------*

      SUBROUTINE RNLD(DRLV,DRRT,DRST,                                   &
     &                RNFLV,RNFRT,RNFST,                                &
     &                RNLDLV,RNLDRT,RNLDST)
!     declarations
      implicit none
!     global
      real(8)  DRLV,DRRT,DRST,RNFLV,RNFRT,RNFST,RNLDLV,RNLDRT,RNLDST

      SAVE        

      RNLDLV= RNFLV * DRLV
      RNLDRT= RNFRT * DRRT
      RNLDST= RNFST * DRST
      
      RETURN
      END

!* --------------------------------------------------------------------*
!*  SUBROUTINE NUTRINIT                                                *
!*  Purpose: Initialization of nutrient parameters                     *
!* --------------------------------------------------------------------*

      SUBROUTINE NUTRINIT(NLOSSL,NLOSSR,NLOSSS,                         &
     &                    NUPTT,RNLV,RNST,RNRT,RNSO,                    &
     &                    RNLDLV,RNLDST,RNLDRT,                         &
     &                    NFIXTT,NNI,NLOSSLDeceasedLvToSoil)
!     declarations
      implicit none
!     global
      real(8)  NLOSSL,NLOSSR,NLOSSS,NUPTT,RNLV,RNST,RNRT,RNSO
      real(8)  RNLDLV,RNLDST,RNLDRT,NFIXTT,NNI,NLOSSLDeceasedLvToSoil


      SAVE        
      
      NLOSSL = 0.0d0
      NLOSSR = 0.0d0
      NLOSSS = 0.0d0
      NUPTT  = 0.0d0
      NFIXTT = 0.0d0

      RNLV   = 0.0d0
      RNST   = 0.0d0
      RNRT   = 0.0d0
      RNSO   = 0.0d0
        
      NNI = 1.0d0

      RNLDST = 0.0d0
      RNLDLV = 0.0d0
      RNLDRT = 0.0d0  
        
      NLOSSLDeceasedLvToSoil = 0.0d0
    
      RETURN
      END
      
      
!* --------------------------------------------------------------------*
!*  SUBROUTINE NUTREMRG                                                *
!*  Purpose: Initialization of nutrient parameters at emergence        *
!* --------------------------------------------------------------------*

      SUBROUTINE NUTREMRG(NMXLV,LSNR,LRNR,WLV,WST,WRT,                  &
     &      ANLV,ANST,ANRT,ANSO,ANLVI,ANSTI,ANRTI,ANSOI,DVS)

      implicit none
!     global
!      integer  ILNMXL  
      real(8)  NMXLV(30),LSNR,LRNR,WLV,WST,WRT,DVS
      real(8)  ANLV,ANST,ANRT,ANSO,ANLVI,ANSTI,ANRTI,ANSOI
!     local
      real(8)  NMAXLVI,NMAXSTI,NMAXRTI,afgen
      SAVE        

      NMAXLVI= afgen (NMXLV, 30, DVS)
      NMAXSTI= LSNR * NMAXLVI
      NMAXRTI= LRNR * NMAXLVI
            
!     initial maximum N concentration in plant organs [kg N ]            
      ANLVI= NMAXLVI * WLV
      ANSTI= NMAXSTI * WST
      ANRTI= NMAXRTI * WRT
      ANSOI= 0.0d0
            
      ANLV = ANLVI
      ANST = ANSTI
      ANRT = ANRTI
      ANSO = ANSOI         
                    
      RETURN
      END
      
      
!* ---------------------------------------------------------------*
!*  SUBROUTINE NUTRSOW                                            *
!*  Purpose: Initialization of nutrient parameters at sowing      *
!* ---------------------------------------------------------------*      
      SUBROUTINE NUTRSOW(ANLV,ANST,ANRT,ANSO)
     
      implicit none
! global
      real(8) ANLV,ANST,ANRT,ANSO
! local
        
      ANLV= 0.0d0
      ANST= 0.0d0
      ANRT= 0.0d0
      ANSO= 0.0d0
            
            
      RETURN
      END
