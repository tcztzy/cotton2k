! File VersionID:
!   $Id: wofost_soil_watern.f90 298 2016-07-25 20:07:31Z kroes006 $
! ----------------------------------------------------------------------
      Subroutine Wofost_SoilWaterN(dz, dt, WFrac_t, WFrac_t0,           &
     &                 Wflux_out, Wflux_transp,Wflux_inBot, Wflux_inTop,&
     &                 Wflux_inLat,TCSF, Ratecon, ProducPot, ProducAct, &
     &                 DryBD,SorpCoef,Cseep,Ctop, Clat, C_t0, C_t, C_av)

!0    Declarations
      implicit none
!0.1  inputs
      real(8) ::  WFrac_t, WFrac_t0, Wflux_out, Wflux_transp,           &
     &            Wflux_inBot, Wflux_inTop, Wflux_inLat
      real(8) ::  C_t0, Cseep, Ctop, Clat
      real(8) ::  dz, dt, TCSF, Ratecon, ProducPot, DryBD, SorpCoef 

!0.2  Outputs
      real(8) :: C_t, C_av, ProducAct

!0.3  intermediate local variables
      Integer :: Iflsol
      real(8) :: A1, A2, B1, B2
      real(8) :: WFrac_av, C_av1
      real(8) :: Hv, Hv1, Hv2, Hv3, Ttry
      real(8) :: Small, Vsmall, Vvsmall, Zero, Half, One
      Logical :: Hv1nil,Hvnil

      Data       Zero, Small, Vsmall, Vvsmall, Half, One                &
     &           /0.0d+00, 1.0d-08, 1.0d-12, 1.0d-20, 0.5d+00, 1.0d+00/
      

!        help variables in analytical solution and test for small values
!     HV       rate of change in moisture content [d-1]
!     HV1      rate of disappearance (outflow+decomposition) [d-1]
!     HV2      total nett inflow of substance per timestep [kg m-1 d-1}

      ProducAct = ProducPot

      WFrac_av = Half * (WFrac_t + WFrac_t0)
      Hv = (WFrac_t - WFrac_t0)/dt
      Hvnil = .False.
      If( abs(Hv) .Lt. Small )Then
        Hv    = Zero
        Hvnil = .True.
      End If

      Hv1 = HV + Wflux_out/dz + TCSF*Wflux_transp/dz + Ratecon*WFrac_av
      Hv1nil = .False.
      If( abs(Hv1) .Lt. Small )Then
        Hv1 = Zero
        Hv1nil = .True.
      End If

      Hv2 = Wflux_inBot * Cseep / dz + Wflux_inTop * Ctop / dz +        &
     &      Wflux_inLat * Clat / dz + ProducPot

!                                   test for type of analytical solution

      If(.Not.Hvnil.And. .Not.Hv1nil)Then
         If(abs(Hv-Hv1).Gt.Small)Then
           Iflsol = 1
         Else
           Iflsol = 2
         End If
      Else If(Hvnil .And. .Not. Hv1nil)Then
         Iflsol = 3
      Else If(.Not.Hvnil .And. Hv1nil)Then
         Iflsol = 4
      Else If(Hvnil .And. Hv1nil)Then
         Iflsol = 5
      End If

      If(Iflsol.Eq.1 .Or. Iflsol.Eq.2)Then

         If( WFrac_t*dz.Le.1.0d-6 .And. WFrac_t0*dz.Gt.1.0d-6)Then
            C_t  = 0.0
            C_av = (Hv2 - Hv * C_t0)/(Hv1 - Hv)
            Return
         Else If(WFrac_t*dz.Gt.1.d-6 .And. WFrac_t0*dz.Le.1.d-6)Then
            C_t  = Hv2/hv1
            C_av = Hv2/hv1
            Return
         Endif

      Else If (Iflsol.Eq.4) Then

         If(WFrac_t*dz.Le.1.0d-6 .And. WFrac_t0*dz.Gt.1.0d-6)Then
            C_t  = 0.0
            C_av = C_t0 - Hv2/hv
            Return
         Else If(WFrac_t*dz.Gt.1.d-6 .And. WFrac_t0*dz.Le.1.d-6)Then
            C_t  = Hv2/hv
            C_av = Hv2/hv
            Return
         Endif

      End If

!     solution when the amount of disappearance equals the amount of 
!     nett inflow (HV2 equals C_t0*HV1):
      If (Iflsol .Le. 3 .And.  abs(Hv1*C_t0-Hv2) .Le. Vvsmall) Then
        C_t  = C_t0
        C_av = C_t0
        Return
      End If
!                                                 determine coefficients
      Call Detcoef(Iflsol,A1,A2,SorpCoef,B1,B2,Hv,Hv1,WFrac_t0,DryBD,dt)

!                                                        compute results
      C_t  = A1 * C_t0 + A2 * Hv2
      C_av = B1 * C_t0 + B2 * Hv2

!     if negative concentrations are calculated reduce REKO and the
!     concentrations 
      If ( C_t .Lt. Zero )Then

!       the long term steady state concentration (HV2/HV1) is greater 
!       than VSMALL  while the amount disappearance exceeds the amount 
!       of nett inflow;
!         or 
!       VSMALL exceeds the long term steady state concentration while
!       the amount of nett inflow exceeds the amount of disappearance:
        If (C_t+Vsmall.Gt.Zero  .Or.  Iflsol.Le.3 .And.                 &
     &      (Hv1*Vsmall-Hv2)/(Hv1*C_t0-Hv2).Lt.Vsmall) Then
            C_t   = C_t0
            C_av  = C_t0
            Return
        End If
 
!       calculate help-variable:
        Hv3 = WFrac_t0 + DryBD*SorpCoef

!       reduce TTRY (analytical solutions):
        If (Iflsol .Le. 2) Then
            Ttry = Hv3/hv*(((Hv1*Vsmall-Hv2)/(Hv1*C_t0-Hv2))**          &
     &                               (-Hv/hv1)-One)
        Else If (Iflsol .Eq. 3) Then
            Ttry = -Hv3/hv1*log((Hv1*Vsmall-Hv2)/(Hv1*C_t0-Hv2))
        Else If (Iflsol .Eq. 4) Then
            Ttry = Hv3/hv*(exp(Hv/hv2*(Vsmall-C_t0))-One)
        Else If (Iflsol .Eq. 5) Then
            Ttry = Hv3/hv2*(Vsmall-C_t0)
        End If     
!                                                 determine coefficients
        C_t  = Vsmall
!ckro - 20160504 solution to prevent Ttry=0 and then divide by zero
       if(Ttry.lt.Vsmall) then
            C_av1  = Vsmall
            C_av  = Vsmall
            ProducAct =  Vsmall
        else
            Call Detcoef(Iflsol,A1,A2,SorpCoef,B1,B2,Hv,Hv1,WFrac_t0,   &
     &               DryBD,Ttry)
            C_av1  = B1*C_t0 + B2*Hv2
            C_av  = (Ttry * (WFrac_t0+Half*Ttry*Hv) * C_av1+(dt-Ttry)*  &
     &          (WFrac_t-Half*(dt-Ttry)*Hv)*Vsmall)/(dt*WFrac_av)
            ProducAct =  ProducPot * Ttry / dt
        End If
      End If

      Return

      End

      Subroutine Detcoef                                                &
     &(Iflsol,A1,A2,Avsocf,B1,B2,Hv,Hv1,Mto,Rhbd,T)
!.......................................................................
!  Subprogram:
!     Subroutine Detcoef - 
!
!  Description:
!     Determines the coefficients A1,A2,B1,B2 of the 
!     analytical solution of the conservation equation. 
!     Description of common variables:
!     Input:
!       IFLSOL        indicates the type of equation.
!       T             length of time interval (d)
!       HV            change of moisture content (d-1)
!       HV1           coefficient A within conservation equation (d-1)
!       MTO           moisture content at start of timestep (-)
!       RHBD          dry bulkdensity (kg.m-3)
!       AVSOCF        average value of adsorption coefficient (kg-1.m3)
!     Output:
!       A1,A2,B1,B2   coefficients (-) in the (average) concentration 
!                     relation according to:
!                            c(t) = A1 * c(t0) + A2 * HV2
!                            cav  = B1 * c(t0) + B2 * HV2   
!  History:
!    2003    Groenendijk, Renaud, Hendriks
!            Release of ANIMO4.0 
!.......................................................................
!---- Declarations
      Implicit None
!                                                       Formal Variables
      Integer :: Iflsol
      real(8) :: A1,A2,B1,B2
      real(8) :: Avsocf,Hv,Hv1,Mto,Rhbd,T
!                                                        Local Variables
      real(8) :: One

      Data    One /1.0d+00/

      select case (Iflsol)
      case (1)
         A1 = ((Mto+Rhbd*Avsocf+Hv*T)/(Mto+Rhbd*Avsocf))**(-Hv1/hv)
         A2 = (One-A1)/hv1
         B1 = ((Mto+Rhbd*Avsocf+Hv*T)/                                  &
     &        (Mto+Rhbd*Avsocf))**((Hv-Hv1)/hv)-One
         B1 = B1*(Mto+Rhbd*Avsocf)/(T* (Hv-Hv1))
         B2 = (One-B1)/hv1
      Return

      case (2)
         A1 = ((Mto+Rhbd*Avsocf+Hv*T)/(Mto+Rhbd*Avsocf))**(-Hv1/hv)
         A2 = (One-A1)/hv1
         B1 = (Mto+Rhbd*Avsocf)/hv/t*log((Mto+Rhbd*Avsocf+Hv*T)/        &
     &        (Mto+Rhbd*Avsocf))
         B2 = (One-B1)/hv1
      Return

      case (3)
         A1 =  exp(-Hv1/(Mto+Rhbd*Avsocf)*T)
         A2 = (One-A1)/hv1
         B1 = (One-A1)/(Hv1/(Mto+Rhbd*Avsocf)*T)
         B2 = (One-B1)/hv1
      Return

      case (4)
         A1 = One
         A2 = One/hv*log((Mto+Rhbd*Avsocf+Hv*T)/(Mto+Rhbd*Avsocf))
         B1 = One
         B2 = A2*(Mto+Rhbd*Avsocf+Hv*T)/(Hv*T) - One/hv
      Return

      case (5)
         A1 = One
         A2 = T/(Mto+Rhbd*Avsocf)
         B1 = One
         B2 = 0.5d+00*T/(Mto+Rhbd*Avsocf)

      case default
         call fatalerr ('DetCoef', 'Illegal value for TASK')
      end select

      Return
      End 