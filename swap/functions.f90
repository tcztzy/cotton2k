! File VersionID:
!   $Id: functions.f90 298 2016-07-25 20:07:31Z kroes006 $
! ----------------------------------------------------------------------
      real(8) function hcomean(swkmean,kup,klow,dzup,dzlow)
! ----------------------------------------------------------------------
!     UpDate             : 20090630
!     Purpose            : calc. mean hydraulic conductivity
!     Subroutines called : -
!     Functions called   : -                                           
!     File usage         : -
! ----------------------------------------------------------------------
      implicit none
      integer swkmean
      real(8) kup,klow,dzup,dzlow, a1, a2

!     unweighted arithmic mean
      if(swkmean.eq.1)then
         hcomean = 0.5d0 * (kup + klow)
!     weighted arithmic mean
      else if(swkmean.eq.2)then
         hcomean = (dzup * kup + dzlow * klow) / (dzup + dzlow)
!     unweighted geometric mean
      else if(swkmean.eq.3)then
         hcomean = dsqrt( kup * klow)
!     weighted geometric mean
      else if(swkmean.eq.4)then
         a1 = dzup / (dzup + dzlow)
         a2 = 1.0d0 - a1
         hcomean =  (kup ** a1) * (klow ** a2)
      end if
      return
      end
!-----------------------------------------------------------------------
      real(8) function watcon (node,head,cofgen,swsophy,numtab,sptab,    &
     &                        ientrytab)
! ----------------------------------------------------------------------
!     UpDate             : 20090630
!     Date               : 20060208
!     Purpose            : calc. water content from pressure head HEAD
!     Subroutines called : -
!     Functions called   : -                                           
!     File usage         : -
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'
! --- global 
      integer swsophy,numtab(macp),node
      real(8) sptab(5,macp,matab),head,cofgen(12,macp)
      integer   ientrytab(macp,0:matabentries)    

! --- local
      real(8) h_crit, h_enpr, help, m, n, s_enpr
      real(8) h105, C105, a, b
      real(8) t105
      real(8) alfamg, thetar, thetas, moiscap
! ----------------------------------------------------------------------
!      sptab(1,node,i):    Theta
!      sptab(2,node,i):    h
!      sptab(3,node,i):    dTheta / dh
!      sptab(4,node,i):    K
!      sptab(5,node,i):    dK / dh

!     Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if(swsophy.eq.0)then

         h_crit = -1.0d-2
         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         n = cofgen(6,node)
         m = cofgen(7,node)
         h_enpr = cofgen(9,node)
! 
         if(h_enpr .gt. h_crit)then

            if (head.ge.0.0d0) then
! ---   saturated moisture content
               watcon = thetas
            else if(head .gt. h_crit)then
               help    = (abs(alfamg*h_crit))** n
               help    = (1.0d0 + help) ** m
               help    = thetar + (thetas-thetar)/help    
               watcon = help + (thetas-help)/(-h_crit)*(head-h_crit) 
               watcon = min(watcon,thetas)  
            else 
! ---   first compute |alpha * h| ** n
               help = (abs(alfamg*head)) ** n

! ---   add 1 and raise to the power m
               help = (1.0d0 + help) ** m

! ---   now compute theta
               watcon = thetar+(thetas-thetar)/help 
            end if
         else

            h105 = 1.05d0 * h_enpr
!         h95 = 0.95d0 * h_enpr

            if (head .ge. h105) then
               t105 = thetar+(thetas-thetar) *                          &
     &             ((1.0d0 + (abs(alfamg*h_enpr)) ** n ) ** m) /        &
     &             ((1.0d0 + (abs(alfamg*h105)) ** n ) ** m)
               C105 = (ThetaS - ThetaR) * alfamg * M * N *              &
     &             (Abs(alfamg * h105) ** (N - 1)) *                    &
     &             ((1 + Abs(alfamg * h_enpr) ** N) ** M) /             &
     &             ((1 + Abs(alfamg * h105) ** N) ** (M + 1))
               a =(t105-thetas-C105*h105)/(C105*h105**2)
               b =(t105**2-2*t105*thetas+thetas**2)/                    &
     &                                           (t105-thetas-C105*h105)
               watcon = thetas + b*a*head/(1.0d0+a*head)

            else 
! ---   first compute |alpha * h| ** n
               help = (abs(alfamg*head)) ** n

! ---   add 1 and raise to the power m
               help = (1.0d0 + help) ** m

! --- For modified VanGenuchten model:
!     - S_enpr: relative saturation at Entry Pressure h_enpr 
               s_enpr = (1.0d0 + (abs(alfamg*h_enpr)) ** n ) ** m

! ---   now compute theta
               watcon = thetar+(thetas-thetar)/help * s_enpr
            end if
         end if

!     Use tabulated function. 
      else if(swsophy.eq.1)then
         if(head.ge.-1.0d-9)then
            watcon = sptab(2,node,numtab(node))
         elseif(head.lt.sptab(1,node,1)) then
            watcon = sptab(2,node,1)
         else
            call EvalTabulatedFunction(0,numtab(node),1,2,4,node,sptab, &
     &                                 ientrytab,head,watcon,moiscap)
         end if
      end if

      return
      end

! ----------------------------------------------------------------------
      real(8) function moiscap (node,head,cofgen,dt,swsophy,numtab,     &
     &                          sptab,ientrytab)
! ----------------------------------------------------------------------
!     UpDate             : 20090630
!     Date               : 20060208
!     Purpose            : calculate the differential moisture         
!                          capacity (as a function of pressure head)
!     Subroutines called : -                                           
!     Functions called   : watcon                                     
!     File usage         : -                                           
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'
! --- global 
      integer swsophy,numtab(macp),node
      real(8) head,dt,sptab(5,macp,matab),cofgen(12,macp)
      integer   ientrytab(macp,0:matabentries)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment

! --- local
      real(8) h_crit,h_enpr,m,n,alphah,s_enpr,watcon
      real(8) alfamg,thetar,thetas
      real(8) h105, C105, a, b, term1, term2
      real(8) t105
      real(8) dummy
! ----------------------------------------------------------------------

!     Use analytical expression. 
      if(swsophy.eq.0)then
         h_crit = -1.0d-2
         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         n      = cofgen(6,node)
         m      = cofgen(7,node)
         h_enpr = cofgen(9,node)

         if(h_enpr .gt. h_crit)then

            if (head.ge. 0.0d0) then

               moiscap = dt * 1.0d-7

            else if(head .gt. h_crit)then

               term1 = watcon(node,h_crit,cofgen,swsophy,numtab,sptab,  &
     &                                                      ientrytab)
               moiscap = (thetas-term1) / (-h_crit)
            else

! ---   use analytical evaluation of capacity
               alphah = abs (alfamg*head)

! ---   compute |alpha * h| to the power n-1
               term1 = alphah ** (n -1.0d0)

! ---   compute |alpha*h| to the power n
               term2 = term1 * alphah

! ---   add one and raise to the power m+1
               term2 = (1.0d0 + term2) ** (m + 1.0d0)

! ---   divide theta-s minus theta-r by term2
               term2 = (thetas-thetar) / term2

! ---   calculate the differential moisture capacity
               moiscap = abs(-1.0d0 * n * m * alfamg * term2 * term1)
            end if
         else

            h105 = 1.05d0 * h_enpr
!         h95 = 0.95d0 * h_enpr

            if (head .ge. h105) then
               t105 = thetar+(thetas-thetar) *                          &
     &             ((1.0d0 + (abs(alfamg*h_enpr)) ** n ) ** m) /        &
     &             ((1.0d0 + (abs(alfamg*h105)) ** n ) ** m)
               C105 = (ThetaS - ThetaR) * alfamg * M * N *              &
     &            (Abs(alfamg * h105) ** (N - 1)) *                     &
     &             ((1 + Abs(alfamg * h_enpr) ** N) ** M) /             &
     &             ((1 + Abs(alfamg * h105) ** N) ** (M + 1))
               a =(t105-thetas-C105*h105)/(C105*h105**2)
               b =(t105**2-2*t105*thetas+thetas**2)/                    &
     &            (t105-thetas-C105*h105)
               moiscap = b*a/((1.0d0+a*head)**2)

            else

               alphah = abs (alfamg*head)
               term1 = alphah ** (n -1.0d0)
               term2 = term1 * alphah
               term2 = (1.0d0 + term2) ** (m + 1.0d0)
               term2 = (thetas-thetar) / term2

! --- For modified VanGenuchten model:
!     - S_enpr: relative saturation at Entry Pressure h_enpr 
               s_enpr = (1.0d0 + (abs(alfamg*h_enpr)) ** n ) ** m

               moiscap = abs(-1.0d0 * n * m * alfamg*term2*term1)*s_enpr
  
            end if

         end if   
      
         if (head.gt.-1.0d0.and.moiscap.lt.(dt * 1.0d-7))               &
     &                       moiscap = dt * 1.0d-7

!     Use tabulated function. 
      else if(swsophy.eq.1)then
         if(head.ge.-1.0d-9)then
            moiscap = 1.0d-08
         elseif(head.lt.sptab(1,node,1)) then
            moiscap = 0.0d0
         else
            call EvalTabulatedFunction(0,numtab(node),1,2,4,node,sptab, &
     &                                 ientrytab,head,dummy,moiscap)
         end if
      end if

      return
      end

      function dhconduc (node,head,theta,dimocap,cofgen,swfrost,rfcp,   &
     &                   swsophy,numtab,sptab,ientrytab)
! ----------------------------------------------------------------------
!     Update             : 20090630
!     Date               : 2005
!     Purpose            : calculate derivative of hydraulic conductivity (as a     
!                          function of THETA)
!     Subroutines called : TTutil
!     Functions called   : - 
!     File usage         : -
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'
! --- global 
      integer swfrost,swsophy,numtab(macp),node
      real(8) sptab(5,macp,matab),dhconduc,dimocap,theta,rfcp,head
      real(8) cofgen(12,macp)
      integer   ientrytab(macp,0:matabentries) 

! --- local
      real(8) term0,term1,term2,term3,term4,relsat, dummy
      real(8) m,n,ksatfit,lambda,h_enpr,s_enpr,thetar,thetas,alfamg
      character(len=200) messag
! ----------------------------------------------------------------------

!     Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if(swsophy.eq.0)then
         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         ksatfit = cofgen(3,node)
         lambda = cofgen(5,node)
         n = cofgen(6,node)
         m = cofgen(7,node)
         h_enpr = cofgen(9,node)

! --- For modified VanGenuchten model:
!     - S_enpr: relative saturation at Entry Pressure h_enpr 
         s_enpr = (abs(alfamg*h_enpr)) ** n
         s_enpr = (1.0d0 + s_enpr) ** (-m)

         relsat = (theta-thetar)/(thetas-thetar)

         if(cofgen(10,node).gt.0.0d0.and. relsat.gt.cofgen(11,node))then
            messag = 'Linear interpoltion option for examined ksat not '&
     &      //'yet implemented for implicit hydraulic conductivity '    &
     &      //'(swKimpl=1) in iteration scheme' 
            call fatalerr ('dhconduc',messag)
         end if

         if (relsat.lt.0.001d0) then
            dhconduc = 0.0d0
         else if(relsat.gt.s_enpr) then
            dhconduc = 1.0d-12
         else
            dhconduc = dimocap / (thetas-thetar)
            term0    = (s_enpr*relsat)**(1.0d0/m)
            term1    = 1.0d0 - term0
            term2    = (2.0d0+lambda)*term0 - lambda
            term3    = term1 ** (m-1.0d0)
            term4    = 1.d0 - (1.d0 - s_enpr**(1.d0/m)) ** m
            dhconduc = dhconduc * ksatfit * relsat ** (lambda-1.0d0)
            dhconduc = dhconduc * (1.0d0 - term1 ** m)
            dhconduc = dhconduc * (lambda + term2 * term3)
            dhconduc = dhconduc / (term4**2)
         endif

!     Use tabulated function. "dhconduc" is calclated as a function of "head"
      else if(swsophy.eq.1)then
         if(theta.ge.sptab(2,node,numtab(node))-1.0d-9)then
            dhconduc = 1.0d+08
         else
            call EvalTabulatedFunction(0,numtab(node),1,3,5,node,sptab, &
     &                                 ientrytab,head,dummy,dhconduc)
         end if
      end if

! --- in case of frost conditions
      if (swfrost.eq.1)then
         dhconduc = dhconduc * rfcp
      endif  

      return
      end
     
      real(8) function hconduc (node,head,theta,cofgen,swfrost,rfcp,    &
     &                         swsophy,numtab,sptab,ientrytab)
! ----------------------------------------------------------------------
!     Update             : 20090630
!     Date               : 19990929
!     Purpose            : calculate hydraulic conductivity (as a     
!                          function of THETA)
!     Subroutines called : -
!     Functions called   : -
!     File usage         : -
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'
! --- global 
      integer swfrost,swsophy,numtab(macp),node
      real(8) sptab(5,macp,matab),theta,rfcp,cofgen(12,macp)
      integer   ientrytab(macp,0:matabentries)    

! --- local
      real(8) term1,relsat,hconode_vsmall,m,ksatfit,lambda,dummy
      real(8) relsatm,relsat1,alfamg,thetar,thetas,head
      real(8) h_enpr,n,term2,h_crit,thetam,relsatthr,ksatthr,ksatexm
      logical flksatexm
      data    hconode_vsmall  /1.0d-10/
! ----------------------------------------------------------------------

!     Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if(swsophy.eq.0)then
         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         ksatfit = cofgen(3,node)
         lambda = cofgen(5,node)
         n = cofgen(6,node)
         m = cofgen(7,node)
         h_enpr = cofgen(9,node)

         h_crit = -1.0d-2

         if(cofgen(10,node).gt.0.0d0)then 
            flksatexm = .true.
            ksatexm   = cofgen(10,node)
            relsatthr = cofgen(11,node)
            ksatthr   = cofgen(12,node)
         else
            flksatexm = .false.
         end if


         relsat = (theta-thetar)/(thetas-thetar)
      
         if(flksatexm .and. relsat.gt.relsatthr)then

            term1   = (relsat-relsatthr)/(1.0d0-relsatthr)
            hconduc = term1 * ksatexm + (1.0d0-term1) * ksatthr

         else
            if(h_enpr .gt. h_crit)then

               if (relsat.lt.0.001d0) then
                  hconduc = hconode_vsmall
               else if(relsat.gt.(1.0d0-1.0d-6))then
                  hconduc = ksatfit
               else
                  term1   = ( 1.0d0-relsat**(1.0d0/m) ) ** m
                  hconduc = ksatfit*(relsat**lambda)*                   &
     &                            (1.0d0-term1)*(1.0d0-term1)
               endif

            else
! --- For modified VanGenuchten model:

               thetam = thetar+(thetas-thetar) *                        &
     &               ((1.0d0 + (abs(alfamg*h_enpr)) ** n ) ** m)

               if (relsat.lt.0.001d0) then
                  hconduc = hconode_vsmall
               else 
                  if (theta .ge. thetam) then
                     hconduc = ksatfit
                  else
                     relsatm = (theta-thetar)/(thetam-thetar)
                     relsat1 = (thetas - thetar)/(thetam-thetar)
                     term1   = ( 1.0d0 - (relsatm)** (1.0d0/m) ) ** m
                     term2   = ( 1.0d0 - (relsat1)** (1.0d0/m) ) ** m
                     hconduc = ksatfit*(relsat**lambda) *               &
     &                        ((1.0d0-term1)/(1.0d0-term2)) ** 2
                  end if
               end if
            end if
            hconduc = min(hconduc,ksatfit)
         end if

!     Use tabulated function. "hconduc" is calclated as a function of "head"
      else if(swsophy.eq.1)then
         if(theta.ge.sptab(2,node,numtab(node))-1.0d-9)then
            hconduc = sptab(3,node,numtab(node))
         else if(theta.le.sptab(2,node,1)+1.0d-9)then
            hconduc = sptab(3,node,1)            
         else
            call EvalTabulatedFunction(0,numtab(node),1,3,5,node,sptab, &
     &                                 ientrytab,head,hconduc,dummy)
         end if
      end if

! --- in case of frost conditions
      if (swfrost.eq.1)then
         hconduc = hconduc * rfcp + hconode_vsmall * (1.0d0 - rfcp)
      endif  

      return
      end
     

! ----------------------------------------------------------------------
      real(8) function afgen(table,iltab,x)
! ----------------------------------------------------------------------
!     source             : kees rappoldt, 1/86
!     purpose            : linear interpolation in table with
!                          length iltab for a given value of the
!                          independent variable x.
! ----------------------------------------------------------------------
      implicit none

      integer i,iltab
      real(8) table(iltab),slope,x
! ----------------------------------------------------------------------
      if (table(1).ge.x)  goto 40
      do 10 i = 3,iltab-1,2
        if (table(i).ge.x) goto 30
        if (table(i).lt.table(i-2)) goto 20
 10   continue
! --- table fully filled, argument larger then last x in table
      afgen = table(iltab)
      return
! --- table partly filled, argument larger then last x in table
 20   afgen = table(i-1)
      return
! --- argument between first and last x in table, interpolation
 30   slope = (table(i+1)-table(i-1))/(table(i)-table(i-2))
      afgen = table(i-1) + (x-table(i-2))*slope
      return
! --- argument less or equal to first x in table
 40   afgen = table(2)
      return
      end

! ----------------------------------------------------------------------
      real(8) function prhead (node,disnod,cofgen,wcon,h,               &
     &                        swsophy,numtab,sptab,ientrytab)
! ----------------------------------------------------------------------
!     update             : 20090630
!     date               : 19990912
!     purpose            : calculate pressure head at nodal      
!                          point node from water content
!     description of coefficient for van Genuchten relation
!        cofgen(1,node) = ores
!        cofgen(2,node) = osat
!        cofgen(3,node) = ksatfit
!        cofgen(4,node) = alfa
!        cofgen(5,node) = lexp
!        cofgen(6,node) = npar
!        cofgen(7,node) = 1.d0 - (1.d0 / npar)
!        cofgen(8,node) = dummy
!        cofgen(9,node) = h_enpr
!        cofgen(10,node)= ksatexm
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'
! --- global
      integer node,swsophy,numtab(macp)
      real(8) disnod,wcon,h(macp),cofgen(12,macp),sptab(5,macp,matab)
      integer   ientrytab(macp,0:matabentries)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment

! --- local
      real(8) alfamg,thetar,thetas,h_enpr,s_enpr,npar,mpar
      real(8) help,dummy, prh
! ----------------------------------------------------------------------


      if(swsophy.eq.0)then
         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         npar   = cofgen(6,node)
         mpar   = cofgen(7,node)
         h_enpr = cofgen(9,node)

         if (thetas-wcon.lt.1.0d-6) then

! --- saturated pressure head
            if (node.eq.1) then
               prhead = disnod
            else
               prhead = h(node-1) + disnod
            endif
            prhead = dmax1(prhead,h_enpr)
         else
            if (wcon-thetar.lt.1.0d-6) then
               prhead = -1.0d12
            else

! --- For modified VanGenuchten model:
!     - S_enpr: relative saturation at Entry Pressure h_enpr 
               s_enpr = (abs(alfamg*h_enpr)) ** npar
               s_enpr = (1.0d0 + s_enpr) ** (-mpar)

! ---     first calculate the inverse of the sorptivity
               help = (thetas - thetar) / (wcon - thetar) / s_enpr    
! ---     raise to the power 1/m
               help = help ** (1.0d0 / mpar)
! ---     subtract one and raise to the power 1/n 
               help = (help - 1.0d0) ** (1.0d0 / npar)
! ---     divide by alpha.
               prhead = -1.0d0 * abs(help/alfamg)
            endif
         endif

      else if(swsophy.eq.1)then
         if (sptab(2,node,numtab(node))-wcon.lt.1.0d-6) then

! --- saturated pressure head
            if (node.eq.1) then
              prhead = disnod
            else
              prhead = h(node-1) + disnod
            endif
            prhead = dmax1(prhead,0.0d0)
         else

            call EvalTabulatedFunction(1,numtab(node),1,2,4,node,sptab, &
     &                                 ientrytab,prh,wcon,dummy)
            prhead = prh
         end if
      end if

      return
      end

! ----------------------------------------------------------------------
      integer function stepnr(array,length,x)
! ----------------------------------------------------------------------
! --- find stepnumber 
! ----------------------------------------------------------------------
      implicit none

      integer i,length
      real(8) x,array(length)
! ----------------------------------------------------------------------
      if (array(1).ge.x) goto 30

      do 10 i = 2,length
        if (array(i).gt.x)          goto 20 
        if (array(i).lt.array(i-1)) goto 20
10    continue  
! --- array fully filled, argument larger than last x in array
      stepnr = length
      return
! --- array partly filled, argument larger than last x in array or
! --- argument between first and last x in table
20    stepnr = i-1
      return
! --- argument less or equal to first x in array
30    stepnr = 1
      return
      end

! ----------------------------------------------------------------------
      real(8) function WLEVST (swstor,sttab)
! ----------------------------------------------------------------------
!     Date               : 23/08/99
!     Purpose            : calculate surface water level from    
!                          surface water storage, using table
!     Subroutines called : -
!     Functions called   : -
!     File usage         : -
!     Differences SWAP/SWAPS: None
! ----------------------------------------------------------------------
      implicit none
! --- global
      real(8) swstor, sttab(22,2)
! --- local
      integer i
      real(8) dswst
      character(len=80) messag      

! ----------------------------------------------------------------------
      if (swstor .lt. sttab(22,2)) then
        messag = 'Surface water storage below bottom of table'
        call fatalerr ('Wlevst',messag)
      endif
      if (swstor .gt. sttab(1,2)) then
        messag = 'Surface water storage above top of table'
        call fatalerr ('Wlevst',messag)
      endif
      i = 0
 100  i = i+1
      if (.not.(swstor.ge.sttab(i+1,2).and.                             &
     &          swstor.le.sttab(i,2))) goto 100
      dswst = (swstor - sttab(i+1,2)) / (sttab(i,2) - sttab(i+1,2))
      wlevst = sttab(i+1,1) + dswst * (sttab(i,1) - sttab(i+1,1))

      return
      end

! ----------------------------------------------------------------------
      real(8) function swstLEV (wlev,sttab)
!-----------------------------------------------------------------------
!     Date               : 23/08/99
!     Purpose            : calculate surface water storage from  
!                          surface water level
!     Subroutines called : -
!     Functions called   : -
!     File usage         : -
!     Differences SWAP/SWAPS: None
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER i

      real(8) wlev, sttab(22,2),dwl
      character(len=200) messag
! ----------------------------------------------------------------------

      if (wlev .lt. sttab(22,1)) then
        messag = 'Surface water storage below bottom of table'
        call fatalerr ('swstlev',messag)
      endif
      if (wlev .gt. sttab(1,1)) then
        messag = 'Surface water storage above top of table'
        call fatalerr ('swstlev',messag)
      endif

      i = 0
 100  i = i+1
      if (.not.(wlev.ge.sttab(i+1,1).and.wlev.le.sttab(i,1))) goto 100
      dwl = (wlev - sttab(i+1,1)) / (sttab(i,1) - sttab(i+1,1))
      swstlev = sttab(i+1,2) + dwl * (sttab(i,2) - sttab(i+1,2))

      return
      end

! ----------------------------------------------------------------------
      real(8) function qhtab (wlev,imper,hqhtab,qqhtab)
!-----------------------------------------------------------------------
!     date               : 23/08/99
!     purpose            : calculate surface water discharge from  
!                          surface water level, using table
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

      integer imper,itab

      real(8) wlev,hqhtab(mamp,mamte),qqhtab(mamp,mamte),dwl
! ----------------------------------------------------------------------

      itab = 2
      do while (wlev.lt.hqhtab(imper,itab))
        itab = itab+1
      enddo
      dwl = (wlev-hqhtab(imper,itab))/                                  &
     &       (hqhtab(imper,itab-1)-hqhtab(imper,itab))
      qhtab = qqhtab(imper,itab) + dwl *                                &
     &        (qqhtab(imper,itab-1)-qqhtab(imper,itab))  

      return
      end

      real(8) function runoff(swdra,pond,pondmx,rsro,rsroexp,wls,       &
     &                       swst,sttab,dt)
      implicit none
      integer swdra
      real(8) pond,pondmx,rsro,rsroexp,wls,swst,sttab(22,2),dt,inun_max,&
     &        swstlev

      runoff = 0.0
      if (pond-pondmx.gt.0.0d0 .and. swdra.ne.2) then
         if(rsro.lt.1.0d-3)then 
            runoff = pond-pondmx
         else         
            runoff = dt/rsro* (pond-pondmx)**rsroexp
         end if
      else if (swdra.eq.2) then
         if(pond.gt.pondmx .and. pond.gt.wls)then
            runoff = dt/rsro* (pond-max(pondmx,wls))**rsroexp
         else if(pond.lt.wls)then
            inun_max = swst - swstlev(pond,sttab)
            runoff = -min(inun_max,wls-max(pond,pondmx))
         end if
      end if
      return
      end
! ----------------------------------------------------------------------
      real(8) function insw (x1,x2,x3)
!-----------------------------------------------------------------------
!     date               : 20150317
!     purpose            : switch routine taken from TTUTIL 
! ----------------------------------------------------------------------
      implicit none
      real(8)  x1,x2,x3
      if(x1.lt.0.0d0) then
         insw = x2
      else
         insw = x3
      endif
      return
      end
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
      real(8) function interpol (mn,mx,x)
!-----------------------------------------------------------------------
!     date               : 20160219
!     purpose            : interpolation routine taken from TTUTIL 
! ----------------------------------------------------------------------
      implicit none
!*     FORMAL_PARAMETERS:
      real(8)  mn,mx,x
!**    local parameters
      character(len=52) messag
      save

      if (mx.lt.mn) then
!*        minimum is larger than maximum, should generate an error
         write (messag, '(2(a,g12.5))')                                 &
     &   'argument error, min = ',mn,', max = ',mx
         call fatalerr ('limit',messag)
      end if

      if (x.lt.mn) then
!*        x below allowed range ; return lower bound
         interpol = mn

      else if (x.le.mx) then
!*        x in range ; return x
         interpol = x

      else
!*        x above allowed range ; return upper bound
         interpol = mx

      end if

      return
      end
! ----------------------------------------------------------------------
