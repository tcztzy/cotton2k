! File VersionID:
!   $Id: rootextraction.f90 303 2016-09-07 11:55:11Z kroes006 $
! ----------------------------------------------------------------------
      subroutine RootExtraction 
! ----------------------------------------------------------------------
!     update    : August 2016: microscopic uptake according to JongvanLier(2013)
!     update    : August 2012: O2-stress according to Bartholomeus(2008)
!     update    : February 2011: macrosopic uptake extended with 
!                                compensation according to Jarvis (1989)
!     date      : August 2004
!     purpose   : Calculate the root water extraction rate as function of soil 
!                 water pressure head and salinity concentration for each node
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local variables
      integer node
      real(8) top,bot,afgen,hlim3,hlim2,qpotrot,qred
      real(8) alpdry,alpwet,alpsol,alpfrs,alptot,vsmall
      real(8) qredwet(macp),qreddry(macp),qredsol(macp),qredfrs(macp)
      real(8) factor
      logical flnus
      
      parameter (vsmall = 1.0d-7)

! ----------------------------------------------------------------------
! --- reset root water extraction array
      do 10 node = 1,numnod
        qrot(node) = 0.0d0
10    continue
      Tactual = qrosum
      qrosum = 0.0d0
      qreddrysum = 0.0d0
      qredwetsum = 0.0d0
      qredsolsum = 0.0d0
      qredfrssum = 0.0d0

! --- skip routine if the crop has not emerged
!      if (.not.flcropEmergence) return

! --- skip routine if there are no roots
      if (rd .lt. vsmall) return

! --- determine lowest compartment containing roots
      node = 1
      do while ((z(node) - dz(node)*0.5d0) .gt. (-rd + 1.d-8))
        node = node + 1
      end do
      noddrz = node

! --- skip routine if transpiration rate is zero                        
      if (ptra .lt. 1.d-10) then                                        
        return                                                          
      endif                                                             
                                                                        
! --- DROUGHT REDUCTION ACCORDING TO FEDDES ET AL. (1978)
      if (swdrought .eq. 1) then                                   

! --- calculate potential root extraction of the compartments
      do node = 1,noddrz
        top = abs((z(node) + 0.5*dz(node)) / rd)
        bot = abs((z(node) - 0.5*dz(node)) / rd)
        qrot(node)=(afgen(cumdens,202,bot)-afgen(cumdens,202,top))* ptra
      enddo

! --- calculating critical point hlim3 according to feddes
        hlim3=hlim3h
        if (atmdem.ge.adcrl) then
          if (atmdem.le.adcrh) hlim3=hlim3h+((adcrh-atmdem)/(adcrh-adcrl))&
     &    *(hlim3l-hlim3h)
        else
          hlim3 = hlim3l
        endif
      endif

! --- DROUGHT REDUCTION ACCORDING TO DE JONG VAN lIER ET AL. (2012)
      if (swdrought .eq. 2) then
        call JongvanLier
      endif

! === COMBINATION OF OXYGEN, DROUGHT, SALT AND FROST STRESS ====

      qrosum = 0.0d0

      do 200 node = 1,noddrz
        alpdry = 1.0d0
        alpwet = 1.0d0
        alpsol = 1.0d0
        alpfrs = 1.0d0

! ---   reduction due to oxygen stress
        if (swoxygen .eq. 1) then
! ---     Feddes linear reduction based on pressure head
          if (node.gt.botcom(1)) then
            hlim2=hlim2l
          else
            hlim2 = hlim2u
          endif
          if (h(node).le.hlim1.and.h(node).gt.hlim2) then
            alpwet = (hlim1-h(node))/(hlim1-hlim2)
          endif
          if (h(node).gt.hlim1) then
            alpwet = 0.0d0
          endif
        else
! ---     Bartholomeus non-linear reduction based on gas filled porosity
          if (swoxygentype .eq. 1) then
! ---       use physical processes
            call OxygenStress(node,alpwet,ResultsOxygenStress) 
          else
! ---       use reproduction functions
            call OxygenReproFunction (OxygenSlope,OxygenIntercept,      &
     &           theta,thetas,tsoil,node,z,dz,alpwet)
          endif
        endif
! ---   reduction due to drought stress
        if (swdrought .eq. 1) then
! ---     Feddes linear reduction based on pressure head
          if (h(node).le.hlim3.and.h(node).ge.hlim4) then
            alpdry = (hlim4-h(node))/(hlim4-hlim3)
          endif
! ---     no transpiration 
          if (h(node).lt.hlim4) then
            alpdry = 0.0d0
          endif
        endif
        if (swdrought .eq. 2) then
! ---     JongvanLier microscopic concept for drought
          alpdry = alpJvLier
        endif

! ---   reduction due to salt stress
        if (swsalinity .eq. 1) then
! ---     reduction according to Maas and Hoffman linear reduction function
          if (cml(node) .gt. saltmax) then
            alpsol = 1.0d0 - (cml(node) - saltmax) * saltslope
            alpsol = max(0.0d0,alpsol)
          endif
        endif
! ---   mind: in case of salt stress with osmotic head, microscopic root water extraction
! ---         according to JongvanLier (2013) should be used (swsalinity = 2); in that case
! ---         salinity stress is included in drought stress and not separately specified
! ---         in output file *.STR

! ----  reduction due to frost conditions
        if (swfrost.eq.1.and.tsoil(node).lt.0.0d0) then
          alpfrs = 0.0d0
        endif

! ----  overall reduction
        qpotrot = qrot(node)
        qrot(node) = qrot(node)*alpwet*alpdry*alpsol*alpfrs
        qrosum = qrot(node) + qrosum

! ----  apportionment to different types stresses (cm)
        qred = qpotrot - qrot(node)
        if(qred.lt.vsmall)then
          qredwet(node) = 0.d0
          qreddry(node) = 0.d0
          qredsol(node) = 0.d0
          qredfrs(node) = 0.d0
        else
          alptot = 0.0d0
!       set flag for non unique stress
          flnus = .true.
          if(alpwet .gt. vsmall) then 
             alptot = alptot + dlog(alpwet)
          else
            qredwet(node) = qred
            qreddry(node) = 0.d0
            qredsol(node) = 0.d0
            qredfrs(node) = 0.d0
            flnus = .false.
          end if
          if(flnus)then
            if(alpdry .gt. vsmall) then 
               alptot = alptot + dlog(alpdry)
            else
              qredwet(node) = 0.d0
              qreddry(node) = qred
              qredsol(node) = 0.d0
              qredfrs(node) = 0.d0
              flnus = .false.
            end if
          end if
          if(flnus)then
            if(alpsol .gt. vsmall) then 
               alptot = alptot + dlog(alpsol)
            else
              qredwet(node) = 0.d0
              qreddry(node) = 0.d0
              qredsol(node) = qred
              qredfrs(node) = 0.d0
              flnus = .false.
            end if
          end if
          if(flnus)then
            if(alpfrs .gt. vsmall) then 
               alptot = alptot + dlog(alpfrs)
            else
              qredwet(node) = 0.d0
              qreddry(node) = 0.d0
              qredsol(node) = 0.d0
              qredfrs(node) = qred
              flnus = .false.
            end if
          end if

          if(flnus)then
            qredwet(node) = dlog(alpwet) / alptot * qred
            qreddry(node) = dlog(alpdry) / alptot * qred
            qredsol(node) = dlog(alpsol) / alptot * qred
            qredfrs(node) = dlog(alpfrs) / alptot * qred
          end if
        end if

        qredwetsum = qredwetsum + qredwet(node)
        qreddrysum = qreddrysum + qreddry(node)
        qredsolsum = qredsolsum + qredsol(node)
        qredfrssum = qredfrssum + qredfrs(node)

200   continue

      if (swdrought .eq. 1) then
! --    Compensated root water uptake according to Jarvis (1989)
        if (abs(alphacrit-1.0d0).ge.vsmall .and.                        &
     &                           abs(qrosum-ptra).ge.vsmall) then
          alptot = qrosum / ptra
          if (alptot .gt. alphacrit) then
            do node = 1,noddrz
              qrot(node) = qrot(node) / alptot            
            enddo
            qredwetsum = 0.0d0
            qreddrysum = 0.0d0
            qredsolsum = 0.0d0
            qredfrssum = 0.0d0
            qrosum = ptra
          else
            do node = 1,noddrz
              qrot(node) = qrot(node) / alphacrit
            enddo
            factor = (ptra - qrosum / alphacrit) / (ptra - qrosum)
            qredwetsum = qredwetsum * factor
            qreddrysum = qreddrysum * factor
            qredsolsum = qredsolsum * factor
            qredfrssum = qredfrssum * factor
            qrosum = qrosum / alphacrit
          endif
        endif
      endif

      return

      end

! ----------------------------------------------------------------------
      subroutine JongvanLier
! ----------------------------------------------------------------------
!     date      : August 2016
!     purpose   : Calculate the root water extraction rate according to 
!                 De Jong van Lier et al. (2013)
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local variables
      integer node,counter
      real(8) afgen,Fy3
      real(8) y1,y2,y3,Fy1,Fy2,hwet,ratio
      real(8) hleafm1,hrootm1(macp)
      real(8) reldepth,rdensity,phi,meandepth
      real(8) maxstep,dHPLant,dHPlantMax
      logical flconverg,flstress
      character(len=200) messag
  
! not used:      include 'params.fi'

! --- initialization
      phi = 3.1415926d0
      counter = 0
      flstress = .false.
      hleafm1 = hleaf
      flconverg = .false.

! --- reset values below root zone to zero
      do node = noddrz+1,numnod
        mflux(node) = 0.d0
        mroot(node) = 0.d0
        hroot(node) = 0.d0
        rootrho(node) = 0.d0
        rootphi(node) = 0.d0
      enddo

! --- give hroot a value when root zone becomes larger and store previous values
      do node = 1, noddrz
        if (abs(hroot(node)) .lt. 1.0d-12) hroot(node) = h(node)
        hrootm1(node) = hroot(node)
      enddo
      
! --- initialization of rootrho and rootphi
      do node = 1,noddrz-1
        reldepth = -z(node)/rd
        rdensity = afgen(rdctb,22,reldepth)
        rmax(node) = 1.d0/sqrt(phi*rdensity)
        rootrho(node) = 4.d0/(rootradius*rootradius-rootcoefa*rmax(node)&
     &      *rootcoefa*rmax(node) + 2.d0 * (rmax(node)*rmax(node) +     &
     &      rootradius*rootradius)*log(rootcoefa*rmax(node)/rootradius))
        rootphi(node) = rootrho(node) * rmax(node)*rmax(node) *         &
     &     log(rootradius/rxylem) / 2.d0 / kroot
      enddo
! --- last node, partly filled with roots
      node = noddrz
      meandepth = (z(node)+0.5*dz(node)-rd)/2.d0
      reldepth = meandepth/(-rd)
      rdensity = afgen(rdctb,22,reldepth)
      rmax(node) = 1.d0/sqrt(phi*rdensity)
      rootrho(node) = 4.d0/(rootradius*rootradius-rootcoefa*rmax(node)  &
     &      *rootcoefa*rmax(node) + 2.d0 * (rmax(node)*rmax(node) +     &
     &      rootradius*rootradius)*log(rootcoefa*rmax(node)/rootradius))
      rootphi(node) = rootrho(node) * rmax(node)*rmax(node) *           &
     &     log(rootradius/rxylem) / 2.d0 / kroot

! --- calculate current matric flux potential in soil water
      do node = 1,noddrz
        h(node) = min(1.0d3,max(h(node),1.0d-8))
        call MatricFlux(2,h(node),node,mflux(node))
      enddo

! --- interpolate matric flux potential in lowest compartment that is partly filled with roots
      node = noddrz
      if (node.gt.1) then
         mflux(node) = mflux(node) + (mflux(node-1)-mflux(node))/       &
     &              (z(node-1)-z(node)) * (meandepth-z(node))
      endif

! --- determine highest pressure head in root zone
      hwet = h(1)
      do node = 2,noddrz
        if (h(node) .gt. hwet) hwet = h(node)
      enddo

! --- determine whether qrosum < ptra (flstress = true)
! --- calculate root water extraction qrosum at Hleaf = wiltpoint and Tactual = Ptra
      Hleaf = wiltpoint
      dHPlant = ptra/Kstem
      if (dHPlant .gt. (hwet - wiltpoint)) flstress = .true.

      HXylem = Hleaf + dHPLant
      call JongvanLierLoop
      if (qrosum .lt. ptra) flstress = .true.

! --- set hroot to previous values
      do node = 1,noddrz
        hroot(node) = hrootm1(node)
      enddo

! --- FLSTRESS = FALSE: DETERMINE HLEAF WITH TACTUAL = PTRA
      if (.not. flstress) then

! --- calculate root extraction at previous value of hleaf
      Hleaf = Hleafm1
      HXylem = Hleaf + dHPLant
      call JongvanLierLoop     

! --- check convergence
      if (abs(qrosum-ptra) .lt. taccur) flconverg = .true.

! --- set initial values y3 and Fy3
      y3 = hleaf
      Fy3 = ptra - qrosum

! --- start search on correct value of pressure head in leaves
      do while (.not. flconverg)
      
! --- reset y1 and Fy1
      y1 = y3
      Fy1 = Fy3

      if (qrosum .lt. ptra) then 
! ---   root water extraction too small, decrease Hleaf
          y2 = y1 - StepHr*log10(max(abs(y1),1.d0))
      else
! ---   root water extraction too large, increase Hleaf
          y2 = y1 + StepHr*log10(max(abs(y1),1.d0))
      endif

! --- calculate root extraction at modified pressure head y2 in leaf
      hleaf = y2
      HXylem = Hleaf + dHPLant
      call JongvanLierLoop     
      Fy2 = ptra - qrosum
      
! --- determine new value of leaf pressure head with Newton Raphson algorithm
      if (abs(Fy1-Fy2).lt. 1.d-10) then
! ---   take maximum step
        if (qrosum .gt. ptra) then
          maxstep = -0.05d0 * y1
          maxstep = max(100.d0,maxstep)
          y3 = y1 + maxstep
        else
          maxstep = 0.05d0 * (y1 - wiltpoint)
          maxstep = max(100.d0,maxstep)
          y3 = y1 - maxstep
        endif
      else
! ---   take Newton Raphson step
        y3 = y1 - (y2-y1)*Fy1 / (Fy2-Fy1)
        if ((y3-y1) .gt. 0.d0) then
          maxstep = -0.05d0 * y1
          maxstep = max(100.d0,maxstep)
          y3 = min((y1 + maxstep),y3)
        else
          maxstep = 0.05d0 * (y1 - wiltpoint)
          maxstep = max(100.d0,maxstep)
          y3 = max((y1 - maxstep),y3)
        endif
      endif
      y3 = min(y3,hwet)
      y3 = max(y3,wiltpoint)
        
! --- new estimated value of leaf pressure head is equal to y3!
 300  Hleaf = y3
      Hxylem = Hleaf + dHPlant

! --- calculate root extraction at y3
      call JongvanLierLoop     
      Fy3 = ptra - qrosum
    
! --- check convergence
      if (abs(qrosum-ptra) .lt. taccur .or.                             &
     &       abs(y3-y1) .lt. 1.0d0) flconverg = .true.

! --- fatal error if too many iterations
      counter = counter + 1
      if (counter .gt. 1000) then
         messag = '4 Too many iterations for microscopic root'          &
     &             //' water uptake. Please adapt input!'
         call warn ('rootextraction',messag,logf,swscre)
         call fatalerr ('rootextraction',messag)
      endif

! --- apply linear interpolation when Fy1 and Fy3 have opposite sign
      if (.not. flconverg) then
        if (((Fy1 .gt. 0.d0 .and. Fy3 .lt. 0.d0) .or.                   & 
     &                     (Fy1 .lt. 0.d0 .and. Fy3 .gt. 0.d0)) .and.   &
     &                      abs(y1-y3) .gt. 1.d0) then

          y3 = y1 + (y3 - y1) * Fy1 / (Fy1 - Fy3)
          y3 = min(y3,hwet)
          y3 = max(y3,wiltpoint)
          goto 300
        endif
      endif

! --- WHILE-DO LOOP TO DETERMINE HLEAF WITHOUT DROUGHT STRESS
      enddo

      qrosum = ptra
      Hleaf = y3
      dHplant = qrosum/kstem
      Hxylem = Hleaf + dHPlant

      else
! --- FLSTRESS = TRUE: DETERMINE QROSUM WITH HLEAF = WILTPOINT
      counter = 0

! --- calculate root extraction at previous value of qrosum
      Hleaf = wiltpoint
      dHplant = tactual/kstem
      if (dHplant .gt. (hwet - wiltpoint)) then
        dHplant = 0.98 * (hwet - wiltpoint)
        tactual = dHplant * kstem
      endif
      HXylem = Hleaf + dHPLant
      call JongvanLierLoop     

! --- check convergence
      if (abs(tactual - qrosum) .lt. taccur) flconverg = .true.

! --- set initial values y3 and Fy3
      y3 = tactual
      Fy3 = tactual - qrosum

! --- start search on correct value of qrosum
      do while (.not. flconverg)
      
! --- reset y1 and Fy1
      y1 = y3
      Fy1 = Fy3

      if (tactual .gt. qrosum) then 
! ---   root water extraction less than adopted, decrease tactual
          y2 = y1 - 0.001d0
      else
! ---   root water extraction larger than adopted, increase tactual
          y2 = y1 + 0.001d0
      endif

! --- calculate root extraction at modified pressure head y2 in leaf
      tactual = y2
      dHplant = tactual/kstem
      HXylem = Hleaf + dHPLant
      call JongvanLierLoop     
      Fy2 = tactual - qrosum
      
! --- determine new value of tactual with Newton Raphson algorithm
        y3 = y1 - (y2-y1)*Fy1 / (Fy2-Fy1)
        maxstep = 0.05d0
        if ((y3-y1) .gt. 0.d0) then
          y3 = min((y1 + maxstep),y3)
        else
          y3 = max((y1 - maxstep),y3)
        endif
      dHPlantMax = hwet - wiltpoint
      y3 = min(y3,(dHPlantMax * Kstem))
      y3 = max(y3,0.d0)
        
! --- new estimated value of tactual is equal to y3!
 400  tactual = y3
      dHPlant = tactual/kstem
      Hxylem = Hleaf + dHPlant

! --- calculate root extraction at y3
      call JongvanLierLoop     
      Fy3 = tactual - qrosum
    
! --- check convergence
      if (abs(tactual - qrosum) .lt. taccur) flconverg = .true.

! --- fatal error if too many iterations
      counter = counter + 1
      if (counter .gt. 1000) then
         messag = '5 Too many iterations for microscopic root'          &
     &             //' water uptake. Please adapt input!'
         call warn ('rootextraction',messag,logf,swscre)
         call fatalerr ('rootextraction',messag)
      endif

! --- apply linear interpolation when Fy1 and Fy3 have opposite sign
      if (.not. flconverg) then
        if ((Fy1 .gt. 0.d0 .and. Fy3 .lt. 0.d0) .or.                    &
     &                     (Fy1 .lt. 0.d0 .and. Fy3 .gt. 0.d0)) then

          y3 = y1 + (y3 - y1) * Fy1 / (Fy1 - Fy3)
          goto 400
        endif
      endif

! --- WHILE-DO LOOP TO DETERMINE HLEAF WITH DROUGHT STRESS
      enddo

      dHplant = qrosum/kstem
      Hxylem = Hleaf + dHPlant

      endif

! --- qrot(node) has been calculated based on Jong van Lier (2013)

      qrosum = 0.0d0
      do node = 1,noddrz
        qrosum = qrot(node) + qrosum
      enddo

      if ( (ptra - qrosum) .lt. taccur) then
! ---   compensate convergence error
        ratio = ptra / qrosum
        do node = 1,noddrz
          qrot(node) = qrot(node) * ratio
        enddo
        qrosum = ptra
        alpJvLier = 1.d0
      else 
! ---   drought reduction factor
        alpJvLier = qrosum / ptra
! ---   calculate potential qrot
        do node = 1,noddrz
          qrot(node) = qrot(node)/alpJvLier
        enddo
      endif

      return

      end

! ----------------------------------------------------------------------
      subroutine JongvanLierLoop 
! ----------------------------------------------------------------------
!     date      : August 2016
!     purpose   : Calculate microscopic root water uptake using hleaf
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local variables
      integer counter,node,lay
      real(8) x1,x2,x3,Fx1,Fx2,qmax,conducsoil,conducroot
      real(8) step,mflux1,mflux2,depth
      character(len=200) messag

! --  initialisatie
      qrosum = 0.d0
      x1 = 999.d0
      counter = 0
      ConducRoot = KRoot / rootradius / log(rootradius/rxylem)

      do node = 1,noddrz
! ---   determine h and matricflux potential at root-soil interface 
        if (h(node) .gt. -1.d0) then
! ---      very wet conditions
           lay = layer(node)
           ConducSoil = ksatfit(lay) / (rootcoefa*rmax(node)) /        &
     &                  log(rootcoefa*rmax(node)/rootradius)
           hroot(node) = (ConducSoil*h(node) + ConducRoot*hxylem) /    &
     &                   (ConducSoil + ConducRoot)
        else
! ---      common conditions
          x3 = hroot(node)
          do while (abs(x3-x1) .gt. CriterHr*log10(max(abs(x3),1.d0)))
            Step = min(abs(x3-x1),StepHr*log10(max(abs(x3),1.d0)))
            x1 = x3
            x2 = x1 - Step
            call MatricFlux(2,x1,node,mflux1)
            call MatricFlux(2,x2,node,mflux2)
            Fx1 = Hxylem -x1 + rootphi(node)*(mflux(node)-mflux1)
            Fx2 = Hxylem -x2 + rootphi(node)*(mflux(node)-mflux2)
            if (abs(Fx2-Fx1).lt.1.d-10) then
              x3 = x2
            else
              x3 = x1 - (x2-x1)*Fx1 / (Fx2-Fx1)
            endif
            x3 = max(x3,min(Hxylem,h(node)))
            x3 = min(x3,max(Hxylem,h(node)))
! ---       fatal error if too many iterations
            counter = counter + 1
            if (counter .gt. 50) then
              messag = '1 Too many iterations for microscopic root'    &
     &               //' water uptake. Please adapt input!'
              call warn ('rootextraction',messag,logf,swscre)
              call fatalerr ('rootextraction',messag)
            endif
          enddo
          x1 = 999.d0
          counter = 0
          hroot(node) = x3
        endif
        call MatricFlux(2,hroot(node),node,MRoot(node))
! ---   calculate root water extraction flux
        if (node .lt. noddrz) then
          qrot(node) = rooteff * rootrho(node) *                        &
     &                 (mflux(node)-mroot(node)) * dz(node)
          if (mflux(node) .gt. mroot(node) ) then
! ---       water extraction, set maximum flux to 10% of available soil water
            qmax = (theta(node) - twilt(node)) * dz(node) * 0.1d0 / dt
            qmax = min(qmax,ptra)
            qrot(node) = min(qrot(node),qmax)
          else
! ---       possible hydraulic lift, set maximum flux to 0.1% change water content
            if (flhydrlift) then
              qmax = -0.001d0 * dz(node) / dt
              qrot(node) = max(qrot(node),qmax)
            else
! ---         no hydraulic lift allowed
              qrot(node) = 0.d0
              hroot(node) = h(node)
              Mroot(node) = Mflux(node)
            endif
          endif
        else
! ---     last node, partly filled with roots
          depth = z(node) + 0.5*dz(node) + rd
          qrot(node) = rooteff * rootrho(node) *                        &
     &                 (mflux(node)-mroot(node)) * depth
          if (mflux(node) .gt. mroot(node) ) then
! ---       water extraction, set maximum flux to 10% of available soil water
            qmax = (theta(node) - twilt(node)) * depth * 0.1d0 / dt
            qmax = min(qmax,ptra)
            qrot(node) = min(qrot(node),qmax)
          else
! ---       possible hydraulic lift, set maximum flux to 0.1% change water content
            if (flhydrlift) then
              qmax = -0.001d0 * depth / dt
              qrot(node) = max(qrot(node),qmax)
            else
! ---         no hydraulic lift allowed
              qrot(node) = 0.d0
              hroot(node) = h(node)
              Mroot(node) = Mflux(node)
            endif
          endif
        endif
        qrosum = qrot(node) + qrosum
      enddo

      return
      end

! ----------------------------------------------------------------------
      subroutine MatricFlux(task,phead,node,outcome) 
! ----------------------------------------------------------------------
!     Date               : February 2010   
!     Purpose            : Initialize and calculate matric flux potential
! ----------------------------------------------------------------------

      use Variables
      implicit none

! --- local variables
      integer task,lay,count,start,node,i
      real(8) phead1,phead2,wcontent,conduc1,conduc2,watcon,hconduc
      real(8) logphead,hosm,hsalt,mfluxsalt,phead,outcome

      select case (task)
      case (1)

! === initialization =========================================================

      do lay = 1,numlay
        do count = 1,801
          mfluxtable(lay,count) = 0.0d0
        enddo
      enddo

      start = int(100.d0*log10(-wiltpoint))
      do lay = 1,numlay
        phead1 = -10.d0**(dble(start)/100.d0)

!       find first Node of the Layer
        i = nod1lay(lay)

        wcontent = watcon(i,phead1,cofgen,swsophy,numtab,sptab,         &
     &                    ientrytab)
        conduc1 = hconduc (i,phead1,wcontent,cofgen,swfrost,10.d0,      &
     &                     swsophy,numtab,sptab,ientrytab)

        do count = start-1,1,-1
          phead2 = -10.d0**(dble(count)/100.d0)
          wcontent = watcon(i,phead2,cofgen,swsophy,numtab,sptab,       &
     &               ientrytab)
          conduc2 = hconduc (i,phead2,wcontent,cofgen,swfrost,10.d0,    &
     &              swsophy,numtab,sptab,ientrytab)
          mfluxtable(lay,count) = mfluxtable(lay,count+1) +             &
     &                 0.5d0 * (conduc1 + conduc2) * (phead2 - phead1) 
          phead1 = phead2
          conduc1 = conduc2
        enddo
      enddo

      return

      case (2)

! === calculation of matric flux potential ===================================

! --- matric flux potential based on soil water pressure head
        lay = layer(node)
        if (phead .lt. wiltpoint) then
! ---     very dry range 
          outcome = 0.0d0
        elseif (phead .gt. -1.023293d0) then
! ---     very wet range (> -10^0.01)
          outcome = mfluxtable(lay,1) + (phead+1.023293d0)*ksatfit(lay)
        else  
! ---     direct access table, with linear interpolation
          logphead = 100.d0*log10(-phead)
          count = int(logphead)
          outcome = (logphead-dble(count))*mfluxtable(lay,count+1) +    &
     &                (dble(count+1)-logphead)*mfluxtable(lay,count)
        endif

! --- correction matric flux potential for osmotic head due to salinity
      if (swsalinity .eq. 2) then
          lay = layer(node)
!         osmotic head in cm
          hosm = salthead * cml(node)
          hsalt = wiltpoint + hosm
          if (hosm .lt. 1.d-3) then
! ---       very dry range 
            mfluxsalt = 0.0d0
          elseif (hsalt .gt. -1.023293d0) then
! ---       very wet range (> -10^0.01)
            mfluxsalt = mfluxtable(lay,1)
          else  
! ---       direct access table, with linear interpolation
            logphead = 100.d0*log10(-hsalt)
            count = int(logphead)
            mfluxsalt = (logphead-dble(count))*mfluxtable(lay,count+1) +&
     &                (dble(count+1)-logphead)*mfluxtable(lay,count)
          endif
          outcome = outcome - mfluxsalt
          outcome = max(0.d0,outcome)
      endif

      case default
         call fatalerr ('MatricFlux', 'Illegal value for TASK')
      end select

      return
      end 
