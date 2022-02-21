! File VersionID:
!   $Id: snow.f90 312 2016-12-22 21:18:18Z kroes006 $
! ----------------------------------------------------------------------
      subroutine snow(task)
! ----------------------------------------------------------------------
!     date               : December 2004
!     purpose            : Simulation snow accumulation and melt
! ----------------------------------------------------------------------
      use Variables
      implicit none

      integer   task        
      real(8)   smelt      ! snowmelt by temperature [cm swe]
      real(8)   smeltr     ! snowmelt by rain [cm swe]
      real(8)   SnDefit
      real(8)   SnLoss
      
! --- constants
      real(8)   cwat       ! specific heat of water [j/kg/k]
      real(8)   lm         ! latent heat of melting [j/kg] 
      real(8)   ts         ! snow temperature [0 oc] 
      real(8)   slw_max    ! max storage of liquid water snow [cm/d]
      real(8)   qlw        ! storage of drained flux from snow pack [cm/d]

 
      parameter(cwat=4180.0d0,lm=333580d0,ts=0.0d0)

! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initialization ===================================================

      if (swinco .eq. 3) then
        snowinco = ssnow
      else
        ssnow = snowinco
      endif

      return

      case (2)

! === snow pack rate and state variables ===============================

! --- reset intermediate snow states
      if (flzerointr) then
        igsnow = 0.0d0
        isubl = 0.0d0
        isnrai = 0.0d0
        ISsnowBeg = Ssnow
      endif

! --- reset cumulative snow states
      if (flzerocumu) then
        cgsnow = 0.0d0
        csubl = 0.0d0
        csnrai = 0.0d0
        cmelt = 0.0d0
        snowinco = ssnow
      endif

! --- when there is snowpack calculate the amount of sublimation
      subl = 0.0d0
      if (swsublim.eq.0) then
         if (ssnow.gt.0.0d0) then
           subl = peva
           if(swetsine.eq.1) subl = pevaday
           empreva = 0.0d0
           peva = 0.0d0
         endif
      endif

! --- when the soil surface is above the freezing point there will be
! --- no accumulation of fresh snow. 
      if (tsoil(1).gt.0.5d0.and.ssnow.lt.1.0d-6.and.gsnow.gt.0.0d0) then
        ssnow = 0.0d0
        melt = gsnow
        subl = 0.d0
      else   
        
! ---   amount of snowmelt [cm swe] negative values of smelt: see 'melt = ' 
        smelt = snowcoef * (tav-ts) 
        
! ---   extra snowmelt when there falls rain on the snowpack [cm swe]
        if (snrai.gt.0.0d0) then
          smeltr = snrai * cwat * (tav-ts) / lm  
        else
          smeltr = 0.0d0
        endif

! ---   total snowmelt [cm swe]; negative values of smelt can partly compensate smeltr
        melt = max(0.0d0,(smelt + smeltr))
      
! ---   amount of snow left [cm swe] without storage of liquid water slw
        ssnow = ssnow + gsnow - subl - melt- slw

! ---   potential amount of liquid water storage
        slw = slw + snrai
! ---   maximum retention of liquid water in snow is fraction 0.07 of total water storage
        slw_max = 0.07 * (slw + ssnow)
! ---   drainage of liquid water from snow
        qlw = max(0.0d0,slw-slw_max)
! ---   remaining storage of liquid water in snow
        slw = slw - qlw
! ---   reset total snow storage and total melt
        ssnow = ssnow + slw
        melt = melt + qlw

! ---   in case of snow deficit: adapt snow loss terms melt and sublimation
        if (ssnow.lt.0.0d0) then
          SnDefit = - Ssnow
          SnLoss  = melt + subl
          melt    = (1.d0 - SnDefit/SnLoss) * melt
          subl    = (1.d0 - SnDefit/SnLoss) * subl
          Ssnow   = 0.d0
          slw     = 0.d0
        endif
      endif

! --- set cumulative amounts
      igsnow = igsnow + gsnow
      isubl = isubl + subl
      isnrai = isnrai + snrai
      cgsnow = cgsnow + gsnow
      csubl = csubl + subl
      cmelt = cmelt + melt
      csnrai = csnrai + snrai

      case default
         call fatalerr ('Snow', 'Illegal value for TASK')
      end select

      return
      end

