! File VersionID:
!   $Id: meteodt.f90 298 2016-07-25 20:07:31Z kroes006 $
!
!     This file contains the following subroutines, in order of calling:
!     1. MeteoDT           : main routine                       ; called in SWAP and ReadMeteo (optional)                          
!     2. ProcessRainEvents : processes input data on rain events; called in MeteoDT (optional)  
!     3. ProcessMeteoTsteps: processes meteo input data per dt  ; called in MeteoDT (optional) 
!     4. ETSine            : distributes potential transpiration &
!                            evaporationaccording to sine wave  ; called in MeteoDT (optional)

! SUBROUTINE 1.
! ----------------------------------------------------------------------
      subroutine MeteoDT
! ----------------------------------------------------------------------
!     Last modified      : February 2014              
!     Purpose            : returns meteorological fluxes of current day
!                        : or of parts of a day (detailed meteo input)
! ----------------------------------------------------------------------
      use variables
      implicit none

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! --- meteo input handling on yearly and daily basis -------------------
! ----------------------------------------------------------------------

! --- beginning of year: process rain events
      if (flYearStart .and. flRainIntens) then
         call ProcessRainEvents(swrain,yearmeteo,dtmin,raintab,tcum,    &
     &               tend,tstart,wet,nmrain,timjan1,rainamount,rainrec, &
     &               arai,rainfluxarray,raintimearray)
!
        flYearStart = .false.
      endif


! ----------------------------------------------------------------------
! --- calculations of meteo variables on time step basis ---------------
! ----------------------------------------------------------------------

! --- update actual rain record and set precipitation fluxes per time step
!     or update actual meteo record and set meteo fluxes per time step
      if (flMeteoDT) then
         call ProcessMeteoTsteps
      endif

! --- distribute potential transpiration and evaporation according to sine wave
      if (flETSine) then
         call ETSine
      endif

      return     
      end subroutine MeteoDT


! SUBROUTINE 2.
! ----------------------------------------------------------------------
      subroutine ProcessRainEvents(swrain,yearmeteo,dtmin,raintab,tcum, &
     &               tend,tstart,wet,nmrain,timjan1,rainamount,rainrec, &
     &               arai,rainfluxarray,raintimearray)
! ----------------------------------------------------------------------
!     Last modified      : February 2014
!     Purpose            : process rain events of one calendar year 
!     Interface:
!       I   - swrain,daynrfirst,daynrlast,yearmeteo,dtmin,raintab,tcum,
!             tend,tstart,wet,nmrain,timjan1,rainamount,rainrec
!       O   - arai,rainfluxarray,raintimearray
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer   nmrain,rainrec,swrain,yearmeteo
      real(8)   arai(366),dtmin,raintab(60),timjan1,wet(366)
      real(8)   tcum                 ! Time cumulative since start simulation (T)
      real(8)   tend                 ! End date of simulation run
      real(8)   tstart               ! End date of simulation run
      real(8)   rainamount(mrain)    ! Array with short duration rainfall (L)
      real(8)   rainfluxarray(mrain) ! Array with short duration rainfall intensities (L/T)
      real(8)   raintimearray(mrain) ! Array with times (T) at which rainfall intensity changes

! --- local
      integer   i,iendyear,j,l,nlack,nn,rday,rdaya(367),rdayold
      real(8)   afgen,araihlp(367),day(mrain),rainam(mrain),rainflux
      real(8)   raintime,ratimar(mrain),tendyear,vsmall,wght,wwet(368)
      vsmall    = 1.0d-8

! ----------------------------------------------------------------------
! --- process rain events on yearly basis ------------------------------
! ----------------------------------------------------------------------

! --- for rain options 1 and 2: convert daily rain quantities and intensities or durations
!     into rain events by creating raintime and rainflux arrays conform rain option 3 
      if (swrain.eq.1 .or. swrain.eq.2) then
         rainrec = 1
         do i = 1, nmrain
           if (raintimearray(i+1).gt.tstart-vsmall) then
             rainrec = rainrec + 1
!   - beginning (00:00) of days of current year within simulation period
             day(rainrec)    = raintimearray(i+1) - timjan1 + 1.d0 
             rainam(rainrec) = 0.1d0 * rainamount(i)  ! convert from mm to cm
             wwet(rainrec)   = wet(i) 
           endif
         enddo

!   - set first record of raintime and rainflux (= 0)
         raintimearray(1) = tcum + dtmin
         rainam(1)        = 0.d0
         rainfluxarray(1) = 0.d0
         
!   - set rest of records of raimtime and rainflux (only when rainam[ount] > 0)     
         nmrain   = rainrec + 1 
         rainrec = 0
         do i = 2, nmrain
           if (rainam(i).gt.vsmall) then
             if (swrain.eq.1) then  
!          + mean rainfall intensities are specified
               rainflux = afgen(raintab,60,day(i))
               raintime = dmin1(0.99d0,rainam(i)/rainflux)

             elseif (swrain.eq.2) then
!          + rainfall durations are specified
               raintime = wwet(i)
             endif
             
             if (i.eq.2) then
               rainrec = rainrec + 1
             else
!          + first raintime of a day: closure of last period of former day with rain = 0
               rainrec = rainrec + 2
               raintimearray(rainrec) = dble(i-2) + tcum
               rainfluxarray(rainrec) = 0.d0
             endif
!          + second raintime of a day: closure of first period of the day, rain = rainam
             raintimearray(rainrec+1)  = dble(i-2) + tcum + raintime
             rainfluxarray(rainrec+1)  = rainam(i) / raintime    
   
           endif
         enddo

! ---    extend array with records at end of current year
         tendyear = 365.d0
         if (mod(yearmeteo,4).eq.0) tendyear = 366.d0
         raintimearray(rainrec+2) = tcum + tendyear + dtmin
         rainfluxarray(rainrec+2) = 0.d0

! --- in case of rain events: 1 calculate daily values
!                             2 fill raintimearray and rainfluxarray 
      elseif (swrain.eq.3) then

! --- total amount of rain per meteo day arai
!   - initialize array with sum of rain
         do i = 1, 366
           araihlp(i)  = 0.d0
         enddo

!   - less rain days than meteo days? Fill gap with dummies
         rdayold = int(raintimearray(1)-timjan1) + 1  ! first day with rain record of the year
         nlack = rdayold - 1
         do j = 1, nlack
           rdaya(j) = j
           araihlp(j)  = 0.d0
         enddo
         rdaya(j) = rdayold
         araihlp(j)  = 0.d0
!   - fill array of daily sums of rain with real values         
         do i = 1, nmrain
           rday = int(raintimearray(i)-timjan1) + 1
           if (rday.gt.rdayold) then
             do l = 1, rday-rdayold-1
               j = j + 1
               rdaya(j) = rdaya(j-1) + 1
               araihlp(j)  = 0.d0
             enddo
!   - in case of rain event exceding current day, calculate weights for assigning parts to current and next day              
             wght =                                                     &
     &        (1.d0-(raintimearray(i-1)-dble(int(raintimearray(i-1)))))/&
     &        (raintimearray(i)-raintimearray(i-1))
             araihlp(j) = araihlp(j) + rainamount(i) * wght
             rdayold  = rday 
             j = j + 1
             rdaya(j) = rday
             araihlp(j) = araihlp(j) + rainamount(i) * (1.d0 - wght)
           else
             if (i.gt.1) then
               araihlp(j) = araihlp(j) + rainamount(i)
             endif
           endif
         enddo
         
!   - rain days missing at the end of the year? Fill gap with dummies     
         iendyear = 365
         if (mod(yearmeteo,4).eq.0) iendyear = 366   
         nlack = iendyear - j
         do i = 1, nlack
           rdaya(j+i) = j + i
           araihlp(j+i)  = 0.d0
         enddo

!   - save help array araihlp into arai
         do i =1,366
            arai(i) = araihlp(i)
!            write(117,*) i, arai(i),araihlp(i+1)
         enddo

! --- assign values to raintimearray and rainam array for calculating rainfluxarray
!   - find time gap without rain events at the beginning of teh year
         rainrec = 1
         ratimar(1) = raintimearray(1) - tstart
         i = 1
         do while (ratimar(i).lt.vsmall)
           i = i + 1
           ratimar(i) = raintimearray(i) - tstart
         enddo
         do while (rainamount(i).lt.vsmall .and.                        &
     &             rainamount(i+1).lt.vsmall)
           i = i + 1
           ratimar(i) = raintimearray(i) - tstart
         enddo

! --- fill arrays with real values of rain events
         nn = i
         do i = nn, nmrain
           ratimar(i+1) = raintimearray(i+1) - tstart
           rainrec = rainrec + 1
           raintimearray(rainrec) = ratimar(i)
           rainam(rainrec) = 0.1d0 * rainamount(i)  ! convert from mm to cm    
         enddo
         nmrain = rainrec

!   - set first record of arrays
         raintimearray(1) = tcum + dtmin
         rainfluxarray(1) = 0.0d0

! --- calculate rainfluxes (cm/d) and fill rainfluxarray 
!     Flx(t1) = P(t1) / (T(t1)-T(t0))  = counts for time interval T(t0) -> T(t1)
!     Flx = flux, P = quantity of rain, T = time
         do i = 2, nmrain
            rainfluxarray(i) = rainam(i) /                              &
     &                         (raintimearray(i)-raintimearray(i-1))
         end do
         
! --- set final values of raintimearray and corresponding rainfluxarray         
         raintimearray(nmrain+1) = dmax1(tend+1.1d0-tstart,             &
     &                                   raintimearray(nmrain)+1.d0)
         rainfluxarray(nmrain+1) = 0.d0
         nmrain = nmrain + 1

      endif

! --- for swrain = 1-3: determine start rain record
      rainrec = 1

      return
      end subroutine ProcessRainEvents  
          

! SUBROUTINE 3.
! ----------------------------------------------------------------------
      subroutine ProcessMeteoTsteps
! ----------------------------------------------------------------------
!     Last modified      : February 2014
!     Purpose            : calculations of meteo variables on time step 
!                          basis (in case of precipitation intensities 
!                          [swrain 1-3] or detailed meteo input) 
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local
!      real(8)   dtCrit

!     critical time-interval 
!      dtCrit = 1.d-6
! ----------------------------------------------------------------------
!
! === Precipitation intensities ========================================

      if (flrainintens) then
! --- per time step: set precipitation fluxes for current time step
         graidt  = fprecnosnow * rainfluxarray(rainrec)
         nraidt  = finterception * graidt
         aintcdt = (1.d0-finterception) * graidt

! --- calculate minimum time step length for occurence of next rain event 
!     (tcum + dt = time at end of current timestep)
         dtEventRain = raintimearray(rainrec) - (tcum + dt)

!
! === Detailed meteo ===================================================

      elseif (flmetdetail) then

        if (flUpdMetDet) then
! --- per meteo time interval: update actual meteo record and set fluxes  
!                              for current time of detailed meteo input
          wrecord = wrecord + 1
          ptra    = tpot(wrecord)
          peva    = epot(wrecord)
          graidt  = grain(wrecord)
          nraidt  = nrain(wrecord)
          aintcdt = graidt - nraidt
!
          flUpdMetDet = .false.
        endif
!
! --- per time step: calculate soil evaporation rate of current time step
        call reduceva (2,swredu,fldaystart,cofred,dt,empreva,           &
     &               ldwet,nird,nraida,peva,pond,rsigni,spev,saev)

      endif
       
      return     
      end subroutine ProcessMeteoTsteps


! SUBROUTINE 4.
! ----------------------------------------------------------------------
      subroutine ETSine
! ----------------------------------------------------------------------
!     Last modified      : october 2008
!     Purpose            : distributes potential transpiration and evaporation
!                          according to sine wave during photoperiodic daylight
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local
      real(8)   daytime,pi,dayl,sinld,cosld,fraction
      real(8), save  :: tsunrise, tsunset
      data      pi/3.14159265d0/    ! number pi [-]

      if (fldaystart) then
! ---   determine duration photoperiodic daylight in hours
        call astro(daynr,lat,rad,dayl,daylp,sinld,cosld,difpp,          &
     &             atmtr,dsinbe)
! ---   determine tsunrise, tsunset and daytime
        tsunrise = 0.5d0 - daylp / 48.d0
        tsunset = 0.5d0 + daylp / 48.d0
      endif

! --- set time as fraction of the day
      daytime = t1900 + dt - int(t1900)

! --- determine fraction of fluxes according to sine wave during this time step
      if (daytime.lt.tsunrise) then
        fraction = 0.d0
      elseif (daytime.gt.tsunrise .and. (daytime-dt).lt.tsunrise) then
        fraction = 0.5d0 * (cos(pi/2.d0 + (tsunrise - 0.5d0)/           &
     &     (tsunset-tsunrise)*pi) - cos(pi/2.d0 + (daytime - 0.5d0)/    &
     &     (tsunset-tsunrise)*pi))
      elseif ((daytime-dt).gt.tsunrise .and. (daytime).lt.tsunset) then
        fraction = 0.5d0 * (cos(pi/2.d0 + (daytime - dt - 0.5d0)/       &
     &     (tsunset-tsunrise)*pi) - cos(pi/2.d0 + (daytime - 0.5d0)/    &
     &     (tsunset-tsunrise)*pi))
      elseif (daytime.gt.tsunset .and. (daytime-dt).lt.tsunset) then
        fraction = 0.5d0 * (cos(pi/2.d0 + (daytime - dt - 0.5d0)/       &
     &     (tsunset-tsunrise)*pi) - cos(pi/2.d0 + (tsunset - 0.5d0)/    &
     &     (tsunset-tsunrise)*pi))
      else
        fraction = 0.d0
      endif

! --- set E and T fluxes
      peva = pevaday * fraction / dt
      ptra = ptraday * fraction / dt

! --- actual soil evaporation rate of current moment 
      call reduceva (2,swredu,fldaystart,cofred,dt,empreva,             &
     &               ldwet,nird,nraida,peva,pond,rsigni,spev,saev)   

      return
      end subroutine ETSine
