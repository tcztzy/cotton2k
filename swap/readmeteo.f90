! File VersionID:
!   $Id: readmeteo.f90 328 2017-06-03 21:05:50Z kroes006 $
!
!     This file contains the following subroutines, in order of calling:
!     1. ReadMeteo      : reads meteorological input data; called in SWAP
!     2. ReadRainEvents : reads input data on rain events; called in ReadMeteo (optional)

! SUBROUTINE 1.
! ----------------------------------------------------------------------
      subroutine ReadMeteo
! ----------------------------------------------------------------------
!     Last modified      : March 2014
!     Purpose            : read meteorological data of one calendar year
! ----------------------------------------------------------------------
      use variables
      implicit none
    
! --- global
!   - general     

! --- local
      character(len=11)  datedum
      character(len=3)   ext
      character(len=200) filnam
      character(len=800) messag
      character(len=80)  station(366)
      integer   ad(mrain),am(mrain),datea(6),daynumber,i,ifnd
      integer   getun2,wth      
      real(4)   fsec
      real(8)   etrmax,etrmin,hummax,hummin,radmax,radmin,raimax,raimin
      real(8)   tmeteo,tmnmax,tmnmin,tmxmax,tmxmin,winmax,winmin
      
!========================= Read Meteo file =============================

! --- detection & handling of missing values ---------------------------
!       default ranges of parameter values
      radmin = -1.0d5
      radmax = 5.0d6
      tmnmin = -1.0d5
      tmnmax = 1.0d5
      tmxmin = -1.0d5
      tmxmax = 1.0d5
      hummin = -1.0d5
      hummax = 1.0d5
      winmin = -1.0d5
      winmax = 1.0d5
      raimin = 0.0d0
      raimax = 1.0d3
      etrmin = -1.0d5
      etrmax = 1.0d5
! --- no missing values allowed if penmon must be executed
      if (swetr.eq.0) then
        radmin = 0.0d0
        tmnmin = -50.0d0
        tmnmax = 35.0d0
        tmxmin = -50.0d0
        tmxmax = 60.0d0
        hummin = 0.0d0
        hummax = 10.0d0
        winmin = 0.0d0
        winmax = 150.0d0
      endif
! --- error in case etref missing
      if (swetr.eq.1) then
        etrmin = -0.00001d0
        etrmax = 1.0d2
      endif  
! --- no missing values for tmn and tmx allowed if crop development or
! ---   numerical soil temperatures must be simulated
      if (swmeteo.eq.2 .or. swcalt.eq.2) then
        tmnmin = -60.0d0
        tmnmax = 50.0d0
        tmxmin = -50.0d0
        tmxmax = 60.0d0
      endif
! --- no missing value for rad allowed in case the detailed crop model
! ---   or the grass routine is active 
      if (swmeteo .eq. 2) then
        radmin = 0.0d0
        radmax = 5.0d6
      endif
! --- end of handling missing values -----------------------------------

! --- compose filename meteorological file
      write (ext,'(i3.3)') mod(yearmeteo,1000)
      filnam = trim(pathatm)//trim(metfil)//'.'//trim(ext)

! --- initialise and start reading
      wth = getun2 (10,90,2)
      call rdinit(wth,logf,filnam)

! --- get values from file; two options:
!   - 1 swmetdetail = 0; daily input
!   - 2 swmetdetail = 1; detailed input for nmetdetail time intervals per day
      if (swmetdetail.eq.0) then
         call rdacha ('station',station,366,ifnd)
         call rdfinr ('dd',1,31,ad,366,ifnd)
         call rdfinr ('mm',1,12,am,366,ifnd)
         call rdfdor ('rad',radmin,radmax,arad,366,ifnd)
         call rdfdor ('tmin',tmnmin,tmnmax,atmn,366,ifnd)
         call rdfdor ('tmax',tmxmin,tmxmax,atmx,366,ifnd)
         call rdfdor ('hum',hummin,hummax,ahum,366,ifnd)
         call rdfdor ('wind',winmin,winmax,awin,366,ifnd)
         if (swrain.lt.3) call rdfdor('rain',raimin,raimax,arai,366,    &
     &                                ifnd)
         call rdfdor ('etref',etrmin,etrmax,aetr,366,ifnd)
         if (swrain.eq.2) call rdfdor('wet',0.d0,1.d0,wet,366,ifnd)
! -      Convert radiation from kj/m2/d to j/m2/d
         arad = 1000.0d0 * arad
!
      elseif (swmetdetail.eq.1) then
         call rdatim ('date',dettime,nmetfile,ifnd)
         call rdfinr ('record',1,nmetdetail,detrecord,nmetfile,ifnd)
         call rdfdor ('rad',radmin,radmax,detrad,nmetfile,ifnd)   
         call rdfdor ('temp',tmxmin,tmxmax,dettav,nmetfile,ifnd)
         call rdfdor ('hum',hummin,hummax,dethum,nmetfile,ifnd)
         call rdfdor ('wind',winmin,winmax,detwind,nmetfile,ifnd)
         call rdfdor ('rain',raimin,raimax,detrain,nmetfile,ifnd)
! -      Convert radiation from kj/m2/d to j/m2/d
         detrad = 1000.0d0 * detrad
      endif
      
! --- close meteorological file
      close (wth)

 
!========================= tests and initialization ====================

! --- perform some reliability tests and some initialization
      if (swmetdetail.eq.0) then

! --- determine first and last day numbers        
         datea(1) = yearmeteo
         datea(2) = 1
         datea(3) = 1
         datea(4) = 0
         datea(5) = 0
         datea(6) = 0
         fsec = 0.0
         call dtardp (datea,fsec,timjan1)
         datea(2) = am(1)
         datea(3) = ad(1)
         call dtardp (datea,fsec,tmeteo)
         daynrfirst = nint ( tmeteo - timjan1 + 1.0d0 )
         datea(2) = am(ifnd)
         datea(3) = ad(ifnd)

         call dtardp (datea,fsec,tmeteo)
         daynrlast = nint ( tmeteo - timjan1 + 1.0d0 )

! --- check date 
         do i = 2, ifnd-1
            datea(2) = am(i)
            datea(3) = ad(i)
            call dtardp (datea,fsec,tmeteo)
            daynumber = nint ( tmeteo - timjan1 + 1.0d0 )
            if (daynumber .ne. daynrfirst+i-1) then 
!           wrong date after daynumber i-1
              datea(2) = am(i-1)
              datea(3) = ad(i-1)
              call dtardp (datea,fsec,tmeteo)
              call dtdpst ('year-month-day',tmeteo,datedum)           
              messag ='In meteo file '//trim(filnam)//' the date after '&
     &        //datedum//' is not correct! First adapt meteo file!'
              call fatalerr ('meteo',messag)
            endif
         end do
         
! --- snow and frost calculation conditions require realistic air temperatures
         do i = 1, ifnd
            if ( (swsnow.eq.1.or.swfrost.eq.1) .and.                    &
     &         (atmn(i).lt.-98.9d0.or.atmx(i).lt.-98.9d0) ) then    
             messag ='In meteo file '//trim(filnam)//' temperatures'//  &
     &       ' must be input to calculate snow conditions (SWSNOW=1)'// &
     &       ' and/or frost conditions (SWFROST=1)! adapt meteo file.'
             call fatalerr ('meteo',messag)
           endif
         enddo
! --- realistic air temperatures are required when soil temp. is simulated 
!     using the numerical model  (not if data from other file 
!! requires additional verification)
!      if (swhea.eq.1 .and. swcalt.eq.2) then
!        if (atmn(i).lt.-98.9d0 .or. atmx(i).lt.-98.9d0) then
!            messag ='In meteo file realistic temperatures are '//
!     &        ' required for simulation of temperature profiles (SWHEA=1) using'//
!     &        'num.model (SWCALT=2)!  Adapt meteo data! '
!            call fatalerr ('meteo',messag)
!        endif
!      endif

! --- in case of swrain=2 then make sure that Rain and Wet correspond
         if (swrain.eq.2) then
            do i = 1, ifnd
              if ((arai(i).gt.1.d-10 .and. wet(i).lt.1.d-10) .or.       &
     &          (arai(i).lt.1.d-10 .and. wet(i).gt.1.d-10) ) then
                write(messag,1001) trim(filnam), am(i), ad(i)
 1001           format(' In meteo file ',a,';  month =',i3,';  day =',  &
     &             i3,'; SwRain= 2      ',                              &
     &             ' Rain and Wet donot correspond.  Adapt meteo data!')
               call fatalerr ('ReadMeteo',messag)
              endif  
            enddo
         endif
!
! --- in case of swrain = 1 or 2 then store daily precipitation (rain) amount 
!     in rainamount and precipitation time in raintimearray
         if (swrain.eq.1 .or. swrain.eq.2) then
            nmrain = 0
            do i = 1, ifnd
               nmrain   = nmrain + 1
               datea(2) = am(i)
               datea(3) = ad(i)
               call dtardp (datea,fsec,tmeteo)
               raintimearray(i+1) = tmeteo
               rainamount(i)      = arai(i)
            enddo
         endif
! --- end of reliability tests and initialization of daily meteo   

      elseif (swmetdetail.eq.1) then
! --- initialization of detailed meteo
  
! ---   initialize total record number for new weather file
        irectotal = int(t1900-dettime(1)+0.1d0)*nmetdetail

! ---   initialize number of days for running average Tmin
        nofd = 0
      
      endif
! --- end of initialization of detailed meteo   
      
! --- end of reading meteo data file **********************************

! --- close present year for further reading      
      flYearStart = .false.

!========================= Read rain file =============================

      if (swrain .eq. 3) then
! ---   rainfall events are specified
        call ReadRainEvents(logf,yearmeteo,pathatm,rainfil,             &
     &               nmrain,rainamount,raintimearray)
      endif
! --- end of reading rain file ****************************************     

! --- reopen present year for processing rain intensity data at beginning MeteoDt      
      if (swrain .gt. 0) then
         flYearStart = .true.
         call MeteoDt
      endif

      return
      end subroutine ReadMeteo


! SUBROUTINE 2.
! ----------------------------------------------------------------------
      subroutine ReadRainEvents(logf,yearmeteo,pathatm,rainfil,         &
     &               nmrain,rainamount,raintimearray)
! ----------------------------------------------------------------------
!     Last modified      : February 2014
!     Purpose            : read rainfall data (events) of one calendar year 
!     Interface:
!       I   - logf,yearmeteo,pathatm,rainfil
!       O   - nmrain,rainamount,raintimearray
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      character(len=*) rainfil,pathatm
      integer   logf,nmrain,yearmeteo
      real(8)   rainamount(mrain)    ! Array with short duration rainfall amounts (L)
      real(8)   raintimearray(mrain) ! Array with times (T) at which rainfall intensity changes

! --- local
      character(len=2)   chday,chmonth
      character(len=5)   chtime
      character(len=3)   ext
      character(len=200) filnam
      character(len=300) messag
      integer   ad(mrain),am(mrain),ay(mrain)
      integer   datea(6),getun2,i,ic,ifnd,pre
      real(4)   fsec,time4(mrain)
      real(8)   tday,tdayold,vsmall
      logical   flrnx
      vsmall = 1.0d-8
!========================= Read rain file =============================

! --- in case of detailed rainfall, open file once a year
      write (ext,'(i3.3)') mod(yearmeteo,1000)
      filnam = trim(pathatm)//trim(rainfil)//'.'//trim(ext)

! --- initialise and start reading
      pre = getun2 (10,90,2)
      call rdinit(pre,logf,filnam)

! --- get values from file
      call rdainr ('day',1,31,ad,mrain,ifnd)
      call rdfinr ('month',1,12,am,mrain,ifnd)
      call rdfinr ('year',1,3000,ay,mrain,ifnd)
      call rdfdor ('time',0.0d0,1.0d0,raintimearray,mrain,ifnd)
      call rdfdor ('amount',0.0d0,10000.0d0,rainamount,mrain,ifnd)

      close (pre)

! --- check whether first rain event starts at t = 0; if not: repair
      if (raintimearray(1).gt.vsmall) then
         do i = ifnd, 1, -1
           ay(i+1) = ay(i)
           am(i+1) = am(i)
           ad(i+1) = ad(i)
           raintimearray(i+1) = raintimearray(i) 
           rainamount(i+1)    = rainamount(i) 
         enddo
         ifnd = ifnd + 1
         raintimearray(i+1) = 0.0d0
         rainamount(i+1)    = 0.0d0 
      endif

! --- convert year, day, month and fraction of day to absolute time
      fsec = 0.0
      ic = 0
      do i = 1,ifnd
         datea(1) = ay(i)
         datea(2) = am(i)
         datea(3) = ad(i)
         call dtardp (datea,fsec,tday)
         time4(i)         = real(raintimearray(i)) 
         raintimearray(i) = tday + raintimearray(i) 
!   - check whether end of previous day is equal to start of this day (original raintimearray(i) = 0.00)
         if (ic.gt.0) then
            if ( (raintimearray(i)-raintimearray(ic)) .lt. vsmall .and. &
     &           tday.gt.tdayold) then
               flrnx = .false.
            else
               flrnx = .true.
            endif
         else
            flrnx = .true.
         endif
         if (flrnx) then
            ic = ic + 1
            raintimearray(ic) = raintimearray(i)
            rainamount(ic)    = rainamount(i) 
         endif
         tdayold = tday
      end do

! --- save number of rain event records
      nmrain = ic   !  ifnd

! --- check whether dates are in ascending order
      do i = 2,ic
         if ( (raintimearray(i)-raintimearray(i-1)) .lt. vsmall) then
            write(chmonth,'(i2)')  am(i)
            write(chday,'(i2)')    ad(i)
            write(chtime,'(f6.3)') time4(i)
           messag ='Error-message from module ReadRainEvents'//         &
     &        '. In detailed rain file '//trim(filnam)//                &
     &        ' the time of a rainrecord should be greater than the'//  &
     &        ' time of its preceding rainrecord! Adapt rain file.'//   &
     &        '      Month: '//chmonth//'; Day: '//chday//'; Time: '//  &
     &        chtime//'!'
            call fatalerr ('ReadRain',messag)
         endif
      enddo

      return
      end subroutine ReadRainEvents