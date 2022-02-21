! File VersionID:
!   $Id: timecontrol.f90 330 2017-06-12 10:22:55Z kroes006 $
! ----------------------------------------------------------------------
      subroutine TimeControl(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : Handles time variables, switches and flags
! ----------------------------------------------------------------------

      use variables
      implicit none
! ----------------------------------------------------------------------
!     DAYNR  : = daynumber relative to start of calendar year
!     DAYCUM : = daynumber relative to start of simulation
!     DAYCROP: = daynumber relative to emergence of crop
!     DAYMETEO: = day number for which meteorological data should be read
!     IMONTH : = actual month number
!     IYEAR  : = actual year
!
!     T      : = time relative to start of calendar year
!     TCUM   : = time relative to start of simulation
!     TEND   : = time (relative to 1900) at which simulation ends
!     TSTART : = time (relative to 1900) at which simulation starts
!     t1900  : = time relative to 1900
!
!     flCultivate : = true if crops are grown on the soil profile
!     flDayEnd : = true if time level equals end of day
!     flRunEnd : = true if time level equals TEND
!     flRunStart  : = true if time level equals TSTART
! ----------------------------------------------------------------------
!     local
      character(len=200) messag
      character(len=11)  tmp
      character(len=80)  filtext
      integer       task, itask
      integer, save :: datea(6),nextyear
      integer, save :: flprevious  ! (status of previous timestep interval)
!                   1 : dt as a result of reduction due to numerical reasons
!                   2 : dt as a result of reduction due to input / output control
      logical, save :: flTnext
      real(4), save :: fsec
      real(8), save :: dtCrit,tchange,dtEvent,tEvent,dtfletsine
      real(8), save :: tcumold
      real(8), save :: dtprevious  ! (length of previous timestep)
      real(8)       :: dtRestDay


!     critical time-interval 
      dtCrit = 1.d-8
      dtfletsine = 0.05d0


      itask = task
      if (itask.eq.2 .and. (fldecdt .or. fldecdtmin)) itask = 3

      select case (task)
      case (1)

! === initialization ===================================================

! --- initialize flags ----------------------------
      fldecdt = .false.
      fldecdtmin = .false.
      flRunEnd = .false.
!      flRunStart = .true.
      flDayStart = .true.
      fldtmin = .false.
      flZeroIntr = .true.
      flZeroCumu = .true.
      floutput = .false.
      flbaloutput = .false.
      flheader = .false.
      flheadirg = .false.
      flIrg1Start = .true.
      flUpdMetDet = .true.
      flYearStart = .true.


      if (nprintday .gt. 1 .or. flprintdt) then
        flprintshort = .true.
        period = 1
      else
        flprintshort = .false.
      endif
      if (swmetdetail .eq. 1) then
        flmetdetail = .true.
        metperiod = 1.0d0 / dble(nmetdetail)
        wrecord = 0
      else
        flmetdetail = .false.
      endif
      if (swrain.gt.0) then    
         flrainintens = .true.
      else
         flrainintens = .false.
      endif   
      if (flmetdetail .or. flrainintens) then 
         flmeteodt = .true.
      else
         flmeteodt = .false.
      endif   
      fletsine = .false.
      if (swetsine .eq. 1) fletsine = .true.
      flIrrigate = .false.
      if (swirfix.eq.1) flIrrigate = .true.
      flMacroPore = .false.
      if (swmacro.eq.1) flMacroPore = .true.
      flDrain = .false.
      if (swdra .eq. 1) flDrain = .true.
      flSurfaceWater = .false.
      if (swdra .eq. 2) flSurfaceWater = .true.
      flTemperature = .false.
      if (swhea .eq. 1) flTemperature = .true.
      flSnow = .false.
      if (swsnow .eq. 1) flSnow = .true.
      flSolute = .false.
      if (swsolu .eq. 1) flSolute = .true.
        
! --- initialize counters ----------------------------
      nirri = 1
      isteps = 0
      ioutdat = 1
      ioutdatint = 1
      cntper = 0
      outper = 0.0d0
      tcumold = 0.0d0
      nprintcount = 1

! --- set main time variable of SWAP ------------------------
      t1900 = tstart

! --- determine time from beginning of calendar year
      call dtdpar (t1900, datea, fsec)
      datea(1) = iyear
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      t = tstart - timjan1
      tEvent = 0.0d0
      tcum = 0.d0
      daynr = nint(t)
      daycum = 0
      daymeteo = daynr + 1

! --- determine year,month and day number of current day
      call dtdpar (t1900+0.1d0, datea, fsec)
      iyear = datea(1)
      imonth = datea(2)
      yearmeteo = iyear

! --- determine date of current day
      call dtdpst ('year-month-day',t1900+0.1d0,date)

! --- output to screen
      if (swscre .eq. 2) then
        filtext = 'Screen output of daynumbers'
        call writehead (5,1,'screen',filtext,project)
        call dtdpst ('year-month-day',tstart,date)
        write (*,'(2x,2a)') 'First day of simulation:  ',date
        call dtdpst ('year-month-day',tend,date)
        write (*,'(2x,2a)') 'Last day of simulation:   ',date
        write (*,'(/,a,/)') '            date  daynr  daycum'
      endif

! -   set crop number
      icrop = 1
      do while (.not. flCropCalendar) 
          
        if (cropstart(icrop) .lt. 1.d0) exit
          
        if (tstart - cropstart(icrop) .gt. -1.d-3                       &
     &                   .and. tstart - cropend(icrop) .lt. 1.d-3) then
          flCropCalendar = .true.
          if (flCropCalendar) then
            if (tstart - cropstart(icrop) .lt. -1.d-3 .and.             &
     &                                              swinco .ne. 3) then
              messag = 'The start of simulation (tstart) begins in '//  &
     &        'crop growing season with swinco 1 or 2'
              call fatalerr ('readswap',messag)
            endif
          end if
        else
          icrop = icrop + 1
        end if
      enddo
      
! --- detailed meteo data needed for crop growth?
      swmeteo = 1
      if (flCropCalendar) then
        if (icrop .gt. 0) then
          if(croptype(icrop).ge.2) then
             swmeteo = 2
          endif
        endif
      endif

! --- initialize
      dtEvent = 1.0d0

! --- limit initial dt in case of short time interval
      if (flprintshort .and. .not.flprintdt) then
        if (dtEvent .gt. dble(nprintcount)/dble(nprintday))then
           dtEvent = dble(nprintcount)/dble(nprintday)
         end if
      endif

! --- intial timestep for rainfall intensities
      if (swmetdetail.eq.0 .and. swrain.gt.0) then
         dtEvent = min(dtEvent, dtmin)
      endif

! --- limit initial dt in case of detailed meteorological input
       if (flmetdetail) then
        if (dtEvent .gt. metperiod)then 
           dtEvent = min(dtEvent, metperiod)
         end if
      endif
      tEvent = dtEvent

! --- initial time step
      if (swinco.eq.3) then
        if (dt.lt.dtmin) then
          messag = 'Initial dt read from file (SWINCO=3)'//             &
     &    ' if absent, then default is assumed'
          call warn ('soilwater',messag,logf,swscre)
          dt = sqrt(dtmin*dtmax)
        endif
      else
        dt = sqrt(dtmin*dtmax)
      endif
      dtprevious = dt
      if(dt+dtCrit .gt. dtEvent)then
         flprevious = 2
         dt = dtEvent
         flTnext = .true.
      else
         flprevious = 1
         flTnext = .false.
      end if

! --- in case of sine wave of ET, limit dt and dtmax
      if (fletsine) then
        dt = min(dt, dtfletsine)
        dtmax = min(dtmax, dtfletsine)
      endif

! --- set value for dtold, used in headcalc for macropore-iteration
      dtold = dt

      return

      case (2)

! === next time step ===================================================

! 2.1  check maximum number of time steps during this day
      isteps = isteps + 1
      if (isteps .gt. msteps) then
        write(tmp,'(i11)') daynr
        tmp = adjustl(tmp)
        messag ='The maximum number of time steps for a day is exceeded'&
     &    //' at daynumber '//trim(tmp)//'. Check input for numerical'  &
     &    //' solution of Richards equation'
        call fatalerr ('timer',messag)
      endif


! 2.2 update time variables

      t = t + dt             ! relative to yyyy0101:00:00:00
      tcum = tcum + dt
      t1900 = tstart + tcum

! 2.3  flag assignments after first time step of a day

      if (flDayStart) then

! ---   set crop conditions
        if (flCropCalendar) then
          flCropOutput = .true.
          if (flCropHarvest) then
            flCropOutput = .false.
          endif
        endif
        
! --    length crop
        if (flCropCalendar) then
          daycrop = daycrop + 1
        endif

        
! ---   set flags for reset intermediate and cumulative fluxes
        flZeroCumu = .false.
        if (flZeroIntr) then
          outper = 0.0d0
          flZeroIntr = .false.
        endif

! ---   set flags for output
        if (floutput) then
         floutput = .false.
         if (swheader .eq. 1) flheader = .false.
        endif
        if (flbaloutput) then
         flbaloutput = .false.
         if (swheader .eq. 1) then
           flheader = .true.
           flheadirg = .true.
         endif
        endif

! 2.4  determine year,month and day number (only during first time step of a day)
        call dtdpar (t1900, datea, fsec)
        iyearm1 = iyear
        iyear = datea(1)
        imonth = datea(2)

! ---   determine date of current day
        call dtdpst ('year-month-day',t1900,date)

! ---   update day numbers
        daynr = daynr+1
        daycum = daycum + 1
        cntper = cntper + 1

! 2.5  in case of detailed meteorological input, reset weather record

        if (iyear .ne. iyearm1) then
! ---     reset daynumber and time because new calender year has started
          daynr = 1
          t = dt

! ---     in case SWRES = 1 reset counter for periodic output to 1
          if (swres.eq.1 .and. period.ne.0) cntper = 1
        endif


      endif


! 2.6  update  logicals for indication of start / end of day
      if ( dble(daycum) - tcum .lt. dtCrit) then 
        flDayEnd = .true.
        flDayStart = .true.
        flTnext = .true.
        if (flmetdetail) wrecord = 1
      else
        flDayEnd = .false.
        flDayStart = .false.
      endif


! 2.7  maximum size of time interval (dtEvent)

      if(flTnext)then 

! 2.7.1  remaining part of a day
        dtEvent = dble(int(tcum+1.0d0+dtCrit)) - tcum

! 2.7.2  printing more than one times a day may limit timestep
         if (flprintshort  .and. .not.flprintdt) then
           if (dble(nprintcount)/dble(nprintday)-tcum .lt. dtCrit) then
              dtEvent = min(dtEvent,1.0d0/dble(nprintday))
           else
              dtEvent = min(dtEvent,                                    &
     &                   (dble(nprintcount)/dble(nprintday)) - tcum)
           endif
         endif

! 2.7.3  input of detailed meteo may limit timestep
         if (flmetdetail) then
            tchange = dble(int(t + dtCrit)) + dble(wrecord) * metperiod
            if ((t+dtEvent) .gt. tchange) dtEvent = tchange - t
         end if

! 2.7.4  precipitation event may limit timestep
         if (swmetdetail.eq.0 .and. swrain.gt.0) then

!        next rainevent! Set new values
           if (raintimearray(rainrec).lt.tcum+dtCrit) then
!        new rain event valid
             rainrec = rainrec + 1
           endif

           dtEvent = min(dtevent,dtEventRain)
           dtEvent = max(dtEvent,dtmin)
         endif

! 2.8  set end of time interval (determined by I/O)
         tEvent = tEvent + dtEvent
         flTnext = .false.

      endif

! 2.9 determine next time step, based on numerical performance
      if (flDayStart) then
        if(flprevious.eq.2)then
           dt = max(dt, sqrt(dtmin*dtmax),dtprevious)
        else
           dt = max(dt, sqrt(dtmin*dtmax))
        end if
        dtprevious = dt
      else
        if(flprevious .eq. 2)then
           dt = dtprevious
        else
           if (numbit.le.3)     dt = min(dt*2.0d0,DtMax)
           if (numbit.ge.MaxIt) dt = max(dt*0.5d0,DtMin)
           dtprevious = dt
         endif
      endif
      flprevious = 1

      if ( tcum + dt - tEvent .gt. dtCrit) then
         dt = tEvent - tcum
         flprevious = 2
         flTnext = .true.
      endif

! 2.10  test last time step of the day: limit dt if it exceeds end of day
      dtRestDay = dble(int(tcum+1.0d0+dtCrit)) - tcum
      if (dtRestDay-dt.lt.1.d-6) then
        dt = dtRestDay
      endif

! JK20131230: when dtmin is large then the timestep-closure may not be correct and errors may occur
! resulting in water balance errors see email Paul v Walsum 210131225. Elimination of next statement
! is not the correct solution
      dt = max(dt,dtmin) 
      dt = min(dt,dtmax)

! 2.11 set flags and variables

! --- in case of output during a day
      if (flprintshort) then
        floutputshort = .false.
        flzerointr = .false.
! ---   determine whether output is required
        if (flprintdt) then 
           outper = tcum - tcumold
           tcumold = tcum
           if (abs(outdatint(ioutdatint) - t1900 + 1.d0).lt.1.d-3) then
              floutputshort = .true.
              flzerointr = .true.
           endif
        else
           if (tcum+dtCrit .gt. dble(nprintcount)/dble(nprintday)) then
              floutputshort = .true.
              flzerointr = .true.
              outper = tcum - tcumold
              tcumold = tcum
           endif
! ---      update counter nprintcount for printing
           do while (tcum+dtCrit .gt. dble(nprintcount)/dble(nprintday)) 
              nprintcount = nprintcount + 1
           end do
         endif
      endif

! --- in case of detailed meteorological input
      if (flmetdetail) then
        tchange = dble(int(t + dtCrit)) + dble(wrecord) * metperiod
        if ((tchange - t) .lt. dtCrit) then
! ---     update actual weather record and fluxes

          flUpdMetDet = .true.
        endif
      endif


! --- update fldtmin
      if (dt .gt. (1.0d0+dtCrit)*dtmin) fldtmin = .false.


! --- procedure when day is finished ---------------------------------------------- 
      if (flDayEnd) then

!! --    length crop
!        if (flCropCalendar) then
!          daycrop = daycrop + 1
!        endif
         
! ---   length output period
        if (.not.flprintshort) then
          outper = outper + 1.0d0
        endif

! ---   write daynumber to screen
        if (swscre .eq. 2) then
          write(*,'("+ ",4x,a11,i6,i8)') date,daynr,daycum
        endif

! ---   end of run?
        if ((tend - t1900 + 1.d0) .lt. 1.d-3) then
          flRunEnd = .true.
          floutput = .true.
          flbaloutput = .true.
          ioutdat = ioutdat + 1
          return
        endif

! ---   in case no end of run, determine whether today output should be written 
        if (cntper .eq. period .and. .not.flprintdt) then
          cntper = 0
          floutput = .true.
          flzerointr = .true.
        endif

        if(flprintdt) then
          if (abs(outdatint(ioutdatint) - t1900 + dt) .lt. 1.d-3) then
            floutput = .true.
            flzerointr = .true.
            ioutdatint = ioutdatint + 1
          endif
        else
          if (abs(outdatint(ioutdatint) - t1900 + 1.d0) .lt. 1.d-3) then
            floutput = .true.
            flzerointr = .true.
            ioutdatint = ioutdatint + 1
          endif
        endif
        if (abs(outdat(ioutdat) - t1900 + 1.d0) .lt. 1.d-3) then
! ---     output of water and solute balances
          floutput = .true.
          flbaloutput = .true.
          flzerointr = .true.
          flzerocumu = .true.
          ioutdat = ioutdat + 1
        endif

! ---   set flharvestday for crop output
        if (abs(t1900 - 1.d0 - cropend(icrop)) .lt. 1.d-3 .and.         &
     &          (.not. flCropHarvest)) flharvestday = .true.

! ---   reset flags for next day
        fldtmin = .false.

! ---   reset counters for next day
        isteps = 0

! ---   determine daynumber and switch for reading meteorological data
        call dtdpar (t1900 + 0.1d0, datea, fsec)
!        call dtdpar (t1900 , datea, fsec)
        nextyear = datea(1)
        if (nextyear .eq. iyear) then
          daymeteo = daynr + 1
        else
          yearmeteo = nextyear
! ---     set flag for new meteo year
          flYearStart = .true.

          daymeteo = 1
! ---     detailed meteo data needed for crop growth?
          swmeteo = 1
          if (flCropCalendar) then
            if (icrop .gt. 0) then
              if (croptype(icrop).ge.2) then
                swmeteo = 2
              endif
            endif
          endif
        endif

      endif

      return

      case (3)

! === reduce time step ===================================================


! --- decrease time step in case of no convergence in headcalc
      if (fldecdt) then
        if (dt .gt. 3.0*dtmin) then
          dt = dt / 3.0
!         force dt to equal multiple dtmin to prevent very small dt-values at end of day 
!          dt = dtmin * dble(max(1,int(dt/dtmin)))
        else
          dt = dtmin
          fldtmin = .true.
        endif
        fldecdt = .false.
        flprevious = 1
        dtprevious = dt
        flTnext = .false.

        return
      endif

! --- decrease time step to dtmin if required by boundtop
      if (fldecdtmin) then
        dt = dtmin
        fldtmin = .true.
        fldecdtmin = .false.
        flprevious = 1
        dtprevious = dt
        return
      endif

! --- decrease in case of Macropores
      if (flMacroPore .and. FlDecMpRat) then
        dt = dsqrt(dtmin*dtmax)
        dtprevious = dt
        return
      endif

      case default
         call fatalerr ('TimeControl', 'Illegal value for TASK')
      end select

      return
      end 


      subroutine IterTime(task)
! ----------------------------------------------------------------------
!     date    : 20080303
!     update  : 20170223: intro part2 - to interrupt (near) endless simulations
!     purpose : statistics of timing and numerical iterations
! ----------------------------------------------------------------------
! --- global variables
      use variables
      implicit none

! --- local variables
      integer task, i, j, timediff
      character(len=400) messag
      real(4)       ::   tmptimeinterrupt
      real(4), save ::   tmptimestart,tmptimeend

      select case (task)

      case (1)
! --- part1 - initial values
      call cpu_time(tmptimestart)
      return

      case (2)
! --- part2 - calculate intermediate time to be able to interrupt (near) endless simulations
      call cpu_time(tmptimeinterrupt)
      timediff = int(tmptimeinterrupt)-MaxIterTime
!     fatal error if cpu time exceeds input value MaxIterTime
      if(timediff.gt.0) then
        write(messag,'(a,i10,3a)')                                      &
     &     'The maximum cpu time of ',MaxIterTime,' (secs)',            &
     &     ' was exceeded.  Therefore simulation was interrupted'
        call fatalerr ('IterTime',messag)
      endif
      return

      case (3)
! --- part3 - write statistics
      write(logf,'(/,a20)')      'Iteration statistics'
      write(logf,'(/,a29,i4)')   'Maximum number of iterations:',MaxIt
      write(logf,'(/,a35/,a35)') 'It Numb  No of Hits  Tot BTr cycles', &
     &                           '-------  ----------  --------------'
      do i=1,100                                       
        if(itnumb(i,1).gt.0)                                            &
     &     write(logf,'(i7,2x,i10,4x,i10)')i,(itnumb(i,j),j=1,2)
      end do

      call cpu_time(tmptimeend)
      write(logf,'(/,a12,f12.2,a4)')                                    &
     &           ' Run-time: ',tmptimeend-tmptimestart,' sec'
     

      case default
         call fatalerr ('IterTime', 'Illegal value for TASK')
      end select

      return
      end
