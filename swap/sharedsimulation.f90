! File VersionID:
!   $Id: sharedsimulation.f90 298 2016-07-25 20:07:31Z kroes006 $
! ----------------------------------------------------------------------
      subroutine SharedSimulation(task) 
! ----------------------------------------------------------------------
!     Date               : July 2009
!     Purpose            : open and write data to shared files
! ----------------------------------------------------------------------

!      use Variables
! Preprocessor directive WINXPIVF for MS-WINdows XP and Intel Visual Fortran compiler
!D#ifdef WINXPIVF
!D     use ifport
!D#endif
! Preprocessor directive WINXPCVF for MS-WINdows XP and Compaq Visual Fortran compiler
!D#ifdef WINXPCVF
!D     use dflib
!D#endif
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer PosArg,status
!D      integer(4) delay
!D      logical    flhold
      integer    unss, ID_Shared, getun       !, IDread, NumChar
      character(len=3)   strIDss
      character(len=80)  strFINA
      character(len=200) messag

      save  unss, ID_Shared
! ----------------------------------------------------------------------
!D      data  delay/100/                       ! delay in milisecs


      select case (task)
      case (1)

      PosArg         = 2
!      Call GetArg (PosArg,strIDss,status)
      Call Get_Command_Argument (number=PosArg,value=strIDss,status=status)

      if(status.ne.2) then
        messag = 'Argument of executable-call is not correct for '//    &
     &                    ' Shared simulation !'
        call fatalerr ('readswap',messag)
      endif

      read(strIDss,'(i3.3)')ID_Shared

      PosArg         = 3
!      Call GetArg (PosArg,strFINA,status)
      Call Get_Command_Argument (number=PosArg,value=strFINA,status=status)

      unss  = getun (20,99)
      open(unit=unss,file=strFINA,status='unknown',                     &
     &     action='READWRITE') !,share='DENYNONE')


!     open shared data file
      call FromSwap(task)
      call ToSwap(task)


      return

      case (2)

!D     flhold = .true.
!D     do while (flhold)
!D#ifdef WINXPCVF
!D        call sleepqq(delay)             ! delay in milisecs
!D#endif
!D#ifdef WINXPIVF
!D        call sleepqq(delay)             ! delay in milisecs
!D#endif
!D        rewind(unss)
!D        read(unss,'(i4)')IDread
!D        if(IDread.eq.ID_Shared) flhold=.false.     
!D      end do

!     read New data

      call ToSwap(2)


      return

      case (3)

      rewind(unss)
      write(unss,'(i4)')-1*ID_Shared

!     write New data

      call FromSwap(2)

      return

      case (4)

! === close Shared Directive file ===========================

      call ToSwap(3)
      call FromSwap(3)

      close (unss)

      case default
         call fatalerr ('SharedSimulation', 'Illegal value for TASK')
      end select

      return
      end
