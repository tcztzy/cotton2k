MODULE ttutilPrefs
   PRIVATE

!  settings
   INTEGER, PUBLIC, PARAMETER :: TTutilNameLEN = 31

!  defaults
   LOGICAL, PUBLIC, PARAMETER :: TOSCR_default = .true.
   LOGICAL, PUBLIC, PARAMETER :: TOLOG_default = .false.
   INTEGER, PUBLIC, PARAMETER :: UNLOG_default = 0

!  message settings
   LOGICAL, PUBLIC, SAVE :: TOSCR = TOSCR_default
   LOGICAL, PUBLIC, SAVE :: TOLOG = TOLOG_default
   INTEGER, PUBLIC, SAVE :: UNLOG = UNLOG_default

!  status
   LOGICAL, PUBLIC, SAVE :: MessagePrefsSet = .false.

!  error message modes
   INTEGER, PUBLIC, PARAMETER :: FatalErrorDefault    = 0
   INTEGER, PUBLIC, PARAMETER :: FatalErrorERRfile    = 1
   INTEGER, PUBLIC, PARAMETER :: FatalErrorInternal   = 2
   INTEGER, PUBLIC, PARAMETER :: FatalErrorExceptFile = 3

!  actual error mode
   INTEGER, PUBLIC, SAVE :: FatalErrorMode = FatalErrorDefault

!  external message
   INTEGER, SAVE      :: MessLineCnt = 0
   INTEGER, SAVE      :: MessActualLEN
   INTEGER, PARAMETER :: MessLEN = 80
   CHARACTER(LEN=MessLEN), DIMENSION(:), ALLOCATABLE, SAVE :: ExternalMessage

!  screen message before STOP
   CHARACTER(LEN=MessLEN), PUBLIC, SAVE :: FatalErrorScreenTextSTOP = ' '

!  supplied error file
   INTEGER, PARAMETER :: FnamLEN = 80
   CHARACTER(LEN=*), PARAMETER :: ErrorFileNameDefault = 'MODEL_ERRORS.TXT'
   CHARACTER(LEN=FnamLEN), PUBLIC, SAVE :: FatalErrorFileName = ErrorFileNameDefault

!  supplied unit number of open exceptions file
   INTEGER, PARAMETER    :: ExceptionUnitDefault = -99
   INTEGER, PUBLIC, SAVE :: ExceptionUnit = ExceptionUnitDefault

   PUBLIC FillExternalMessage, WriteExternalMessage, SetErrorFileName, SetExceptionFile
CONTAINS
! =================================================================
SUBROUTINE FillExternalMessage (CharArray)
! =================================================================
!  allocates and fills special external error message with input text
   IMPLICIT NONE
!  administration
   CHARACTER (LEN=*), PARAMETER :: SubProgram = 'FillExternalMessage'
!  formal
   CHARACTER(LEN=*), DIMENSION(:) :: CharArray
!  local
   INTEGER :: iw

   MessActualLEN = len(CharArray)
   if (MessActualLEN > MessLEN) then
      FatalErrorMode = 0
      call FatalERR ('FillMessage','Message lines too long')
   end if

   MessLineCnt = size(CharArray) ; iw = 0
   if (allocated(ExternalMessage))  deallocate (ExternalMessage,  stat=iw)
   if (iw/=0) call FatalERR (SubProgram,'cannot de-allocate ExternalMessage')
   allocate (ExternalMessage(MessLineCnt),  stat=iw)
   if (iw/=0) call FatalERR (SubProgram,'cannot allocate ExternalMessage')

   ExternalMessage = CharArray
Return
END SUBROUTINE FillExternalMessage
! =================================================================
SUBROUTINE WriteExternalMessage
! =================================================================
!  writes special external error message
   IMPLICIT NONE
!  local
   INTEGER :: i

   if (TOSCR .and. MessLineCnt>0) then
      write (*,'(1x,a)') Repeat('-',MessActualLEN)
      do i=1,MessLineCnt
         write (*,'(1x,a)') TRIM(ExternalMessage(i))
      end do
      write (*,'(1x,a,/)') Repeat('-',MessActualLEN)
   end if

   if (TOLOG .and. MessLineCnt>0) then
      write (UNLOG,'(1x,a)') Repeat('-',MessActualLEN)
      do i=1,MessLineCnt
         write (UNLOG,'(1x,a)') TRIM(ExternalMessage(i))
      end do
      write (UNLOG,'(1x,a,/)') Repeat('-',MessActualLEN)
   end if

   if (FatalErrorMode == FatalErrorExceptFile .and. MessLineCnt>0 .and. ExceptionUnit > 0) then
      write (ExceptionUnit,'(a)') Repeat('-',MessActualLEN)
      do i=1,MessLineCnt
         write (ExceptionUnit,'(a)') TRIM(ExternalMessage(i))
      end do
      write (ExceptionUnit,'(a,/)') Repeat('-',MessActualLEN)
   end if
Return
END SUBROUTINE WriteExternalMessage
! =================================================================
SUBROUTINE SetErrorFileName (Filename)
! =================================================================
!  sets the name of the error file name used in mode FatalErrorUSERfile
   IMPLICIT NONE
!  administration
   CHARACTER (LEN=*), PARAMETER :: SubProgram = 'SetErrorFileName'
!  formal
   CHARACTER(LEN=*) :: Filename

   FatalErrorFileName = Filename
Return
END SUBROUTINE SetErrorFileName
! =================================================================
SUBROUTINE SetExceptionFile (Unit, Filename)
! =================================================================
!  sets the name of the error file name used in mode FatalErrorUSERfile
   IMPLICIT NONE
!  administration
   CHARACTER (LEN=*), PARAMETER :: SubProgram = 'SetErrorFileName'
!  formal
   INTEGER, INTENT(IN)          :: Unit
   CHARACTER(LEN=*), INTENT(IN) :: Filename

!  store in module data
   ExceptionUnit      = Unit
   FatalErrorFileName = Filename
Return
END SUBROUTINE SetExceptionFile
END MODULE ttutilPrefs
