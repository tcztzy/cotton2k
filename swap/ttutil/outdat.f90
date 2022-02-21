MODULE Module_OutDat
! storage of variable values and output table construction for simulation
! oct 2011 (DK) - addition of optional memory storage
! jan 2012 (KR) - restructured source code, functionality for all output formats unchanged
!                 also variable names kept the same, use of NAME(36:36) by OUTAR2 removed
  IMPLICIT NONE

! default
  PRIVATE

! administration
  CHARACTER (LEN=*), PARAMETER :: ModuleName = 'OUTDAT'

! NAMES_MXN - maximum number of names of dependent variables, can be increased without problems
! VARL_M    - maximum length of a variable name is 36 characters (31 for name, 5 for possible index e.g. '(134)')
  INTEGER, PARAMETER :: NAMES_MXN      = 1000
  INTEGER, PARAMETER :: VARL_M         = 36             ! DO NOT CHANGE
  INTEGER, PARAMETER :: NAMES_IB_MXN_4 = 8              ! Maximum number of names in block in itask=4
  INTEGER, PARAMETER :: NAMES_IB_MXN_5 = 100            ! Maximum number of names in block in itask=5
  INTEGER, PARAMETER :: NAMES_IB_MXN_7 = 8              ! Maximum number of names in block in itask=7
  INTEGER, PARAMETER :: NAMES_IB_MXN_8 = NAMES_IB_MXN_5 ! Maximum number of names in block in itask=8
  INTEGER, PARAMETER :: NAMES_IB_MXN_9 = NAMES_IB_MXN_5 ! Maximum number of names in block in itask=9

! memory storage replacing the use of a binary file
! the record structure is the same as the one on file
  TYPE :: TMemBin
     CHARACTER(LEN=1)      :: RunType
     INTEGER               :: InfoID
     CHARACTER(LEN=VARL_M) :: VarName
     REAL                  :: VarValue
  END TYPE TMemBin
! array allocated in MemWrite
  TYPE(TMemBin), ALLOCATABLE, DIMENSION(:), SAVE :: MemBin
  INTEGER, SAVE                                  :: MemBinHigh = 0
  LOGICAL, SAVE                                  :: UseMemoryStorage
  LOGICAL, SAVE, PUBLIC                          :: OutDat_UseMemory = .false.

! variable name storage
  CHARACTER(LEN=VARL_M), DIMENSION(NAMES_MXN), SAVE :: AVN, ASELN ! Array of (selected) variable names
  INTEGER, DIMENSION(NAMES_MXN), SAVE               :: FND        ! Array with counts how many times a variable is seen at itask=2
  LOGICAL, DIMENSION(NAMES_MXN), SAVE               :: MGIVEN     ! Array of flags whether message about repeated input is given
  CHARACTER(LEN=VARL_M), SAVE                       :: LXN        ! Name of independent variable

! control variables
  INTEGER, SAVE ::  ITOLD = -1       ! previous task
  INTEGER, SAVE ::  IRUN1 =  0       ! number of initializations (runs) sofar (NOT reset at 99 call)
  INTEGER, SAVE ::  INSEL =  0       ! number of selected variables
  LOGICAL, SAVE :: RECOVR = .false.  ! flags recovery just taking place
  LOGICAL, SAVE :: FIRST8 = .true.   ! comment header for end of run output has been written
  LOGICAL, SAVE :: FIRST9 = .true.   ! similar flag for TAB-delimited output

! other saved parameters
  LOGICAL, SAVE :: TOSCR     ! flags screen use for messages
  LOGICAL, SAVE :: TOLOG     ! flags logfile use (currently not used)
  INTEGER, SAVE :: UNLOG     ! logfile unit number (currently not used)
  INTEGER, SAVE :: ILU1      ! unit number temporary, binary storage file
  INTEGER, SAVE :: ILU2      ! unit number formatted output file
  INTEGER, SAVE :: IRUN3     ! run number being stored (current runID)
  INTEGER, SAVE :: IRUN4     ! number of runs in BinaryFile
  INTEGER, SAVE :: IREC      ! first record to write to
  INTEGER, SAVE :: ISREC     ! start record of set
  INTEGER, SAVE :: IFND      ! number of different variable names that have been found in a particular run
  INTEGER, SAVE :: ISAVE     ! number of previous variable name

  CHARACTER(LEN=1) , SAVE :: RUNTYP     ! run type ('N' = for just counting, 'R' = obtained from RDFROM)
  CHARACTER(LEN=80), SAVE :: BinaryFile ! name of binary storage file, which is derived from name of formatted file

! declare public SUBROUTINEs and functions
  PUBLIC :: ProcessOutputData

CONTAINS

! =================================================================
SUBROUTINE ProcessOutputData (ITASK, IUNIT, RN, R)
! =================================================================
! write output variables to temporary storage and construct output table
  IMPLICIT NONE

! formal parameters
  INTEGER, INTENT(IN)           :: ITASK, IUNIT
  REAL, INTENT(IN)              :: R
  CHARACTER(LEN=*), INTENT(IN)  :: RN

! non-saved variables
  INTEGER               :: INAME   ! variable number
  CHARACTER(LEN=512)    :: String  ! used to derive BinaryFile from name of formatted file
  CHARACTER(LEN=VARL_M) :: LN      ! local variable name

! help variables
  LOGICAL :: OPEND, OK
  INTEGER :: il, ip, IRUN2, I1

! function without Fortran-90 interface
  INTEGER :: IFINDC

  if (ITASK == 1) then
!    --------------
!    INITIALIZATION
!    --------------
!    unit number and status check:
!    unit number must be > 0 at first call, may be zero or equal to value at first call
     if (ITOLD == -1) then
!       first initialization, or initialization after reset call (itask=99)
        UseMemoryStorage = OutDat_UseMemory

!       valid initialization call ; desired message output
        call MESSINQ (TOSCR, TOLOG, UNLOG)

!       outdat writes messages to its formatted output file and no additional logfile output is generated.
!       so messages to screen require toscr to tbe true tolog and unlog are ignored.

        if (IUNIT == 0) call FatalERR (ModuleName,'no unit number supplied')
        ILU2 = IUNIT

     else
!       initialization after previous activity ; check
        if (IUNIT /= 0 .and. IUNIT /= ILU2) call FatalERR (ModuleName,'change of unit number not allowed')

        if (ITOLD == 1) then
!          repeated initialization is taking place
           if (TOSCR) write (*,'(A)') ' WARNING from OUTDAT: ignoring repeated initialization'
                   write (ILU2,'(A)') ' WARNING from OUTDAT: ignoring repeated initialization'
           RETURN

        else if (ITOLD == 2) then
           CONTINUE

        else if (ITOLD == 3) then
!          during a previous call, one or more variables were selected, this selection is discarded after initialization
           if (TOSCR) write (*,'(A)') ' WARNING from OUTDAT: selected variables discarded'
                   write (ILU2,'(A)') ' WARNING from OUTDAT: selected variables discarded'

        else if (ITOLD >= 4) then
           if (RECOVR) call FatalERR (ModuleName,'normal outdat call impossible after recovery')
        end if
     end if

!    make unit number available also to outplt
     call AMBUSY (1, 'OUTDAT', ILU2)

!    open file check: see if units ilu2 is open. if not open them here
     INQUIRE (UNIT=ILU2, OPENED=OPEND)
     if (.not.OPEND) call FOPENG (ILU2, 'RES.DAT', 'NEW', 'SF', 0, 'DEL')

!    make unit number for temporary file local
     ILU1 = ILU2+1
!    initialize storage
     if (.not. UseMemoryStorage) then
!       see if unit ilu1 is open, if open use it, if not, open using default output file name
        INQUIRE (UNIT=ILU1, OPENED=OPEND)
        if (.not.OPEND) then
!          derive the name of binary storage file
           INQUIRE (UNIT=ILU2, NAME=String)  !  this is name of formatted output file
!          search directory separator '/', ':' or '\ '; remove directory part from name
           il = len_trim(String)
           ip = il
           do
              if (INDEX ('/\:', String(ip:ip)) > 0) then
!                directory separator found
                 BinaryFile = String(ip+1:il)
                 EXIT
              end if
              ip = ip - 1
              if (ip==0) then
!                no directory separator found
                 BinaryFile = String
                 EXIT
              end if
           end do

           call EXTENS (BinaryFile, 'bin', 1, BinaryFile)
           call FOPENG (ILU1, BinaryFile, 'NEW', 'UD', 48, 'DEL')

!          reset number of runs in BinaryFile, first record to write to and start record of set
           IRUN4 = 0
           IREC  = 1
           ISREC = 0
        else
!          temporary file is open
           if (ITOLD == -1) call FatalERR (ModuleName, 'temporary file may not be opened outside outdat')
        end if

     else
!       for ITOLD=-1 (first call or after reset call (99)), no memory should be allocated
        OPEND = Allocated(MemBin)

        if (.not.OPEND) then
!          reset number of runs in storage, first record to write to and start record of set
           IRUN4 = 0
           IREC  = 1
           ISREC = 0
        else
           if (ITOLD == -1) call FatalERR (ModuleName, 'memory storage allocation error (INTERNAL error 1)')
        end if
     end if

!    initialize routine that writes comment lines to output file
     call OUTCOM ('<INIT$$$>')

!    reset arrays with names
     AVN    = ' '
     FND    = 0
     MGIVEN = .false.

!    find out if initialization if generated by reruns
     call AMBUSY (2, 'RDFROM', IRUN2)
     if (IRUN2 == 0) then
!       run number not obtained from rdfrom or first run
        IRUN3  = IRUN1
        RUNTYP = 'N'
     else
!       run number obtained from rdfrom
        IRUN3  = IRUN2
        RUNTYP = 'R'
     end if

!    increase run number and number of runs in file
     IRUN1 = IRUN1+1
     IRUN4 = IRUN4+1

!    write total number of runs in BinaryFile file to BinaryFile file
     if (.not. UseMemoryStorage) then
        write (ILU1,REC=1) '-', IRUN4, '....................................', 0.
     else
        call Memwrite (1, '-', IRUN4, '....................................', 0.)
     end if

!    update pointer record from previous set only if it is not
!    the first initialization to the same BinaryFile
     if (IRUN4 >= 2) then
        if (.not. UseMemoryStorage) then
           write (ILU1,REC=ISREC) ' ', IREC, '....................................', 0.
        else
           call Memwrite (ISREC, ' ', IREC, '....................................', 0.)
        end if
     end if
     IREC  = IREC+1
     ISREC = IREC

!    check on length of name
     if (LEN (RN) <= VARL_M) then
!       length of name is ok
        LXN = RN
     else
!       length of name is not ok
        LXN = RN(1:VARL_M)
        if (TOSCR) write (*,'(2A,/,2A)') ' WARNING from OUTDAT: Variable ',trim(RN),' is truncated to: ', LXN
                write (ILU2,'(2A,/,2A)') ' WARNING from OUTDAT: Variable ',trim(RN),' is truncated to: ', LXN
     end if
     call UPPERC (LXN)

     if (.not. UseMemoryStorage) then
        write (ILU1, REC=IREC) ' ', -99, '....................................', 0.
     else
        call Memwrite (IREC, ' ', -99, '....................................', 0.)
     end if

     IREC  = IREC+1

     if (.not. UseMemoryStorage) then
        write (ILU1, REC=IREC) RUNTYP, IRUN3, LXN, R
     else
        call Memwrite (IREC, RUNTYP, IRUN3, LXN, R)
     end if

!    initialize other control
     IFND  = 0
     INSEL = 0
     ISAVE = 0

!    save ITASK
     ITOLD = ITASK

  else if (ITASK == 2) then
!    ---------------------
!    DUMP VARIABLE TO FILE
!    ---------------------
!    check status first
     if (ITOLD >= 3) call FatalERR (ModuleName, 'Initialization not done')

!    check on length of name
     if (LEN (RN) <= VARL_M) then
!       length of name is ok
        LN = RN
     else
!       length of name is not ok
        LN = RN(1:VARL_M)
        if (TOSCR) write (*,'(2A,/,2A)') ' WARNING from OUTDAT: Variable ',trim(RN),' is truncated to: ',LN
                write (ILU2,'(2A,/,2A)') ' WARNING from OUTDAT: Variable ',trim(RN),' is truncated to: ',LN
     end if

!    values and name are written to file ; initialize
     call UPPERC (LN) ; OK = .true. ; INAME = 0

     if (LN == LXN) then
!       variable is independent variable ; reset Variable_Found_Counters
        FND(1:IFND) = 0

     else
!       variable is not the independent variable, look in list
!       search to end of list starting at next variable relative to previous call
        SearchEnd: do I1=ISAVE+1,IFND
           if (LN == AVN(I1)) then
              INAME = I1
              EXIT SearchEnd
           end if
        end do SearchEnd

        if (INAME == 0) then
!          match not found to end of list, try beginning
           SearchBegin: do I1=1,ISAVE
              if (LN == AVN(I1)) then
                 INAME = I1
                 EXIT SearchBegin
              end if
           end do SearchBegin
        end if

        if (INAME == 0) then
!          name not found in list, add to list
           IFND = IFND+1
           if (IFND > NAMES_MXN) call FatalERR (ModuleName, 'too many variables for output')
           AVN(IFND)  = LN
           INAME      = IFND
           FND(INAME) = 1
        else
!          name found in list
           FND(INAME) = FND(INAME)+1
           if (FND(INAME) >= 2) then
!             variable supplied more than once, prevent writing to file
              OK = .false.
              if (.not.MGIVEN(INAME)) then
!                give warning only once per run
                 MGIVEN(INAME) = .true.
                 if (TOSCR) write (*,'(4A)') &
                   ' WARNING from OUTDAT: variable ',trim(LN),' supplied without new independent ', trim(LXN)
                 write (ILU2,'(4A)') &
                   ' WARNING from OUTDAT: variable ',trim(LN),' supplied without new independent ', trim(LXN)
              end if
           end if
        end if
     end if

     if (OK) then
        IREC = IREC+1

        if (.not. UseMemoryStorage) then
           write (ILU1,REC=IREC) ' ', INAME, LN, R
        else
           call Memwrite (IREC, ' ', INAME, LN, R)
        end if
     end if

!    save previous name
     ISAVE = INAME

!    save ITASK
     ITOLD = ITASK

  else if (ITASK == 3) then
!    ---------------------------------------
!    SELECTION OF OUTPUT VARIABLES FOR TABLE
!    ---------------------------------------
     if (ITOLD == 1 .or. ITOLD == 99) then
        if (TOSCR) write (*,'(A)') ' WARNING from OUTDAT: wrong task preceeds selection of output'
                write (ILU2,'(A)') ' WARNING from OUTDAT: wrong task preceeds selection of output'
        RETURN
     else if (ITOLD == 2 .or. (ITOLD >= 4 .and. ITOLD <= 16)) then
!       FIRST TIME
        INSEL = 0
     end if

     LN = RN ; call UPPERC (LN)

!    lookup name in list
     if (INSEL == 0) then
        I1 = 0
     else
        I1 = IFINDC (ASELN,INSEL,1,INSEL,LN)
     end if

     if (I1 == 0) then
!       name not found in list, add
        INSEL = INSEL+1
        if (INSEL > NAMES_MXN) call FatalERR (ModuleName, 'too many variables selected')
        ASELN(INSEL) = LN
     end if

!    save ITASK
     ITOLD = ITASK

  else if ((ITASK >= 4 .and. ITASK <= 9) .or. (ITASK >= 14 .and. ITASK <= 19)) then
!    --------------------------
!    WRITING OF OUTPUT TABLE(S)
!    --------------------------
     call OutdatConstructTable (ITASK, IUNIT, RN)

!    save ITASK
     ITOLD = ITASK

  else if (ITASK == 99) then
!    --------------------------------------------------------------
!    THIS OPTION DELETES THE TEMPORARY FILE AND RESETS ALL CONTROLS
!    --------------------------------------------------------------
     if (ITOLD == 99) call FatalERR (ModuleName, 'temporary file already deleted')

     if (.not. UseMemoryStorage) then
        close (ILU1, STATUS='DELETE')
     else
        call MemDeallocate
     end if

!    reset initialized control variables except IRUN1
     ITOLD  = -1       ! previous task
     INSEL  =  0       ! number of selected variables
     RECOVR = .false.  ! flags recovery just taking place
     FIRST8 = .true.   ! comment header for end of run output has been written
     FIRST9 = .true.   ! similar flag for TAB-delimited output

  else if (ITASK == 0) then
!    --------------------------------------------------------------------
!    ALLOW ZERO TASK, (OCCURS WITH IPFORM=0 IN RKDRIV AND EUDRIV DRIVERS)
!    --------------------------------------------------------------------
     CONTINUE
  else
     call FatalERR (ModuleName, 'WRONG ITASK')
  end if
Return
END SUBROUTINE ProcessOutputData
! =================================================================
SUBROUTINE Memwrite (RecNo, RunType, InfoID, VarName, VarValue)
! =================================================================
! writes record to memory ; same datastructure as BIN file
  IMPLICIT NONE

! formal parameters
  integer, intent(in)            :: RecNo
  integer, intent(in)            :: InfoID
  character*1, intent(in)        :: RunType
  character*(VARL_M), intent(in) :: VarName
  real, intent(in)               :: VarValue

! local variables
  integer, parameter                       :: callocAmount = 50000 ! amount by which MemBin is increased when it is too small
  type(TMemBin), allocatable, dimension(:) :: MemBinSave
  integer :: NewMembinHigh

! see if memory needs to be increased
  if (Recno > MemBinHigh) then
!   add integer times callocAmount elements
    NewMembinHigh = MemBinHigh + callocAmount * (1 + (Recno - MemBinHigh - 1) / callocAmount)

    if (allocated(MemBin)) then
!     array already allocated, increase but retain contents
      allocate(MemBinSave(NewMembinHigh))
      MemBinSave(1:MemBinHigh) = MemBin
      call move_alloc (MemBinSave, MemBin)  ! this de-allocates MemBinSave
    else
!     array not allocated
      allocate(MemBin(NewMembinHigh))
    end if
    MembinHigh = NewMembinHigh ! this is the current size of MemBin
  end if

  MemBin(Recno)%RunType  = RunType
  MemBin(Recno)%InfoID   = InfoID
  MemBin(Recno)%VarName  = VarName
  MemBin(Recno)%VarValue = VarValue
END SUBROUTINE Memwrite
! =================================================================
SUBROUTINE MemRead (RecNo, RunType, InfoID, VarName, VarValue, EOM)
! =================================================================
! reads record from memory ; same datastructure as BIN file
  IMPLICIT NONE

! formal parameters
  integer, intent(in)             :: RecNo
  integer, intent(out)            :: InfoID
  character*1, intent(out)        :: RunType
  character*(VARL_M), intent(out) :: VarName
  real, intent(inout)             :: VarValue
  integer, intent(out)            :: EOM ! ' end of memory :-)'

  if ((RecNo >= 1) .and. (RecNo <= MemBinHigh)) then
    RunType  = MemBin(Recno)%RunType
    InfoID   = MemBin(Recno)%InfoID
    VarName  = MemBin(Recno)%VarName
    VarValue = MemBin(Recno)%VarValue
    EOM      = 0
  else
    EOM      = -1
  end if
END SUBROUTINE MemRead
! =================================================================
SUBROUTINE MemDeallocate
! =================================================================
! deallocate memory storage
  IMPLICIT NONE
  integer :: iw
  Deallocate (MemBin, stat=iw)
  if (iw/=0) call FatalERR (ModuleName,'cannot de-allocate storage')
! reset current size
  MembinHigh = 0
END SUBROUTINE MemDeallocate
! =================================================================
SUBROUTINE OutdatConstructTable (ITASK, IUNIT, RN)
! =================================================================
! construct output table
  IMPLICIT NONE

! formal parameters
  INTEGER, INTENT(IN)           :: ITASK, IUNIT
  CHARACTER(LEN=*), INTENT(IN)  :: RN

! non-saved array variables
  INTEGER, DIMENSION(NAMES_MXN)      :: IAVNL    ! Array of printed lengths of AVN (last paranthesis of array variable is omitted)
  INTEGER, DIMENSION(NAMES_MXN)      :: BLK      ! Assigned block number, 0 if not assigned
  INTEGER, DIMENSION(NAMES_MXN)      :: SEQ2     ! Array with numbers that link to selected arrays
  REAL,    DIMENSION(NAMES_IB_MXN_5) :: AVV      ! Array with values of AVN
  LOGICAL, DIMENSION(NAMES_IB_MXN_5) :: FNDA     ! Value found vlag used during table construction
  INTEGER, DIMENSION(NAMES_IB_MXN_5) :: EXTRA_SP ! Extra width of columns
  INTEGER, DIMENSION(NAMES_IB_MXN_5) :: SEQ      ! SEQ(i) = element of AVN that is first, second, third in current block,

! output line length depending on output table format
  INTEGER, PARAMETER :: COL_WIDTH_MNN = 12   !  maximum width of a Y column, excluding spaces
  INTEGER, PARAMETER :: CENTRE = 8
  INTEGER, PARAMETER :: LINE_LEN_4 = COL_WIDTH_MNN + 2 + NAMES_IB_MXN_4 * (COL_WIDTH_MNN+1)
  INTEGER, PARAMETER :: LINE_LEN_5 = COL_WIDTH_MNN + 2 + NAMES_IB_MXN_5 * (COL_WIDTH_MNN+1)
  INTEGER, PARAMETER :: LINE_LEN_7 = COL_WIDTH_MNN + 2 + NAMES_IB_MXN_7 * (COL_WIDTH_MNN+1)
  INTEGER, PARAMETER :: LINE_LEN_8 = LINE_LEN_5
  INTEGER, PARAMETER :: LINE_LEN_9 = LINE_LEN_5

! output table format description
  CHARACTER(LEN=18), DIMENSION(4:9), PARAMETER :: TEXT = (/'Table output      ', 'Spreadsheet output', '2 column output   ', &
                                                           'ICASA output      ', 'End of run output ', 'Greenery output   '/)
! non-saved variables
  INTEGER :: INAME       ! name number found
  INTEGER :: LINE_LEN    ! length of output record
  REAL    :: DumVarValue ! dummy value read from memory storage
  INTEGER :: DumIOSValue ! dummy value read from memory storage
  INTEGER :: LINE_L      ! length of output line
  INTEGER :: IFND2       ! the true number of variables to be printed in a particular block
  LOGICAL :: FULL_BLOCK  ! flags block is full with columns
  LOGICAL :: EXTRA_F     ! flags whether extra spaces between columns are necessary
  INTEGER :: ILTASK      ! task modulo 10
  INTEGER :: EXTRA_SPX   ! Extra width for X columns
  INTEGER :: ILXN        ! Length of independent variable name

  CHARACTER(LEN=1)          :: RUNDUM     ! dummy value of RUNTYP
  CHARACTER(LEN=1)          :: COMMCHR    ! comment character
  CHARACTER(LEN=1)          :: CHR        ! collumn separation character
  CHARACTER(LEN=VARL_M)     :: LN         ! local variable name
  CHARACTER(LEN=VARL_M)     :: DumVarName ! dummy value read from memory storage
  CHARACTER(LEN=LINE_LEN_5) :: LINE       ! declaration of longest line, others should fit into this one
  CHARACTER(LEN=80)         :: SPACE      ! just space

! help variables
  LOGICAL :: OPEND, YFND, SELECTED
  INTEGER :: IRUN5, ISEL, IOS, ICHECK, IBLOK, IR
  INTEGER :: IEREC, IB, I1, I2, I3, I4, I5, I6
  REAL    :: LV, LVO

! function without Fortran-90 interface
  INTEGER :: IFINDC

! error checks on old task
  if (ITOLD == -1) then
!    recovery of output data in BinaryFile after crash ; recovery from a file called RES.BIN
     UseMemoryStorage = .false.

!    recovery of default output file, open RES.BIN file
     BinaryFile = 'RES.BIN'
     TOSCR  = .true.
     RECOVR = .true.
     if (IUNIT <= 10) call FatalERR (ModuleName, 'illegal unit number')
     ILU2 = IUNIT
     ILU1 = ILU2+1
     call FOPENG (ILU1, BinaryFile, 'OLD', 'UD', 48, ' ')

!    check if output file is open, make new one if it is not
     INQUIRE (UNIT=ILU2, OPENED=OPEND)
     if (.not.OPEND) call FOPENG (ILU2,'RECOVER.DAT','NEW','SF',0,'DEL')

!    find out if end record of last set is in the file
     read (ILU1, REC=1) RUNTYP, IRUN5

     if (RUNTYP == '-') then
!       write end record of last set in the file

!       go to start of last set first
        IEREC = 1
        do I1=1,IRUN5-1
           ISREC = IEREC+1
           read (ILU1, REC=ISREC) RUNDUM, IEREC
        end do

!       read rest until end_of_file
        IOS = 0
        IR = IEREC+2
        do WHILE (IOS == 0)
           read (ILU1, REC=IR, IOSTAT=IOS) RUNDUM
           if (IOS == 0) IR = IR+1
        end do

        write (ILU1, REC=IEREC+1) ' ', IR-1 , '....................................', 0.
        write (ILU1, REC=1)       '+', IRUN5, '....................................', 0.
     end if

  else if (ITOLD == ITASK) then

     if (TOSCR) write (*,'(A)') ' WARNING from OUTDAT: new output and previous output have same format'
             write (ILU2,'(A)') ' WARNING from OUTDAT: new output and previous output have same format'

  else if (ITOLD == 99) then

     if (TOSCR) write (*,'(A)') ' WARNING from OUTDAT: output table cannot be created, no .BIN file'
             write (ILU2,'(A)') ' WARNING from OUTDAT: output table cannot be created, no .BIN file'
     RETURN

  else if (ITOLD == 3.OR.ITOLD == 2.OR.ITOLD == 1) then

!    previous call was with itask=3, 2, or 1
     IEREC = IREC
     if (.not. UseMemoryStorage) then
        write (ILU1, REC=ISREC) ' ', IREC , '....................................', 0.
        read  (ILU1, REC=1)     RUNDUM,IRUN5
        write (ILU1, REC=1)     '+', IRUN5, '....................................', 0.
     else
        call Memwrite (ISREC, ' ', IREC, '....................................', 0.)
        call MemRead (1, RUNDUM, IRUN5, DumVarName, DumVarValue, DumIOSValue)
        call Memwrite (1,     '+', IRUN5, '....................................', 0.)
     end if
  end if

  ILTASK = MOD (ITASK,10)

  if (ITASK >= 14) then
!    whole BinaryFile should be converted
     I6 = 1
  else
!    only last set in BinaryFile should be converted
     I6 = IRUN5
  end if

! fill space variable
  SPACE  = ' '

! rerun loop
  RerunLoop: do I5=I6,IRUN5

  if (ITASK >= 14 .or. ITOLD == -1) then
!    read names from file into avn array of the set which is converted into an output table. this is not necessary if
!    only the last set has to be converted and this last set was just finished with an itask=2 call (name of independent
!    variable and array avn are already known !!!).

!    search for first record of set first
     IEREC = 1
     do I1=1,I5
        ISREC = IEREC+1
        if (.not. UseMemoryStorage) then
           read (ILU1, REC=ISREC) RUNDUM, IEREC
        else
           call MemRead (ISREC, RUNDUM, IEREC, DumVarName, DumVarValue, DumIOSValue)
        end if
     end do

!    read names and determine total number of names
!    read name of independent variable
!    determine length of name of independent variable

     if (.not. UseMemoryStorage) then
        read (ILU1, REC=ISREC+1) RUNTYP, IRUN3, LXN
     else
        call MemRead ((ISREC+1), RUNTYP, IRUN3, LXN, DumVarValue, DumIOSValue)
     end if

     IFND = 0

     do IR=ISREC+2,IEREC
        if (.not. UseMemoryStorage) then
           read (ILU1, REC=IR) RUNDUM, INAME, LN
        else
           call MemRead (IR, RUNDUM, INAME, LN, DumVarValue, DumIOSValue)
        end if
        IFND = MAX (INAME, IFND)
        if (INAME > 0) AVN(INAME) = LN
     end do
  end if

! length of independent variable
  ILXN = LEN_TRIM (LXN)

! assign each selected variable a block number, this is based on the length of the name of the selected variable,
! when a new name does not fit anymore on the record, a new block is generated
  FULL_BLOCK = .true.
  IBLOK      = 0
  I2         = 1 + MAX (COL_WIDTH_MNN+1,ILXN)

  if (ILTASK == 4) then
     LINE_LEN = LINE_LEN_4
  else if (ILTASK == 5) then
     LINE_LEN = LINE_LEN_5
  else if (ILTASK == 6) then
     LINE_LEN = LINE_LEN_4
  else if (ILTASK == 7) then
     LINE_LEN = LINE_LEN_7
  else if (ILTASK == 8) then
     LINE_LEN = LINE_LEN_8
  else if (ILTASK == 9) then
     LINE_LEN = LINE_LEN_9
  end if

! assign block numbers to variables
  do I1=1,IFND

     if (INSEL == 0) then
!       no variables were selected, select each one
        SELECTED = .true.
     else
        ISEL = IFINDC (ASELN,INSEL,1,INSEL,AVN(I1))
        if (ISEL > 0) then
!          variable found in list of selected variables
           SELECTED = .true.
        else
!          variable not found in list of selected variables
           SELECTED = .false.
        end if
     end if

     if (SELECTED) then

        if (FULL_BLOCK) then
!          previous block full
           IBLOK      = IBLOK+1
           FULL_BLOCK = .false.
        end if

!       variable selected
        IAVNL(I1) = LEN_TRIM (AVN(I1))
        if (IAVNL(I1) == 0) then
           if (TOSCR) write (*,'(A)') ' WARNING from OUTDAT: zero length variable name'
                   write (ILU2,'(A)') ' WARNING from OUTDAT: zero length variable name'
           IAVNL(I1) = 1
        end if

        if (ILTASK == 6) then
           FULL_BLOCK   = .true.
        else
           I3 = 1+MAX (COL_WIDTH_MNN,IAVNL(I1))
           if (I2+I3 > LINE_LEN) then
!             new variable exceeds output file record, increase block number etc.
              IBLOK = IBLOK+1
              I2    = 1+MAX (COL_WIDTH_MNN+1,ILXN)
           end if
           I2 = I2+I3
           FULL_BLOCK = .false.
        end if
        BLK(I1) = IBLOK
     else
        BLK(I1) = 0
     end if
  end do

! check number of found blocks
  if (IBLOK == 0) then
     if (TOSCR) write (*,'(A)') ' WARNING from OUTDAT: no output values given'
           write (ILU2,'(/,A)') ' WARNING from OUTDAT: no output values given'
  else if (IBLOK > 1.AND.ILTASK == 8) then
     call FatalERR (ModuleName, 'more than one block with end of run option')
  end if

! write comment header stuff
! --------------------------
  if (ILTASK == 4.OR.ILTASK == 5.OR.ILTASK == 7) then

!    display comment header for normal, tab-delimited output and ICASA output
     COMMCHR = '*'
     if (ILTASK == 7) COMMCHR = '!'

!    write line to mark start of new run
     write (ILU2,'(A,76A1)') COMMCHR,('-',I1=1,76)

!    write header to output file and possible extra comment lines
     if (IRUN3 == 0) then
        write (ILU2, '(2A)') COMMCHR, ' Output table number  :  0 (=first output table)'
     else if (RUNTYP == 'R') then
        write (ILU2, '(2A,I5)') COMMCHR, ' Output from rerun set:', IRUN3
     else if (RUNTYP == 'N') then
        write (ILU2, '(2A,I5)') COMMCHR,' Output table number  :', IRUN3
     else
        call FatalERR (ModuleName,'unknown run type')
     end if

     write (ILU2,'(3A,/,3A)') COMMCHR,' Output table format  : ',TEXT(ILTASK), &
                              COMMCHR,' ',RN

     if (ILTASK /= 7) then
!       instruct outcom to write comment lines to output file
        call OUTCOM ('<PRINT$$$>')
     end if

!    do not write outcom lines with other output formats, because they will be written
!    with an asterisk instead of with an exclamation mark

  else if (ILTASK == 8) then

!    write comment header for end of run output only once
     if (FIRST8) then
!       write line to mark start of run
        write (ILU2,'(A,76A1)') '*',('-',I1=1,76)

        write (ILU2,'(2A)') '* Output table format  : ', TEXT(ILTASK)
!       instruct outcom to write comment lines to output file
        call OUTCOM ('<PRINT$$$>')
     end if
  else if (ILTASK == 9) then
!    do not write a comment header for greenery output
     CONTINUE
  end if

  do IB=1,IBLOK
!    search stuff for block number ib
     IFND2 = 0
     do I1=1,IFND
        if (BLK(I1) == IB) then
!          variable is in current block, add to list, put pointer in long list
           IFND2      = IFND2+1
           SEQ(IFND2) = I1
           SEQ2(I1)   = IFND2
        end if
     end do

!    header with variable names is written dependent on itask
     if (ILTASK == 4.OR.ILTASK == 7) then

        write (ILU2,'(A)') ' '
        CHR = ' '

!       write name of independent variable
!       ----------------------------------
!       initialize and add leading space
        LINE    = ' '
        LINE_L  = 1

        if (ILTASK == 4) then
!          normal table format
           if (ILXN <= CENTRE+1) then
!             variable name fits left of centre point
              I3          = LINE_L + 1 + CENTRE + 1-ILXN
              I4          = LINE_L + CENTRE + 1
              LINE(I3:I4) = LXN(1:ILXN)
              LINE_L      = LINE_L + COL_WIDTH_MNN + 1
              EXTRA_SPX   = 0
           else
!             variable name extends beyond the centre point
              LINE(LINE_L + 1:LINE_L + ILXN) = LXN(1:ILXN)
              EXTRA_SPX                      = MAX (ILXN - COL_WIDTH_MNN + 1,0)
              LINE_L                         = LINE_L + COL_WIDTH_MNN + 1 + EXTRA_SPX
           end if
        else if (ILTASK == 7) then
!          ICASA format
           LINE   = '@'
           LINE_L = 1
           call ADDSTR (LINE, LINE_L, LXN)
           call ADDSTR (LINE, LINE_L, '>')
           LINE_L    = 14
           EXTRA_SPX = 0
        end if

!       set flag whether extra space is required while writing
        EXTRA_F = (EXTRA_SPX > 0)

!       write names of dependent variables centered above column
        do I1=1,IFND2
           I2 = SEQ(I1)

!          add separating space between variable names
           LINE_L = LINE_L+1

!          variable should appear in this block
           if (IAVNL(I2) <= CENTRE) then
!             variable name fits to the left of centre point
              I3 = LINE_L + 1 + CENTRE - IAVNL(I2)
              I4 = LINE_L + CENTRE

              LINE(I3:I4)  = AVN(I2)(1:IAVNL(I2))
              LINE_L       = LINE_L + COL_WIDTH_MNN
              EXTRA_SP(I1) = 0
           else
!             variable name extends beyond centre point
              LINE(LINE_L + 1:LINE_L + IAVNL(I2)) = AVN(I2)(1:IAVNL(I2))
              EXTRA_SP(I1)                        = MAX (IAVNL(I2)-COL_WIDTH_MNN,0)
              LINE_L                              = LINE_L+COL_WIDTH_MNN+EXTRA_SP(I1)
           end if

           if (EXTRA_SP(I1) > 0) EXTRA_F = .true.

        end do

!       write line to file
        write (ILU2,'(A)') LINE(1:LINE_L)
        write (ILU2,'(A)') ' '

!       add a single space to all columns
        EXTRA_SPX = EXTRA_SPX+1
        do I1=1,IFND2
           EXTRA_SP(I1) = EXTRA_SP(I1)+1
        end do

     else if (ILTASK == 5 .or. (ILTASK == 8 .and. FIRST8)) then

!       write tab delimited variable header
        CHR = CHAR (9)

!       write tabs between the names, allow no spaces !!
        LINE   = ' '
        LINE_L = 0
        call ADDSTR (LINE,LINE_L,LXN)
        do I1=1,IFND2
           call ADDSTR (LINE, LINE_L, CHR)
           call ADDSTR (LINE, LINE_L, AVN(SEQ(I1))(1:IAVNL(SEQ(I1))))
        end do

        EXTRA_F = .false.

!       write line to file
        write (ILU2,'(A)') ' '
        write (ILU2,'(A)') LINE(1:LINE_L)
        write (ILU2,'(A)') ' '

        if (ILTASK == 8.AND.FIRST8) FIRST8 = .false.

     else if (ILTASK == 6) then

!       two column output

        CHR = ' '
        write (ILU2,'(4A,/,A,/,1X,4A)') &
          '* ', LXN(1:ILXN), CHR,AVN(SEQ(1))(1:IAVNL(SEQ(1))), &
          ' 1 1 1', &
          AVN(SEQ(1))(1:IAVNL(SEQ(1))),'(',LXN(1:ILXN),')'

        EXTRA_F = .false.

      else if (ILTASK == 9.AND.FIRST9) then

!       write comma delimited variable header
        CHR = ','

!       write tabs between the names, allow no spaces !!
        LINE   = ' '
        LINE_L = 0
        call ADDSTR (LINE,LINE_L,LXN)
        do I1=1,IFND2
           call ADDSTR (LINE, LINE_L, CHR)
           call ADDSTR (LINE, LINE_L, AVN(SEQ(I1))(1:IAVNL(SEQ(I1))))
        end do

        EXTRA_F = .false.

!       write LINE TO FILE
        write (ILU2,'(A)') LINE(1:LINE_L)

        if (ILTASK == 9.AND.FIRST9) FIRST9 = .false.

     end if

!    initialize output
     YFND = .false.

!    ONLY END OF RUN VALUES, SEARCH FOR LAST OUTPUT SET
!    (THIS PART WAS TAKEN FROM CRETTP ROUTINE OF TTSELECT)

     if (ILTASK == 8) then
        LN = ' '
        IR = IEREC
        do WHILE (LN /= LXN.AND.IREC > ISREC)
           IR = IR-1
           if (.not. UseMemoryStorage) then
              read (ILU1,REC=IR) RUNDUM, INAME, LN
           else
              call MemRead (IR, RUNDUM, INAME, LN, DumVarValue, DumIOSValue)
           end if
        end do
        ISREC = IR-1
     end if

     do IR=ISREC+1,IEREC
!       read next record
        if (.not. UseMemoryStorage) then
           read (ILU1,REC=IR) RUNDUM, INAME, LN, LV
        else
           call MemRead (IR, RUNDUM, INAME, LN, LV, DumIOSValue)
        end if
!       see if variable name is the independent variable
        if (LN == LXN) then
           if (YFND) then
              if (ICHECK == IFND2) then
!                no missing variables, write line directly to file
                 if (ILTASK >= 4.AND.ILTASK <= 8) then
                    if (.not.EXTRA_F) then
!                      no extra spaces required, no missing variables
                       write (ILU2,'(1X,1P,G13.6,255(A,G12.5))') LVO,(CHR,AVV(I1),I1=1,IFND2)
                    else
!                      extra spaces required
                       write (ILU2,'(A,1P,G13.6,255(A,G12.5))') &
                          SPACE(1:EXTRA_SPX), LVO, (SPACE(1:EXTRA_SP(I1)),AVV(I1),I1=1,IFND2)
                    end if
                 else if (ILTASK == 9) then
                    LINE   = ' '
                    LINE_L = 0
                    call ADDINT (LINE, LINE_L, INT (LVO))
                    call ADDSTR (LINE, LINE_L, ',')
                    call ADDINT (LINE, LINE_L, INT (AVV(1)))
                    call ADDSTR (LINE, LINE_L, ',')
                    call ADDINT (LINE, LINE_L, INT (AVV(2)))
                    call ADDSTR (LINE, LINE_L, ',')
                    call ADDINT (LINE, LINE_L, INT (AVV(3)))
                    call ADDSTR (LINE, LINE_L, ',')
                    do I1=4,IFND2
                       call ADDREA (LINE, LINE_L, AVV(I1), '1P,G12.5')
                       call ADDSTR (LINE, LINE_L, ',')
                    end do
                    write (ILU2,'(A)') LINE(1:LINE_L)
                 else
                    call FatalERR (ModuleName, 'INTERNAL ERROR 2')
                 end if
              else
!                one or more missing values found write to string line first
                 if (ILTASK >= 4.AND.ILTASK <= 8) then
                    if (.not.EXTRA_F) then
!                      no extra spaces required between columns
                       LINE   = ' '
                       LINE_L = 1
                       call ADDREF (LINE, LINE_L,LVO, '1P,G13.6')
                       do I1=1,IFND2
                          call ADDSTF (LINE,LINE_L,CHR)
                          if (FNDA(I1)) then
!                            value was found
                             call ADDREF (LINE, LINE_L, AVV(I1), '1P,G12.5')
                          else
!                            value was not found
                             LINE(LINE_L + CENTRE:LINE_L + CENTRE) = '-'
                             LINE_L                                = LINE_L + COL_WIDTH_MNN
                          end if
                       end do
                    else
!                      extra spaces required between columns
                       LINE   = ' '
                       LINE_L = EXTRA_SPX
                       call ADDREF (LINE, LINE_L, LVO, '1P,G13.6')
                       do I1=1,IFND2
                          LINE_L = LINE_L+EXTRA_SP(I1)
                          if (FNDA(I1)) then
!                            value was found
                             call ADDREF (LINE, LINE_L, AVV(I1), '1P,G12.5')
                          else
!                            value was not found
                             LINE(LINE_L + CENTRE:LINE_L + CENTRE) = '-'
                             LINE_L                                = LINE_L + COL_WIDTH_MNN
                          end if
                       end do
                    end if
                 else if (ILTASK == 9) then
                    LINE   = ' '
                    LINE_L = 0
                    call ADDINT (LINE, LINE_L, INT (LVO))
                    call ADDSTR (LINE, LINE_L, ',')
                    call ADDINT (LINE, LINE_L, INT (AVV(1)))
                    call ADDSTR (LINE, LINE_L, ',')
                    call ADDINT (LINE, LINE_L, INT (AVV(2)))
                    call ADDSTR (LINE, LINE_L, ',')
                    call ADDINT (LINE, LINE_L, INT (AVV(3)))
                    do I1=4,IFND2
                       call ADDSTR (LINE, LINE_L, ',')
                       if (FNDA(I1)) then
!                         value was found
                          call ADDREA (LINE, LINE_L, AVV(I1), '1P,G12.5')
                       else
!                         value was not found
                          call ADDSTR (LINE, LINE_L, '-')
                       end if
                    end do
                 else
                    call FatalERR (ModuleName, 'INTERNAL ERROR 3')
                 end if

                 write (ILU2,'(A)') LINE(1:LINE_L)
              end if
!             initialize search for 'y' values
              YFND = .false.
           end if

!          reinitialize things
           do I1=1,IFND2
              FNDA(I1) = .false.
              AVV(I1)  = -99.
           end do
           ICHECK = 0
           LVO = LV
        else
!          record contains 'y' ; check names i2...i3
           if (BLK(INAME) == IB) then
              AVV(SEQ2(INAME))  = LV
              YFND              = .true.
              ICHECK            = ICHECK+1
              FNDA(SEQ2(INAME)) = .true.
           end if
        end if
     end do

!    write last line at E_O_F
!    unfortunately this section must be identical to the previous one
     if (YFND) then
        if (ICHECK == IFND2) then
!          no missing variables, write line directly to file
           if (ILTASK >= 4.AND.ILTASK <= 8) then
              if (.not.EXTRA_F) then
!                no extra spaces required, no missing variables
                 write (ILU2,'(1X,1P,G13.6,255(A,G12.5))') LVO, (CHR,AVV(I1),I1=1,IFND2)
              else
!                extra spaces required
                 write (ILU2,'(A,1P,G13.6,255(A,G12.5))') &
                        SPACE(1:EXTRA_SPX), LVO, (SPACE(1:EXTRA_SP(I1)),AVV(I1),I1=1,IFND2)
              end if
           else if (ILTASK == 9) then
              LINE   = ' '
              LINE_L = 0
              call ADDINT (LINE, LINE_L, INT (LVO))
              call ADDSTR (LINE, LINE_L, ',')
              call ADDINT (LINE, LINE_L, INT (AVV(1)))
              call ADDSTR (LINE, LINE_L, ',')
              call ADDINT (LINE, LINE_L, INT (AVV(2)))
              call ADDSTR (LINE, LINE_L, ',')
              call ADDINT (LINE, LINE_L, INT (AVV(3)))
              call ADDSTR (LINE, LINE_L, ',')
              do I1=4,IFND2
                 call ADDREA (LINE, LINE_L, AVV(I1), '1P,G12.5')
                 call ADDSTR (LINE, LINE_L, ',')
              end do
              write (ILU2,'(A)') LINE(1:LINE_L)
           else
              call FatalERR (ModuleName, 'INTERNAL ERROR 4')
           end if
        else
!          one or more missing values found write to string line first
           if (ILTASK >= 4.AND.ILTASK <= 8) then
              if (.not.EXTRA_F) then
!                no extra spaces required between columns
                 LINE   = ' '
                 LINE_L = 1
                 call ADDREF (LINE, LINE_L, LVO, '1P,G13.6')
                 do I1=1,IFND2
                    call ADDSTF (LINE, LINE_L,CHR)
                    if (FNDA(I1)) then
!                      value was found
                       call ADDREF (LINE, LINE_L, AVV(I1), '1P,G12.5')
                    else
!                      value was not found
                       LINE(LINE_L + CENTRE:LINE_L + CENTRE) = '-'
                       LINE_L                                = LINE_L + COL_WIDTH_MNN
                    end if
                 end do
              else
!                extra spaces required between columns
                 LINE   = ' '
                 LINE_L = EXTRA_SPX
                 call ADDREF (LINE, LINE_L, LVO, '1P,G13.6')
                 do I1=1,IFND2
                    LINE_L = LINE_L + EXTRA_SP(I1)
                    if (FNDA(I1)) then
!                      value was found
                       call ADDREF (LINE, LINE_L, AVV(I1), '1P,G12.5')
                    else
!                      value was not found
                       LINE(LINE_L + CENTRE:LINE_L + CENTRE) = '-'
                       LINE_L                                = LINE_L + COL_WIDTH_MNN
                    end if
                 end do
              end if
           else if (ILTASK == 9) then
              LINE   = ' '
              LINE_L = 0
              call ADDINT (LINE, LINE_L, INT (LVO))
              call ADDSTR (LINE, LINE_L, ',')
              call ADDINT (LINE, LINE_L, INT (AVV(1)))
              call ADDSTR (LINE, LINE_L, ',')
              call ADDINT (LINE, LINE_L, INT (AVV(2)))
              call ADDSTR (LINE, LINE_L, ',')
              call ADDINT (LINE, LINE_L, INT (AVV(3)))
              do I1=4,IFND2
                 call ADDSTR (LINE,LINE_L,',')
                 if (FNDA(I1)) then
!                   value was found
                    call ADDREA (LINE, LINE_L, AVV(I1), '1P,G12.5')
                 else
!                   value was not found
                    call ADDSTR (LINE, LINE_L, '-')
                 end if
              end do
           else
              call FatalERR (ModuleName, 'INTERNAL ERROR 5')
           end if

           write (ILU2,'(A)') LINE(1:LINE_L)
        end if
     end if
  end do

! add some blank lines if table or spreadsheet output was chosen
  if (ILTASK == 4 .or. ILTASK == 5) write (ILU2,'(/,/,/,1X)')

  end do RerunLoop
Return
END SUBROUTINE OutdatConstructTable
END MODULE Module_OUTDAT

SUBROUTINE OUTDAT (ITASK, IUNIT, RN, R)
! Functionality of old OUTDAT, making use of Module_OUTDAT
  USE Module_OUTDAT, ONLY: ProcessOutputData

  IMPLICIT NONE
! formal parameters
  INTEGER, INTENT(IN)           :: ITASK, IUNIT
  REAL, INTENT(IN)              :: R
  CHARACTER(LEN=*), INTENT(IN)  :: RN

  call ProcessOutputData (ITASK, IUNIT, RN, R)
Return
END SUBROUTINE OUTDAT

SUBROUTINE OUTDAT_Memory_Storage (InMemory)
! Sets storage mode of Module_OUTDAT
  USE Module_OUTDAT, ONLY: OutDat_UseMemory

  IMPLICIT NONE
! formal parameter
  LOGICAL, INTENT(IN) :: InMemory

  OutDat_UseMemory = InMemory
Return
END SUBROUTINE OUTDAT_Memory_Storage
