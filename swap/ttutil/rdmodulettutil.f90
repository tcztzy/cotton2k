MODULE RDmoduleTTutil
!  module RD which can be extended until containing the entire RD set
   USE ttutilPrefs, ONLY: TTutilNameLEN

!  default
   PRIVATE

!  administration
   CHARACTER (LEN=*), PARAMETER :: ModuleName    = 'RDmoduleTTutil'
   CHARACTER (LEN=*), PARAMETER :: ModuleAuthor  = 'Kees Rappoldt'
   CHARACTER (LEN=*), PARAMETER :: ModuleVersion = '0.01'

!  Data file inquiry
   TYPE, PUBLIC :: RDdataVarInfoType
      CHARACTER(LEN=TTutilNameLEN) :: Name
      CHARACTER(LEN=1)             :: Type
      LOGICAL                      :: Array
      INTEGER                      :: Size
   END TYPE RDdataVarInfoType
   INTEGER, PUBLIC                                            :: RDdataVariableCnt
   TYPE(RDdataVarInfoType), DIMENSION(:), ALLOCATABLE, PUBLIC :: RDdataVariable

!  declare public subroutines and functions
   PUBLIC RDinqDataFile
CONTAINS
!  =================================================================
SUBROUTINE RDinqDataFile
!  =================================================================
!  writes variable info into Fortran-90 allocatable info records
   IMPLICIT NONE

!  administration
   CHARACTER (LEN=*), PARAMETER :: SubProgram = 'RDinqDataFile'

!  read rddata common blocks
   INCLUDE 'rdmachin.inc'
   INCLUDE 'rddata.inc'

!  local
   INTEGER :: i, iw

!  make avaiable the variable count
   RDdataVariableCnt = INFDAT

!  allocate variable inquiry records
   select case (RDdataVariableCnt)
   case (0)
      call Warning (SubProgram,'no RDdata file active')
      iw = 0
      if (allocated(RDdataVariable))  deallocate (RDdataVariable,  stat=iw)
      if (iw/=0) call FatalERR (SubProgram,'cannot de-allocate RDdataVariable')

   case default
      iw = 0
      if (allocated(RDdataVariable))  deallocate (RDdataVariable,  stat=iw)
      if (iw/=0) call FatalERR (SubProgram,'cannot de-allocate RDdataVariable')
      allocate (RDdataVariable(RDdataVariableCnt),  stat=iw)
      if (iw/=0) call FatalERR (SubProgram,'cannot allocate RDdataVariableCnt')

!     fill variable inquiry records
      do i=1, RDdataVariableCnt
         RDdataVariable(i)%Name  = DATLIS(i)
         RDdataVariable(i)%Type  = DATTYP(i)
         RDdataVariable(i)%Array = DATARR(i) == 'a'
         RDdataVariable(i)%Size  = IARDAT(i)
      end do
   end select
Return
END SUBROUTINE RDinqDataFile
END MODULE RDmoduleTTutil
