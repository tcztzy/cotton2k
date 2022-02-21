FUNCTION NewLine ()
! returns machine dependent newline code
  IMPLICIT NONE
! function name
  CHARACTER(LEN=2) NewLine

! standard fortran intrinsic function call
  NewLine = new_line (' ')
Return
END FUNCTION NewLine
