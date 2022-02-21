SUBROUTINE OUTAR2 (NAME, ARRAY, LDEC, UDEC, I1, I2)
! array elements to outdat
  IMPLICIT NONE

! formal
  CHARACTER(LEN=*)           :: NAME
  INTEGER                    :: LDEC, UDEC, I1, I2
  REAL, DIMENSION(LDEC:UDEC) :: ARRAY

! local variables
  CHARACTER(LEN=5), DIMENSION(0:999), SAVE :: IND

  LOGICAL, SAVE           :: INIT     = .false. ! generate index strings
  INTEGER, PARAMETER      :: VARL_M   = 36      ! DO NOT CHANGE
  INTEGER, PARAMETER      :: IndexLen =  5      ! maximum leng of index string, including parentheses
  CHARACTER(LEN=IndexLen) :: IndexString
  CHARACTER(LEN=VARL_M)   :: TEMP
  INTEGER                 :: I3, IL, STEP

! initialize character array of subscripts
  if (.not.INIT) then
     do I3=0,999
        if (I3.GE.0 .and. I3.LE.9  ) write (IND(I3),'(A,I1,A)') '(',I3,')'
        if (I3.GT.9 .and. I3.LE.99 ) write (IND(I3),'(A,I2,A)') '(',I3,')'
        if (I3.GT.99.and. I3.LE.999) write (IND(I3),'(A,I3,A)') '(',I3,')'
     end do
     INIT = .true.
  end if

  if (I1 < LDEC .or. I1 > UDEC .or. I2 < LDEC .or. I2 > UDEC) then
     call FatalERR ('OUTAR2','Output range of '//trim(NAME)//' outside declared bounds')
  end if

  STEP = 1 ; if (I2 < I1) STEP = -1

  IL = LEN_TRIM (NAME)
  if (IL > VARL_M - IndexLen) call FatalERR ('OUTAR2','Variable name '//trim(NAME)//' too long')

! dump names and values from start element to finish element
  do I3=I1,I2,STEP
!    send element to OUTDAT
     if (I3 >= 0 .and. I3 <= 999) then
!       index was created in local character array, concatenation with name
        IndexString = IND(I3)
     else
!       index string was not created in local character array,
        if (I3 >= -99 .and. I3 <= -10) then
           write (IndexString,'(A,I3,A)') '(',I3,')'
        else if (I3 >= -9 .and. I3 <= -1) then
           write (IndexString,'(A,I2,A)') '(',I3,')'
        ELSE
           call FATALERR ('OUTAR2', 'subscript out of range')
        end if
     end if

!    send array element to OUTDAT
     TEMP = NAME(1:IL)//trim(IndexString)
     call OUTDAT (2, 0, TEMP, ARRAY(I3))
  end do
Return
END SUBROUTINE OUTAR2
