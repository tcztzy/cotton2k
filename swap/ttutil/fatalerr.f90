SUBROUTINE FatalERR (MODULE,MESSAG)
  USE ttutilPrefs, ONLY: FatalErrorMode, FatalErrorDefault, FatalErrorInternal, FatalErrorERRfile, FatalErrorExceptFile, &
                         WriteExternalMessage, FatalErrorFileName, ExceptionUnit, FatalErrorScreenTextSTOP

  IMPLICIT NONE
! formal
  CHARACTER(LEN=*) :: MODULE, MESSAG

! local variables
  INTEGER          :: il1,il2,UNLOG, ErrUnit, UnitRD
  LOGICAL          :: TOSCR, TOLOG
  CHARACTER(LEN=1) :: DUMMY

! function
  INTEGER :: Getun2

  il1 = LEN_TRIM (MODULE)
  il2 = LEN_TRIM (MESSAG)

! desired output type
  call MESSINQ (TOSCR, TOLOG, UNLOG)

! fudge construction to fool the compiler about the return statement
  if (FatalErrorMode == FatalErrorDefault .or. FatalErrorMode == FatalErrorInternal) then
!    special message ?
     if (FatalErrorMode == FatalErrorInternal) call WriteExternalMessage

     if (il1==0 .and. il2==0) then
        if (TOSCR) write (*,'(A)')     ' Fatal execution error, press <Enter>'
        if (TOLOG) write (UNLOG,'(A)') ' Fatal execution error'
     else if (il1>0 .and. il2==0) then
        if (TOSCR) write (*,'(3A)')     ' Fatal execution error in ',MODULE(1:il1),', press <Enter>'
        if (TOLOG) write (UNLOG,'(2A)') ' Fatal execution error in ',MODULE(1:il1)
     else
        if (TOSCR) write (*,'(4A,/,A)') ' ERROR in ',MODULE(1:il1),': ',MESSAG(1:il2),' Press <Enter>'
        if (TOLOG) write (UNLOG,'(4A)') ' ERROR in ',MODULE(1:il1),': ',MESSAG(1:il2)
     end if

     if (len_trim(FatalErrorScreenTextSTOP) > 0) then
        write (*,'(a)') trim(FatalErrorScreenTextSTOP)
     end if
!    ================================================
!    replace by any other EXIT procedure if necessary
!    ================================================
!    delete temporaries
     UnitRD = Getun2 (10,99,1)
     call RDDTMP (UnitRD)
     if (TOSCR) READ (*,'(A)') DUMMY
     STOP

  else if (FatalErrorMode == FatalErrorERRfile) then
!    special error file to be created
     if (il1==0 .and. il2==0) then
        if (TOSCR) write (*,'(A)')     ' Fatal execution error'
        if (TOLOG) write (UNLOG,'(A)') ' Fatal execution error'
     else if (il1>0 .and. il2==0) then
        if (TOSCR) write (*,'(2A)')     ' Fatal execution error in ',MODULE(1:il1)
        if (TOLOG) write (UNLOG,'(2A)') ' Fatal execution error in ',MODULE(1:il1)
     else
        if (TOSCR) write (*,'(4A)')     ' ERROR in ',MODULE(1:il1),': ',MESSAG(1:il2)
        if (TOLOG) write (UNLOG,'(4A)') ' ERROR in ',MODULE(1:il1),': ',MESSAG(1:il2)
     end if

     ErrUnit = Getun2 (10,99,1)
     call FOPENS (ErrUnit, FatalErrorFileName, 'NEW','DEL')
     if (il1==0 .and. il2==0) then
        write (ErrUnit,'(A)') ' Fatal execution error'
     else if (il1>0 .and. il2==0) then
        write (ErrUnit,'(2A)') ' Fatal execution error in ',MODULE(1:il1)
     else
        write (ErrUnit,'(4A)') ' ERROR in ',MODULE(1:il1),': ',MESSAG(1:il2)
     end if
     close (ErrUnit)

     if (len_trim(FatalErrorScreenTextSTOP) > 0) then
        write (*,'(a)') trim(FatalErrorScreenTextSTOP)
     end if
!    ================================================
!    replace by any other EXIT procedure if necessary
!    ================================================
!    delete temporaries
     UnitRD = Getun2 (10,99,1)
     call RDDTMP (UnitRD)
     STOP

  else if (FatalErrorMode == FatalErrorExceptFile) then
!    error is written to open exceptions file
     if (il1==0 .and. il2==0) then
        if (TOSCR) write (*,'(A)')     ' Fatal execution error'
        if (TOLOG) write (UNLOG,'(A)') ' Fatal execution error'
     else if (il1>0 .and. il2==0) then
        if (TOSCR) write (*,'(2A)')     ' Fatal execution error in ',MODULE(1:il1)
        if (TOLOG) write (UNLOG,'(2A)') ' Fatal execution error in ',MODULE(1:il1)
     else
        if (TOSCR) write (*,'(4A)')     ' ERROR in ',MODULE(1:il1),': ',MESSAG(1:il2)
        if (TOLOG) write (UNLOG,'(4A)') ' ERROR in ',MODULE(1:il1),': ',MESSAG(1:il2)
     end if

     if (il1==0 .and. il2==0) then
        write (ExceptionUnit,'(A)') ' Fatal execution error'
     else if (il1>0 .and. il2==0) then
        write (ExceptionUnit,'(2A)') ' Fatal execution error in ',MODULE(1:il1)
     else
        write (ExceptionUnit,'(4A)') ' ERROR in ',MODULE(1:il1),': ',MESSAG(1:il2)
     end if
!    close exceptions file
     close (ExceptionUnit)

     if (len_trim(FatalErrorScreenTextSTOP) > 0) then
        write (*,'(a)') trim(FatalErrorScreenTextSTOP)
     end if
!    ================================================
!    replace by any other EXIT procedure if necessary
!    ================================================
!    delete temporaries
     UnitRD = Getun2 (10,99,1)
     call RDDTMP (UnitRD)
     STOP
  end if
Return
END SUBROUTINE FatalERR
