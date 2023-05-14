SUBROUTINE WWM_ABORT(string)
  IMPLICIT NONE
  character(*), intent(in) :: string
  Print *, 'We have to abort. Reason:'
  Print *, TRIM(string)
  STOP 'WWM_ABORT'
END SUBROUTINE WWM_ABORT

SUBROUTINE TEST_FILE_EXIST_DIE(string1, string2)
  CHARACTER(LEN=*), intent(in) :: string1
  CHARACTER(LEN=*), intent(in) :: string2
  CHARACTER(LEN=512) :: ErrMsg
  LOGICAL :: LFLIVE
  INQUIRE( FILE = TRIM(string2), EXIST = LFLIVE )
  IF ( .NOT. LFLIVE ) THEN
     Print *, 'TRIM(string1)=', TRIM(string1)
     Print *, 'TRIM(string2)=', TRIM(string2)
     WRITE(ErrMsg,10) TRIM(string1), TRIM(string2)
10   FORMAT(a, ' ', a)
     CALL WWM_ABORT(TRIM(ErrMsg))
  END IF
END SUBROUTINE TEST_FILE_EXIST_DIE
