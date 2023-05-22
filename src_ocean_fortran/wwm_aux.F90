!**********************************************************************
!*                                                                    *
!**********************************************************************
SUBROUTINE WWM_ABORT(string)
  IMPLICIT NONE
  character(*), intent(in) :: string
  Print *, 'We have to abort. Reason:'
  Print *, TRIM(string)
  STOP 'WWM_ABORT'
END SUBROUTINE WWM_ABORT
!**********************************************************************
!*                                                                    *
!**********************************************************************
SUBROUTINE TEST_FILE_EXIST_DIE(string1, string2)
  CHARACTER(LEN=*), intent(in) :: string1
  CHARACTER(LEN=*), intent(in) :: string2
  CHARACTER(LEN=512) :: ErrMsg
  LOGICAL :: LFLIVE
  INQUIRE( FILE = TRIM(string2), EXIST = LFLIVE )
  IF ( .NOT. LFLIVE ) THEN
     Print *, "Error in TEST_FILE_EXIST_DIE"
     Print *, 'Context  : TRIM(string1)=', TRIM(string1)
     Print *, 'FileName : TRIM(string2)=', TRIM(string2)
     WRITE(ErrMsg,10) TRIM(string1), TRIM(string2)
10   FORMAT(a, ' ', a)
     CALL WWM_ABORT(TRIM(ErrMsg))
  END IF
END SUBROUTINE TEST_FILE_EXIST_DIE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT_COEF(X,Y,XP,YP,WI)
      USE DATAPOOL, ONLY : RKIND
      IMPLICIT NONE
      REAL(rkind),    INTENT(IN)  :: X(3), Y(3)
      REAL(rkind),    INTENT(IN)  :: XP, YP
      REAL(rkind),  INTENT(OUT) :: WI(3)

      REAL(rkind) :: y1,y2,y3,x1,x2,x3
      REAL(rkind) :: n1, d1, n2, d2, n3, d3
      x1 = X(1)
      x2 = X(2)
      x3 = X(3)
      y1 = Y(1)
      y2 = Y(2)
      y3 = Y(3)
      n1=(XP-x2)*(y3-y2) - (YP-y2)*(x3-x2)
      d1=(x1-x2)*(y3-y2) - (y1-y2)*(x3-x2)
      n2=(XP-x1)*(y3-y1) - (YP-y1)*(x3-x1)
      d2=(x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
      n3=(XP-x1)*(y2-y1) - (YP-y1)*(x2-x1)
      d3=(x3-x1)*(y2-y1) - (y3-y1)*(x2-x1)
      Wi(1)=n1/d1
      Wi(2)=n2/d2
      Wi(3)=n3/d3
      END SUBROUTINE INTELEMENT_COEF
