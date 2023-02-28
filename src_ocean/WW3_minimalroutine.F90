!**********************************************************************
!*   This is the minimal code for writing a WAVEWATCH III             *
!*   Header.                                                          *
!**********************************************************************
        SUBROUTINE WRITE_WAVEWATCH_HEADER(ChoiceFile, NX, NY, GTYPE, IsFormatted)
        IMPLICIT NONE
        INTEGER, intent(in) :: ChoiceFile
        INTEGER, intent(in) :: NX, NY, GTYPE, IsFormatted
        !
        CHARACTER(LEN=20) :: FileName
        CHARACTER(LEN=3) :: IDFLD
        INTEGER :: FILLER(3)
        CHARACTER(LEN=13) :: IDSTR = 'WAVEWATCH III'
        INTEGER :: TIDEFLAG = 0
        INTEGER :: TheOut = 10
        IF (ChoiceFile .eq. 1) THEN
          FileName = 'wind.ww3'
          IDFLD = 'WND'
        END IF
        IF (ChoiceFile .eq. 2) THEN
          FileName = 'current.ww3'
          IDFLD = 'CUR'
        END IF
        IF (ChoiceFile .eq. 3) THEN
          FileName = 'level.ww3'
          IDFLD = 'LEV'
        END IF
        Print *, 'WRITE_WAVEWATCH_HEADER, step 0'
        FILLER(:)=0
        Print *, 'WRITE_WAVEWATCH_HEADER, step 1'
        Print *, 'ChoiceFile=', ChoiceFile
        Print *, 'NX=', NX
        Print *, 'NY=', NY
        Print *, 'GTYPE=', GTYPE
        Print *, 'FileName=', FileName
        Print *, 'IDFLD   =', IDFLD

        IF (IsFormatted .eq. 1) THEN
           OPEN(TheOut, FILE=TRIM(FileName), FORM='FORMATTED', status='replace', action='write')
           Print *, 'WRITE_WAVEWATCH_HEADER, step 2'
           WRITE (TheOut) IDSTR, IDFLD, NX, NY, GTYPE, FILLER(1:2), TIDEFLAG
           Print *, 'WRITE_WAVEWATCH_HEADER, step 3'
           CLOSE(TheOut)
           Print *, 'WRITE_WAVEWATCH_HEADER, step 4'
        ELSE
           OPEN(TheOut, FILE=TRIM(FileName), FORM='UNFORMATTED', status='replace', action='write')
           Print *, 'WRITE_WAVEWATCH_HEADER, step 5'
           WRITE (TheOut) IDSTR, IDFLD, NX, NY, GTYPE, FILLER(1:2), TIDEFLAG
           Print *, 'WRITE_WAVEWATCH_HEADER, step 6'
           CLOSE(TheOut)
           Print *, 'WRITE_WAVEWATCH_HEADER, step 7'
        END IF
        RETURN

900     FORMAT (1X,A13,1X,A3,6I12)
        END SUBROUTINE WRITE_WAVEWATCH_HEADER
!**********************************************************************
!*   This is the minimal code for writing a WAVEWATCH III             *
!*   Entry. Suitable for wind or currents                             *
!**********************************************************************
        SUBROUTINE WRITE_WAVEWATCH_ENTRY_TWO_FIELD(FileName, TFN, NX, NY, U, V, IsFormatted)
        IMPLICIT NONE
        CHARACTER(LEN=20), intent(in) :: FileName
        INTEGER, intent(in) :: TFN(2)
        INTEGER, intent(in) :: NX, NY, IsFormatted
        REAL, intent(in) :: U(NX*NY), V(NX*NY)
        REAL Uprov(NX,NY), Vprov(NX,NY)
        INTEGER IX,IY
        INTEGER :: TheOut = 10
        Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 1'
        DO IX=1,NX
           DO IY=1,NY
              Uprov(IX,IY) = U(IX + NX*(IY-1))
              Vprov(IX,IY) = V(IX + NX*(IY-1))
           END DO
        END DO
        IF (IsFormatted .eq. 1) THEN
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 2'
           OPEN(TheOut, FILE=TRIM(FileName),FORM='FORMATTED',status='old',position='append',action='write')
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 3'
           WRITE(TheOut) TFN
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 4'
           WRITE(TheOut) ((Vprov(IX,IY),IX=1,NX),IY=1,NY)
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 5'
           WRITE(TheOut) ((Uprov(IX,IY),IX=1,NX),IY=1,NY)
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 6'
           CLOSE(TheOut)
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 7'
        ELSE
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 8'
           OPEN(TheOut, FILE=TRIM(FileName),FORM='UNFORMATTED',status='old',position='append',action='write')
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 9'
           WRITE(TheOut) TFN
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 10'
           WRITE(TheOut) ((Vprov(IX,IY),IX=1,NX),IY=1,NY)
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 11'
           WRITE(TheOut) ((Uprov(IX,IY),IX=1,NX),IY=1,NY)
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 12'
           CLOSE(TheOut)
           Print *, 'WRITE_WAVEWATCH_ENTRY_TWO_FIELD, step 13'
        END IF
        END SUBROUTINE
!**********************************************************************
!*   This is the minimal code for writing a WAVEWATCH III             *
!*   Entry. Suitable for water levels                                 *
!**********************************************************************
        SUBROUTINE WRITE_WAVEWATCH_ENTRY_ONE_FIELD(FileName, TFN, NX, NY, F, IsFormatted)
        IMPLICIT NONE
        CHARACTER(LEN=20), intent(in) :: FileName
        INTEGER, intent(in) :: TFN(2)
        INTEGER, intent(in) :: NX, NY, IsFormatted
        REAL, intent(in) :: F(NX*NY)
        REAL Fprov(NX,NY)
        INTEGER IX,IY
        INTEGER :: TheOut = 10
        DO IX=1,NX
           DO IY=1,NY
              Fprov(IX,IY) = F(IX + NX*(IY-1))
           END DO
        END DO
        IF (IsFormatted .eq. 1) THEN
           OPEN(TheOut, FILE=TRIM(FileName),FORM='FORMATTED',status='old',position='append',action='write')
           WRITE(TheOut) TFN
           WRITE(TheOut) ((Fprov(IX,IY),IX=1,NX),IY=1,NY)
           CLOSE(TheOut)
        ELSE
           OPEN(TheOut, FILE=TRIM(FileName),FORM='UNFORMATTED',status='old',position='append',action='write')
           WRITE(TheOut) TFN
           WRITE(TheOut) ((Fprov(IX,IY),IX=1,NX),IY=1,NY)
           CLOSE(TheOut)
        END IF
        END SUBROUTINE
