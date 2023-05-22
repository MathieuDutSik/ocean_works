

PROGRAM INTERPOLATE_SPECTRUM
  USE DATAPOOL
  IMPLICIT NONE
  character(len=140) :: FileInput
  character(len=140) :: FileInterpolatePoint
  INTEGER I
  integer nbArg
  NAMELIST /GRID/ NUMSIG, NUMDIR, FRLOW, FRHIGH, EXTRAPOLATION_ALLOWED_BOUC
  NAMELIST /PROC/ FileInput, PrefixOutput, FileInterpolatePoint
  EXTRAPOLATION_ALLOWED_BOUC = .FALSE.
  nbArg=command_argument_count()
  IF (nbArg .ne. 1) THEN
     CALL WWM_ABORT('Number of argument is 0 or 1')
  ENDIF
  CALL GET_COMMAND_ARGUMENT(1, INP%FNAME)
  CALL TEST_FILE_EXIST_DIE("Input fileName", INP%FNAME)
  OPEN(INP%FHNDL, FILE = INP%FNAME)
  READ(INP%FHNDL, NML = GRID)
  READ(INP%FHNDL, NML = PROC)
  CLOSE(INP%FHNDL)

  WAV % FNAME = FileInput

  ! Reading the output points
  CALL TEST_FILE_EXIST_DIE("FileInterpolatePoint", FileInterpolatePoint)

  NODES % FNAME = FileInterpolatePoint
  OPEN(NODES%FHNDL, FILE = NODES%FNAME, STATUS='OLD')
  READ(NODES%FHNDL,*) IWBMNP
  allocate(XP(IWBMNP), YP(IWBMNP), stat=istat)
  DO I=1,IWBMNP
     READ(NODES%FHNDL,*) XP(I), YP(I)
  END DO
  CLOSE(NODES%FHNDL)
  Print *, "IWBMNP=", IWBMNP

  !  Reading the data

  CALL INIT_GRIB_WAM_BOUNDARY
  CALL OUTPUT_SPECTRUM
END PROGRAM INTERPOLATE_SPECTRUM