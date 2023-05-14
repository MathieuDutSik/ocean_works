

PROGRAM INTERPOLATE_SPECTRUM

  character(len=200) :: INP_FNAME
  character(len=200) :: FileInput
  character(len=200) :: PrefixOutput
  character(len=200) :: FileInterpolatePoint
  integer Npoint_ouput
  REAL(8), allocatable :: ListLonOutput(:), ListLatOutput(:)
  integer nbArg

  NAMELIST /GRID/ NUMSIG, NUMDIR, FRLOW, FRHIGH
  NAMELIST /PROC/ FileInput, PrefixOutput, FileInterpolatePoint, Npoint_ouput, ListLonOutput, ListLatOutput
  nbArg=command_argument_count()
  IF (nbArg .ne. 1) THEN
     CALL WWM_ABORT('Number of argument is 0 or 1')
  ENDIF
  CALL GET_COMMAND_ARGUMENT(1, INP_FNAME)
  CALL TEST_FILE_EXIST_DIE(INP_FNAME)
  OPEN(INP_FHNDL, FILE = INP_FNAME)
  READ(INP_FHNDL, NML = GRID)
  READ(INP_FHNDL, NML = PROC)
  CLOSE(INP_FHNDL)

  WAV % FNAME = FileInput

  ! Reading the output points
  NODES % FNAME = FileInterpolatePoint
  OPEN(NODES%FHNDL, FILE = WAV%FNAME, STATUS='OLD')
  READ(NODES%FHDNL,*) IWBMNP
  allocate(XP(IWBMNP), YP(IWBMNP), stat=istat)
  DO I=1,IWBMNP
     READ(NODES%FHDNL,*) XP(I), YP(I)
  END DO
  CLOSE(NODES%FHNDL)

  !  Reading the data

  CALL INIT_GRIB_WAM_BOUNDARY
END PROGRAM INTERPOLATE_SPECTRUM
