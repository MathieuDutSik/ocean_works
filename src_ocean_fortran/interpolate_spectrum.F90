

PROGRAM INTERPOLATE_SPECTRUM
  USE DATAPOOL
  IMPLICIT NONE
  character(len=400) :: FILE_NAME
  character(len=140) :: FileInput
  character(len=140) :: FileInterpolatePoint
  character(len=140) :: PrefixBounc
  INTEGER I
  integer nbArg
  NAMELIST /GRID/ NUMSIG, NUMDIR, FRLOW, FRHIGH, EXTRAPOLATION_ALLOWED_BOUC
  NAMELIST /PROC/ FileInput, PrefixOutput, FileInterpolatePoint, PrefixBounc
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
  CALL INIT_SPECTRAL_GRID

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
  OPEN(OUT%FHNDL, FILE= "ww3_bounc.inp", status='unknown')
  WRITE(OUT%FHNDL,'(A)') "$ -------------------------------------------------------------------- $"
  WRITE(OUT%FHNDL,'(A)') "$ WAVEWATCH III NetCDF boundary input processing                       $"
  WRITE(OUT%FHNDL,'(A)') "$--------------------------------------------------------------------- $"
  WRITE(OUT%FHNDL,'(A)') "$"
  WRITE(OUT%FHNDL,'(A)') "$ Boundary option: READ or WRITE"
  WRITE(OUT%FHNDL,'(A)') "$"
  WRITE(OUT%FHNDL,'(A)') " WRITE"
  WRITE(OUT%FHNDL,'(A)') "$"
  WRITE(OUT%FHNDL,'(A)') "$ Interpolation method: 1: nearest"
  WRITE(OUT%FHNDL,'(A)') "$                       2: linear interpolation"
  WRITE(OUT%FHNDL,'(A)') " 1"
  WRITE(OUT%FHNDL,'(A)') "$ Verbose (0, 1, 2)"
  WRITE(OUT%FHNDL,'(A)') "1"
  WRITE(OUT%FHNDL,'(A)') "$"
  WRITE(OUT%FHNDL,'(A)') "$ List of spectra files. These NetCDF files use the WAVEWATCH III"
  WRITE(OUT%FHNDL,'(A)') "$ format as described in the ww3_ounp.inp file. The files are"
  WRITE(OUT%FHNDL,'(A)') "$ defined relative to the directory in which the program is run."
  WRITE(OUT%FHNDL,'(A)') "$"
  DO I=1,IWBMNP
     WRITE (FILE_NAME,40) TRIM(PrefixBounc),I
40   FORMAT (a,'_',i4.4,'.nc')
     WRITE(OUT%FHNDL,'(A)') TRIM(FILE_NAME)
  END DO
  WRITE(OUT%FHNDL,'(A)') "'STOPSTRING'"
  WRITE(OUT%FHNDL,'(A)') "$"
  WRITE(OUT%FHNDL,'(A)') "$ -------------------------------------------------------------------- $"
  WRITE(OUT%FHNDL,'(A)') "$ End of input file                                                    $"
  WRITE(OUT%FHNDL,'(A)') "$ -------------------------------------------------------------------- $"
  CLOSE(OUT%FHNDL)


  
END PROGRAM INTERPOLATE_SPECTRUM
