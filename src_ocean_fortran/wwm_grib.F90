# define MyREAL(xinp) DBLE(xinp)
!****************************************************************************
!* Reading grid information from a GRIB file (case 1)                       *
!****************************************************************************
      SUBROUTINE READ_GRID_INFO_FROM_GRIB_TYPE1(TheInfo, eGrib)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(out) :: TheInfo
      integer, intent(in) :: eGrib
      !
      REAL(rkind) ::longitudeOfFirstPointInDegrees, latitudeOfFirstPointInDegrees, longitudeOfLastPointInDegrees, latitudeOfLastPointInDegrees
      REAL(rkind) :: deltaLAT, deltaLON
      REAL(rkind) :: iDirectionIncrement, jDirectionIncrement
      integer nx_dim, ny_dim, iX, iY
      call grib_get(eGrib,"numberOfPointsAlongAParallel", nx_dim)
      call grib_get(eGrib,"numberOfPointsAlongAMeridian", ny_dim)
      WRITE(STAT%FHNDL, *) 'nx_dim=', nx_dim
      WRITE(STAT%FHNDL, *) 'ny_dim=', ny_dim
      TheInfo % nx_dim = nx_dim
      TheInfo % ny_dim = ny_dim
      allocate(TheInfo % LON(nx_dim, ny_dim), TheInfo % LAT(nx_dim, ny_dim), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')
      call grib_get(eGrib, 'longitudeOfFirstGridPointInDegrees', longitudeOfFirstPointInDegrees)
      call grib_get(eGrib, 'latitudeOfFirstGridPointInDegrees', latitudeOfFirstPointInDegrees)
      call grib_get(eGrib, 'longitudeOfLastGridPointInDegrees', longitudeOfLastPointInDegrees)
      call grib_get(eGrib, 'latitudeOfLastGridPointInDegrees', latitudeOfLastPointInDegrees)
      call grib_get(eGrib, 'iDirectionIncrementInDegrees', iDirectionIncrement)
      call grib_get(eGrib, 'jDirectionIncrementInDegrees', jDirectionIncrement)
      deltaLON=(longitudeOfLastPointInDegrees - longitudeOfFirstPointInDegrees)/(nx_dim - 1)
      deltaLAT=(latitudeOfLastPointInDegrees - latitudeOfFirstPointInDegrees)/(ny_dim - 1)
      DO iX=1,nx_dim
        DO iY=1,ny_dim
          TheInfo % LON(iX,iY)=longitudeOfFirstPointInDegrees + (iX-1)*deltaLON
          TheInfo % LAT(iX,iY)=latitudeOfFirstPointInDegrees + (iY-1)*deltaLAT
        END DO
      END DO
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 2)                       *
!****************************************************************************
      SUBROUTINE phirot2phi(eLatOut, phirot, rlarot, polphi, pollam, polgam)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(out) :: eLatOut
      real(rkind), intent(in) :: phirot, rlarot, polphi, pollam, polgam
      real(rkind) :: zsinpol
      real(rkind) :: zcospol
      real(rkind) :: zphis
      real(rkind) :: zrlas
      real(rkind) :: zarg
      real(rkind) :: zgam
      zsinpol = sin(DEGRAD * polphi)
      zcospol = cos(DEGRAD * polphi)
      zphis  = DEGRAD * phirot

      if (rlarot .gt. 180.) THEN
        zrlas = rlarot - 360.
      ELSE
        zrlas = rlarot
      END IF
      zrlas = DEGRAD * zrlas;
      if (ABS(polgam) .gt. 0) THEN
        zgam  = DEGRAD * polgam;
        zarg = zsinpol*sin(zphis) + zcospol*cos(zphis) * ( cos(zrlas)*cos(zgam) - sin(zgam)*sin(zrlas))
      ELSE
        zarg = zcospol * cos(zphis) * cos(zrlas) + zsinpol * sin(zphis)
      END IF
      eLatOut = RADDEG * asin(zarg);
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 2)                       *
!****************************************************************************
      SUBROUTINE rlarot2rla(eLonOut, phirot, rlarot, polphi, pollam, polgam)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(out) :: eLonOut
      real(rkind), intent(in) :: phirot, rlarot, polphi, pollam, polgam
      !
      real(rkind) :: zsinpol
      real(rkind) :: zcospol
      real(rkind) :: zphis
      real(rkind) :: zrlas
      real(rkind) :: zgam, zlampol
      real(rkind) :: zarg1, zarg2
      zsinpol = sin(DEGRAD * polphi)
      zcospol = cos(DEGRAD * polphi)
      zphis  = DEGRAD * phirot


      if (rlarot .gt. 180.) THEN
        zrlas = rlarot - 360.
      ELSE
        zrlas = rlarot
      END IF
      zrlas = DEGRAD * zrlas
      zlampol = DEGRAD * pollam
      IF (ABS(polgam) .gt. 0) THEN
        zgam    = DEGRAD * polgam;
        zarg1   = sin (zlampol) *                                                       &
   &     (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))          &
   &     + zcospol * sin(zphis))                                                        &
   &     - cos (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam))
        zarg2   = cos (zlampol) *                                                       &
   &     (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))          &
   &     + zcospol * sin(zphis))                                                        &
   &     + sin (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam))
      ELSE
        zarg1   = sin (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis)  +                &
   &     zcospol *              sin(zphis)) -                                           &
   &     cos (zlampol) *             sin(zrlas) * cos(zphis)
        zarg2   = cos (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis)  +                &
   &     zcospol *              sin(zphis)) +                                           &
   &     sin (zlampol) *             sin(zrlas) * cos(zphis)
      END IF
      if (zarg2 .eq. 0) zarg2=1.0e-20
      eLonOut = RADDEG * atan2(zarg1,zarg2);
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 2)                       *
!****************************************************************************
      SUBROUTINE READ_GRID_INFO_FROM_GRIB_TYPE2(TheInfo, eGrib)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(out) :: TheInfo
      integer, intent(in) :: eGrib
      real(rkind) :: latitudeOfSouthernPoleInDegrees, longitudeOfSouthernPoleInDegrees
      real(rkind) :: angleOfRotationInDegrees
      real(rkind) :: latitudeOfFirstGridPointInDegrees, longitudeOfFirstGridPointInDegrees
      real(rkind) :: latitudeOfLastGridPointInDegrees, longitudeOfLastGridPointInDegrees
      real(rkind) :: iDirectionIncrementInDegrees, jDirectionIncrementInDegrees
      real(rkind) :: pollat_sp, pollon_sp, polgam, zstartlon_tot, zstartlat_tot
      real(rkind) :: zendlon_tot, zendlat_tot, dlon, dlat
      real(rkind) :: eLonR, eLatR, eLonOut, eLatOut
      real(rkind) :: pollat, pollon
      real(rkind) :: startlon_tot, startlat_tot
      integer :: nx_dim, ny_dim, iX, iY
      !
      call grib_get(eGrib, 'latitudeOfSouthernPoleInDegrees',latitudeOfSouthernPoleInDegrees)
      call grib_get(eGrib, 'longitudeOfSouthernPoleInDegrees',longitudeOfSouthernPoleInDegrees)
      call grib_get(eGrib, 'angleOfRotationInDegrees',angleOfRotationInDegrees)
      call grib_get(eGrib, 'latitudeOfFirstGridPointInDegrees',latitudeOfFirstGridPointInDegrees)
      call grib_get(eGrib, 'longitudeOfFirstGridPointInDegrees',longitudeOfFirstGridPointInDegrees)
      call grib_get(eGrib, 'latitudeOfLastGridPointInDegrees',latitudeOfLastGridPointInDegrees)
      call grib_get(eGrib, 'longitudeOfLastGridPointInDegrees',longitudeOfLastGridPointInDegrees)
      call grib_get(eGrib, 'iDirectionIncrementInDegrees',iDirectionIncrementInDegrees)
      call grib_get(eGrib, 'jDirectionIncrementInDegrees',jDirectionIncrementInDegrees)
      pollat_sp=latitudeOfSouthernPoleInDegrees
      pollon_sp=longitudeOfSouthernPoleInDegrees
      polgam=angleOfRotationInDegrees
      zstartlon_tot=longitudeOfFirstGridPointInDegrees
      zstartlat_tot=latitudeOfFirstGridPointInDegrees
      zendlon_tot=longitudeOfLastGridPointInDegrees
      zendlat_tot=latitudeOfLastGridPointInDegrees
      dlon=iDirectionIncrementInDegrees
      dlat=jDirectionIncrementInDegrees
      !
      ! Now reading the mapped grid
      !
      call grib_get(eGrib,"Nx", nx_dim)
      call grib_get(eGrib,"Ny", ny_dim)
!      WRITE(STAT%FHNDL, *) 'nx_dim = ', nx_dim
!      WRITE(STAT%FHNDL, *) 'ny_dim = ', ny_dim
      TheInfo % nx_dim = nx_dim
      TheInfo % ny_dim = ny_dim
      pollat= - pollat_sp
      pollon= pollon_sp - 180.
      startlon_tot=zstartlon_tot
      startlat_tot=zstartlat_tot

      allocate(TheInfo % LON(nx_dim, ny_dim), TheInfo % LAT(nx_dim, ny_dim), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')
      DO iX=1,nx_dim
        DO iY=1,ny_dim
          eLonR = startlon_tot + MyREAL(iX-1)*dlon
          eLatR = startlat_tot + MyREAL(iY-1)*dlat
          CALL phirot2phi(eLatOut, eLatR, eLonR, pollat, pollon, polgam)
          CALL rlarot2rla(eLonOut, eLatR, eLonR, pollat, pollon, polgam)
          TheInfo % LON(iX,iY) = eLonOut
          TheInfo % LAT(iX,iY) = eLatOut
        END DO
      END DO
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 1)                       *
!****************************************************************************
      SUBROUTINE READ_GRID_INFO_FROM_GRIB_TYPE3(TheInfo, eGrib)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(out) :: TheInfo
      integer, intent(in) :: eGrib
      REAL(rkind), allocatable :: LON_serial(:), LAT_serial(:), DATA_Serial(:)
      integer nx_dim, ny_dim, idx, eProd
      integer status, iX, iY
      !
      call grib_get(eGrib,"Nx", nx_dim)
      call grib_get(eGrib,"Ny", ny_dim)
      WRITE(STAT%FHNDL, *) 'nx_dim = ', nx_dim
      WRITE(STAT%FHNDL, *) 'ny_dim = ', ny_dim
      TheInfo % nx_dim = nx_dim
      TheInfo % ny_dim = ny_dim
      allocate(TheInfo % LON(nx_dim, ny_dim), TheInfo % LAT(nx_dim, ny_dim), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')
      eProd=nx_dim*ny_dim
      allocate(LON_serial(eProd), LAT_serial(eProd), DATA_serial(eProd), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 48')
      call grib_get_data(eGrib, LAT_serial, LON_serial, DATA_serial, status)
      idx=0
      DO iY=1,ny_dim
        DO iX=1,nx_dim
          idx=idx+1
          TheInfo % LON(iX,iY)=LON_serial(idx)
          TheInfo % LAT(iX,iY)=LAT_serial(idx)
        END DO
      END DO
      DEALLOCATE(LON_serial, LAT_serial, DATA_serial)
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file                                *
!****************************************************************************
      SUBROUTINE READ_GRID_INFO_FROM_GRIB(TheInfo, TheFile, shortName, GRIB_TYPE)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(out) :: TheInfo
      character(len=*), intent(in) :: TheFile
      character(len=20), intent(in) :: shortName
      integer, intent(in) :: GRIB_TYPE
      !
      integer ifile, i, n
      logical WeFound
      integer, allocatable :: igrib(:)
      character(len=100) eShortName
      !
      WRITE(STAT%FHNDL,*) 'TheFile=', TheFile
      FLUSH(STAT%FHNDL)
      CALL TEST_FILE_EXIST_DIE("Missing wind grib file 2: ", TheFile)
      CALL GRIB_OPEN_FILE(ifile, TheFile, 'r')
      call grib_count_in_file(ifile,n)
      allocate(igrib(n))
      !
      WRITE(STAT%FHNDL,*) 'n=', n
      WRITE(STAT%FHNDL,*) 'eShortName=', eShortName
      WRITE(STAT%FHNDL,*) 'GRIB_TYPE=', GRIB_TYPE
      FLUSH(STAT%FHNDL)
      WeFound=.FALSE.;
      DO i=1,n
        call grib_new_from_file(ifile, igrib(i))
        call grib_get(igrib(i), 'shortName', eShortName)
        WRITE(STAT%FHNDL,*) 'i=', i, ' WeFound=', WeFound
        FLUSH(STAT%FHNDL)
        IF ((TRIM(eShortName) .eq. shortName).and.(WeFound .eqv. .FALSE.)) THEN
          IF (GRIB_TYPE .eq. 1) THEN
            CALL READ_GRID_INFO_FROM_GRIB_TYPE1(TheInfo, igrib(i))
          END IF
          IF (GRIB_TYPE .eq. 2) THEN
            CALL READ_GRID_INFO_FROM_GRIB_TYPE2(TheInfo, igrib(i))
          END IF
          IF (GRIB_TYPE .eq. 3) THEN
            CALL READ_GRID_INFO_FROM_GRIB_TYPE3(TheInfo, igrib(i))
          END IF
          WeFound=.TRUE.
        END IF
        call grib_release(igrib(i))
      END DO
      WRITE(STAT%FHNDL,*) 'After loop'
      FLUSH(STAT%FHNDL)
      IF (WeFound .eqv. .FALSE.) THEN
        Print *, 'Failed to find the wind variable in the grib file'
        CALL WWM_ABORT("Wind has not been found in grib file")
      END IF
      WRITE(STAT%FHNDL, *) 'WeFound=', WeFound
      FLUSH(STAT%FHNDL)
      CALL GRIB_CLOSE_FILE(ifile)
      deallocate(igrib)
      END SUBROUTINE












