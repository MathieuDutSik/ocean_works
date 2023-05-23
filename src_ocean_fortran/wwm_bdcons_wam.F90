# define MyREAL(xinp) DBLE(xinp)
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_BND_INTERPOLATION_ARRAY(TheInfo)
      USE DATAPOOL
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(in) :: TheInfo
      integer IP
      real(rkind) eX, eY
      integer eCF_IX, eCF_IY
      real(rkind) eCF_COEFF(4)
      LOGICAL EXTRAPO_OUT
      integer nbExtrapolation
      Print *, 'Begin COMPUTE_BND_INTERPOLATION_ARRAY'
      Print *, 'EXTRAPOLATION_ALLOWED_BOUC=', EXTRAPOLATION_ALLOWED_BOUC
      allocate(CF_IX_BOUC(IWBMNP), CF_IY_BOUC(IWBMNP), CF_COEFF_BOUC(4,IWBMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('CF_*_BOUC allocation error')
      nbExtrapolation = 0
      DO IP=1,IWBMNP
        eX=XP(IP)
        eY=YP(IP)
        CALL COMPUTE_SINGLE_INTERPOLATION_INFO(TheInfo, EXTRAPOLATION_ALLOWED_BOUC, eX, eY, eCF_IX, eCF_IY, eCF_COEFF, EXTRAPO_OUT)
        Print *, "IP=", IP, " eCF_IX/IY=", eCF_IX, eCF_IY, " sum(eCF_COEFF)=", sum(eCF_COEFF)
        CF_IX_BOUC(IP) = eCF_IX
        CF_IY_BOUC(IP) = eCF_IY
        CF_COEFF_BOUC(:,IP) = eCF_COEFF
        IF (EXTRAPO_OUT .eqv. .TRUE.) THEN
          nbExtrapolation=nbExtrapolation + 1
        END IF
      END DO
      IF (EXTRAPOLATION_ALLOWED_BOUC) THEN
         Print *, "nbExtrapolation=", nbExtrapolation
      END IF
      END SUBROUTINE COMPUTE_BND_INTERPOLATION_ARRAY
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_GRIB_WAM_BOUNDARY
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      INTEGER ifile, IFILE_IN
      LOGICAL STEPRANGE_IN
      type(FD_FORCING_GRID) :: TheInfo
      character(len=20) shortName
      integer GRIB_TYPE
      integer, allocatable :: ListDir_i(:), ListFreq_i(:)
      integer nbTotalNumberEntry
      LOGICAL IsFirst
      character(len=281) FULL_FILE
      integer i, idir, ifreq, n
      integer nbdir_wam_read, nbfreq_wam_read
      integer freqScal, dirScal
      integer, allocatable :: igrib(:)
      real(rkind) eDIR, eFR, eFreq
      real(rkind) eDiff, eDiff1, eDiff2
      real(rkind) eWD1, eWD2
      logical IsAssigned
      integer IS, ID, idx
      integer ID1, ID2
      real(rkind) eTimeOut, DeltaDiff
      real(rkind) DELTH_WAM, CO1, WETAIL_WAM
      integer M
      CALL TEST_FILE_EXIST_DIE("Missing list of WAM files: ", TRIM(WAV%FNAME))
      OPEN(WAV%FHNDL,FILE=WAV%FNAME,STATUS='OLD')
      Print *, "WAV%FHNDL=", WAV%FHNDL, " WAV%FNAME=", WAV%FNAME
      STEPRANGE_IN = .TRUE.
      !
      ! Determining the number of times
      !
      NUM_WAM_SPEC_FILES = 0
      DO
         READ( WAV%FHNDL, *, IOSTAT = ISTAT )
         IF ( ISTAT /= 0 ) EXIT
         NUM_WAM_SPEC_FILES = NUM_WAM_SPEC_FILES + 1
      END DO
      REWIND(WAV%FHNDL)
      IF (NUM_WAM_SPEC_FILES .eq. 0) THEN
         CALL WWM_ABORT('We need at least one file in order for this to work')
      END IF
      !
      ! Reading the file names
      !
      ALLOCATE(WAM_SPEC_FILE_NAMES_BND(NUM_WAM_SPEC_FILES), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons_wam, allocate error 9')
      DO IFILE_IN = 1, NUM_WAM_SPEC_FILES
         READ( WAV%FHNDL, *) WAM_SPEC_FILE_NAMES_BND(IFILE_IN)
      END DO
      CLOSE (WAV%FHNDL)
      !
      ! Determining total number of entries
      ! and also the number of directions and frequencies
      !
      nbTotalNumberEntry=0
      IsFirst=.TRUE.
      DO IFILE_IN = 1, NUM_WAM_SPEC_FILES
         FULL_FILE=TRIM(PREFIX_WAVE_FILE) // TRIM(WAM_SPEC_FILE_NAMES_BND(IFILE_IN))
         CALL TEST_FILE_EXIST_DIE("Missing wam grib file 1: ", TRIM(FULL_FILE))
         CALL GRIB_OPEN_FILE(ifile, TRIM(FULL_FILE), 'r')
         call grib_count_in_file(ifile,n)
         allocate(igrib(n))
         DO i=1,n
            call grib_new_from_file(ifile, igrib(i))
            call grib_get(igrib(i), 'numberOfDirections', nbdir_wam_read)
            call grib_get(igrib(i), 'numberOfFrequencies', nbfreq_wam_read)
            IF (IsFirst .eqv. .TRUE.) THEN
               nbdir_wam = nbdir_wam_read
               nbfreq_wam = nbfreq_wam_read
               call grib_get(igrib(i), 'directionScalingFactor', dirScal)
               call grib_get(igrib(i), 'frequencyScalingFactor', freqScal)
               allocate(ListDir_i(nbdir_wam), ListFreq_i(nbfreq_wam), ListDir_wam(nbdir_wam), ListFreq_wam(nbfreq_wam), DFIM_wam(nbFreq_wam), stat=istat)
               call grib_get(igrib(i), 'scaledDirections', ListDir_i)
               call grib_get(igrib(i), 'scaledFrequencies', ListFreq_i)
               DO idir=1,nbdir_wam
                  eDir = MyREAL(ListDir_i(idir)) / MyREAL(dirScal)
                  eDir = 270 - eDir
                  IF (eDir .le. ZERO) THEN
                     eDir = eDir+ 360
                  END IF
                  IF (eDir .ge. 360) THEN
                     eDir = eDir - 360
                  END IF
                  ListDir_wam(idir) = eDir
               END DO
               DO ifreq=1,nbfreq_wam
                  eFreq = MyREAL(ListFreq_i(ifreq)) / MyREAL(freqScal)
                  ListFreq_wam(ifreq) = eFreq
               END DO
               FRATIO = ListFreq_wam(2) / ListFreq_wam(1)
               DELTH_WAM = PI2 / MyREAL(nbdir_wam)
               CO1 = 0.5*(FRATIO-1.)*DELTH_WAM
               DFIM_WAM(1) = CO1 * ListFreq_wam(1)
               DO M=2,nbFreq_wam-1
                  DFIM_WAM(M) = CO1 * (ListFreq_wam(M) + ListFreq_wam(M-1))
               ENDDO
               DFIM_WAM(nbFreq_wam) = CO1 * ListFreq_wam(nbFreq_wam-1)
               WETAIL_WAM = 0.25
               DELT25_WAM = WETAIL_WAM*ListFreq_wam(nbFreq_wam)*DELTH_WAM
               deallocate(ListDir_i, ListFreq_i)
            ELSE
               IF ((nbdir_wam .ne. nbdir_wam_read).or.(nbfreq_wam .ne. nbfreq_wam_read)) THEN
                  Print *, 'nbdir_wam =', nbdir_wam,  'nbdir_wam_read =', nbdir_wam_read
                  Print *, 'nbfreq_wam=', nbfreq_wam, 'nbfreq_wam_read=', nbfreq_wam_read
                  CALL WWM_ABORT('number of frequencies/directions is inconsistent')
               END IF
            END IF
            IsFirst=.FALSE.
            call grib_release(igrib(i))
         END DO
         deallocate(igrib)
         CALL GRIB_CLOSE_FILE(ifile)
         nbTotalNumberEntry = nbTotalNumberEntry + n / (nbdir_wam * nbfreq_wam)
      END DO
      ALLOCATE(eVAR_BOUC_WAM % ListTime(nbTotalNumberEntry), ListIFileWAM(nbTotalNumberEntry), stat=istat)
      eVAR_BOUC_WAM % nbTime = nbTotalNumberEntry
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons_wam, allocate error 9')
      idx=0
      DO IFILE_IN = 1, NUM_WAM_SPEC_FILES
         FULL_FILE = TRIM(PREFIX_WAVE_FILE) // TRIM(WAM_SPEC_FILE_NAMES_BND(IFILE_IN))
         CALL TEST_FILE_EXIST_DIE("Missing wam grib file 2: ", TRIM(FULL_FILE))
         CALL GRIB_OPEN_FILE(ifile, TRIM(FULL_FILE), 'r')
         call grib_count_in_file(ifile,n)
         allocate(igrib(n))
         DO i=1,n
            call grib_new_from_file(ifile, igrib(i))
            call grib_get(igrib(i), 'directionNumber', idir)
            call grib_get(igrib(i), 'frequencyNumber', ifreq)
            IF ((idir .eq. 1).and.(ifreq .eq. 1)) THEN
               CALL RAW_READ_TIME_OF_GRIB_FILE(igrib(i), STEPRANGE_IN, eTimeOut)
               !
               idx=idx+1
               eVAR_BOUC_WAM % ListTime(idx) = eTimeOut
               ListIFileWAM(idx) = IFILE_IN
            END IF
            call grib_release(igrib(i))
         END DO
         deallocate(igrib)
         CALL GRIB_CLOSE_FILE(ifile)
      END DO
      !
      ! reading the grid
      !
      shortName='2dfd'
      GRIB_TYPE=1 ! 1 for ECMWF
      IFILE_IN = 1
      FULL_FILE = TRIM(PREFIX_WAVE_FILE) // TRIM(WAM_SPEC_FILE_NAMES_BND(IFILE_IN))
      CALL READ_GRID_INFO_FROM_GRIB(TheInfo, TRIM(FULL_FILE), shortName, GRIB_TYPE)
      CALL COMPUTE_BND_INTERPOLATION_ARRAY(TheInfo)
      deallocate(TheInfo % LON, TheInfo % LAT)
      nx_wam = TheInfo % nx_dim
      ny_wam = TheInfo % ny_dim
      !
      ! Now the spectral interpolation arrays
      !
      allocate(WAM_ID1(NUMDIR), WAM_ID2(NUMDIR), WAM_WD1(NUMDIR), WAM_WD2(NUMDIR), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('CF_*_BOUC allocation error')
      WAM_ID1=0
      WAM_ID2=0
      DO ID=1,NUMDIR
        eDIR=SPDIR(ID) * RADDEG
        IsAssigned=.false.
        DO ID1=1,nbdir_wam
          IF (ID1 .lt. nbdir_wam) THEN
            ID2=ID1+1
          ELSE
            ID2=1
          END IF
          IF (IsAssigned .eqv. .false.) THEN
            eDiff = ListDir_wam(ID2) - ListDir_wam(ID1)
            IF (eDiff .gt. 180) THEN
              eDiff = eDiff - 360.0
            END IF
            IF (eDiff .lt. -180) THEN
              eDiff = eDiff + 360.0
            END IF
            !
            eDiff1=eDIR - ListDir_wam(ID1)
            IF (eDiff1 .gt. 180) THEN
              eDiff1 = eDiff1 - 360.0
            END IF
            IF (eDiff1 .lt. -180) THEN
              eDiff1 = eDiff1 + 360.0
            END IF
            !
            eDiff2=ListDir_wam(ID2) - eDir
            IF (eDiff2 .gt. 180) THEN
              eDiff2 = eDiff2 - 360.0
            END IF
            IF (eDiff2 .lt. -180) THEN
              eDiff2 = eDiff2 + 360.0
            END IF
            DeltaDiff = abs(eDiff) - abs(eDiff1) - abs(eDiff2)
            IF (abs(DeltaDiff) .le. 1.0) THEN
              eWD1 = eDiff2 / eDiff
              eWD2 = eDiff1 / eDiff
              IsAssigned=.TRUE.
              WAM_ID1(ID) = ID1
              WAM_ID2(ID) = ID2
              WAM_WD1(ID) = eWD1
              WAM_WD2(ID) = eWD2
            END IF
          END IF
        END DO
        IF (IsAssigned .eqv. .FALSE.) THEN
          CALL WWM_ABORT('Error in the interpolation direction')
        END IF
      END DO
      allocate(WAM_IS1(NUMSIG), WAM_IS2(NUMSIG), WAM_WS1(NUMSIG), WAM_WS2(NUMSIG), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('CF_*_BOUC allocation error')
      WAM_IS1=0
      WAM_IS2=0
      DO IS=1,NUMSIG
        IsAssigned=.FALSE.
        eFR=FR(IS)
        DO iFreq=1,nbfreq_wam-1
          IF (IsAssigned .eqv. .FALSE.) THEN
            eDiff=ListFreq_wam(iFreq+1) - ListFreq_wam(iFreq)
            eDiff1=eFR - ListFreq_wam(iFreq)
            eDiff2=ListFreq_wam(iFreq+1) - eFR
            IF ((eDiff1 .ge. 0).and.(eDiff2 .ge.0)) THEN
              IsAssigned=.TRUE.
              WAM_IS1(IS)=iFreq
              WAM_IS2(IS)=iFreq+1
              WAM_WS1(IS)=eDiff2 / eDiff
              WAM_WS2(IS)=eDiff1 / eDiff
            END IF
          END IF
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL_NAKED(WBAC_WAM, IFILE_IN, eTimeSearch)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      real(rkind), intent(out) :: WBAC_WAM(nbdir_wam, nbfreq_wam, nx_wam, ny_wam)
      integer, intent(in) :: IFILE_IN
      real(rkind), intent(in) :: eTimeSearch
      !
      real(rkind) :: values(nx_wam*ny_wam)
      integer :: DirFreqStatus(nbdir_wam, nbfreq_wam)
      character(len=281) FULL_FILE
      integer i, n
      integer eDiff
      integer idx, idir, ifreq
      real(rkind) DeltaDiff
      integer, allocatable :: igrib(:)
      LOGICAL :: STEPRANGE_IN = .TRUE.
      real(rkind) eTimeOut
      character(len=140) eShortName
      integer iX, iY, ifile
      real(rkind) eVal
      DirFreqStatus=0
      FULL_FILE=TRIM(PREFIX_WAVE_FILE) // TRIM(WAM_SPEC_FILE_NAMES_BND(IFILE_IN))
      CALL TEST_FILE_EXIST_DIE("Missing wam grib file 3: ", TRIM(FULL_FILE))
      CALL GRIB_OPEN_FILE(ifile, TRIM(FULL_FILE), 'r')
      call grib_count_in_file(ifile,n)
      allocate(igrib(n))
      WBAC_WAM=0
      DO i=1,n
        call grib_new_from_file(ifile, igrib(i))
        CALL RAW_READ_TIME_OF_GRIB_FILE(igrib(i), STEPRANGE_IN, eTimeOut)
        DeltaDiff = abs(eTimeOut - eTimeSearch)
        IF (DeltaDiff .le. 1.0E-8) THEN
          call grib_get(igrib(i), 'shortName', eShortName)
          IF (TRIM(eShortName) .eq. '2dfd') THEN
            call grib_get(igrib(i), 'directionNumber', idir)
            call grib_get(igrib(i), 'frequencyNumber', ifreq)
            CALL grib_get(igrib(i), 'values', values)
            DirFreqStatus(idir, ifreq) = 1
            idx=0
            DO iY=1,ny_wam
              DO iX=1,nx_wam
                idx=idx+1
                IF (values(idx) .lt. MyREAL(9990)) THEN
                  eVal=EXP(values(idx)*LOG(10.))
                  WBAC_WAM(idir, ifreq, iX,iY) = eVal
                ENDIF
              END DO
            END DO
          END IF
        END IF
        CALL grib_release(igrib(i))
      END DO
      deallocate(igrib)
      CALL GRIB_CLOSE_FILE(ifile)
      eDiff= sum(DirFreqStatus) - nbdir_wam * nbfreq_wam
      if (eDiff .ne. 0) THEN
        CALL WWM_ABORT('Error reading WAM file. Some direction/frequencies not assigned')
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL(WBACOUT, IFILE, eTimeSearch)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: WBACOUT(NUMSIG,NUMDIR,IWBMNP)
      integer, intent(in) :: IFILE
      real(rkind), intent(in) :: eTimeSearch
      !
      real(rkind) :: WBAC_WAM    (nbdir_wam, nbfreq_wam, nx_wam, ny_wam)
      real(rkind) :: WBAC_WAM_LOC(nbdir_wam, nbfreq_wam)
      integer ID1, ID2, IS1, IS2
      integer ID, IS, J, IP
      real(rkind) WD1, WD2, WS1, WS2
      real(rkind) WALOC(NUMSIG,NUMDIR)
      integer IX, IY
      real(rkind) eAC_1, eAC_2, eAC
      real(rkind) EM, HS_WAM, eSum, quot
      integer M, K
      LOGICAL :: DoHSchecks = .TRUE.
      real(rkind) ETOT, tmp(NUMSIG), DS, ETAIL, HS_WWM, EMwork
      INTEGER SHIFTXY(4,2)
      SHIFTXY(1,1)=0
      SHIFTXY(1,2)=0
      SHIFTXY(2,1)=1
      SHIFTXY(2,2)=0
      SHIFTXY(3,1)=0
      SHIFTXY(3,2)=1
      SHIFTXY(4,1)=1
      SHIFTXY(4,2)=1
      CALL READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL_NAKED(WBAC_WAM, IFILE, eTimeSearch)
      Print *, "READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL iFile=", iFile, " sum=", sum(WBAC_WAM)
      DO IP=1,IWBMNP
        IX=CF_IX_BOUC(IP)
        IY=CF_IY_BOUC(IP)
        WBAC_WAM_LOC=0
        DO J=1,4
          WBAC_WAM_LOC(:,:) = WBAC_WAM_LOC(:,:) + CF_COEFF_BOUC(J,IP)*WBAC_WAM(:,:,IX+SHIFTXY(J,1),IY+SHIFTXY(J,2))
        END DO
        Print *, "IP=", IP, " sum(WBAC_WAM_LOC)=", sum(WBAC_WAM_LOC)
        !
        IF (DoHSchecks) THEN
          EM=0
          DO M=1,nbfreq_wam
            eSum=0
            DO K=1,nbdir_wam
              eSum = eSum + WBAC_WAM_LOC(K,M)
            END DO
            EM = EM + DFIM_WAM(M)*eSum
          END DO
          EM = EM + DELT25_WAM*eSum
          EMwork=MAX(ZERO, EM)
          HS_WAM = 4.*SQRT(EMwork)
        END IF
        WALOC=0
        DO IS=1,NUMSIG
          DO ID=1,NUMDIR
            ID1=WAM_ID1(ID)
            ID2=WAM_ID2(ID)
            WD1=WAM_WD1(ID)
            WD2=WAM_WD2(ID)
            !
            IS1=WAM_IS1(IS)
            IS2=WAM_IS2(IS)
            WS1=WAM_WS1(IS)
            WS2=WAM_WS2(IS)
            !
            IF (IS1 .gt. 0) THEN
              eAC_1=WD1 * WBAC_WAM_LOC(ID1, IS1) + WD2 * WBAC_WAM_LOC(ID2, IS1)
              eAC_2=WD1 * WBAC_WAM_LOC(ID1, IS2) + WD2 * WBAC_WAM_LOC(ID2, IS2)
              eAC=WS1 * eAC_1 + WS2 * eAC_2
              WALOC(IS,ID)=eAC / (SPSIG(IS) * PI2)
            END IF
          END DO
        END DO
        IF (DoHSchecks) THEN
          ETOT=0
          DO ID=1,NUMDIR
            tmp(:) = WALOC(:,id) * spsig
            ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
            do is = 2, NUMSIG
              ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
            end do
            ETOT = ETOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
          END DO
          DS    = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)
          ETAIL = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,2) * DDIR * DS
          ETOT  = ETOT + TAIL_ARR(6) * ETAIL
          HS_WWM = 4*SQRT(MAX(0.0, ETOT))
          IF (ETOT .gt. 0) THEN
            quot = EM/ETOT
          ELSE
            quot = -1
          END IF
        END IF
        WBACOUT(:,:,IP)=WALOC
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_GRIB_WAM_BOUNDARY_WBAC(WBACOUT)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: WBACOUT(NUMSIG, NUMDIR, IWBMNP, eVAR_BOUC_WAM % nbTime)
      REAL(rkind)   :: WBACREAD(NUMSIG, NUMDIR, IWBMNP)
      !
      integer iTime
      real(rkind) eTimeDay
      integer iFile
      DO iTime=1, eVAR_BOUC_WAM % nbTime
        eTimeDay=eVAR_BOUC_WAM % ListTime(iTime)
        iFile=ListIFileWAM(iTime)
        CALL READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL(WBACREAD, iFile, eTimeDay)
        Print *, "READ_GRIB_WAM_BOUNDARY_WBAC iTime=", iTime, " sum=", sum(WBACREAD)
        WBACOUT(:,:,:,iTime) = WBACREAD
      END DO
      END SUBROUTINE READ_GRIB_WAM_BOUNDARY_WBAC
!****************************************************************************
!* Raw reading of time entry for GRIB                                       *
!****************************************************************************
      SUBROUTINE RAW_READ_TIME_OF_GRIB_FILE(eGrib, STEPRANGE_IN, eTimeOut)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      integer, intent(in) :: eGrib
      LOGICAL, intent(in) :: STEPRANGE_IN
      real(rkind), intent(out) :: eTimeOut
      !
      LOGICAL :: USE_DATATIME = .TRUE.
      integer eYear, eMonth, eDay, resYear, resMonth
      integer eHour, eMin, eSec
      integer dataDate, stepRange, dataTime
      character (len=15) :: eStrTime
      REAL(rkind) :: eTimeBase
      call grib_get(eGrib, 'dataDate', dataDate)
      eYear=(dataDate - mod(dataDate,10000))/10000
      resYear=dataDate - 10000*eYear
      eMonth=(resYear - mod(resYear,100))/100
      resMonth=resYear - 100*eMonth;
      eDay=resMonth
      IF (STEPRANGE_IN) THEN
        call grib_get(eGrib, 'stepRange', stepRange)
      ELSE
        stepRange=0
      END IF
      IF (USE_DATATIME) THEN
        call grib_get(eGrib, 'dataTime', dataTime)
        eHour=(dataTime - mod(dataTime,100))/100
        eMin=dataTime - 100*eHour
        eSec=0
      ELSE
        eHour=0
        eMin=0
        eSec=0
      END IF
      WRITE(eStrTime,10) eYear, eMonth, eDay, eHour, eMin, eSec
 10   FORMAT(i4.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)
      CALL CT2MJD(eStrTime, eTimeBase)
      eTimeOut=eTimeBase + MyREAL(stepRange)/24.0_rkind
      END SUBROUTINE
