!**********************************************************************
!*                                                                    *
!**********************************************************************
SUBROUTINE GENERIC_NETCDF_ERROR_WWM(CallFct, idx, iret)
  USE NETCDF
  USE DATAPOOL
  implicit none
  integer, intent(in) :: iret, idx
  character(*), intent(in) :: CallFct
  character(len=500) :: CHRERR
  character(len=1000)  :: wwmerr
  IF (iret .NE. nf90_noerr) THEN
     CHRERR = nf90_strerror(iret)
     WRITE(wwmerr,*) 'NETCDF error in routine ', TRIM(CallFct), ' Error Message: ', TRIM(CHRERR), ' Position in the routine :', idx
     CALL WWM_ABORT(wwmerr)
  ENDIF
END SUBROUTINE GENERIC_NETCDF_ERROR_WWM
!**********************************************************************
!*                                                                    *
!**********************************************************************
SUBROUTINE OUTPUT_SPECTRUM
  USE NETCDF
  USE DATAPOOL
  IMPLICIT NONE
  character(len=400) FILE_NAME
  REAL(8), allocatable :: ListTimeWrite(:)
  REAL(4), allocatable :: longitude(:), latitude(:)
  REAL(4), allocatable :: frequency(:), frequency1(:), frequency2(:)
  REAL(4), allocatable :: direction(:)
  REAL(4), allocatable :: efth_write(:, :, :, :)
  REAL(rkind), allocatable :: WBACOUT(:, :, :, :)
  character (len = *), parameter :: CallFct = "OUTPUT_SPECTRUM"
  integer nbTime, iTime, iFreq, iDir, IB
  integer var_id, ncid, iret
  integer station_dims, time_dims, frequency_dims, string16_dims
  integer direction_dims
  integer ARR_I1(1)
  REAL(rkind) FRATIO_sqrt, eFR
  nbTime = eVAR_BOUC_WAM % nbTime
  allocate(ListTimeWrite(nbTime))
  allocate(longitude(nbTime), latitude(nbTime))
  allocate(frequency(NUMSIG), frequency1(NUMSIG), frequency2(NUMSIG))
  allocate(direction(NUMDIR))
  allocate(efth_write(nbTime, 1, NUMSIG, NUMDIR))
  FRATIO_sqrt = SQRT(FRATIO)
  DO iFreq=1,NUMSIG
     eFR = SPSIG(iFreq ) / PI2
     frequency(iFreq) = REAL(eFR)
     frequency1(iFreq) = REAL(eFR / FRATIO_sqrt)
     frequency2(iFreq) = REAL(eFR * FRATIO_sqrt)
  END DO
  DO iDir=1,NUMDIR
     direction(iDir) = REAL(290.0_rkind - SPDIR(iDir))
  END DO
  allocate(WBACOUT(NUMSIG, NUMDIR, IWBMNP, nbTime))
  CALL READ_GRIB_WAM_BOUNDARY_WBAC(WBACOUT)
  DO IB=1,IWBMNP
     WRITE (FILE_NAME,40) TRIM(PrefixOutput),IB
40   FORMAT (a,'_',i4.4,'.nc')
     !
     ! Creating the file
     !
     iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
     !
     iret = nf90_def_dim(ncid, 'time', 0, time_dims)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
     !
     iret = nf90_def_dim(ncid, 'station', 1, station_dims)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
     !
     iret = nf90_def_dim(ncid, 'string16', 1, string16_dims)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
     !
     iret = nf90_def_dim(ncid, 'frequency', NUMSIG, frequency_dims)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
     !
     iret = nf90_def_dim(ncid, 'direction', NUMDIR, direction_dims)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
     !
     iret=nf90_def_var(ncid,"time",NF90_DOUBLE,(/ time_dims/),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
     iret=nf90_put_att(ncid,var_id,"units","days since 1990-01-01T00:00:00Z")
     iret=nf90_put_att(ncid,var_id,"long_name","julian day (UT)")
     iret=nf90_put_att(ncid,var_id,"standard_name","time")
     iret=nf90_put_att(ncid,var_id,"conventions","Relative julian days with decimal part (as parts of the day)")
     iret=nf90_put_att(ncid,var_id,"axis","T")
     !
     iret=nf90_def_var(ncid,"station",NF90_INT,(/ station_dims/),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
     iret=nf90_put_att(ncid,var_id,"long_name","station id")
     iret=nf90_put_att(ncid,var_id,"axis","X")
     !
     iret=nf90_def_var(ncid,"string16",NF90_INT,(/ string16_dims/),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
     iret=nf90_put_att(ncid,var_id,"long_name","station_name number of characters")
     iret=nf90_put_att(ncid,var_id,"axis","W")
     !
     iret=nf90_def_var(ncid,"station_name",NF90_CHAR,(/ station_dims, string16_dims/),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
     iret=nf90_put_att(ncid,var_id,"long_name","station name")
     iret=nf90_put_att(ncid,var_id,"content","XW")
     iret=nf90_put_att(ncid,var_id,"associates","station string16")
     !
     iret=nf90_def_var(ncid,"longitude",NF90_REAL,(/ station_dims, time_dims/),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
     iret=nf90_put_att(ncid,var_id,"units","degree_east")
     iret=nf90_put_att(ncid,var_id,"long_name","longitude")
     iret=nf90_put_att(ncid,var_id,"standard_name","longitude")
     iret=nf90_put_att(ncid,var_id,"valid_min",-180)
     iret=nf90_put_att(ncid,var_id,"valid_max", 180)
     iret=nf90_put_att(ncid,var_id,"content","TX")
     iret=nf90_put_att(ncid,var_id,"associates","time station")
     !
     iret=nf90_def_var(ncid,"latitude",NF90_REAL,(/ station_dims, time_dims /),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
     iret=nf90_put_att(ncid,var_id,"units","degree_north")
     iret=nf90_put_att(ncid,var_id,"long_name","latitude")
     iret=nf90_put_att(ncid,var_id,"standard_name","latitude")
     iret=nf90_put_att(ncid,var_id,"valid_min",-90)
     iret=nf90_put_att(ncid,var_id,"valid_max", 90)
     iret=nf90_put_att(ncid,var_id,"content","TX")
     iret=nf90_put_att(ncid,var_id,"associates","time station")
     !
     iret=nf90_def_var(ncid,"frequency",NF90_REAL,(/ frequency_dims /),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
     iret=nf90_put_att(ncid,var_id,"units","s-1")
     iret=nf90_put_att(ncid,var_id,"long_name","frequency of center band")
     iret=nf90_put_att(ncid,var_id,"standard_name","sea_surface_wave_frequency")
     iret=nf90_put_att(ncid,var_id,"globwave_name","frequency")
     iret=nf90_put_att(ncid,var_id,"valid_min", 0)
     iret=nf90_put_att(ncid,var_id,"valid_max",10)
     iret=nf90_put_att(ncid,var_id,"axis","Y")
     !
     iret=nf90_def_var(ncid,"frequency1",NF90_REAL,(/ frequency_dims /),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
     iret=nf90_put_att(ncid,var_id,"units","s-1")
     iret=nf90_put_att(ncid,var_id,"long_name","frequency of lower band")
     iret=nf90_put_att(ncid,var_id,"standard_name","frequency_of_lower_band")
     iret=nf90_put_att(ncid,var_id,"globwave_name","frequency_lower_band")
     iret=nf90_put_att(ncid,var_id,"valid_min", 0)
     iret=nf90_put_att(ncid,var_id,"valid_max",10)
     iret=nf90_put_att(ncid,var_id,"axis","Y")
     iret=nf90_put_att(ncid,var_id,"associates","frequency")
     !
     iret=nf90_def_var(ncid,"frequency2",NF90_REAL,(/ frequency_dims /),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
     iret=nf90_put_att(ncid,var_id,"units","s-1")
     iret=nf90_put_att(ncid,var_id,"long_name","frequency of upper band")
     iret=nf90_put_att(ncid,var_id,"standard_name","frequency_of_upper_band")
     iret=nf90_put_att(ncid,var_id,"globwave_name","frequency_upper_band")
     iret=nf90_put_att(ncid,var_id,"valid_min", 0)
     iret=nf90_put_att(ncid,var_id,"valid_max",10)
     iret=nf90_put_att(ncid,var_id,"axis","Y")
     iret=nf90_put_att(ncid,var_id,"associates","frequency")
     !
     iret=nf90_def_var(ncid,"direction",NF90_REAL,(/ direction_dims /),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
     iret=nf90_put_att(ncid,var_id,"units","degree")
     iret=nf90_put_att(ncid,var_id,"long_name","sea surface wave to direction")
     iret=nf90_put_att(ncid,var_id,"standard_name","sea_surface_wave_to_direction")
     iret=nf90_put_att(ncid,var_id,"globwave_name","direction")
     iret=nf90_put_att(ncid,var_id,"valid_min", 0)
     iret=nf90_put_att(ncid,var_id,"valid_max",360)
     iret=nf90_put_att(ncid,var_id,"axis","Z")
     !
     iret=nf90_def_var(ncid,"efth",NF90_REAL,(/ station_dims, frequency_dims, direction_dims, time_dims /),var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
     iret=nf90_put_att(ncid,var_id,"units","m2 s rad-1")
     iret=nf90_put_att(ncid,var_id,"long_name","sea surface wave directional variance spectral density")
     iret=nf90_put_att(ncid,var_id,"standard_name","sea_surface_wave_directional_variance_spectral_density")
     iret=nf90_put_att(ncid,var_id,"globwave_name","directional_variance_spectral_density")
     iret=nf90_put_att(ncid,var_id,"scale_factor", 1.0)
     iret=nf90_put_att(ncid,var_id,"add_offset", 0.0)
     iret=nf90_put_att(ncid,var_id,"valid_min", 0)
     iret=nf90_put_att(ncid,var_id,"valid_max",1e+20)
     iret=nf90_put_att(ncid,var_id,"content","TXYZ")
     iret=nf90_put_att(ncid,var_id,"associates","time station frequency direction")
     !
     iret = nf90_close(ncid)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)
     !
     ! Writing the data
     !
     iret = nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 19, iret)
     !
     iret=nf90_inq_varid(ncid, "time", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, iret)
     iret=nf90_put_var(ncid,var_id,ListTimeWrite,start=(/ 1 /), count=(/ nbTime /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 21, iret)
     !
     iret=nf90_inq_varid(ncid, "station", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 22, iret)
     iret=nf90_put_var(ncid,var_id,ARR_I1,start=(/ 1 /), count=(/ 1 /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 23, iret)
     !
     longitude(:) = REAL(XP(IB))
     iret=nf90_inq_varid(ncid, "longitude", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 24, iret)
     iret=nf90_put_var(ncid,var_id,longitude,start=(/ 1, 1 /), count=(/ 1, nbTime /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 25, iret)
     !
     latitude(:) = REAL(YP(IB))
     iret=nf90_inq_varid(ncid, "latitude", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 26, iret)
     iret=nf90_put_var(ncid,var_id,latitude,start=(/ 1, 1 /), count=(/ 1, nbTime /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 27, iret)
     !
     iret=nf90_inq_varid(ncid, "frequency", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 28, iret)
     iret=nf90_put_var(ncid,var_id,frequency,start=(/ 1 /), count=(/ NUMSIG /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 29, iret)
     !
     iret=nf90_inq_varid(ncid, "frequency1", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 30, iret)
     iret=nf90_put_var(ncid,var_id,frequency1,start=(/ 1 /), count=(/ NUMSIG /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 31, iret)
     !
     iret=nf90_inq_varid(ncid, "frequency2", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 32, iret)
     iret=nf90_put_var(ncid,var_id,frequency2,start=(/ 1 /), count=(/ NUMSIG /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 33, iret)
     !
     iret=nf90_inq_varid(ncid, "direction", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 34, iret)
     iret=nf90_put_var(ncid,var_id,direction,start=(/ 1 /), count=(/ NUMDIR /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 35, iret)
     !
     DO iTime=1,nbTime
        DO iFreq=1,NUMSIG
           DO iDir=1,NUMDIR
              efth_write(iTime, 1, iFreq, iDir) = REAL(WBACOUT(iFreq, iDir, IB, iTime))
           END DO
        END DO
     END DO
     iret=nf90_inq_varid(ncid, "efth", var_id)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 36, iret)
     iret=nf90_put_var(ncid,var_id,efth_write,start=(/ 1,1,1,1 /), count=(/ 1, NUMSIG, NUMDIR, nbTime /))
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 37, iret)
     !
     iret = nf90_close(ncid)
     CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 38, iret)
  END DO
  deallocate(ListTimeWrite)
  deallocate(longitude, latitude)
  deallocate(frequency, frequency1, frequency2)
  deallocate(direction)
  deallocate(WBACOUT)
END SUBROUTINE OUTPUT_SPECTRUM
