      SUBROUTINE INIT_SPECTRAL_GRID()
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER :: IS, ID
      REAL(rkind)    :: SGLOW, SGHIGH, FRINTF, FRINTH, SFAC
      REAL(rkind)    :: MAXDIR, MINDIR

      ALLOCATE( SPSIG(NUMSIG), SPDIR(NUMDIR), FR(NUMSIG), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 5')
      SPSIG = zero
      SPDIR = zero
      FR    = zero

      SGLOW  = PI2*FRLOW
      SGHIGH = PI2*FRHIGH

      FRINTF = LOG(SGHIGH/SGLOW)/DBLE(NUMSIG-1)
      SFAC   = EXP(FRINTF)


      FRATIO = SFAC
      FRINTH = SQRT(SFAC)
      FR(1)  = FRLOW

      DO IS = 2, NUMSIG
        FR(IS) = FR(IS-1) * SFAC
      END DO

      SPSIG = FR * PI2

      MAXDIR = PI2
      MINDIR = 0
      DDIR = ABS(MAXDIR-MINDIR)/DBLE(NUMDIR)
      DO ID = 1, NUMDIR
         SPDIR(ID) = MINDIR + DDIR * DBLE(ID-1)
      END DO

      ALLOCATE(SIGPOW(NUMSIG,6), DS_BAND(0:NUMSIG+1), DS_INCR(0:NUMSIG+1), stat=istat)
      DS_BAND(0)     = SPSIG(2)- SPSIG(1)
      DS_BAND(1)     = DS_BAND(0)
      DS_BAND(NUMSIG)   = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)
      DS_BAND(NUMSIG+1) = DS_BAND(NUMSIG)
      DS_INCR(0)     = DS_BAND(0)
      DS_INCR(1)     = DS_BAND(0)
      DS_INCR(NUMSIG)   = DS_BAND(NUMSIG)
      DS_INCR(NUMSIG+1) = DS_INCR(NUMSIG)
      DO IS = 2, NUMSIG-1 ! Bandwith at gridpoints
         DS_BAND(IS) = (SPSIG(IS)-SPSIG(IS-1))/2. + (SPSIG(IS+1)-SPSIG(IS))/2.
      END DO
      DO IS = 2, NUMSIG ! Stepwidth between gridpoints K and K-1
         DS_INCR(IS) = SPSIG(IS) - SPSIG(IS-1)
      END DO
      !
      ! The sigma powers
      SIGPOW(:,1) = SPSIG(:)
      SIGPOW(:,2) = SPSIG(:)**2
      SIGPOW(:,3) = SPSIG(:) * SIGPOW(:,2)
      SIGPOW(:,4) = SPSIG(:) * SIGPOW(:,3)
      SIGPOW(:,5) = SPSIG(:) * SIGPOW(:,4)
      SIGPOW(:,6) = SPSIG(:) * SIGPOW(:,5)

      END SUBROUTINE
