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
      
      END SUBROUTINE
