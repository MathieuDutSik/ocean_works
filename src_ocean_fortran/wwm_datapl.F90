!**********************************************************************
!*                                                                    *
!**********************************************************************
      MODULE DATAPOOL
        IMPLICIT NONE
        integer,parameter :: rkind = 8      ! Default real datatype
        integer NUMSIG, NUMDIR
        integer istat
        REAL(rkind), PARAMETER             :: ONEHALF   = 0.5_rkind
        REAL(rkind), PARAMETER             :: ZERO   = 0.0_rkind
        REAL(rkind), PARAMETER             :: TWO   = 2.0_rkind
        REAL(rkind), PARAMETER             :: PI        = 3.141592653589793_rkind
        REAL(rkind), PARAMETER             :: PIHALF    = PI*ONEHALF
        REAL(rkind), PARAMETER             :: PI2       = TWO*PI
        REAL(rkind), PARAMETER             :: DEGRAD    = PI/180._rkind
        REAL(rkind), PARAMETER             :: RADDEG    = 180._rkind/PI
        REAL(rkind), PARAMETER             :: SMALL     = 10E-7
        REAL(rkind), PARAMETER             :: LARGE     = 1./SMALL
        REAL(rkind), PARAMETER             :: THR       = TINY(1.)
        REAL(rkind),  PARAMETER            :: DAY2SEC  = 86400.d0
        REAL(rkind),  PARAMETER            :: SEC2DAY  = 1.d0/DAY2SEC
        REAL(rkind)  :: FRLOW, FRHIGH, FRATIO
        REAL(rkind)  :: ddir
        real(rkind) DELT25_WAM
        REAL(rkind), ALLOCATABLE      :: SPSIG(:)
        REAL(rkind), ALLOCATABLE      :: SIGPOW(:,:)
        REAL(rkind), ALLOCATABLE      :: SPDIR(:)
        REAL(rkind), ALLOCATABLE      :: FR(:)
        character(len=200) :: PrefixOutput
        LOGICAL EXTRAPOLATION_ALLOWED_BOUC
        !
        CHARACTER(LEN=140)     :: PREFIX_WAVE_FILE = ''
        INTEGER                :: NUM_WAM_SPEC_FILES
        real(rkind), allocatable :: WAM_SPEC_ListTime(:)
        character(len=140), allocatable :: WAM_SPEC_FILE_NAMES_BND(:)
        integer, allocatable :: ListIFileWAM(:)
        REAL(rkind)                   :: TAIL_ARR(8)
        REAL(rkind), ALLOCATABLE      :: DS_INCR(:)
        REAL(rkind), ALLOCATABLE      :: DS_BAND(:)
        integer, allocatable :: WAM_ID1(:), WAM_ID2(:), WAM_IS1(:), WAM_IS2(:)
        real(rkind), allocatable :: ListDir_wam(:), ListFreq_wam(:)
        real(rkind), allocatable :: DFIM_wam(:)
        real(rkind), allocatable :: WAM_WD1(:), WAM_WD2(:), WAM_WS1(:), WAM_WS2(:)
        integer, allocatable :: CF_IX_BOUC(:)
        integer, allocatable :: CF_IY_BOUC(:)
        real(rkind), allocatable :: CF_COEFF_BOUC(:,:)
        integer IWBMNP
        real(rkind), allocatable :: XP(:), YP(:)
        integer nbdir_wam, nbfreq_wam, nx_wam, ny_wam

        TYPE FILEDEF
           CHARACTER(LEN=140)  :: FNAME
           INTEGER             :: FHNDL
        END TYPE FILEDEF
        TYPE (FILEDEF)         :: STAT, WAV, INP, NODES, OUT

        TYPE FD_FORCING_GRID
           integer nx_dim, ny_dim
           real(rkind), dimension(:,:), pointer :: LON
           real(rkind), dimension(:,:), pointer :: LAT
        END TYPE FD_FORCING_GRID

        TYPE VAR_TIME
           integer nbTime
           real(rkind), allocatable :: ListTime(:)
        END TYPE VAR_TIME
        TYPE(VAR_TIME) :: eVAR_BOUC_WAM

      END MODULE DATAPOOL
