PROGRAM	= AltimeterComparison \
	  ComparisonAltimetrySource \
	  PLOT_diff_results \
	  CREATE_sflux \
	  MERGE_field \
	  ScaleXYgrid \
	  COMP_GridTimeStepCFL \
	  MERGE_RiverFiles \
	  PLOT_grid \
	  PLOT_PointOutput \
	  PLOT_transect \
	  ComputeWWMboundary \
	  ROMS_SetTracerInitial \
	  ROMS_CopyROMSfields \
	  ROMS_CreateDefaultSetting \
	  SetStraightLine \
	  IDEAL_current \
	  ROMS_SetTime \
	  SET_TIME_RomsInitial \
	  SHP_GetShapeType \
	  SHP_PrintShapeType \
	  COAST_ConvertCoastlineFormat \
	  COAST_ConvertMeshToCoast \
	  DBF_PrintDatabase \
	  GRIB_FindDisturbanceSequence \
	  ShowGridInformation \
	  ConvertGrid \
	  CHECK_HYCOM \
	  ConvertGridBound \
	  ComputeWindRot2grid \
	  GRIB_PrintSequenceFiles \
	  SetBathymetry \
	  AquaOperation \
	  MakeRiverFile \
	  PLOT_results \
	  PLOT_buoy \
	  SHOW_river \
	  PLOT_river \
	  PLOT_bound


#PROGRAM = COMP_drifters CREATE_LTransInput CREATE_RomsOfflineInput PLOT_drifters



#	  CREATE_TracerSourceTerm.cpp \ Deleted file, maybe interesting.

FILE_F90 = WW3_minimalroutine


OBJECTS		= $(PROGRAM:%=%.o)

OBJECTS_F90     = $(FILE_F90:%=%.o)


#CC		= icpc -std=c++11 -Wall -g -O2 
#CC              = clang++ -Weverything -std=c++11
#CC              = clang++ -Weverything -std=c++11 -fsanitize=address
#CC              = clang++ -Weverything -std=c++11 -fsanitize=memory
#CC		= g++ -std=c++14 -Wall -g -O2 
#CC		= g++ -std=c++14 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fstack-protector-all
#CC		= g++ -std=c++11 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fstack-protector-all
CC		= g++ -std=c++11 -Wall -g -O2
#FC		= gfortran -Wall -g -O1
#FC		= gfortran -Wall -g -O1 -fno-underscoring
FC		= $(STANDARD_FORTRAN_FC)

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_matrix -I../src_comb -I../src_graph -I../src_svg


CFLAGS 		= -I$(EIGEN_PATH) $(NETCDF_CXX_ALLINC) $(GRIB_INCDIR) $(SHAPELIB_INCLUDE) $(GLPK_INCLUDE) -I$(BOOST_INCDIR) $(OPENCV_KDTREE_INCLUDE)

LIBS		= -lm -lpthread $(STANDARD_FORTRAN_LIBS) $(NETCDF_CXX_LINK) $(GRIB_C_LINK) $(SHAPELIB_LINK) $(GLPK_LINK) $(OPENCV_KDTREE_LINK) $(FORTRAN_LINK) $(GMP_CXX_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

%.o: %.F90
	@echo
	@echo compile $<
	$(FC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS) $(OBJECTS_F90)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(OBJECTS_F90) $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
