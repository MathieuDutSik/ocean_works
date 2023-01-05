PROGRAM	= MEAS_AltimeterComparison \
	  MEAS_ComparisonAltimetrySource \
	  MEAS_AquaOperation \
	  MEAS_SST_Comparison \
	  MEAS_CTD_Comparison \
	  GRID_ScaleXYgrid \
	  GRID_CompWWMboundary \
	  GRID_SetStraightLine \
	  GRID_ShowGridInformation \
	  GRIB_FindDisturbanceSequence \
	  GRID_ConvertGrid \
	  GRID_ConvertGridBound \
	  DATA_ComputeRomsDensity \
	  GRID_ComputeWindRot2grid \
	  GRIB_PrintSequenceFiles \
	  GRID_CompGridTimeStepCFL \
	  GRID_SetBathymetry \
	  RIVER_MakeRiverFile \
	  RIVER_ShowRiver \
	  RIVER_PlotRiverOnGrid \
	  RIVER_MergeRiverFiles \
	  PLOT_results \
	  PLOT_bound \
	  PLOT_buoy \
	  PLOT_float \
	  PLOT_diff_results \
	  PLOT_grid \
	  PLOT_PointOutput \
	  PLOT_transect \
	  ROMS_VerticalStratification_Diagnostics \
	  ROMS_SetTracerInitial \
	  ROMS_CreateDefaultSetting \
	  ROMS_SetTime \
	  ROMS_SetTimeRomsInitial \
	  DATA_CreateIdealCurrent \
	  DATA_CreateSflux \
	  DATA_MergeField \
	  DATA_AverageField \
	  CHECK_HYCOM


#PROGRAM = COMP_drifters DATA_CreateLTransInput DATA_CreateRomsOfflineInput PLOT_drifters



#	  CREATE_TracerSourceTerm.cpp \ Deleted file, maybe interesting.

FILE_F90 = WW3_minimalroutine


OBJECTS		= $(PROGRAM:%=%.o)

OBJECTS_F90     = $(FILE_F90:%=%.o)


#CC		= icpc -std=c++17 -Wall -g -O2 
#CC              = clang++ -Weverything -std=c++17
#CC              = clang++ -Weverything -std=c++17 -fsanitize=address
#CC              = clang++ -Weverything -std=c++17 -fsanitize=memory
#CC		= g++ -std=c++17 -Wall -g -O2 
#CC		= g++ -std=c++17 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fstack-protector-all
#CC		= g++ -std=c++17 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fstack-protector-all
CC		= g++ -std=c++17 -Wall -g -O2
#FC		= gfortran -Wall -g -O1
#FC		= gfortran -Wall -g -O1 -fno-underscoring
FC		= $(STANDARD_FORTRAN_FC)

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_svg


CFLAGS 		= -I$(EIGEN_PATH) $(NETCDF_CXX_ALLINC) $(GRIB_INCDIR) $(SHAPELIB_INCLUDE) -I$(BOOST_INCDIR) $(OPENCV_KDTREE_INCLUDE)

LIBS		= -lm -lpthread $(STANDARD_FORTRAN_LIBS) $(NETCDF_CXX_LINK) $(GRIB_C_LINK) $(SHAPELIB_LINK) $(OPENCV_KDTREE_LINK) $(FORTRAN_LINK) $(GMP_CXX_LINK)

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
