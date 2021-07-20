PROGRAM	= Compute_Seiche_Information

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

CLOCAL          = -I../src_basic -I../src_matrix -I../src_comb -I../src_graph -I../src_svg -I../src_ocean -I../src_poly


CFLAGS 		= -I$(EIGEN_PATH) -I$(SPECTRA_PATH) $(NETCDF_CXX_ALLINC) $(GRIB_INCDIR) $(SHAPELIB_INCLUDE) $(GLPK_INCLUDE) -I$(BOOST_INCDIR) $(OPENCV_KDTREE_INCLUDE)

LIBS		= -lm -lpthread $(STANDARD_FORTRAN_LIBS) $(NETCDF_CXX_LINK) $(GRIB_C_LINK) $(SHAPELIB_LINK) $(GLPK_LINK) $(OPENCV_KDTREE_LINK) $(GMP_CXX_LINK)

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
