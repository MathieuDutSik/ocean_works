#include ${PETSC_DIR}/lib/petsc/conf/variables
PROG = interpolate_spectrum

FFLAGS = -DST42 -DWWM_SETUP


#LIBS += -L/opt/mpi/openmpi-icc/4.0.3.3/lib -Wl,-rpath,/opt/mpi/openmpi-icc/4.0.3.3/lib
#LIBS += -L/opt/software/intel/parallel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin -Wl,-rpath,/opt/software/intel/parallel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin

F90 = gfortran

#GFORTRAN
F90OPTI = -g -O3 -ffree-line-length-0
#F90OPTI = -g -O3 -ffree-line-length-0 -fbacktrace -Wall
#F90OPTI = -g -O3 -ffree-line-length-0 -fbacktrace -Wall -fcheck=all
#F90OPTI = -O3 -ffree-line-length-0
#F90OPTI    = -Wall -O0 -g -ggdb -ffree-line-length-none -fcheck=all -fsanitize=thread -fbacktrace -fstack-protector-all -ffpe-trap=invalid,zero,overflow
#F90OPTI    = -Wall -O0 -g -ggdb -ffree-line-length-none -fcheck=all -fsanitize=address -fsanitize=undefined   -fbacktrace -fstack-protector-all -ffpe-trap=invalid,zero,overflow
#F90OPTI    = -Wall -O3 -g -ggdb -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -DDEBUG
#F90OPTI    = -Wall -O3 -g -ggdb -ffree-line-length-none -ffpe-trap=invalid,zero,overflow
#F90OPTI    = -Wall -O9 -ffree-line-length-none -g
#F90OPTI    = -g -Wall -O9 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow


OBJS = wwm_datapl.o wwm_aux.o wwm_mjdv2.o wwm_grib.o wwm_bdcons_wam.o interpolate_spectrum.o

FFLAGS += -DWWM_SOLVER
all: $(PROG)

LIBS += $(GRIB_FORTRAN_LINK) $(NETCDF_FORTRAN_LINK)
FFLAGS += $(GRIB_FORTRAN_INCDIR) -I$(NETCDF_INCDIR)



F90FLAGS = ${FFLAGS} ${F90OPTI}

NETC = NETC

$(NETC):
	cp -f $(NETCDF_INCDIR)/netcdf.mod .
	cp -f $(NETCDF_INCDIR)/typesizes.mod .

$(PROG): $(OBJS)
	$(F90) -o $(PROG) $(F90FLAGS)  $(OBJS) $(LIBS)

clean:
	rm -f *.o *.oo *.obj *.ipo *.mod *.map *__genmod.f90 *.ilk *.pdb $(PROG)

cleanall:
	rm -f *.*~ *.spg *.smb *.o *.oo *.obj *.ipo *.mod *.map *.ilk *.pdb *genmod* $(PROG)

.SUFFIXES: $(SUFFIXES) .F90 .f .ftn

.F90.o:
	$(F90) $(F90FLAGS)   -c $<

.ftn.o:
	$(F90) $(F90FLAGS) -c $<
