# to do anything with Docker:
#  docker build -t ocean_works_environment:build-env .

# Set the base image
FROM ubuntu


# Dockerfile maintainer (original author: Mathieu Dutour Sikiric)
MAINTAINER Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>


# Install Ubuntu packages
RUN apt-get update -y
RUN apt-get install -y libgdbm-dev libsnappy-dev make pkg-config rlwrap screen software-properties-common sudo unzip zlib1g-dev zsh
RUN apt-get install -y wget
RUN apt-get install -y g++ gcc gfortran
RUN apt-get install -y git cmake


# BOOST
RUN apt-get install -y libboost-dev
ENV BOOST_INCDIR "/usr/include"


# GMP
RUN apt-get install -y libgmp-dev
ENV GMP_INCDIR "/usr/include"
ENV GMP_CXX_LINK "-lm -lgmp -lgmpxx"


# downloading eigen (no compiling needed)
RUN mkdir -p /opt
RUN cd /opt && git clone https://gitlab.com/libeigen/eigen.git
ENV EIGEN_PATH /opt/eigen


# downloading grib and compiling it
RUN wget https://confluence.ecmwf.int/download/attachments/45757960/eccodes-2.9.2-Source.tar.gz
RUN tar -vxzf eccodes-2.9.2-Source.tar.gz && cd eccodes-2.9.2-Source && mkdir -p build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/opt/eccodes-2.9.2 && make && make install
ENV GRIB_PATH /opt/eccodes-2.9.2
ENV GRIB_FORTRAN_LINK "-L$GRIB_PATH/lib -Wl,-rpath,$GRIB_PATH/lib -leccodes_f90 -leccodes"
ENV GRIB_FORTRAN_INCDIR "-I$GRIB_PATH/include"
ENV GRIB_C_LINK "-L$GRIB_PATH/lib -leccodes -Wl,-rpath,$GRIB_PATH/lib"
ENV GRIB_INCDIR "-I$GRIB_PATH/include"


# Threading Building Blocks
RUN apt-get install -y libtbb-dev
ENV TBB_INCDIR=/usr/include
ENV TBB_LINK="-ltbb"


# Gnu Scientific library
RUN apt-get install -y libgsl-dev
ENV GSL_INCDIR=/usr/include
ENV GSL_C_LINK="-lgsl -lgslcblas"


# GLPK
RUN wget https://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz
RUN tar -vxzf glpk-4.65.tar.gz && cd glpk-4.65 && ./configure --prefix=/opt/glpk-4.65 && make && make install
ENV GLPK_PATH /opt/glpk-4.65
ENV GLPK_INCLUDE "-I$GLPK_PATH/include"
ENV GLPK_LINK="-L$GLPK_PATH/lib -lglpk -Wl,-rpath,$GLPK_PATH/lib"
RUN ls -la /opt/glpk-4.65/lib


# Downloading hdf5 and compiling it
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/src/hdf5-1.8.21.tar.gz
RUN tar -vxzf hdf5-1.8.21.tar.gz && cd hdf5-1.8.21 && ./configure --prefix=/opt/hdf5-1.8.21 && make && make install


# downloading netcdf functionality and compiling it
RUN apt-get install -y libcurl4-openssl-dev m4
RUN wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-c-4.6.2.tar.gz
RUN wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz
RUN wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-cxx4-4.3.0.tar.gz
RUN tar -vxzf netcdf-c-4.6.2.tar.gz && cd netcdf-c-4.6.2 && export LDFLAGS="-L/opt/hdf5-1.8.21/lib -lhdf5_hl -lhdf5 -Wl,-rpath,/opt/hdf5-1.8.21/lib" && export CPPFLAGS="-I/opt/hdf5-1.8.21/include" && ./configure --prefix=/opt/netcdf-4.6.2 && make && make install
RUN tar -vxzf netcdf-fortran-4.4.4.tar.gz && cd netcdf-fortran-4.4.4 && export LDFLAGS="-L/opt/netcdf-4.6.2/lib -lnetcdf -Wl,-rpath,/opt/netcdf-4.6.2/lib" && export CPPFLAGS="-I/opt/netcdf-4.6.2/include" && ./configure --prefix=/opt/netcdf-fortran-4.4.4 && make && make install
RUN tar -vxzf netcdf-cxx4-4.3.0.tar.gz && cd netcdf-cxx4-4.3.0 && export LDFLAGS="-L/opt/netcdf-4.6.2/lib -lnetcdf -Wl,-rpath,/opt/netcdf-4.6.2/lib" && export CPPFLAGS="-I/opt/netcdf-4.6.2/include" && ./configure --prefix=/opt/netcdf-cxx4-4.3.0 && make && make install
ENV NETCDF_CXX_PATH /opt/netcdf-cxx4-4.3.0
ENV NETCDF_C_PATH /opt/netcdf-4.6.2
ENV NETCDF_CXX_ALLINC "-I$NETCDF_CXX_PATH/include -I$NETCDF_C_PATH/include"
ENV NETCDF_CXX_LINK "-L$NETCDF_CXX_PATH/lib -lnetcdf_c++4 -Wl,-rpath,$NETCDF_CXX_PATH/lib"


# shapelib
RUN apt-get install -y libproj-dev
RUN wget http://download.osgeo.org/shapelib/shapelib-1.4.1.tar.gz
RUN tar -vxzf shapelib-1.4.1.tar.gz && cd shapelib-1.4.1 && ./configure --prefix=/opt/shapelib-1.4.1 && make && make install
ENV SHAPELIB_PATH "/opt/shapelib-1.4.1"
ENV SHAPELIB_INCLUDE "-I$SHAPELIB_PATH/include"
ENV SHAPELIB_LINK "$SHAPELIB_PATH/lib/libshp.a"


# cloning ocean works and compiling it
ENV STANDARD_FORTRAN_FC "gfortran"
ENV STANDARD_FORTRAN_LIBS "-lgfortran"
RUN mkdir GIT
RUN cd GIT && git clone https://github.com/MathieuDutSik/ocean_works.git --recursive
RUN cd GIT/ocean_works/src_ocean && make

