ifeq ($(COMPILER),intel)
  FC_AUTO_R8 :=  -r8 
  MPIFC :=  mpiifort 
  FFLAGS_NOOPT :=  -O0 
  CXX_LDFLAGS :=  -cxxlib 
  SUPPORTS_CXX := TRUE
  LAPACK_LIBDIR := /usr/lib64
  FFLAGS :=  -qno-opt-dynamic-align -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source 
  FIXEDFLAGS :=  -fixed 
  SCC :=  icc 
  SFC :=  ifort 
  SLIBS := -L/opt/hpc/software/mathlib/netcdf/4.4.1/intel/lib -lnetcdf -lnetcdff -L/public/home/chengxl/lib -llapack -L/public/home/chengxl/lib -lblas -mkl
  MPICC :=  mpiicc 
  MPI_PATH := /public/software/mpi/intelmpi/2017.4.239/intel64
  CFLAGS :=  -qno-opt-dynamic-align -fp-model precise -std=gnu99 
  MPICXX :=  mpiicpc 
  FREEFLAGS :=  -free 
  CXX_LINKER := FORTRAN
  SCXX :=  icpc 
endif
ifeq ($(MODEL),pop)
  CPPDEFS := $(CPPDEFS)  -D_USE_FLOW_CONTROL 
endif
ifeq ($(COMPILER),intel)
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DCPRINTEL
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DCPRINTEL
  ifeq ($(compile_threaded),true)
    FFLAGS := $(FFLAGS)  -qopenmp 
    CFLAGS := $(CFLAGS)  -qopenmp 
  endif
  ifeq ($(DEBUG),TRUE)
    FFLAGS := $(FFLAGS)  -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created 
    CFLAGS := $(CFLAGS)  -O0 -g 
  endif
  ifeq ($(DEBUG),FALSE)
    FFLAGS := $(FFLAGS)  -O2 -debug minimal 
    CFLAGS := $(CFLAGS)  -O2 -debug minimal 
  endif
  ifeq ($(MPILIB),mvapich2)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpich2)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpt)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),openmpi)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpich)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mvapich)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),impi)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpi-serial)
    SLIBS := $(SLIBS)  -mkl 
  endif
  ifeq ($(compile_threaded),true)
    FFLAGS_NOOPT := $(FFLAGS_NOOPT)  -qopenmp 
    LDFLAGS := $(LDFLAGS)  -qopenmp 
    LDFLAGS := $(LDFLAGS)  -qopenmp 
  endif
endif
