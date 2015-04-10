CC?=g++
LD=$(CC)

nthreads ?= 48
config ?= debug
bc ?= mixed
precision ?= single
poisson ?= hypre
multiphase ?= false
bs ?= 32

CPPFLAGS+= -w

ifeq "$(bc)" "periodic"
CPPFLAGS+= -D_PERIODIC_
endif

ifeq "$(precision)" "single"
CPPFLAGS+= -D_SP_COMP_
endif

CPPFLAGS+= -D_BS_=$(bs)

ifeq "$(config)" "production"
	CPPFLAGS += -DNDEBUG -O3
else
	CPPFLAGS = -g
endif

ifeq "$(poisson)" "split-fftw"
	CPPFLAGS+= -D_SPLIT_
endif
ifeq "$(poisson)" "hypre"
	CPPFLAGS+= -D_MULTIGRID_
endif

ifeq "$(multiphase)" "true"
	CPPFLAGS+= -D_MULTIPHASE_
endif

vtk-inc ?= 
vtk-lib ?=
fftw-inc ?= 
fftw-lib ?= 
hypre-inc ?= 
hypre-lib ?=

# needed in order to compile on nodes
include make.brutus

ifneq "$(findstring brutus,$(shell hostname))" ""
include make.brutus
endif

ifneq "$(findstring euler,$(shell hostname))" ""
include make.euler
endif

CPPFLAGS += -std=c++11 -fopenmp -fpermissive
CPPFLAGS += -DNTHREADS=$(nthreads)
CPPFLAGS += -I../source/ -I../Cubism/ -I$(vtk-inc) -I$(fftw-inc)
CPPFLAGS += -D_USE_FPZIP_  -I../tools/fpzip/inc

LDFLAGS = $(CPPFLAGS)

LIBS += -lstdc++ -lm \
	-L$(vtk-lib) \
	-lvtkHybrid \
	-lvtkVolumeRendering \
	-lvtkRendering \
	-lvtkIO \
	-lvtkGenericFiltering \
	-lvtkGraphics \
	-lvtkImaging \
	-lvtkFiltering \
	-lvtkCommon \
	-lvtkftgl \
	-lvtkfreetype \
	-lvtkDICOMParser \
	-lvtkexpat \
	-lvtktiff \
	-lvtkpng \
	-lvtkjpeg \
	-lvtkzlib \
	-lvtksys -lz \
	-L$(fftw-lib) -lfftw3 -lfftw3f  -lfftw3f_threads
LIBS += -L../tools/fpzip/lib -lfpzip

ifeq "$(poisson)" "hypre"
	CC=mpic++
	LD=$(CC)
	CPPFLAGS += -I$(hypre-inc)
	LIBS += -L$(hypre-lib) -lHYPRE
endif



BUILDDIR = .
SRC_DIR = $(BUILDDIR)/../source/
CUBISM_DIR = $(BUILDDIR)/../Cubism/
OBJECTS = PoissonSolverScalarFFTW.o ProcessOperatorsOMP.o Sim_Multiphase.o Sim_FSI_Gravity.o Sim_FSI_Moving.o Sim_FSI_Fixed.o Profiler.o main.o
OBJECTS_TEST = PoissonSolverScalarFFTW.o ProcessOperatorsOMP.o TestRotation.o TestTranslation.o TestPenalization.o TestGravity.o TestVarCoeffPoisson.o TestPressure.o TestAdvection.o TestDiffusion.o Profiler.o mainTest.o

VPATH := $(SRC_DIR):$(CUBISM_DIR)

.DEFAULT: all;

all: simulation test

simulation: $(OBJECTS)
#	echo Linking
#	echo OBJECTS = $(OBJECTS)
#	echo ************************
	$(LD) $(OBJECTS) -o $@ $(LIBS) $(LDFLAGS)
#	echo ************************

test: $(OBJECTS_TEST)
#	echo Linking
#	echo OBJECTS = $(OBJECTS_TEST)
#	echo ************************
	$(LD) $(OBJECTS_TEST) -o $@ $(LIBS) $(LDFLAGS)
#	echo ************************

%.o: $(SRC_DIR)/%.cpp
#	echo PRODUCING $@
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f simulation test
	rm -f $(OBJECTS) $(OBJECTS_TEST)
