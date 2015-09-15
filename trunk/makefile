
CFLAGS	         = 
FFLAGS	         = 
CPPFLAGS         =
FPPFLAGS         =
#HOME_DIR         = /home/tzuchen
LOCDIR           = ${HOME_DIR}/proj/fluid_top/trunk
#setting of these dirs is done in shell script
#PETSC_DIR        = ${HOME_DIR}/petsc/petsc-2.3.2-p10
#SLEPC_DIR        = ${HOME_DIR}/slepc-2.3.2
#FFTW_DIR         = ${HOME_DIR}/fftw-3.1.2/.libs/

FFTW_FLAG        = -L${FFTW_DIR} -lfftw3 
#FFTW_INCLUDE     = -I/home/tzuchen/fftw/include/


MANSEC           = KSP


FLUIDSOLVER_OBJS     = fluidsolver.o PetscReceive.o
FLUIDSOLVER2_OBJS    = fluidsolver2.o PetscStreamline.o PetscReceive.o MarkovMapFunction.o 
FLUIDOPT_OBJS        = fluidopt.o PetscStreamline.o PetscReceive.o MarkovMapFunction.o DataProcess.o 
STREAMLINE_OBJS      = streamline.o PetscStreamline.o PetscReceive.o MarkovMapFunction.o
STREAMLINER_OBJS     = streamliner.o PetscStreamline.o PetscReceive.o MarkovMapFunction.o
TESTMARKOVMAP_OBJS   = testmarkovmap.o PetscStreamline.o PetscReceive.o MarkovMapFunction.o 
TESTMARKOVMAP2_OBJS  = testmarkovmap2.o PetscStreamline.o PetscReceive.o MarkovMapFunction.o 
TASKMANAGER_OBJS     = taskmanager.o PetscStreamline.o PetscReceive.o MarkovMapFunction.o
MAPREFINE_OBJS       = Maprefine.o PetscReceive.o FFTRoutine.o IterationFunction.o LargeVecFunction.o
MAPREFINE2_OBJS      = Maprefine2.o PetscReceive.o FFTRoutine.o IterationFunction.o LargeVecFunction.o
MAPREFINE3_OBJS      = Maprefine3.o PetscReceive.o FFTRoutine.o

PROGS            = fluidsolver fluidsolver2 streamline streamliner testmarkovmap testmarkovmap2 taskmanager Maprefine fluidopt

all: ${PROGS}

#include ${PETSC_DIR}/bmake/common/base
#Slepc part
include ${SLEPC_DIR}/bmake/slepc_common

#include ${FFTW_DIR}
	
#%.o: %.c
#	gcc -c -o $@ $<
#CLINKER =  mpicc -fPIC -Wall -Wwrite-strings -fomit-frame-pointer -Wno-strict-aliasing -g3

 
fluidsolver: ${FLUIDSOLVER_OBJS}
	     ${CLINKER} -o fluidsolver ${FLUIDSOLVER_OBJS} ${PETSC_KSP_LIB} ${SLEPC_LIB}

    
fluidsolver2: ${FLUIDSOLVER2_OBJS}
	     ${CLINKER} -o fluidsolver2 ${FLUIDSOLVER2_OBJS} ${PETSC_KSP_LIB} ${SLEPC_LIB}

fluidopt:    ${FLUIDOPT_OBJS}
	     ${CLINKER} -o fluidopt ${FLUIDOPT_OBJS} ${PETSC_KSP_LIB} ${SLEPC_LIB}

streamline: ${STREAMLINE_OBJS} 
	    ${CLINKER} -o streamline ${STREAMLINE_OBJS} ${PETSC_KSP_LIB} ${SLEPC_LIB} 

streamliner: ${STREAMLINER_OBJS} 
	     ${CLINKER} -o streamliner ${STREAMLINER_OBJS} ${PETSC_KSP_LIB} ${SLEPC_LIB} 

testmarkovmap: ${TESTMARKOVMAP_OBJS} 
	       ${CLINKER} -o testmarkovmap ${TESTMARKOVMAP_OBJS} ${PETSC_KSP_LIB} ${SLEPC_LIB}

testmarkovmap2: ${TESTMARKOVMAP2_OBJS} 
	       ${CLINKER} -o testmarkovmap2 ${TESTMARKOVMAP2_OBJS} ${PETSC_KSP_LIB} ${SLEPC_LIB}

taskmanager: ${TASKMANAGER_OBJS}
	     ${CLINKER} -o taskmanager ${TASKMANAGER_OBJS} ${PETSC_KSP_LIB} ${SLEPC_LIB}

Maprefine: ${MAPREFINE_OBJS}
	   ${CLINKER} -o Maprefine ${MAPREFINE_OBJS} ${PETSC_KSP_LIB} ${FFTW_FLAG} ${SLEPC_LIB}

#Maprefine2: ${MAPREFINE2_OBJS}
#	   ${CLINKER} -o Maprefine2 ${MAPREFINE2_OBJS} ${PETSC_KSP_LIB} ${FFTW_FLAG} ${SLEPC_LIB}  

#Maprefine3: ${MAPREFINE3_OBJS}
#	   ${CLINKER} -o Maprefine3 ${MAPREFINE3_OBJS} ${PETSC_KSP_LIB} ${FFTW_FLAG} ${SLEPC_LIB}

#clean:
#	rm -f *.o ${PROGS}



#----------------------------------------------------------------------------
#runex2:
#	-@${MPIRUN} -np 1 ./ex2 -ksp_smonitor -m 5 -n 5 -ksp_gmres_cgs_refinement_type refine_always > ex2_1.tmp 2>&1; \
#	   if (${DIFF} output/ex2_1.out ex2_1.tmp) then true; \
#	   else echo "Possible problem with ex2_1, diffs above"; fi; \
#	   ${RM} -f ex2_1.tmp
#
#

#include ${PETSC_DIR}/bmake/common/test
