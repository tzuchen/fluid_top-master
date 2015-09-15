
#
#requirespackage  'PETSC_HAVE_MATLAB'
#requiresscalar real 
#requiresprecision  double
#
#   matlabcodes makes the Matlab mex files for receiving PETSc objects
# in Matlab. You cannot use a complex numbers. If mex fails 
# in compiling these files and you are using the gcc compiler then try 
# removing the -fPIC from the CC declaration in ${PETSC_DIR}/bmake/${PETSC_ARCH}/variables
#
SPECIALLIB     = yes
SPECIALFASTLIB = yes
speciallib:
	-@echo "BEGINNING TO COMPILE MATLAB INTERFACE"
	-@${OMAKE}  matlabcodes PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR}
	-@if [ ! -d "${PETSC_DIR}/bin/matlab" ] ; then ${MKDIR}  ${PETSC_DIR}/bin/matlab; fi
	-@${MV} send.mex*  ${PETSC_DIR}/bin/matlab
specialfastlib: speciallib

matlabcodes: send

CFLAGS    =
FFLAGS    =
CPPFLAGS  = 
SOURCEC   = 
SOURCEF   =
SOURCEH   = send.c send2file.c bread.c
OBJSC     = 
OBJSF     =
MANSEC    = Viewer
LIBBASE   = libpetsc
DIRS      = 
#LOCDIR    = src/sys/viewer/impls/socket/matlab/
LOCDIR    = home/tzuchen/proj/fluid_top/trunk/matlab_send/
NOADIC    = true

# The following parameters are added by Tzuchen and may change
PETSC_DIR = /home/tzuchen/petsc/petsc-2.3.2-p10
PETSC_ARCH = linux-gnu-c-debug



include ${PETSC_DIR}/bmake/common/base
include ${PETSC_DIR}/bmake/common/test

MATLAB_MEX = /opt/matlab/bin/mex
MATLAB_CC  = gcc32

matlabcodes 	= send send2file

# MATLAB_MEX and MATLAB_CC are defined in bmake/${PETSC_ARCH}/packages
#
send:  
	-@${MATLAB_MEX} -g CC=${MATLAB_CC} CFLAGS='${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS}' send.c bread.c 
	-@${RM} -f send.o 
	-@if [ ${PETSC_ARCH} = IRIX64 ] ; then mv -f send.mexsg64 send.mexsg ; fi 

send2file:  
	-@${MATLAB_MEX} -g CC=${MATLAB_CC} CFLAGS='${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS}' send2file.c bread.c 
	-@${RM} -f send2file.o 
	-@if [ ${PETSC_ARCH} = IRIX64 ] ; then mv -f send2file.mexsg64 send2file.mexsg ; fi 

