alltree: lib

CPPFLAGS = 	-std=c++11 -O2 \
		-I${DGPLUS_HOME}/include \
		-I${TRIRODDIR}/include \
		-I${GSLPATH}/include

SOURCECXX = 	TriRodsConfigurationNew1.cpp  \
		TriRodsProperties.cpp \
		TriRodsUtils.cpp \
                TriRodsOps.cpp

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

INSTALL_LIB_DIR = ${TRIRODDIR}/lib/${PETSC_ARCH}/
LIBBASE  = libTriRods
LOCDIR   = ${CURDIR}
