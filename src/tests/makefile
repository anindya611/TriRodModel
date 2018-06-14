ALL: 	testPointKinematics \
        newTestPointKinematics \
	testUtils \
	testConfiguration \
	testOps \
	testAssembly

testPointKinematics: testPointKinematics.o
	${CLINKER} ${CFLAGS} ${CPPFLAGS} -o $@ $^ ${LIBS} 
	${RM} $^

newTestPointKinematics: newTestPointKinematics.o
	${CLINKER} ${CFLAGS} ${CPPFLAGS} -o $@ $^ ${LIBS} 
	${RM} $^

Clean:
	-rm -f testPointKinematics newTestPointKinematics 
	-rm -f *.o
	-rm -rf *.dSYM


CPPFLAGS = 	-std=c++11 -O2 \
		-I${DGPLUS_HOME}/include -L${DGPLUS_HOME}/lib/${PETSC_ARCH} \
		-I${TRIRODDIR}/include -L${TRIRODDIR}/lib/${PETSC_ARCH} \
		-I${GSLPATH}/include -L${GSLPATH}/lib \
		-L${PETSC_LIB_DIR}

LIBS =		-lTriRods -lElm -lElOp  \
		${BLASLAPACK_LIB} ${PETSC_KSP_LIB} -lgsl -lgslcblas

LOCDIR   = ${CURDIR}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test
