include Makefile.inc
INCL= ##-I/usr/local/include/ -L/usr/local/lib/
VENDORINCLUDE= -Ivendor/include/
VENDORLIB= -Lvendor/lib/
LIBS= ${VENDORINCLUDE} ${VENDORLIB} ${INCL}
LD= -lcsv_parser -lgsl -lgmp -lgmpxx -lgslcblas
BUILDDIR=build

all: est

est: est.cc problem libs
	${CPPC} ${FLAGS} ${LIBS} ${LD} -o ${BUILDDIR}/est est.cc

problem: problem.c data
	${CPPC} ${FLAGS} ${LIBS} -c -o ${BUILDDIR}/problem.o problem.c ${BUILDDIR}/data.o

data: data.cc
	${CPPC} ${FLAGS} ${LIBS} -c -o ${BUILDDIR}/data.o data.cc

check: test

test: problem data
	${CPPC} ${FLAGS} -g ${LIBS} -lcheck ${LD} -o ${BUILDDIR}/test test/test.c problem.c ${BUILDDIR}/data.o
	${BUILDDIR}/test


libs:
	cd vendor/; make
