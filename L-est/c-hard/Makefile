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

problem: problem.c
	${CC} ${FLAGS} -c -o ${BUILDDIR}/problem.o problem.c

check: test

test: problem
	${CC} ${FLAGS} -g ${LIBS} -lcheck ${LD} -o ${BUILDDIR}/test test/test.c problem.c
	${BUILDDIR}/test


libs:
	cd vendor/; make
