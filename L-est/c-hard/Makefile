include Makefile.inc
INCL= ##-I/usr/local/include/ -L/usr/local/lib/
VENDORINCLUDE= -Ivendor/include/
VENDORLIB= -Lvendor/lib/
LIBS= ${VENDORINCLUDE} ${VENDORLIB} ${INCL}
LD= -lcsv_parser -lgsl -lgmp -lgmpxx -lgslcblas
BUILDDIR=build

all: mcmc.o

clean: 
	rm -r build/*

# est: est.cc problem libs
# 	${CPPC} ${FLAGS} ${LIBS} ${LD} -o build/est est.cc

problem.o: problem.c problem.h data.o
	${CPPC} ${FLAGS} ${LIBS} -c -o build/problem.o problem.c build/data.o

data.o: data.cc data.h
	${CPPC} ${FLAGS} ${LIBS} -c -o build/data.o data.cc

mcmc.o: mcmc.cc mcmc.h data.o problem.o
	${CPPC} ${FLAGS} ${LIBS} -c -o build/mcmc.o mcmc.cc

check: test

test: problem.o data.o mcmc.o
	${CPPC} ${FLAGS} -g ${LIBS} -lcheck ${LD} -o build/test test/test.c build/*.o
	build/test


libs:
	cd vendor/; make
