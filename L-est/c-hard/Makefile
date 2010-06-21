include Makefile.inc
VENDORINCLUDE= -Ivendor/include/
VENDORLIB= -Lvendor/lib/
LIBS= ${VENDORINCLUDE} ${VENDORLIB} -lcsv_parser -lgmp -lgmpxx

all: est

est: est.cc libs
	${CC} ${FLAGS} ${LIBS} -o est est.cc

libs:
	cd vendor/; make