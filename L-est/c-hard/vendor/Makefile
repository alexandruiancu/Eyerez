include ../Makefile.inc

all: csv_parser.o

CSVPARSERDIR=./libcsv_parser++-1.0.0
csv_parser.o: ${CSVPARSERDIR}/include/csv_parser/csv_parser.hpp
	cd ${CSVPARSERDIR}; make