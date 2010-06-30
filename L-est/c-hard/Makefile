
COREDIR= .
VENDORDIR=$(COREDIR)/vendor
TESTDIR=$(COREDIR)/test
BD=$(COREDIR)/build
CLD=$(COREDIR)/classes
CLBD=$(BD)/classes

CC=g++
FLAGS=-O3 -Wall
INCL=#-I/usr/local/include/ -L/usr/local/lib/
LD=-lcsv_parser -lgsl -lgslcblas
VENDOR=-I$(VENDORDIR)/include/ -L$(VENDORDIR)/lib/
LIBS= $(VENDOR) $(INCL)

GTEST_DIR=$(VENDORDIR)/gtest-1.5.0
GTEST_INCL=$(LIBS) -I$(GTEST_DIR)/include
GTEST_FLAGS= $(FLAGS) -g -Wextra
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

TESTS=metropolis_test.o problem_def_test.o

#-----------------------------------------------------------------------------

all: 

clean: clean-tests
	rm -f $(BD)/**/*.a
	rm -f $(BD)/**/*.o

CLASS_NAMES=metropolis.o variance.o random.o problem.o observations.o
CLASSES=$(addprefix $(CLBD)/, $(CLASS_NAMES))

$(CLBD)/metropolis.o: $(CLD)/metropolis.cc $(CLD)/metropolis.hpp
	$(CC) $(FLAGS) -I$(COREDIR) $(LIBS) -c \
	-o $(CLBD)/metropolis.o $(CLD)/metropolis.cc

$(CLBD)/variance.o: $(CLD)/variance.cc $(CLD)/variance.hpp
	$(CC) $(FLAGS) -I$(COREDIR) $(LIBS) -c \
	-o $(CLBD)/variance.o $(CLD)/variance.cc

$(CLBD)/random.o: $(CLD)/random.cc $(CLD)/random.hpp
	$(CC) $(FLAGS) -I$(COREDIR) $(LIBS) -c \
	-o $(CLBD)/random.o $(CLD)/random.cc

$(CLBD)/problem.o: $(CLD)/problem.cc $(CLD)/problem.hpp
	$(CC) $(FLAGS) -I$(COREDIR) $(LIBS) -c \
	-o $(CLBD)/problem.o $(CLD)/problem.cc

$(CLBD)/observations.o: $(CLD)/observations.cc $(CLD)/observations.hpp
	$(CC) $(FLAGS) -I$(COREDIR) $(LIBS) -c \
	-o $(CLBD)/observations.o $(CLD)/observations.cc

$(BD)/classes.a: $(CLASSES)
	$(AR) $(ARFLAGS) $@ $^

#-----------------------------------------------------------------------------

metropolis_test.o: $(TESTDIR)/metropolis_test.cc metropolis.o classes.a
	$(CC) $(GTEST_FLAGS) $(LIBS) -I$(CLD) -c $(TESTDIR)/metropolis_test.cc \
	  -o $(TESTDIR)/$@ \
	  $(CLBD)/metropolis.o

problem_def_test.o: $(TESTDIR)/problem_def_test.cc problem.o mathy.a
	$(CC) $(GTEST_FLAGS) $(LIBS) -I$(CLD) -c $(TESTDIR)/problem_def_test.cc \
	  -o $(TESTDIR)/$@ \
	  $(BD)/problem.o $(BD)/data.o

run_tests: $(TESTS) gtest_main.a
	$(CC) $(GTEST_FLAGS) $(LIBS) $(LD) -lpthread \
	  -o $(TESTDIR)/$@ \
	  $(addprefix $(TESTDIR)/, $(TESTS)) $(TESTDIR)/gtest_main.a \
	  $(BD)/mathy.a $(BD)/classes.a

#-----------------------------------------------------------------------------

check: $(TESTDIR)/test.cc $(BD)/classes.a
	$(CC) $(GTEST_FLAGS) $(LIBS) $(LD) -I$(COREDIR) -o $(TESTDIR)/test_me \
	  $(TESTDIR)/test.cc \
	  $(BD)/classes.a

	$(TESTDIR)/test_me > out.csv

# check: mathy.a classes.a run_tests
# 	$(TESTDIR)/run_tests

clean-tests:
	rm -f $(TESTDIR)/*.a
	rm -f $(TESTDIR)/*.o
	rm -f $(TESTDIR)/run_tests

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

gtest-all.o: $(GTEST_SRCS_)
	$(CC) $(FLAGS) -I$(GTEST_DIR) $(LIBS) -c \
	  $(GTEST_DIR)/src/gtest-all.cc \
	  -o $(TESTDIR)/$@

gtest_main.o: $(GTEST_SRCS_)
	$(CC) $(FLAGS) -I$(GTEST_DIR) $(LIBS) -c \
          $(GTEST_DIR)/src/gtest_main.cc \
	  -o $(TESTDIR)/$@

gtest.a: gtest-all.o
	$(AR) $(ARFLAGS) $(TESTDIR)/$@ $(TESTDIR)/$^

gtest_main.a: gtest_main.o
	$(AR) $(ARFLAGS) $(TESTDIR)/$@ $(addprefix $(TESTDIR)/, $^)

#-----------------------------------------------------------------------------

# recursive, unfortunately
vendor:
	cd vendor/; make
