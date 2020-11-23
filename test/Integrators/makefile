CXX=g++
CFLAGS=-O3

GSLDIR=$(shell pwd)/../util/built/gsl-2.6
CUBADIR=$(shell pwd)/../util/built/cuba-4.2

ifdef GSLDIR
CFLGAS+=-static
COMMON_INC+=-I$(GSLDIR)
COMMON_LIB_DIRS+=-L$(GSLDIR)
COMMON_LIBS+=-lgsl -lgslcblas
endif

ifdef CUBADIR
CFLGAS+=
COMMON_INC+=-I$(CUBADIR)
COMMON_LIBS+=-lcuba
COMMON_LIB_DIRS+=-L$(CUBADIR)
endif

TEST_OBJ=Integrator_Tester.o

all: Test

Test: Integrator_Tester.o
	@echo $(CXX) Integrator_Tester
	@$(CXX) $(CFLAGS) -o Integrator_Tester $(COMMON_INC) $(COMMON_LIB_DIRS) -L. $(INTEGRAND_OBJ) Integrator_Tester.o $(COMMON_LIBS)

%.o:%.cpp
	@echo $(CXX) -c $< 
	@$(CXX) -c $(CFLAGS) $< -o $@ $(COMMON_INC) $(COMMON_LIB_DIRS) $(COMMON_LIBS)

clean:
	@rm -f *.o