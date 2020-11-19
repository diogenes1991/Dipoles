CXX=g++
CFLAGS=-O3

LHAPDF_DIR=/home/diogenes1991/LHAPDF-6.3.0/include/LHAPDF

COMMON_INC=
COMMON_LIBS=
COMMON_LIB_DIRS=	

ifdef LHAPDF_DIR
COMMON_INC+=-I$(LHAPDF_DIR)
COMMON_LIBS+=-lLHAPDF
COMMON_LIB_DIRS+=-L$(LHAPDF_DIR)
endif

TEST_OBJ=PDF_Tester.o

all: Test

Test: $(TEST_OBJ)
	@echo $(CXX) Integrator_Tester
	@$(CXX) $(CFLAGS) -o PDF_Tester $(COMMON_INC) $(COMMON_LIB_DIRS) -L. $(TEST_OBJ) $(COMMON_LIBS)

%.o:%.cpp
	@echo $(CXX) -c $< 
	@$(CXX) -c $(CFLAGS) $< -o $@ $(COMMON_INC) $(COMMON_LIB_DIRS) $(COMMON_LIBS)

clean:
	@rm -f *.o