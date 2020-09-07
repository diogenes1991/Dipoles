CXX=g++
FC=gfortran

#PROCPATH=./
PROCPATH=$(shell pwd)/

CFLAGS= -O3 -Wall -Wno-reorder -DNLOXPROCPATH=\"$(PROCPATH)\"
FFLAGS= -O3

####TREDDIR

####ONELOOPDIR
####QCDLOOPDIR
####QCDLOOP2DIR
####LOOPTOOLSDIR

COMMON_INC=-I$(TREDDIR) \
    -I$(PROCPATH)/code
COMMON_LIB_DIRS=-L$(TREDDIR)
COMMON_LIBS=-Wl,-rpath,$(TREDDIR) -ltred -lquadmath -lgcc -lstdc++ -lm

ifdef QCDLOOPDIR
CFLAGS+=-DQCDLOOP
COMMON_LIB_DIRS+=-L$(QCDLOOPDIR)
COMMON_LIB_DIRS+=-L$(QCDLOOPDIR)/ql
COMMON_LIB_DIRS+=-L$(QCDLOOPDIR)/ff
COMMON_LIB_DIRS+=-L$(QCDLOOPDIR)/lib
COMMON_LIB_DIRS+=-L$(QCDLOOPDIR)/lib64
COMMON_LIB_DIRS+=-L$(QCDLOOPDIR)/lib/NLOX_util
COMMON_LIBS+=-lqcdloop -lff
fortran_needed = yes
endif

ifdef LOOPTOOLSDIR
CFLAGS+=-DLOOPTOOLS
COMMON_INC+=-I$(LOOPTOOLSDIR)
COMMON_INC+=-I$(LOOPTOOLSDIR)/include
COMMON_INC+=-I$(LOOPTOOLSDIR)/include/NLOX_util
COMMON_LIB_DIRS+=-L$(LOOPTOOLSDIR)
COMMON_LIB_DIRS+=-L$(LOOPTOOLSDIR)/lib
COMMON_LIB_DIRS+=-L$(LOOPTOOLSDIR)/lib64
COMMON_LIB_DIRS+=-L$(LOOPTOOLSDIR)/lib/NLOX_util
COMMON_LIBS+=-looptools
fortran_needed = yes
endif

## Currently NLOX/src/tred/interfaces.h looks for the includes of QCDLoop2 in a subfolder qcdloop/.
## The library files of QCDLoop2 are currently located in a qcdloop2/ subfolder in the NLOX_util installation. 
ifdef QCDLOOP2DIR
CFLAGS+=-DQCDLOOP2 -std=gnu++11
COMMON_INC+=-I$(QCDLOOP2DIR)
COMMON_INC+=-I$(QCDLOOP2DIR)/include
COMMON_INC+=-I$(QCDLOOP2DIR)/include/NLOX_util
COMMON_LIB_DIRS+=-L$(QCDLOOP2DIR)
COMMON_LIB_DIRS+=-L$(QCDLOOP2DIR)/lib
COMMON_LIB_DIRS+=-L$(QCDLOOP2DIR)/lib64
COMMON_LIB_DIRS+=-L$(QCDLOOP2DIR)/lib/NLOX_util
COMMON_LIB_DIRS+=-L$(QCDLOOP2DIR)/lib/NLOX_util/qcdloop2
COMMON_LIBS+=-Wl,-Bstatic -lqcdloop2 -Wl,-Bdynamic
endif

ifdef ONELOOPDIR
CFLAGS+=-DONELOOP
COMMON_LIB_DIRS+=-L$(ONELOOPDIR)
COMMON_LIB_DIRS+=-L$(ONELOOPDIR)/lib
COMMON_LIB_DIRS+=-L$(ONELOOPDIR)/lib64
COMMON_LIB_DIRS+=-L$(ONELOOPDIR)/lib/NLOX_util
COMMON_LIBS+=-lavh_olo
fortran_needed = yes
endif

ifdef fortran_needed
COMMON_LIBS+=-lgfortran
endif

TEST_INC=$(COMMON_INC)
TEST_LIB_DIRS=$(COMMON_LIB_DIRS)
TEST_LIBS=$(COMMON_LIBS)

BORN_OBJS=$(patsubst %.cc,%.o,$(wildcard $(PROCPATH)/*/*born*/C/*.cc))
VIRT_OBJS=$(patsubst %.cc,%.o,$(wildcard $(PROCPATH)/*/*virt*/C/*.cc))
CT_OBJS=$(patsubst %.cc,%.o,$(wildcard $(PROCPATH)/*/*CT*/C/*.cc))

SUBPROC_CODE_OBJS=$(patsubst %.cc,%.o,$(wildcard $(PROCPATH)/*/code/*.cc))
PROC_CODE_OBJS=$(patsubst %.cc,%.o,$(wildcard $(PROCPATH)/code/*.cc))

ALL_OBJECTS=$(patsubst %.cc,%.o,$(wildcard *.cc)) $(BORN_OBJS) $(VIRT_OBJS) $(CT_OBJS) $(SUBPROC_CODE_OBJS) $(PROC_CODE_OBJS)

LIB_OBJ=$(filter-out $(PROCPATH)/code/ftest_process.o $(PROCPATH)/code/test_process.o $(PROCPATH)/code/nlox_olp_fortran.o, $(ALL_OBJECTS))
FLIB_OBJ=$(filter-out $(PROCPATH)/code/ftest_process.o $(PROCPATH)/code/test_process.o $(PROCPATH)/code/nlox_olp.o, $(ALL_OBJECTS))

all: test

test: library examples/test_process.cc
	$(shell cp examples/test_process.cc .)
	$(shell cp examples/nlox_parameters.par .)
	$(CXX) $(CFLAGS) $(TEST_INC) -c test_process.cc -o test_process.o
	$(CXX) $(CFLAGS) -o test_process $(TEST_INC) $(TEST_LIB_DIRS) -L. test_process.o -lprocess $(TEST_LIBS)

ftest: flibrary examples/ftest_process.f90
	$(shell cp examples/ftest_process.f90 .)
	$(shell cp examples/nlox_parameters.par .)
	$(FC) $(FFLAGS) $(TEST_INC) -c ftest_process.f90 -o ftest_process.o
	$(FC) $(FFLAGS) -o ftest_process $(TEST_INC) $(TEST_LIB_DIRS) -L. ftest_process.o -lfprocess $(TEST_LIBS) 

library: $(LIB_OBJ)
	ar rvs libprocess.a $(LIB_OBJ)

flibrary: $(FLIB_OBJ)
	ar rvs libfprocess.a $(FLIB_OBJ)

%.o: %.cc
	$(CXX) $(CFLAGS) $(TEST_INC) -c $< -o $@

clean:
	rm -f *# *.tmp *~
	rm -f $(ALL_OBJECTS)
	rm -f test_process.o test_process.cc test_process nlox_parameters.par libprocess.a
	rm -f ftest_process.o ftest_process.f90 ftest_process nlox_parameters.par libfprocess.a

print-%  : ; @echo $* = $($*)


