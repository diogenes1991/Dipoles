CXX=g++

PROCPATH=$(shell pwd)/Matrix_Elements

CFLAGS= -O3 -Wall -Wno-reorder -DNLOXPROCPATH=\"$(PROCPATH)\"

TREDDIR=/home/diogenes1991/NLOX/src/tred
ONELOOPDIR=/home/diogenes1991/NLOX_util/built
QCDLOOPDIR=/home/diogenes1991/NLOX_util/built
QCDLOOP2DIR=/home/diogenes1991/NLOX_util/built
LOOPTOOLSDIR=/home/diogenes1991/NLOX_util/built

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

COMMON_INC += -I$(PROCPATH)/code
COMMON_INC += -I$(PROCPATH)
COMMON_LIB_DIRS += -L$(PROCPATH)

all: test

test:
	$(CXX) -c $(CFLAGS) Test_Code.cpp $(COMMON_INC) $(COMMON_LIB_DIRS) $(COMMON_LIBS)
	$(CXX) $(CFLAGS) -o Test_Code $(COMMON_INC) $(COMMON_LIB_DIRS) -L. Test_Code.o -lprocess $(COMMON_LIBS)
