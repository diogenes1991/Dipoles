CCP	= g++
FCP	= gfortran

CFLAGS 	= -w -fpermissive -static -c -Wno-write-strings -Wall -O3
FFLAGS	= -c

LDIR	= $(HOME)/local
IFLAG 	= -I$(LDIR)/include 

LOGFILE = $(LOGPATH) `date +'%y-%m-%d'`

PROCESS = QCD_Dipoles
NAME 	= $(PROCESS)\$(LOGFILE)

SOURCE	= Test_QCD_Poles
OUT 	= Test_QCD

NEW 	= .NEW_\$(SOURCE)

all: compile clean run


refresh:
# 	@echo "Pulling the newests versions from Dropbox's shared folder"
	@cp -f ~/Dropbox/Projects/$(PROCESS)/$(SOURCE).cpp $(NEW).cpp
	@cmp -s $(NEW).cpp $(SOURCE).cpp; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "The file " $(SOURCE).cpp " is already in its newest version"; \
	else \
		echo "The file " $(SOURCE).cpp " is not in its newest version"; \
		read -p "Do you want to update it? [Y\N]: " ANS; \
		if [ $$ANS = y ] || [ $$ANS = Y ] || [ $$ANS = YES ] || [ $$ANS = yes ]; then \
			echo "The local version of " $(SOURCE).cpp "was updated"; \
			cp $(SOURCE).cpp .OLD_\$(SOURCE).cpp;\
			cp $(NEW).cpp $(SOURCE).cpp; \
		else \
			echo "The local version of " $(SOURCE).cpp " was kept, but it is out of date"; \
		fi; \
	fi
	@rm $(NEW).cpp

	

compile: 
#	$(FCP) $(FFLAGS)   $(FFILE).f
	
	$(CCP) $(CFLAGS) $(IFLAG) $(SOURCE).cpp 
#	$(CCP) $(SOURCE).o $(FFILE).o -o $(OUT) $(LCUBA) $(LIBS)
	$(CCP) $(SOURCE).o -o $(OUT) $(LCUBA) $(LIBS)



clean:
	@rm -f *.o *~ *.mod *# *.str

run: 
	mkdir -p $(OUT)\_Runs
	
	./$(OUT) 
# #	kate $(OUT)\_1.txt &
# 	mv $(OUT)\_1.txt $(OUT)\_Runs
# 	mv $(OUT) $(OUT)\_Runs
	
	
	
backup: 
	@mkdir -p ~/Dropbox/Backups/$(NAME)
	@cp -f *.cpp *.py *.f *.h ~/Dropbox/Backups/$(NAME)
	@cp -f *.cpp *.py *.f *.h ~/Dropbox/Projects/$(PROCESS)
	
	