# Build mode
# 0: development (max safety, no optimisation)
# 1: release (min safety, optimisation)
# 2: fast and furious (no safety, optimisation)
BUILD_MODE?=0

all: main

# Makefile definitions
MAKEFILE_INC=../../PBMake/Makefile.inc
include $(MAKEFILE_INC)

# Rules to make the executable
repo=pbdataanalysis
main: \
		main.o \
		$($(repo)_EXE_DEP) \
		$($(repo)_DEP)
	$(COMPILER) `echo "$($(repo)_EXE_DEP) main.o" | tr ' ' '\n' | sort -u` $(LINK_ARG) $($(repo)_LINK_ARG) -o main 
	
main.o: \
		main.c \
		$($(repo)_INC_H_EXE) \
		$($(repo)_EXE_DEP)
	$(COMPILER) $(BUILD_ARG) $($(repo)_BUILD_ARG) `echo "$($(repo)_INC_DIR)" | tr ' ' '\n' | sort -u` -c main.c
	
input:
	cat Data/*.txt > targets.txt; \
	wc -l targets.txt | awk '{print $$1}' > input.txt; \
	cat targets.txt | awk '{printf "%f %f\n", $$4, $$5}' >> input.txt

run: input
	main 250 250 9 input.txt
