CC=gcc
AR=ar
RANLIB=ranlib

BSSN_DIR=BSSN_Hydro_without_Hydro_Ccodes
# Create output directory
$(shell mkdir -p output)

#Enable compiler debug options?
DEBUG=0

#Default compiler flags:
CFLAGS=-fPIC#-std=gnu99 -Wall -I../common_functions/ -DENABLE_SPHERICAL_TO_CARTESIAN_ETK_LAYER -fPIC

OPTS = -Ofast -march=native -funroll-loops
ifeq ($(CC), icc)
	OPTS=-O3 -ip -xHOST
endif

ifeq ($(DEBUG),1)
# If we ARE debugging, then disable all compiler optimizations and OpenMP
CFLAGS += -O0 -g -Wno-unknown-pragmas
else
# If we are not debugging, then enable OpenMP and compiler optimizations, according to Intel compiler version
#CFLAGS += -O2 -g -DPOISON_GFS
CFLAGS += $(OPTS) -fopenmp # <- slight perf gain
endif

default: libBHaH.a

libBHaH.a: bhah_lib.o
	$(AR)  r libBHaH.a bhah_lib.o
	$(RANLIB) libBHaH.a

bhah_lib.o: $(BSSN_DIR)/bhah_lib.c
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $(DEFINES) $(BSSN_DIR)/bhah_lib.c

clean:
	rm -f bhah_lib_main libBHaH.a bhah_lib.o bhah_lib_main.o
