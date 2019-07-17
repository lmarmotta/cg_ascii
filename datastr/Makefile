# Because of the LAPACK dependency, this code shall be compiled as:
#
# gcc -Ddgetrf=dgetrf_ -Ddgetri=dgetri_ m_inversion.c -lblas -llapack 
#
# Name of the exec file.
TARGET=code

#
# Link flags abreviation.
LIBS=-lm -lcgns -lhdf5 -lblas -llapack

#
# Path to hdf5 parallel compiler (needs compilation first).
CC=/home/leonardo/opt/hdf5/bin/h5pcc

# 
# Optimization or debug flags.
#
# -O0: Reduce compilation time and make debugging produce the expected results. This is the default.
# -std=c99: Use the 1999 standard.
# -ggdb: Produces debugging information specifically intended for gdb
# -Wall: This enables all the warnings about constructions that some users consider questionable
# -Wextra: This enables some extra warning flags that are not enabled by -Wall
# -pedantic: Standard diagnostics.
# -Wuninitialized: Warn if an automatic variable is used without first being initialized.
# -Winit-self: Warn about uninitialized variables that are initialized with themselves.
# -Wstrict-prototypes: Warn if a function is declared or defined without specifying the argument types.
# -Wformat-security: If -Wformat is specified, also warn about uses of format functions that represent possible security problems.
#
CFLAGS=-O0 -std=c99 -ggdb -Wall -Wextra -pedantic -Wuninitialized -Winit-self -Wstrict-prototypes -Wformat -Wformat-security -Warray-bounds -Ddgetrf=dgetrf_ -Ddgetri=dgetri_ 
# CFLAGS=-O3 -std=c99 -Ddgetrf=dgetrf_ -Ddgetri=dgetri_

#
# CGNS include path.
INCLUDE=/home/leonardo/opt/cgns/include

#
# CGNS lib path.
LDFLAGS=/home/leonardo/opt/cgns/lib

#
# ----------- DO NOT MODIFY -----------  
#

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS=$(patsubst %.c, %.o, $(wildcard *.c))
HEADERS=$(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(INCLUDE) -c $< -o $@ -L$(LDFLAGS) $(LIBS)

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -I$(INCLUDE) -o $@ -L$(LDFLAGS) $(LIBS)

clean:
	-rm -f *.o *.dat *.pdf
	-rm -f $(TARGET)
