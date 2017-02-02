COMP=-c

CC=icc
CCFLAGS= -Wall -g -lcgns
CCLOAD= $(CCFLAGS)

EXEC=cgns2sparta

OBJ=\
	lm_converter.o

all:	$(EXEC)


$(EXEC):	$(OBJ)
	$(CC) -o $(EXEC) $(CCLOAD) $(OBJ)

MM.o:	lm_converter.o

clean:
	rm -f *.o *.dat $(EXEC)


.SUFFIXES	:	.o .c

.c.o:
	$(CC) $(COMP) $(CCFLAGS) $*.c 



