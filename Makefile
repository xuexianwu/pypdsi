####################################################################################
#    Description:
#    Author: 
#    Date:
####################################################################################

CC=g++
CFLAGS=-O2 -g -funroll-loops  -Wall

# executable name
TARGET=pdsi

# object files
OBJ=pdsi.o

# assume these are in standard places
LIBS=
INC=-I./

all: ${TARGET}

${TARGET}: ${OBJ} 
	${CC} ${CFLAGS} ${OBJ} -o $@ ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${INC} -c $< 

clean:
	rm -f ${OBJ} ${TARGET}

