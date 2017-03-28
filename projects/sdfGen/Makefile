include ../include.mk

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON} -O3
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/sdfGen 

SOURCES    = main.cpp makelevelset3.cpp \
             ../../src/algebra/MATRIX3.cpp \
             ../../src/geometry/FIELD_2D.cpp \
             ../../src/geometry/FIELD_3D.cpp \
             ../../src/util/SIMPLE_PARSER.cpp \
             ../../src/util/TIMER.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
