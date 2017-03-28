include ../include.mk

CC         = g++
CFLAGS     = ${CFLAGS_COMMON}
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/marchingCubes3D

SOURCES    = marchingCubes3D.cpp \
					../../src/algebra/POLYNOMIAL_4D.cpp \
					../../src/algebra/MATRIX3.cpp \
					../../src/algebra/QUATERNION.cpp \
					../../src/algebra/QUATERNION_T.cpp \
					../../src/geometry/FIELD_2D.cpp \
					../../src/geometry/FIELD_3D.cpp \
					../../src/geometry/TRIANGLE_MESH.cpp \
					../../src/geometry/TRIANGLE.cpp \
					../../src/util/TIMER.cpp \
					../../src/solvers/OPTIMIZE_3D.cpp 

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
