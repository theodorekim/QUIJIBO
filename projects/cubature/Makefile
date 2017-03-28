include ../include.mk

CFLAGS_COMMON = -c -Wall -I./ -I/opt/local/include/ -O3 -fopenmp -msse2 ${MY_INCLUDES}

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON}
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/cubature

SOURCES    = cubature.cpp \
		../../src/geometry/VECTOR3_FIELD_3D.cpp \
		../../src/geometry/FIELD_2D.cpp \
		../../src/geometry/FIELD_3D.cpp \
		../../src/geometry/TRIANGLE.cpp \
		../../src/geometry/TRIANGLE_MESH.cpp \
		../../src/geometry/SPHERE.cpp \
		../../src/solvers/OPTIMIZE_3D.cpp \
		../../src/algebra/POLYNOMIAL_4D.cpp \
		../../src/algebra/MATRIX3.cpp \
		../../src/algebra/QUATERNION.cpp \
		../../src/algebra/QUATERNION_T.cpp \
		../../src/util/TIMER.cpp \
		../../src/util/SIMPLE_PARSER.cpp

OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
