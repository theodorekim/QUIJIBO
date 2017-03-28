#include ../include.mk

# calls:
CC         = g++
CFLAGS     = -c -Wall -I./vcglib/ -I./vcglib/eigenlib -I../../src/solvers -I../../src/TOUPEE -I../../src/integrators -I../../src/algebra -I../../src/geometry -I../../src/util -I../../src/solvers -I./ -I/opt/local/include -I../../ -O3 -fopenmp -msse2
LDFLAGS    = -L/usr/local/lib/x86_64 -L/opt/local/lib -lstdc++ -lpng -lz -ljpeg -msse2
EXECUTABLE = ../../bin/poissonDisk

SOURCES    = vcg_mesh.cpp \
						 utils_sampling.cpp \
						 ../../src/geometry/VECTOR3_FIELD_3D.cpp \
						 ../../src/geometry/FIELD_3D.cpp \
						 ../../src/geometry/FIELD_2D.cpp \
						 ../../src/geometry/TRIANGLE_MESH.cpp \
						 ../../src/geometry/TRIANGLE.cpp \
						 ../../src/algebra/MATRIX3.cpp \
						 ../../src/algebra/QUATERNION.cpp \
						 ../../src/algebra/QUATERNION_T.cpp \
						 ../../src/algebra/POLYNOMIAL_4D.cpp \
						 ../../src/util/SIMPLE_PARSER.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/solvers/OPTIMIZE_3D.cpp \
						 ./vcglib/wrap/ply/plylib.cpp \
						 poissonDisk.cpp
#SOURCES    = example2.cpp ../../src/geometry/FIELD_2D.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o; rm ./vcglib/wrap/ply/plylib.o
