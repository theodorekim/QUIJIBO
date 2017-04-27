include ../include.mk

#CC         = g++
CC         = /opt/local/bin/gcc-mp-5
CFLAGS     = ${CFLAGS_COMMON}
LDFLAGS    = ${LDFLAGS_COMMON} -framework GLUT -framework OpenGL -framework Cocoa
EXECUTABLE = ../../bin/bunnyMarchingCubes3D

SOURCES    = bunnyMarchingCubes3D.cpp \
					../../src/algebra/POLYNOMIAL_4D.cpp \
					../../src/algebra/MATRIX3.cpp \
					../../src/algebra/QUATERNION.cpp \
					../../src/algebra/QUATERNION_T.cpp \
					../../src/geometry/FIELD_2D.cpp \
					../../src/geometry/FIELD_3D.cpp \
					../../src/geometry/TRIANGLE_MESH.cpp \
					../../src/geometry/TRIANGLE.cpp \
					../../src/util/SIMPLE_PARSER.cpp \
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
