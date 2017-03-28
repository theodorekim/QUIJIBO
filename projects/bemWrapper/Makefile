include ../include.mk

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON} -O3
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/bemWrapper

SOURCES    = bemWrapper.cpp \
						../../src/algebra/QUATERNION.cpp \
						../../src/algebra/POLYNOMIAL_4D.cpp \
						../../src/algebra/MATRIX3.cpp \
						../../src/geometry/TRIANGLE.cpp \
						../../src/geometry/TRIANGLE_MESH.cpp \
						../../src/geometry/FIELD_2D.cpp \
						../../src/geometry/VECTOR3_FIELD_3D.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/SIMPLE_PARSER.cpp \
						../../src/geometry/FIELD_3D.cpp
#SOURCES    = example2.cpp ../../src/geometry/FIELD_2D.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
