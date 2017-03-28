include ../include.mk

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON} -O3
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/toothOptimizer

SOURCES    = toothOptimizer.cpp \
						 ../../src/TOUPEE/TAO_SOLVER.cpp \
						 ../../src/TOUPEE/LINE_SEARCH.cpp \
						 ../../src/TOUPEE/BLMVM.cpp \
						 ../../src/TOUPEE/MORE_THUENTE.cpp \
						 ../../src/TOUPEE/GPCG_LINE_SEARCH.cpp \
						 ../../src/TOUPEE/ARMIJO.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/SIMPLE_PARSER.cpp \
						 ../../src/geometry/FIELD_3D.cpp \
						 ../../src/geometry/FIELD_2D.cpp \
						 ../../src/geometry/TRIANGLE_MESH.cpp \
						 ../../src/geometry/TRIANGLE.cpp \
						 ../../src/algebra/POLYNOMIAL_4D.cpp \
						 ../../src/algebra/QUATERNION.cpp \
						 ../../src/algebra/MATRIX3.cpp \
						 ../../src/solvers/OPTIMIZE_3D.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
