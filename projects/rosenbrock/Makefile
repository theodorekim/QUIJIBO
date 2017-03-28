CFLAGS_COMMON = -c -Wall -I../../Eigen/ -I../../ -I../../src/TOUPEE -O3 
LDFLAGS_COMMON = -lstdc++

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON} -O3
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/rosenbrock

SOURCES    = rosenbrock.cpp \
						../../src/TOUPEE/ARMIJO.cpp \
						../../src/TOUPEE/BLMVM.cpp \
						../../src/TOUPEE/GPCG_LINE_SEARCH.cpp \
						../../src/TOUPEE/LINE_SEARCH.cpp \
						../../src/TOUPEE/MORE_THUENTE.cpp \
						../../src/TOUPEE/TAO_SOLVER.cpp

OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
