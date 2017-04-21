include ../include.mk
# common stuff
#LDFLAGS_COMMON = -framework GLUT -framework OpenGL -framework Cocoa -framework Accelerate -lz
#CFLAGS_COMMON = -O3 -c -Wall -I../../src/glvu/ -I../../src/geometry/ -I/opt/local/include -I/opt/X11/include

#CC         = g++
CC         = /opt/local/bin/gcc-mp-5
CFLAGS     = ${CFLAGS_COMMON} -I../../src/glvu/
LDFLAGS    = ${LDFLAGS_COMMON} -framework GLUT -framework OpenGL -framework Cocoa
EXECUTABLE = ../../bin/objCutter

GEOMETRYSRC = ../../src/geometry/OBJ.cpp\
							../../src/geometry/TRIANGLE.cpp \
							../../src/geometry/BOX.cpp
UTILSRC = ../../src/glvu/glvu.cpp \
					../../src/util/TIMER.cpp
SOURCES    = main.cpp $(GEOMETRYSRC) $(UTILSRC)
OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
