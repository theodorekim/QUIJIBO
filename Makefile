SHELL := /bin/bash -e

SRCS = geometry algebra util solvers TOUPEE

all:
	cd projects/bemWrapper; make
	cd projects/cubature; make
	cd projects/marchingCubes3D; make
	cd projects/objToGmsh; make
	cd projects/poissonDisk; make
	cd projects/sdfGen; make
	cd projects/toothOptimizer; make
	cd projects/rosenbrock; make

clean: 
	-for d in $(SRCS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./src/$$d; rm *.o; cd ../..); done
	cd projects/bemWrapper; make clean
	cd projects/cubature; make clean
	cd projects/marchingCubes3D; make clean
	cd projects/objToGmsh; make clean
	cd projects/poissonDisk; make clean
	cd projects/sdfGen; make clean
	cd projects/toothOptimizer; make clean
	cd projects/rosenbrock; make clean

ctags:
	ctags -R *
