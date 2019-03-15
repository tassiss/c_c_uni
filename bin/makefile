all: cond.o
	gfortran -o cond cond.o
cond.o: ../source/cond.f90
	gfortran -c ../source/cond.f90
clean:
	rm -rf  cond *.*
