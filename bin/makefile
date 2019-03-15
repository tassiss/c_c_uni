all: cond.o
	gfortran -o cond cond.o
cond.o: ../source/cond.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/cond.f90
clean:
	rm -rf  cond *.*
