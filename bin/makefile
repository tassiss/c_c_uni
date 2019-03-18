all: mod_cond.o cond.o
	gfortran -o cond cond.o mod_cond.o
cond.o: ../source/cond.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/cond.f90
mod_cond.o: ../source/mod_cond.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/mod_cond.f90

test: all
	valgrind --leak-check=full --show-reachable=yes --track-origins=yes ./cond
clean:
	rm -rf  cond *.*
