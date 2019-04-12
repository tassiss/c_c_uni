all: in_out.o calc.o
	gfortran -o cond calc.o in_out.o
calc.o: ../source/calc.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/calc.f90
in_out.o: ../source/in_out.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/in_out.f90

test: all
	valgrind --leak-check=full --show-reachable=yes --track-origins=yes ./cond
clean:
	rm -rf  cond *.*

