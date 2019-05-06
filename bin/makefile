all:functions.o in_out.o calc.o 
	gfortran -o cond calc.o in_out.o functions.o
calc.o: ../source/calc.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/calc.f90
in_out.o: ../source/in_out.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/in_out.f90
functions.o: ../source/functions.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/functions.f90
output: 
	mkdir ../output
test: all
	valgrind --leak-check=full --show-reachable=yes --track-origins=yes ./cond
clean:
	rm -rf  cond *.*
python:
	cd ../source/ && python out.py
gif: python
	cd ../output/ && convert -delay 10 -loop 0 *.png ../out.gif


out:
	cd ../output/ && rm -rf *.*
