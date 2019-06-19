all:functions.o in_out.o calc.o lapack_parcer.o
	gfortran -o cond *.o
calc.o: ../source/calc.f90 lapack_parcer.o
	gfortran -c -Wall -Werror -fcheck=all ../source/calc.f90
in_out.o: ../source/in_out.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/in_out.f90
functions.o: ../source/functions.f90
	gfortran -c -Wall -Werror -fcheck=all ../source/functions.f90
lapack_tools.o: ../source/pacotes/lapack_tools.f
	gfortran -c -Wall -fcheck=all ../source/pacotes/lapack_tools.f
lapack.o: ../source/pacotes/lapack.f lapack_tools.o
	gfortran -c -Wall -fcheck=all ../source/pacotes/lapack.f
lapack_parcer.o: ../source/pacotes/lapack_parcer.f08 lapack.o
	gfortran -c -Wall -Werror -fcheck=all ../source/pacotes/lapack_parcer.f08
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


