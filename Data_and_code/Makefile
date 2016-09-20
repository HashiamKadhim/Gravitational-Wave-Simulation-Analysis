# Makefile for unit test

CXX = g++
LDLIBS = -lgsl -lgslcblas -lfftw3 -lcblas

hw6: 
	${CXX}  main.cpp ${LDLIBS} -o hw6 -O2

run:
	./hw6

clean: 
	rm -f hw6

