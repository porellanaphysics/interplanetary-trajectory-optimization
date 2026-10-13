# Makefile for Mars Trajectory Simulation & Optimization
FC ?= gfortran
FFLAGS = -O3 -Wall -fdefault-real-8
PYTHON ?= python3

.PHONY: all clean f2py run test

all: marte

marte: marte.f95
	$(FC) $(FFLAGS) -o $@ marte.f95

f2py: simulacion2.f95
	$(PYTHON) -m numpy.f2py -c simulacion2.f95 -m simulacion2

run: marte
	./marte

clean:
	rm -f marte pruebaSim *.o *.mod *.so salida.log
