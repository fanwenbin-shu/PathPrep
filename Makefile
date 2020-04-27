FC = ifort

all: rm pot compile

pot:
	$(FC) -c fh2n5z.f utility.f

compile: 
	$(FC) -c interface.f90
	$(FC) -c main.f90 -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include -L${MKLROOT}/include -llapack
	$(FC) -o PathPrep main.o interface.o fh2n5z.o utility.o -mkl ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a
	./PathPrep

clean:
	rm -f *.o *.mod log

rm: clean
	rm -f a.out opt.xyz spline.xyz umberlla.xyz
	rm -f opt_conv_difference.txt opt_conv_energy.txt opt_conv_xi.txt
	rm -f opt_difference.txt opt_energy.txt opt_ene-xi.txt opt_xi.txt out_ene-xi.txt
	rm -f umbrella_configurations jacobi.txt
