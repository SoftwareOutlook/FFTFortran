
##SCARF
#FCC = gfortran  
#FFLAGS77 = -O3  -fopenmp

#FFLAGS = -O3  -fopenmp  -std=f2003   -fdefault-integer-8 -I${MKLROOT}/include

#FFLAGS = -O3  -fopenmp -fimplicit-none -Wall -std=f2003   -fdefault-integer-8  -I${MKLROOT}/include -g -fbacktrace  

##ARCHER with module swap PrgEnv-cray PrgEnv-gnu
FCC = ftn

#FFLAGS77 = -O3  -fopenmp
#FFLAGS = -O3  -fopenmp  -std=f2003  -I${MKLROOT}/include

FLAGS77 = -O3 -openmp -g -traceback

FFLAGS = -O3 -openmp -g -traceback -warn all -I${MKLROOT}/include -I${HOME}/local/include


#FCC = mpif90
#FFLAGS = -O3 -qopenmp -integer-size 64 -I${MKLROOT}/include 
#FFLAGS77 = -O3  -qopenmp

#MKLROOT=  /apps/intel/2017/compilers_and_libraries_2017.2.174/linux/mkl


#f90 = g95 
#FFLAGS= -fintrinsic-extensions -std=f2003  -Wimplicit-none -ftrace=full

#f90 = f95
#FFLAGS = -O3 -fomit-frame-pointer -fopenmp -I..

#LDFLAGS = -L/opt/intel/composer_xe_2013_sp1.1.106/mkl

f90 = $(FCC) $(FFLAGS)
f77 = $(FCC) $(FFLAGS77)

#LIBS = 	-lfftw3_threads -lfftw3 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5  -ldl

#LIBS =  -lfftw3_threads -lfftw3 -Wl,  -L{$MKLROOT}/lib/intel64/ -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -Wl,--end-group -liomp5 -ldl

LIBS = ${HOME}/local/lib/libp3dfft.a   -lfftw3_mpi -lfftw3_threads -lfftw3 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_cdft_core.a  ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -ldl -liomp5 -lm -lpthread

#LIBS =   -dynamic      -lfftw3_threads -lfftw3 -Wl, -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl


ddeps= 	


ddeps90= 


dtdeps=  


dtdeps90=  

ffteparam= ffte-6.0/param.h


all:  bench1d bench1dc bench2d bench2dc bench3d bench3dc
#all: bench1d_mpi bench1dc_mpi 
#all: bench2dc_mpi bench2d_mpi

#all: bench3dc_mpi bench3d_mpi

bench2d:   mkl_dfti.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o bench2d.o
	 $(f90) mkl_dfti.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o bench2d.o $(LIBS) -o bench2d.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	 ./bench2d.exe  6  3 2 2 2

bench2d_mpi:   mkl_cdft.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o bench2d_mpi.o
	$(f90) mkl_cdft.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o bench2d_mpi.o $(LIBS) -o bench2d_mpi.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)


bench2dc:   mkl_dfti.o zfft2d.o fft235.o factor.o kernel.o bench2dc.o
	 $(f90) mkl_dfti.o zfft2d.o fft235.o factor.o kernel.o bench2dc.o $(LIBS) -o bench2dc.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench2dc.exe  6  3 2 2 2

bench2dc_mod:   mkl_dfti.o zfft2d_mod.o fft235.o factor.o kernel.o bench2dc.o
	$(f90) mkl_dfti.o zfft2d_mod.o fft235.o factor.o kernel.o bench2dc.o $(LIBS) -o bench2dc_mod.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench2dc_mod.exe  6  3 2 2 2

bench2dc_mpi:   mkl_cdft.o zfft2d.o fft235.o factor.o kernel.o bench2dc_mpi.o
	$(f90) mkl_cdft.o zfft2d.o fft235.o factor.o kernel.o bench2dc_mpi.o $(LIBS) -o bench2dc_mpi.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
#        ./bench2dc_mpi.exe  6  3 2 2 2



bench1d:     mkl_dfti.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o bench1d.o
	$(f90) mkl_dfti.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o bench1d.o $(LIBS) -o bench1d.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench1d.exe  6 3 4 2

bench1d_mpi:     mkl_cdft.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o bench1d_mpi.o
	$(f90) mkl_cdft.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o bench1d_mpi.o $(LIBS) -o bench1d_mpi.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
#        ./bench1d.exe  6 3 4 2


bench1dc:   mkl_dfti.o zfft2d.o fft235.o factor.o kernel.o bench1dc.o
	$(f90) mkl_dfti.o zfft2d.o fft235.o factor.o kernel.o bench1dc.o $(LIBS) -o bench1dc.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench1dc.exe  6  3 4 2

bench1dc_mod:   mkl_dfti.o zfft2d_mod.o fft235.o factor.o kernel.o bench1dc.o
	$(f90) mkl_dfti.o zfft2d_mod.o fft235.o factor.o kernel.o bench1dc.o $(LIBS) -o bench1dc_mod.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench1dc_mod.exe  6  3 4 2



bench1dc_mpi: mkl_cdft.o  pzfft3d.o zfft3d.o fft235.o factor.o pfactor.o kernel.o bench1dc_mpi.o
	$(f90) mkl_cdft.o pzfft3d.o zfft3d.o fft235.o factor.o pfactor.o kernel.o bench1dc_mpi.o $(LIBS) -o bench1dc_mpi.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
#	./bench1dc_mpi.exe  6  3 4 2



bench3d:   mkl_dfti.o bench3d.o
	 $(f90) mkl_dfti.o bench3d.o $(LIBS) -o bench3d.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	 ./bench3d.exe  6  3 2 2 2

bench3d_mpi:   mkl_cdft.o bench3d_mpi.o
	$(f90) mkl_cdft.o bench3d_mpi.o $(LIBS) -o bench3d_mpi.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)



bench3dc:   mkl_dfti.o zfft3d.o fft235.o factor.o kernel.o bench3dc.o
	 $(f90) mkl_dfti.o zfft3d.o fft235.o factor.o kernel.o bench3dc.o $(LIBS) -o bench3dc.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench3dc.exe  6  3 2 2 3

bench3dc_mpi:   mkl_cdft.o zfft3d.o fft235.o factor.o kernel.o bench3dc_mpi.o
	$(f90) mkl_cdft.o zfft3d.o fft235.o factor.o kernel.o bench3dc_mpi.o $(LIBS) -o bench3dc_mpi.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)



bench3dc_mod:   mkl_cdft.o zfft3d_mod.o fft235.o factor.o kernel.o bench3dc.o
	$(f90) mkl_cdft.o zfft3d_mod.o fft235.o factor.o kernel.o bench3dc.o $(LIBS) -o bench3dc_mod.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench3dc_mod.exe  6  3 2 2 3


#> temp >&1 < hsl_minresds.data
#	 diff temp hsl_minresds.output


 #> temp >&1 < hsl_minreszs.data
 #	 diff temp hsl_minreszs.output

ddeps90.f90:
	cat  $(ddeps90) >ddeps90.f90
	echo $(ddeps90) >ddeps90
ddeps.f:
	cat  $(ddeps) $(DBLAS) >ddeps.f
	echo $(ddeps) >ddeps
	echo $(DBLAS) >dblas
param.h:
	cat $(ffteparam) >param.h
	echo $(ffteparam) >param


ddeps.o:	ddeps.f
	$(f77) -c ddeps.f
ddeps90.o:	ddeps90.f90
	$(f90) -c ddeps90.f90

bench2dc.o:      bench2dc.f90
	$(f90) -c bench2dc.f90

bench2dc_mpi.o:      bench2dc_mpi.f90
	$(f90) -c bench2dc_mpi.f90

bench2d_mpi.o:      bench2d_mpi.f90
	$(f90) -c bench2d_mpi.f90


bench2d.o:	bench2d.f90
	$(f90) -c bench2d.f90

bench1d.o:      bench1d.f90
	$(f90) -c bench1d.f90

bench1d_mpi.o:      bench1d_mpi.f90
	$(f90) -c bench1d_mpi.f90

bench1dc.o:      bench1dc.f90
	$(f90) -c bench1dc.f90

bench1dc_mpi.o:      bench1dc_mpi.f90
	$(f90) -c bench1dc_mpi.f90

bench3d.o:      bench3d.f90
	$(f90) -c bench3d.f90

bench3d_mpi.o:      bench3d_mpi.f90
	$(f90) -c bench3d_mpi.f90


bench3dc.o:      bench3dc.f90
	$(f90) -c bench3dc.f90

bench3dc_mpi.o:      bench3dc_mpi.f90
	$(f90) -c bench3dc_mpi.f90

mkl_dfti.o : mkl_dfti.f90 
	$(f90) -c mkl_dfti.f90


mkl_cdft.o : mkl_cdft.f90
	$(f90) -c mkl_cdft.f90

zfft1d.o : ffte-6.0/zfft1d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zfft1d.f   -o zfft1d.o

pzfft1d.o : ffte-6.0/mpi/pzfft1d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/mpi/pzfft1d.f   -o pzfft1d.o

pzfft2d.o : ffte-6.0/mpi/pzfft2d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/mpi/pzfft2d.f   -o pzfft2d.o

pzfft3d.o : ffte-6.0/mpi/pzfft3d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/mpi/pzfft3d.f   -o pzfft3d.o

dzfft2d.o : ffte-6.0/dzfft2d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/dzfft2d.f   -o dzfft2d.o


zdfft2d.o : ffte-6.0/zdfft2d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zdfft2d.f   -o zdfft2d.o


zfft2d.o : ffte-6.0/zfft2d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zfft2d.f   -o zfft2d.o

zfft2d_mod.o : ffte-6.0/zfft2d_mod.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zfft2d_mod.f   -o zfft2d_mod.o


dzfft3d.o : ffte-6.0/dzfft3d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/dzfft3d.f   -o dzfft3d.o


zdfft3d.o : ffte-6.0/zdfft3d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zdfft3d.f   -o zdfft3d.o

zdfft3d_mod.o : ffte-6.0/zdfft3d_mod.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zdfft3d_mod.f   -o zdfft3d_mod.o

zfft3d.o : ffte-6.0/zfft3d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zfft3d.f   -o zfft3d.o


fft235.o : ffte-6.0/fft235.f kernel.o  param.h
	$(f77) -c ffte-6.0/fft235.f   -o fft235.o


factor.o : ffte-6.0/factor.f  param.h
	$(f77) -c ffte-6.0/factor.f   -o factor.o
pfactor.o : ffte-6.0/mpi/pfactor.f  param.h
	$(f77) -c ffte-6.0/mpi/pfactor.f   -o pfactor.o


kernel.o : ffte-6.0/kernel.f  param.h
	$(f77) -c ffte-6.0/kernel.f   -o kernel.o

clean:
	rm a.out *.o  temp *deps.f *.mod *deps90.f90 *.exe

ffte=./ffte-6.0/

