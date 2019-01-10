

FCC = gfortran  
FFLAGS77 = -O3  -fopenmp
FFLAGS = -O3  -fopenmp -std=f2003   -fdefault-integer-8  -I${MKLROOT}/include -g -fbacktrace  


#FCC = mpif90
#FFLAGS = -O3 -qopenmp -integer-size 64 -I${MKLROOT}/include 
#FFLAGS77 = -O3  -qopenmp

MKLROOT=  /apps/intel/2017/compilers_and_libraries_2017.2.174/linux/mkl


#f90 = g95 
#FFLAGS= -fintrinsic-extensions -std=f2003  -Wimplicit-none -ftrace=full

#f90 = f95
#FFLAGS = -O3 -fomit-frame-pointer -fopenmp -I..


f90 = $(FCC) $(FFLAGS)
f77 = $(FCC) $(FFLAGS77)

LIBS = 	-lfftw3 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl



ddeps= 	


ddeps90= 


dtdeps=  


dtdeps90=  

ffteparam= ffte-6.0/param.h


all:  bench2dc


bench2d:   mkl_dfti.o bench2d.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o
	 $(f90) bench2d.o mkl_dfti.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o $(LIBS) -o bench2d.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	 ./bench2d.exe  256  256 256 2 2

bench2dc:   mkl_dfti.o bench2dc.o zfft2d.o fft235.o factor.o kernel.o
	 $(f90) bench2dc.o mkl_dfti.o zfft2d.o fft235.o factor.o kernel.o $(LIBS) -o bench2dc.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench2dc.exe  256  256 256 2 2


bench1d:   mkl_dfti.o bench1d.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o
	$(f90) bench1d.o mkl_dfti.o dzfft2d.o zdfft2d.o fft235.o factor.o kernel.o $(LIBS) -o bench1d.exe
	echo $(LIBS)
	echo $(OMP_NUM_THREADS)
	./bench1d.exe  6  2 4 3



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

bench2d.o:	bench2d.f90
	$(f90) -c bench2d.f90

bench1d.o:      bench1d.f90
	$(f90) -c bench1d.f90


mkl_dfti.o : mkl_dfti.f90 
	$(f90) -c mkl_dfti.f90


dzfft2d.o : ffte-6.0/dzfft2d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/dzfft2d.f   -o dzfft2d.o


zdfft2d.o : ffte-6.0/zdfft2d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zdfft2d.f   -o zdfft2d.o


zfft2d.o : ffte-6.0/zfft2d.f fft235.o factor.o param.h
	$(f77) -c ffte-6.0/zfft2d.f   -o zfft2d.o


fft235.o : ffte-6.0/fft235.f kernel.o  param.h
	$(f77) -c ffte-6.0/fft235.f   -o fft235.o


factor.o : ffte-6.0/factor.f  param.h
	$(f77) -c ffte-6.0/factor.f   -o factor.o


kernel.o : ffte-6.0/kernel.f  param.h
	$(f77) -c ffte-6.0/kernel.f   -o kernel.o

clean:
	rm a.out *.o  temp *deps.f *.mod *deps90.f90 

ffte=./ffte-6.0/

