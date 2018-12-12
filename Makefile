#f90= nagfor -dcfuns -C=all -C=undefined -gline -f95 
# ifort not available on fox
#f90= ifort -C -check noarg_temp_created -u
#f90 = g95 -fintrinsic-extensions -std=f2003  -Wimplicit-none -ftrace=full

f90 = gfortran  $(FFLAGS)
FC = $(f90)
FFLAGS = -fimplicit-none  -fbounds-check
F77 = $(f90)

LIBS = 	



ddeps= 	


ddeps90= 


dtdeps=  


dtdeps90=  

ffteparam= ffte-6.0/param.h


all:  bench2d


bench2d:   bench2d.o dzfft2d.o fft235.o factor.o kernel.o
	 $(f90) $(FFLAGS) bench2d.o  dzfft2d.o fft235.o factor.o kernel.o $(LIBS) -o bench2d.exe
	 ./bench2d.exe  5  5 5 4 5
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
	$(f90) $(FFLAGS) -c ddeps.f
ddeps90.o:	ddeps90.f90
	$(f90) $(FFLAGS) -c ddeps90.f90

hsl_minress.o:	hsl_minress.f90
	$(f90) $(FFLAGS) -c hsl_minress.f90
bench2d.o:	bench2d.f90
	$(f90) $(FFLAGS) -c bench2d.f90

dzfft2d.o : ffte-6.0/dzfft2d.f fft235.o factor.o param.h
	$(F77) $(FFLAGS) -c ffte-6.0/dzfft2d.f   -o dzfft2d.o


fft235.o : ffte-6.0/fft235.f kernel.o  param.h
	$(F77) $(FFLAGS) -c ffte-6.0/fft235.f   -o fft235.o


factor.o : ffte-6.0/factor.f  param.h
	$(F77) $(FFLAGS) -c ffte-6.0/factor.f   -o factor.o


kernel.o : ffte-6.0/kernel.f  param.h
	$(F77) $(FFLAGS) -c ffte-6.0/kernel.f   -o kernel.o

clean:
	rm a.out *.o  temp *deps.f *.mod *deps90.f90 

ffte=./ffte-6.0/

