#f90= nagfor -dcfuns -C=all -C=undefined -gline -f95 
# ifort not available on fox
#f90= ifort -C -check noarg_temp_created -u
#f90 = g95 -fintrinsic-extensions -std=f2003  -Wimplicit-none -ftrace=full

f90 = gfortran  -fimplicit-none  -fbounds-check
FC = $(f90)


LIBS = 	



ddeps= 	


ddeps90= 


dtdeps=  


dtdeps90=  


all:  bench2d


bench2d:   bench2d.o 
	 $(f90)  bench2d.o   $(LIBS) -o bench2d.exe
	 ./bench2d.exe  1  2 3 4
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

ddeps.o:	ddeps.f
	$(f90) -c ddeps.f
ddeps90.o:	ddeps90.f90
	$(f90) -c ddeps90.f90

hsl_minress.o:	hsl_minress.f90
	$(f90) -c hsl_minress.f90
bench2d.o:	bench2d.f90
	$(f90) -c bench2d.f90
clean:
	rm a.out *.o  temp *deps.f *.mod *deps90.f90 



