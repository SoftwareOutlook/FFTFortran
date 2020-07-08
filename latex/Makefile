all:
	latex FFTFortran
	latex FFTFortran
	bibtex FFTFortran
	latex FFTFortran
	latex FFTFortran
	dvips -Pcmz -Pamz -Ppdf -G0 FFTFortran -o FFTFortran.ps
	ps2pdf -sPAPERSIZE=a4 FFTFortran.ps FFTFortran.pdf 
#	ps2pdf -dCompatibilityLevel=1.2 FFTFortran.ps FFTFortran.pdf 
#	acroread FFTFortran.pdf
clean: 
	rm -f *~ *.log *.aux *.dvi *.toc *.bbl *.blg *.ps 
