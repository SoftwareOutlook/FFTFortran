	latex FFTFortran
	latex FFTFortran
	bibtex FFTFortran
	latex FFTFortran
	latex FFTFortran
	dvips -Pcmz -Pamz -Ppdf -G0 FFTFortran -o FFTFortran.ps
	ps2pdf -sPAPERSIZE=a4 FFTFortran.ps FFTFortran.pdf 
