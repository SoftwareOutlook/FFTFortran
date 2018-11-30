all:
	latex benchmark_plan
	latex benchmark_plan
#	bibtex benchmark_plan
#	latex benchmark_plan
#	latex benchmark_plan
	dvips -Pcmz -Pamz -Ppdf -G0 benchmark_plan -o benchmark_plan.ps
	ps2pdf -sPAPERSIZE=a4 benchmark_plan.ps benchmark_plan.pdf 
#	ps2pdf -dCompatibilityLevel=1.2 benchmark_plan.ps benchmark_plan.pdf 
#	acroread benchmark_plan.pdf
clean: 
	rm -f *~ *.log *.aux *.dvi *.toc *.bbl *.blg *.ps *.pdf
