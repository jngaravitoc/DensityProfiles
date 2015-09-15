LATEX = pdflatex
MAIN  = DensityProfiles 

main: 	        $(MAIN).tex
		$(LATEX) $(MAIN).tex 
		bibtex $(MAIN).aux
		$(LATEX) $(MAIN).tex 
		$(LATEX) $(MAIN).tex
		#dvips $(MAIN).dvi -o $(MAIN).ps
		#ps2pdf $(MAIN).ps


clean:
	rm -f $(MAIN).{aux,bbl,ps,pdf,div,blg,log}
	rm -f *~

view: 
	evince $(MAIN).pdf &
