# stub makefile, does nothing, but causes R CMD build not to force rebuild pdf file

REVAL=TRUE

all:	accuracy_vignette.pdf

%.tex:	%.Rnw
	echo "Sweave(\"$<\", debug=TRUE, eval=$(REVAL))" | R --slave

%.pdf:	%.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

clean:	
	rm -f accuracy_vignette*.pdf accuracy_vignette*.tex accuracy_vignette*.eps accuracy_vignette.out accuracy_vignette.aux accuracy_vignette.dvi accuracy_vignette.bbl accuracy_vignette.blg Rplots.* accuracy_vignette.log

