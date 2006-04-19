# stub makefile, does nothing, but causes R CMD build not to force rebuild pdf file

REVAL=TRUE

all:	accuracy_vignette.pdf

%.tex:	%.Rnw
	echo "Sweave(\"$<\", debug=TRUE, eval=$(REVAL))" | R --slave

%.pdf:	%.tex
	pdflatex $<
