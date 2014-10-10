TEXFILE=SEQLinkagePosterASHG2014

all: $(TEXFILE).pdf; make clean

.PHONY: all clean

$(TEXFILE).pdf: $(TEXFILE).tex
	pdflatex --enable-write18 $(TEXFILE)
	@if(grep "Citation" $(TEXFILE).log > /dev/null);\
	then \
		bibtex $(TEXFILE);\
		pdflatex $(TEXFILE);\
	fi

	@if(grep "Rerun" $(TEXFILE).log > /dev/null);\
	then \
		pdflatex $(TEXFILE);\
	fi

clean:
	rm -f $(TEXFILE).aux $(TEXFILE).bbl $(TEXFILE).blg $(TEXFILE).log
