# Usually, only this line needs changing
PRESENTATION := fbtsa

# Source files
RFILES := $(wildcard *.R) $(wildcard *.rda)
TEXFILES := $(wildcard *.tex) $(wildcard *.sty) $(wildcard *.cls)

all: $(PRESENTATION).pdf 
.PHONY: all

$(PRESENTATION).pdf : $(PRESENTATION).Rmd $(RFILES) $(TEXFILES)
	Rscript -e "rmarkdown::render('$<')"

clean:
	rm -fv $(PRESENTATION).pdf
	latexmk -c
	rm -fv $(PRESENTATION).tex
	rm -rfv $(PRESENTATION)_cache
	rm -rfv $(PRESENTATION)_files
