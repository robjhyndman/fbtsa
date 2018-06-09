# Usually, only this line needs changing
PRESENTATION := fbtsa

# Source files
RFILES := $(wildcard *.R) $(wildcard *.rda)
TEXFILES := $(wildcard *.tex) $(wildcard *.sty) $(wildcard *.cls)
FIGFILES := $(wildcard figs/*)

all: $(PRESENTATION).pdf
.PHONY: all

$(PRESENTATION).pdf : $(PRESENTATION).Rmd $(RFILES) $(TEXFILES) $(FIGFILES)
	Rscript -e "rmarkdown::render('$<')"

clean:
	rm -fv $(PRESENTATION).pdf
	latexmk -c
	rm -fv $(PRESENTATION).tex
	rm -rf $(PRESENTATION)_cache
	rm -rf $(PRESENTATION)_files
