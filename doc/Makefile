all: \
	Kurzanleitung.html Kurzanleitung.pdf\
	README.html README.pdf

%.html: %.txt
	pandoc $(^) -o $(@)

%.pdf: %.txt
	pandoc --latex-engine=xelatex $(^) -o $(@)

#
# PDF4 format is produced by XeTeX
#
%.pdf: %.lyx
	lyx --export pdf4 $(<)