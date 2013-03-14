all: \
	Kurzanleitung.html Kurzanleitung.pdf\
	README.html README.pdf

%.html: %.txt
	pandoc $(^) -o $(@)

%.pdf: %.txt
	pandoc --latex-engine=xelatex $(^) -o $(@)