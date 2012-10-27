all: Kurzanleitung.html Kurzanleitung.pdf

%.html: %.txt
	pandoc $(^) -o $(@)

%.pdf: %.txt
	pandoc $(^) -o $(@)