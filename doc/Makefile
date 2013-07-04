all: \
	Kurzanleitung.html Kurzanleitung.pdf\
	README.html README.pdf

%.html: %.txt
	pandoc $(^) -o $(@)

%.pdf: %.txt
	pandoc --latex-engine=xelatex $(^) -o $(@)

#
# There are a few figures in this document:
#
free-energy.pdf: \
	free-energy.dir/components.png \
	free-energy.dir/flow.png

#
# PDF4 format is produced by XeTeX
#
%.pdf: %.lyx
	lyx --export pdf4 $(<)

%.png: %.svg
	inkscape $(<) --export-dpi=300 --export-png=$(@)