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
	free-energy.dir/components.pdf \
	free-energy.dir/flow.pdf

#
# PDF4 format is produced by XeTeX
#
%.pdf: %.lyx
	lyx --export pdf4 $(<)

%.pdf: %.svg
	inkscape $(<) --export-dpi=300 --export-pdf=$(@)