all: uo2.pdf

# Auxilary files are put there:
d=./uo2.d

#
# There are a few figures in this document:
#
$d/fig1.svg: $d/mkfig.py $d/gw96_rdf.dat $d/kl1_rdf.dat $d/kl2_rdf.dat $d/pm13_rdf.dat $d/spc_rdf.dat
	(cd $(@D); ./mkfig.py $(@F))

# LyX uses Inkscape to convert SVG -> PDF:
figures: $d/fig1.svg

.PHONY: figures

#
# PDF2 format is produced by pdfTeX
#
%.pdf: %.lyx bibuo2.bib figures
	lyx --export pdf2 $(<)

%.pdf: %.svg
	inkscape $(<) --export-dpi=300 --export-pdf=$(@)
