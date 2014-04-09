all: uo2.pdf

# Auxilary files are put there:
d=./uo2.d

#
# There are a few figures in this document:
#
$d/fig2.svg: $d/mkfig2.py $d/gw96_rdf.dat $d/kl1_rdf.dat $d/kl2_rdf.dat $d/pm13_rdf.dat $d/spc_rdf.dat
	(cd $(@D); ./mkfig2.py $(@F))

#
# PDF2  format  is  produced  by  pdfTeX.  LyX  uses  rsvg-convert  or
# Inkscape to convert SVG -> PDF:
#
%.pdf: %.lyx %.bib $d/fig2.svg
	lyx --export pdf2 $(<)

%.pdf: %.svg
	inkscape $(<) --export-dpi=300 --export-pdf=$(@)
