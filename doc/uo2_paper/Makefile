all: uo2.pdf

#
# There are a few figures in this document:
#

#
# PDF2 format is produced by pdfTeX
#
%.pdf: %.lyx bibuo2.bib
	lyx --export pdf2 $(<)

%.pdf: %.svg
	inkscape $(<) --export-dpi=300 --export-pdf=$(@)
