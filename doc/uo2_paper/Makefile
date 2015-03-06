all: uo2.pdf

# Auxilary files are put there:
d=uo2.d

#
# There are a few figures in this document:
#
$d/fig2.svg: $d/mkfig2.py $d/gw96_rdf.dat $d/kl1_rdf.dat $d/kl2_rdf.dat $d/pm13_rdf.dat $d/spc_rdf.dat
	(cd $(@D); python ./mkfig2.py $(@F))

$d/fig1.svg: $d/mkfig1.py $d/md-oh.dat $d/md-oo.dat $d/md-uh.dat $d/md-uo.dat $d/cspce-pm13_rdf.dat $d/pm13_rdf.dat $d/shspce-pm13_rdf.dat
	(cd $(@D); python ./mkfig1.py $(@F))

#
# PDF2  format  is  produced  by  pdfTeX.  LyX  uses  rsvg-convert  or
# Inkscape to convert SVG -> PDF:
#
%.pdf: %.lyx ../bibliography.bib $d/fig1.svg $d/fig2.svg $d/fig3.pdf
	lyx --export pdf2 $(<)

%.pdf: %.svg
	inkscape $(<) --export-dpi=300 --export-pdf=$(@)

# 300  DPI is a  moderate quality,  just make  sure that  the physical
# dimensions of the figure are reasonable:
$d/fig3.png: DPI=300
$d/fig1.png $d/fig2.png: DPI=300

%.png: %.svg
	inkscape $(<) --export-dpi=$(DPI) --export-png=$(@)

#
# Convert XYZ file to PNG using Jmol.  The flag -n tells Jmol not open
# a graphics  window. The  jmol.script says to  write the  output into
# jmol.png --- we have to rename it here:
#
%.png: %.xyz
	jmol $(<) -s $d/jmol.script -n
	mv jmol.png $(@)

bitmaps = \
	$d/4h2o,mm+rism.png \
	$d/5h2o,mm+rism.png \
	$d/6h2o,mm+rism.png \
	$d/4h2o,mm.png \
	$d/5h2o,mm.png \
	$d/6h2o,mm.png

$(bitmaps): $d/jmol.script

# This bitmap is overlayed with text annotation in fig3.svg:
$d/fig3-bitmap.png: $(bitmaps)
	montage -mode concatenate -tile 3x2 $(bitmaps) $(@)

# Therefore the additional dependence:
$d/fig3.pdf: $d/fig3-bitmap.png
