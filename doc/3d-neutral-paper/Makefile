all: 3d-neutral-paper.pdf

# Auxilary files are put there:
d=3d-neutral.d

#
# There are a few figures in this document:
#
$d/fig2.svg: $d/mkfig2.py $d/dipole.dat
	(cd $(@D); python ./mkfig2.py $(@F))

$d/fig1.svg: $d/mkfig1.py $d/mm.dat $d/fed.dat $d/qm.dat $d/scf-qm.dat
	(cd $(@D); python ./mkfig1.py $(@F))

$d/fig3.svg: $d/mkfig3.py $d/relax.dat
	(cd $(@D); python ./mkfig3.py $(@F))
#
# PDF2  format  is  produced  by  pdfTeX.  LyX  uses  rsvg-convert  or
# Inkscape to convert SVG -> PDF:
#
%.pdf: %.lyx %.bib $d/fig1.svg $d/fig2.svg $d/fig3.svg
	lyx --export pdf2 $(<)

# 300  DPI is a  moderate quality,  just make  sure that  the physical
# dimensions of the figure are reasonable:
$d/fig1.png $d/fig2.png $d/fig3.png: DPI=300

%.png: %.svg
	inkscape $(<) --export-dpi=$(DPI) --export-png=$(@)
