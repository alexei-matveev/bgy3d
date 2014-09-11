#!/bin/bash -x

L=10
N=256
rm -f g0.bin g1.bin
mpiexec="mpiexec"
for n in 0 4 5 6; do
    for clo in KH PSE2 PSE3; do
        rdf=${n}w,$clo,OU.rdf
        if [ -e $rdf ]; then
            continue
        fi
        xyz=${n}w,$clo.xyz
        gx2xyz --no-dummies ${n}w,$clo.gx > $xyz
        $mpiexec ~/darcs/bgy3d-wheezy/guile/runbgy.scm \
            --solvent="water, PR-SPC/E" \
            --solute="uranyl, ${n}w, pcm" \
            --norm-tol=1e-14 \
            --dielectric=78.4 \
            --rho=0.0333295 \
            --beta=1.6889 \
            --L=$L \
            --N=$N \
            --hnc \
            --closure=$clo \
            --solute-geometry=$xyz \
            --save-binary \
            energy
        $mpiexec ~/darcs/bgy3d-wheezy/guile/runbgy.scm \
            --solvent="water, PR-SPC/E" \
            --solute="uranyl, ${n}w, pcm" \
            --norm-tol=1e-14 \
            --dielectric=78.4 \
            --rho=0.0333295 \
            --beta=1.6889 \
            --L=$L \
            --N=$N \
            --hnc \
            --closure=$clo \
            --solute-geometry=$xyz \
            rdf  "U" g0.bin g1.bin > ${n}w,$clo,U.rdf
        $mpiexec ~/darcs/bgy3d-wheezy/guile/runbgy.scm \
            --solvent="water, PR-SPC/E" \
            --solute="uranyl, ${n}w, pcm" \
            --norm-tol=1e-14 \
            --dielectric=78.4 \
            --rho=0.0333295 \
            --beta=1.6889 \
            --L=$L \
            --N=$N \
            --hnc \
            --closure=$clo \
            --solute-geometry=$xyz \
            rdf  "OU" g0.bin g1.bin > $rdf
    done
done
