#!/bin/bash -x

for n in 0 4 5 6; do
    for clo in KH PSE2 PSE3; do
        xyz=${n}w,$clo.xyz
        gx2xyz --no-dummies ${n}w,$clo.gx > $xyz
        mpiexec -n 16 ~/darcs/bgy3d-wheezy/guile/runbgy.scm \
            --solvent="water, PR-SPC/E" \
            --solute="uranyl, ${n}w, pcm" \
            --norm-tol=1e-14 \
            --dielectric=78.4 \
            --rho=0.0333295 \
            --beta=1.6889 \
            --L=10 \
            --N=96 \
            --hnc \
            --closure=$clo \
            --solute-geometry=$xyz \
            energy
        mpiexec -n 16 ~/darcs/bgy3d-wheezy/guile/runbgy.scm \
            --solvent="water, PR-SPC/E" \
            --solute="uranyl, ${n}w, pcm" \
            --norm-tol=1e-14 \
            --dielectric=78.4 \
            --rho=0.0333295 \
            --beta=1.6889 \
            --L=10 \
            --N=96 \
            --hnc \
            --closure=$clo \
            --solute-geometry=$xyz \
            rdf  "U" g0.bin g1.bin > ${n}w,$clo,U.rdf
        mpiexec -n 16 ~/darcs/bgy3d-wheezy/guile/runbgy.scm \
            --solvent="water, PR-SPC/E" \
            --solute="uranyl, ${n}w, pcm" \
            --norm-tol=1e-14 \
            --dielectric=78.4 \
            --rho=0.0333295 \
            --beta=1.6889 \
            --L=10 \
            --N=96 \
            --hnc \
            --closure=$clo \
            --solute-geometry=$xyz \
            rdf  "OU" g0.bin g1.bin > ${n}w,$clo,OU.rdf
    done
done
