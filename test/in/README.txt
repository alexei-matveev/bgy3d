
Information of generating input g2 files
========================================

Hydrogen chloride
-----------------

  g2 files: g2HH, g2HCl and g2ClCl have been exsisted in the original
copy of the code (with original file names: g2H, g2HCl and g2Cl in test
subdirectory). No clue how one could reproduce them.


Carbon disulfide
----------------

  g2 files: g2CC, g2CS and g2SS are obtained from the RDF data by MD
runs with LAMMPS (version: 14 May 2012), in.CS2 and CS2.data are inputs
for MD simulation. In which:

  in.CS2 is the main input script for LAMMPS which controls the
parameter settings, one can see the comments in the script and read the
manual of LAMMPS (http://lammps.sandia.gov/doc/Manual.html) for details.
As cited by Jager's thesis, the force field parameters are from [1]. But
it should be noted that since the harmonic intramolecular potential
expression described in [1] is not compatible in the current version of
LAMMPS, a simple derivation was applied to use the bond parameters. See
the comments in "bond section" of in.CS2.

  CS2.data recored the initial structure of 80 CS2 molecules in a cubic
box with the length of 20 angstroms in each dimension and is read in
when executing LAMMPS with the command "read_data CS2.data" from in.CS2.
Those molecules are created by gencs2.py randomly but without distance
check, so minimization runs or long-enough equilibration runs are
necessary.

  MD calculation runs at temperature T = 360 K with timestep = 1.0 fs
and 125,000 steps for equilibration and 1.2 million steps for production.
  RDF data is sampled every 500,000 steps by averaging the values of
each 500 steps for 1000 times.


Reference
=========
1. Zhu, S. B., Lee, J. and Robinson, G. W. Molecular dynamics simulation
    of liquid carbon disulphide with a harmonic intramolecular potential.
   Molecular Physics 65, 65-75, doi:10.1080/00268978800100851 (1988).
