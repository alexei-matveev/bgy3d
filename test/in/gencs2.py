import numpy as np
import random as rdm

XYZ = 0

def create_latt(a, b, c):
    '''
    Create the lattice matrix
    '''
    latt = np.array([[float(a), 0.0, 0.0],
                     [0.0, float(b), 0.0],
                     [0.0, 0.0, float(c)]])
    
    return latt

def assign_molecule(N, latt, bond):
    '''
    Generate a rigid 3-site molecule at a random place in the lattice
    with random orientation, return the coordinates of all three atoms
    '''
    rdm.seed(118245)

    # two random numbers, R for molecule palce, rvec for orientation
    R = np.random.rand(N, 3);
    rvec = np.random.rand(N, 3);

    # normalization of the orientation random
    rvec2 = rvec[:, 0]**2 + rvec[:, 1]**2 + rvec[:, 2]**2
    rnorm = 1.0 / np.sqrt(rvec2)

    # the norm of rvec should be 1
    rvec = rvec * rnorm.reshape(-1, 1)

    # coordinates for center atom
    R = np.dot(R, np.transpose(latt))

    # coordinates for the rest two atoms
    R1 = rvec * bond + R

    R2 = rvec * bond * (-1.0)  + R

    return R, R1, R2


def output(N, charge, R, R1, R2):
    '''
    output the atom information in XYZ format or LAMMPS data format
    '''
    # output as XYZ format
    if (XYZ):
        print N
        print 
        for i in range(N):
            print "C %12.6f %12.6f %12.6f " % (R[i, 0], R[i, 1], R[i, 2])
            print "S %12.6f %12.6f %12.6f " % (R1[i, 0], R1[i, 1], R1[i, 2])
            print "S %12.6f %12.6f %12.6f " % (R2[i, 0], R2[i, 1], R2[i, 2])
    # output LAMMPS data format
    else:
        print "%10d atoms" % (3 * N)
        print "%10d bonds" % (3 * N)
        print "%10d angles" % (N)
        print "%10d dihedrals" %0
        print "%10d impropers\n" %0

        print "%10d atom types" %2
        print "%10d bond types" %2
        print "%10d angle types" %1
        print "%10d dihedral types" %0
        print "%10d improper types\n" %0

        xlo = 0.0
        ylo = 0.0
        zlo = 0.0
        xhi = 20.0
        yhi = 20.0
        zhi = 20.0

        print "%12.6f  %12.6f xlo xhi" %(xlo, xhi)
        print "%12.6f  %12.6f ylo yhi" %(ylo, yhi)
        print "%12.6f  %12.6f zlo zhi" %(zlo, zhi)

        print "\nMasses\n"
        print "%5d %10.4f # C" %(1, 12.0107)
        print "%5d %10.4f # S" %(2, 32.065)

        print "\nAtoms\n"
        for i in range(N):
            print "%-5d %-5d %-5d %10.4f %12.5f %12.5f %12.5f" % (3 * i + 1, i + 1, 1, charge[0], R[i, 0], R[i, 1], R[i, 2])
            print "%-5d %-5d %-5d %10.4f %12.5f %12.5f %12.5f" % (3 * i + 2, i + 1, 2, charge[1], R1[i, 0], R1[i, 1], R1[i, 2])
            print "%-5d %-5d %-5d %10.4f %12.5f %12.5f %12.5f" % (3 * i + 3, i + 1, 2, charge[2], R2[i, 0], R2[i, 1], R2[i, 2])

        print "\nBonds\n"
        for i in range(N):
            print "%-5d %-5d %-5d %-5d" % (3 * i + 1, 1, 3 * i + 1, 3 * i + 2)
            print "%-5d %-5d %-5d %-5d" % (3 * i + 2, 1, 3 * i + 1, 3 * i + 3)
            print "%-5d %-5d %-5d %-5d" % (3 * i + 3, 2, 3 * i + 2, 3 * i + 3)

        print "\nAngles\n"
        for i in range(N):
            print"%-5d %-5d %-5d %-5d %-5d" % (i+1, 1, 3 * i + 2, 3 * i + 1, 3 * i + 3)

if __name__ == "__main__":


    ma = create_latt(20, 20, 20)

    r, r1, r2 = assign_molecule(80, ma, 1.56)

    charge = np.array([-0.308, 0.154, 0.154])

    output(80, charge, r, r1, r2)

