#!/usr/bin/python

import bgy3d as bgy
import numpy as np

def delta_V(r, dr):
    '''
    Return the volume of layer between shell(r + dr) and shell(r)

          4                  3   3
    dV = --- * PI * [(r + dr) - r ]
          3

    while r is shell radius and dr is bin width between two shells.

    Note the difference in precision:

    >>> r, dr = 1.0e6, 1.0e-6
    >>> 4. / 3. * np.pi * ((r + dr)**3 - r**3)
    12566638.696932279

    >>> delta_V(r, dr)
    12566370.614371737
    '''
    return 4.0 * np.pi * (3.0 * r * dr * (r + dr) + dr**3) / 3.0

def grid_distance(N, center):
    '''
    Return the distance of each grid point to the center of the box

    >>> r = grid_distance(3, (1, 1, 1))

    Middle plane:

    >>> r[1]
    array([[ 1.41421356,  1.        ,  1.41421356],
           [ 1.        ,  0.        ,  1.        ],
           [ 1.41421356,  1.        ,  1.41421356]])
    '''
    # Generate a N x N x N grid mesh
    t = np.mgrid[0:N, 0:N, 0:N]

    # return the square root of (x - centerx)^2 + (y - centery)^2 + (z -centerz)^2
    r = np.sqrt((t[0] - center[0])**2 + (t[1] - center[1])**2 + (t[2] - center[2])**2)

    # shape(r) = (N, N, N)
    return r

def layer_interval(r, dr):
    '''
    Return  number N, if  r locates  between Nth  and N+1th  shell, of
    which layer distance between each layer is dr. FIXME: maybe return
    an integer instead?

    >>> layer_interval(1.1, 1.0)
    1.0
    >>> layer_interval(0.9, 1.0)
    0.0
    >>> layer_interval(0.5, 1.0)
    0.0
    '''
    from math import floor

    return floor(r / dr)


def get_rdf(g, dr):
    '''
    Return shell radius spacing with dr and radial components of distribution functions
    '''

    # Number of grids in each direction
    N = bgy.root3(np.size(g))

    # Get the zero and first moments 
    m = bgy.moments1(1.0 - g)

    # Get the center of the grid
    center = m[1:4] / m[0]

    # Get the coordinate distances from each grid point to the center
    rgrid = grid_distance(N, center)

    # the interval of grid box is [-10 : 10]
    interval = [-10, 10]

    # Size of each grid point
    gsize = (interval[1] - interval[0]) / float(N) 

    # Now convert rgrid to REAL unit
    rgrid *= gsize

    # volume of each grid point
    gvol = gsize*gsize*gsize

    # Set the maximum shell radius equal to half of grid box length
    Rmax = (interval[1] - interval[0]) / 2.0

    # Arrays to store shell radius
    Rshell = np.empty(int(Rmax / dr))

    # Array to store accumulated g in each layer
    acc_g = np.zeros(int(Rmax / dr))

    # Loop over the grid storing g and accumulate the value in each layer
    for i in range(N):
        for j in range(N):
            for k in range(N):
                R = rgrid[i, j, k]
                # Ensure we'are still in the max shell
                if R <= Rmax:
                    # r<= R < r + dr
                    R_interval = layer_interval(R, dr)
                    # Accumulate the value of g in this layer
                    acc_g[R_interval] += g[i, j, k] 

    # First mutiply accumulated g by volume of each grid occupy
    acc_g *= gvol

    # Then divide it by the volume of the layer between two shells
    for ii in range(int(Rmax / dr)):
        # Get the radius of each shell
        Rshell[ii] = ii * dr
        acc_g[ii] /= delta_V(Rshell[ii], dr)


    return Rshell, acc_g

def save_file(r, g):
    '''
    Save Rshell and rdf in a single npz file,
    keyword "R" for Rshell and "rdf" for value of rdf
    '''
    import os
    inputting = True
    if os.path.isfile('./rdf.npz'):
        print "file rdf.npz already exists, overwrite it (y/n)"
        opt = raw_input()
        while inputting:
            if opt == 'y':
                f = file("rdf.npz", "wb")
                inputting = False
            elif opt == 'n':
                nfname = raw_input("Enter a new file name:")
                f = file(nfname, "wb")
                inputting = False
            else:
                opt = raw_input("Unknown options, enter 'y' or 'n'")
    np.savez(f, R = r, rdf = g)
    f.close()


def test(path):

    g = bgy.from_file(path)

    # For N = 32, try dr = 20 / 32
    dr = 0.625
    # For N = 128, try dr = 20 / 128
    # dr = 0.15625
    # For N = 256, try dr = 20 / 256
    # dr = 0.078125
    r, rdf = get_rdf(g, dr)

#    print r
#    print rdf

    save_file(r, rdf)

if __name__ == "__main__":
    
    test("./g2H_32.bin")



