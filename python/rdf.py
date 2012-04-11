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

def GCD(a, b):
    '''
    Return the GCD of two numbers

    Greatest Common Divisor algorithm, get from google
    http://code.activestate.com/recipes/577282-finding-the-gcd-of-a-list-of-numbers-aka-reducing-/

    >>> GCD(6,9)
    3
    >>> GCD(4.0, 8.0)
    4.0
    '''
    return abs(a) if b == 0 else GCD(b, a % b)

def GCD_List(list):
    '''
    Return the GCD of a list of numbers

    >>> a = [2, 4, 6, 8]
    >>> GCD_List(a)
    2
    >>> b = [3.0, 6.0, 9.0]
    >>> GCD_List(b)
    3.0
    '''
    return reduce(GCD, list)


def get_one(g, vec=[1, 1, 1], interval=(-10, 10)):
    '''
    Return radial components of distribution function along with relevant distances 
    in a single direction determined by vec, which defaulted as [1, 1, 1]

    >>> N, R = 20, 2.0
    >>> g3 = bgy.sinc_hole(N, R=R)
    >>> rs, gs = get_one(g3)

    >>> gs
    array([ 0.28111293,  1.19780139,  0.92713167,  1.01024038,  1.02594318,
            0.95473599,  1.05195458,  0.95099858,  1.03919611,  0.97466649])
    >>> 1.0 - np.sinc(rs / R)
    array([ 0.28111293,  1.19780139,  0.92713167,  1.01024038,  1.02594318,
            0.95473599,  1.05195458,  0.95099858,  1.03919611,  0.97466649])
    '''
    # Number of grids in each direction
    N = bgy.root3(np.size(g))

    # Get the zero and first moments 
    m = bgy.moments1(1.0 - g)

    # Get the center of the grid
    center = np.empty(3,dtype='float')
    center = m[1:4] / m[0]
    center = [round(center[i], 1) for i in range(3)]

    # Get the coordinate distances from each grid point to the center
    rgrid = grid_distance(N, center)

    # Size of each grid point
    gsize = (interval[1] - interval[0]) / float(N) 

    # Now convert rgrid to REAL unit
    rgrid *= gsize

    # Get GCD of three component of direction vector
    v_gcd = GCD_List(vec)

    # Then get the reduced direction vector
    vec_new = [vec[i] / v_gcd for i in range(3)]

    # cor_ini: initial grid for iteration in each direction
    # Nmax: max iterations times in each direction
    # if the center looks like:
    #   center = [10, 9.5, 10]
    # and shape(g) = [20, 20, 20]
    # then we have different strategies for vec = [1, 1, 1] and vec = [1, -1, 1]
    # if vec = [1, 1, 1], iteration begins from [10, 10, 10], and ends after 20 - 10 = 10 times
    # if vec = [1, -1, 1], iteration begins from [10, 9, 10], and ends after 9 - 0 + 1 = 10 times
    cor_ini = np.empty(3)
    Nmax = np.empty(3)

    # Get initial grid and max iteration times in each direction
    for i in range(3):
        # center coordinate is fractional
        if center[i] % 1 != 0:
            # begins from the left nearest grid if vec[i] < 0
            if vec_new[i] < 0:
                cor_ini[i] = np.round(center[i] - 0.5, 1)
                Nmax[i] = abs(int((cor_ini[i] + 1) / vec_new[i]))
            # begins from the right nearest grid if vec[i] > 0
            elif vec_new[i] > 0:
                cor_ini[i] = np.round(center[i] + 0.5, 1)
                Nmax[i] = abs(int((N - cor_ini[i]) / vec_new[i]))
            # always begins from the left nearest grid if vec[i] = 0
            else:
                cor_ini[i] = np.round(center[i], 1)
                Nmax[i] = float("inf")
        else:
        # simpler if center coordinates is integers, but iteration times in different directions might vary
            cor_ini[i] = center[i]
            if vec_new[i] < 0:
                Nmax[i] = abs(int((cor_ini[i] + 1) / vec_new[i]))
            elif vec_new[i] > 0:
                Nmax[i] = abs(int((N - cor_ini[i]) / vec_new[i]))
            else:
                Nmax[i] = float("inf")

    # Actual iteration times is the mininum of those in three direction
    # Might improve this if PBS applied
    N_iter = int(min(Nmax))

    # rs: distances
    # gs: value of distribution function
    rs = np.empty(N_iter)
    gs = np.empty(N_iter)

    for ii in range(N_iter):
        # Get grid coordinates in each direction
        cor_xx = cor_ini[0] + ii * vec_new[0]
        cor_yy = cor_ini[1] + ii * vec_new[1]
        cor_zz = cor_ini[2] + ii * vec_new[2]
        rs[ii] = rgrid[cor_xx, cor_yy, cor_zz]
        gs[ii] = g[cor_xx, cor_yy, cor_zz]

    return rs, gs


def get_rdf(g, dr, interval=(-10, 10)):
    '''
    Return  shell radius  spacing  with dr  and  radial components  of
    distribution functions.

    The interval  of grid box is  (-10, 10) by default  which is often
    used in BGY3D code.

    >>> N, R, dR = 20, 2.0, 1.0
    >>> g3 = bgy.sinc_hole(N, R=R)
    >>> rs, gs = get_rdf(g3, dr=dR)

    FIXME: the interface of get_rdf()  is somewhat cryptic, why do the
    following two arrays show no sign of similarity?

    #>> gs
    #>> 1.0 - np.sinc(rs / R)
    '''

    # Number of grids in each direction
    N = bgy.root3(np.size(g))

    # Get the zero and first moments 
    m = bgy.moments1(1.0 - g)

    # Get the center of the grid
    center = m[1:4] / m[0]
    # print "center=", center

    # Get the coordinate distances from each grid point to the center
    rgrid = grid_distance(N, center)

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
    # acc_g = np.zeros(int(Rmax / dr))

    # import time
    # print "A", time.time()

    # FIXME: why stopping with array operations here?
    ind = np.empty(rgrid.shape, dtype='int')

    # fill an integer array with data, converted to integer:
    ind[:, :, :] = np.floor(rgrid / dr)

    acc = np.zeros(1 + np.max(ind))

    for i, gi in zip(ind.flat, g.flat):
        acc[i] += gi

    # print "B", time.time()

    # Loop over the grid storing g and accumulate the value in each layer
    # for i in xrange(N):
    #     for j in xrange(N):
    #         for k in xrange(N):
    #             R = rgrid[i, j, k]
    #             # Ensure we'are still in the max shell
    #             if R <= Rmax:
    #                 # r<= R < r + dr
    #                 R_interval = layer_interval(R, dr)
    #                 # assert R_interval == bins[i, j, k]
    #                 # Accumulate the value of g in this layer
    #                 acc_g[R_interval] += g[i, j, k] 

    # print "C", time.time()

    # print "1=", acc_g
    # print "2=", acc
    # First mutiply accumulated g by volume of each grid occupy
    # acc_g *= gvol
    acc *= gvol

    # Then divide it by the volume of the layer between two shells
    for ii in range(int(Rmax / dr)):
        # Get the radius of each shell
        Rshell[ii] = ii * dr
        # acc_g[ii] /= delta_V(Rshell[ii], dr)
        acc[ii] /= delta_V(Rshell[ii], dr)


    # return Rshell, acc_g
    return Rshell, acc

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
    else:
        f = file("rdf.npz", "wb")

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
    # r, rdf = get_rdf(g, dr)
    r, rdf = get_one(g)


#    print r
#    print rdf

    save_file(r, rdf)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    # test("./g2H_32.bin")



