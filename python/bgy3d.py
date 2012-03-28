import sys
import petsc4py

petsc4py.init(sys.argv)

from petsc4py import PETSc
from numpy import array, size, shape, mgrid, sum

VERBOSE = 0

def root3(n):
    """
    Poor man cubic root. Fails if the root is not integer.

    >>> root3(8)
    2
    """

    i = 0
    while i**3 < n:
        i = i + 1
    assert i**3 == n

    return i

def from_file(path):
    """
    Returns a PETSc.Vec converted to N x N x N numpy array.
    """

    viewer = PETSc.Viewer()
    viewer.createBinary(path, "r")

    # convert to an array immediately (FIXME: copy here):
    vec = array(PETSc.Vec.Load(viewer))

    # reshape it:
    n = root3(size(vec))
    vec.shape = (n, n, n)

    return vec

def moments1(h, r0=(0, 0, 0)):
    """
    Computes the  lowest momenta of  the distribution in  integer grid
    coordinates.

    >>> m = moments1(sinc_hole(64) - 1.0)
    >>> m[1:4] / m[0]
    array([ 31.5,  31.5,  31.5])
    """

    N = root3(size(h))

    #
    # Works much like range(N), but is an (integer) array:
    #
    t = mgrid[0:N]
    # print type(t), type(t[0])

    #
    # These are  arrays with three  axes, comfortable with  each other
    # and any  N x N x  N array but  consuming only a fraction  of the
    # space:
    #
    x = t[:, None, None] - r0[0]
    y = t[None, :, None] - r0[1]
    z = t[None, None, :] - r0[2]
    # print shape(x * y * z)

    #
    # First moments approximated by discrete integrals:
    #
    m0 = sum(1 * h)
    mx = sum(x * h)
    my = sum(y * h)
    mz = sum(z * h)

    if VERBOSE:
        print "<1> =", m0
        print "<x> =", mx / m0
        print "<y> =", my / m0
        print "<z> =", mz / m0

    # Return a flat array of four moments:
    return array((m0, mx, my, mz))

def moments2(h, r0=(0, 0, 0)):
    """
    Computes the  lowest momenta of  the distribution in  integer grid
    coordinates.

    >>> m = moments2(sinc_hole(64) - 1.0)
    >>> m[1:4] / m[0]
    array([ 31.5,  31.5,  31.5])
    """

    N = root3(size(h))

    #
    # FIXME: does it consume 3*N**3 memory?
    #
    r = mgrid[0:N, 0:N, 0:N]

    m0 = sum(h)
    m1 = [sum((r[i] - r0[i]) * h) for i in range(3)]

    if VERBOSE:
        for i, m in enumerate(m1):
            print "<r[%d]> =" % i, m / m0

    # Return a flat array of four moments:
    return array([m0] + m1)

def sinc_hole(N, a=None, R=None):
    """
    Returns an N x N x N distribution g(x) = 1 - sinc(|x - a| / R) for
    use in testing.

    >>> vol = sinc_hole(5)

    Middle plane:

    >>> vol[2]
    array([[ 1.18607792,  1.        ,  0.88411815,  1.        ,  1.18607792],
           [ 1.        ,  0.539657  ,  0.29800204,  0.539657  ,  1.        ],
           [ 0.88411815,  0.29800204,  0.        ,  0.29800204,  0.88411815],
           [ 1.        ,  0.539657  ,  0.29800204,  0.539657  ,  1.        ],
           [ 1.18607792,  1.        ,  0.88411815,  1.        ,  1.18607792]])
    """

    from numpy import sinc, sqrt

    if a is None:
        # fp-division here:
        a = array([N - 1, N - 1, N - 1]) / 2.0

    if R is None:
        # R/N -> 0 for large N:
        R = sqrt(N)

    x = mgrid[0:N, 0:N, 0:N]
    r = sqrt((x[0] - a[0])**2 + (x[1] - a[1])**2 + (x[2] - a[2])**2)

    return 1.0 - sinc(r / R)

def test(path):
    g = from_file(path)

    #
    # Differential redistribution density:
    #
    h = 1.0 - g
    print moments1(h)
    print moments2(h)
    # print moments1(h, (16, 16, 16))
    # print moments2(h, (16, 16, 16))

#
# Run doctests by default, python path.py [-v]:
#
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    # test("../test/g2H.bin")
