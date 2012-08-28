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

    #
    # Python wrappers change from version to version, try to guess:
    #
    try:
        READ_MODE = PETSc.Viewer.Mode.READ
    except:
        READ_MODE = "r"

    viewer.createBinary(path, READ_MODE)

    try:
        vec = PETSc.Vec().load(viewer)
    except:
        vec = PETSc.Vec.Load(viewer)

    # convert to an array immediately (FIXME: copy here):
    vec = array(vec)

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

def moments2nd(h, r0 = (0, 0, 0)):
    """
    Computes the second momenta of the distribution in integer grid
    coordinates.

    >>> m1 = moments1(sinc_hole(64) - 1.0)
    >>> center = m1[1:4] / m1[0]
    >>> m2 = moments2nd(sinc_hole(64) - 1.0, center)
    >>> m2 / m1[0]
    array([ -8.08456013e-13,   3.06434477e-14,   8.46290161e-14,
             5.96045623e+02,   1.29571417e-11,   2.68220530e+03])

    """

    N = root3(size(h))

    t = mgrid[0:N]
 
    x = t[:, None, None] - r0[0]
    y = t[None, :, None] - r0[1]
    z = t[None, None, :] - r0[2]

    # Second moments
    # < x * y >
    mxy = sum(x * y * h)
    # < y * z >
    myz = sum(y * z * h)
    # < z * x >
    mzx = sum(z * x * h)
    # < z^2 - 1/3 * r^2 >
    mz2 = sum((2 * z**2 - x**2 + y**2) / 3.0 * h)
    # < x^2 - y^2 >
    mx2y2 = sum(( x**2 - y**2) * h)
    # < r^2 >
    mr2 = sum( ( x**2 + y**2 + z**2) *h )

    return array((mxy, myz, mzx, mz2, mx2y2, mr2))


def sinc_hole(N, a=None, R=None):
    """
    Returns an N x N x N distribution g(x) = 1 - sinc(|x - a| / R) for
    use in testing.

    >>> from numpy import round
    >>> vol = round(100 * sinc_hole(9, R=2.0))

    Middle plane:

    >>> vol[4]
    array([[  94.,   87.,   90.,   97.,  100.,   97.,   90.,   87.,   94.],
           [  87.,   94.,  110.,  119.,  121.,  119.,  110.,   94.,   87.],
           [  90.,  110.,  122.,  110.,  100.,  110.,  122.,  110.,   90.],
           [  97.,  119.,  110.,   64.,   36.,   64.,  110.,  119.,   97.],
           [ 100.,  121.,  100.,   36.,    0.,   36.,  100.,  121.,  100.],
           [  97.,  119.,  110.,   64.,   36.,   64.,  110.,  119.,   97.],
           [  90.,  110.,  122.,  110.,  100.,  110.,  122.,  110.,   90.],
           [  87.,   94.,  110.,  119.,  121.,  119.,  110.,   94.,   87.],
           [  94.,   87.,   90.,   97.,  100.,   97.,   90.,   87.,   94.]])
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
