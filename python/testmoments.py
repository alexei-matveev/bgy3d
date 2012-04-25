
import bgy3d as bgy
import contf as ctf
import sys
import os

def checkfile(path):
    """
    Check the file suffix provided by path, then 
    use different method to read data.
    """

    filename = path.split(os.sep)[-1]

    filesuffix = filename.split(".")[-1]

    if filesuffix == 'bin':
        g = bgy.from_file(path)
    elif filesuffix == 'm':
        g = ctf.m2dat(path)
    else:
        print 'Unknown file suffix.'
        exit()

    return g

# main function
if __name__ == "__main__":

    path = sys.argv[1]

    g = checkfile(path)

    h = 1.0 - g

    m1 = bgy.moments1(h)

    # Get the center of distribution
    center = m1[1:4] / m1[0]

    # Use center to compute 2nd momenta
    m2 = bgy.moments2nd(h, center)

    print "<1> = ", m1[0]
    print "<x> = ", m1[1] / m1[0]
    print "<y> = ", m1[2] / m1[0]
    print "<z> = ", m1[3] / m1[0]
    print "<xy> = ", m2[0] / m1[0]
    print "<yz> = ", m2[1] / m1[0]
    print "<zx> = ", m2[2] / m1[0]
    print "<z^2 - 1/3 * r^2> = ", m2[3] / m1[0]
    print "<x^2 - y^2> = ", m2[4] / m1[0]
    print "<r^2> = ", m2[5] / m1[0]

