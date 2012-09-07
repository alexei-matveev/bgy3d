#!/usr/bin/env python
import bgy3d
import contf
import sys
import os

def from_file(path):
    """
    Check the file suffix provided  by path, then use different method
    to read data.
    """

    filename = os.path.basename(path)

    base, suffix = os.path.splitext(filename);

    if suffix == '.bin':
        g = bgy3d.from_file(path)
    elif suffix == '.m':
        g = contf.m2dat(path)
    else:
        print 'Unknown file suffix.'
        exit()

    return g

def moments(path):
    """
    Computes and prints moments 0-2 to stdout.
    """

    g = from_file(path)

    h = 1.0 - g

    m1 = bgy3d.moments1(h)

    # Get the center of distribution
    center = m1[1:4] / m1[0]

    # Use center to compute 2nd momenta
    m2 = bgy3d.moments2nd(h, center)

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

# main function
if __name__ == "__main__":
    map(moments, sys.argv[1:])
