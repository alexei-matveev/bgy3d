from evaluate import  get_interval, interplt_fit # poly_fit is not working well
import sys
import numpy as np
from scipy.interpolate import interp1d

def uo2_ni(path):
    """
    Load data file containing RDF (g) and number integral (ni), find first
    local maxmum g(r_max) and local minium g(r_min) in interpolated RDF , then
    find number integral ni(r_min)
    """
    # no way to know n and m from pure data file, for UO2 + H2O, n = 3 (number of
    # solute  sites) and m = 3 (number of solvent sites)
    n = 3
    m = 3

    # path only contains numbers and nothing else
    data = np.loadtxt(path)

    # first column is r
    r = data[:, 0]

    # 1:n*m+1 columns are g(r)
    g = data[:, 1:n*m+1]

    # n*m+2 : 2*n*m+2 columns are ni(r)
    ni = data[:, n*m+2:2*n*m+2]

    # Get a short part from g, which contains the first local maximum and minimum
    rinter, ginter = get_interval(r, g[:, 0], interval = (0.0, 3.5))

    # find the local maximum and minimum. x0 are the approximate x-coordinates,
    # max is flag for whether it's a maximum or a minimum
    xmi, ymi = interplt_fit(rinter, ginter, x0 = [2.4, 3.0], max = [True, False], kind = 'linear')

    # liear interpolat number integral
    ni_f = interp1d(r, ni[:, 0], 'linear')
    print "Number integral to first local minimum  ", xmi[1], " is ",  ni_f(xmi[1])

if __name__ == "__main__":
    map(uo2_ni, sys.argv[1:])
