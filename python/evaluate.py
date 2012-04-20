
import numpy as np
import pylab as plt
from scipy.interpolate import interp1d
from scipy.optimize import fmin

# Set to 1 if need to show plots
PLOT = 1

def loadbin(path):
    '''
    Load binary file which saved by rdf.py
    '''
    f = file(path, "rb")
    fdata = np.load(f)
    r = fdata["R"]
    g = fdata["rdf"]
    f.close()

    return r, g

def loadtxt(path):
    '''
    Load plain text file
    '''
    f = file(path)
    fdata = np.loadtxt(f)
    r = fdata[:, 0]
    g = fdata[:, 1]
    f.close()

    return r, g


def get_interval(r, g, interval):
    '''
    Return "real" interval from r, by comparing the nearest value of dataset to input interval
    >>> r = np.linspace(0, 2.0 * np.pi, 10)
    >>> g = np.sin(r)
    >>> interval = (1, 3)
    >>> realr, realg = get_interval(r, g, interval)
    Real interval: [0.6981, 2.7925]
    '''

    # Convert list to numpy.array
    # so we can use array.argmin()
    r = np.array(r)
    g = np.array(g)

    # Substrating the upper and lower value of interval from r
    sub_lower = abs(r - interval[0])
    sub_upper = abs(r - interval[1])

    # Find the index of nearest r[] to the interval bounds
    idx_lower = sub_lower.argmin()
    idx_upper = sub_upper.argmin()

    # Real interval from data series
    print "Real interval: [%6.4f, %6.4f]" %( r[idx_lower],  r[idx_upper] )

    return r[idx_lower : idx_upper + 1], g[idx_lower : idx_upper + 1]


def interplt_fit(x, y, x0, max, kind = 5):
    '''
    Return local extrema of interpolated function near given x0[],
    array max[] indicates it's max or min.

    >>> x = np.linspace(0, 2.0 * np.pi, 10)
    >>> y = np.sin(x)
    >>> x0 = [1.5, 4.5]
    >>> max = [True, False]
    >>> interplt_fit(x, y, x0, max)
    Interpolating: 
    Optimization terminated successfully.
             Current function value: -0.999934
             Iterations: 11
             Function evaluations: 22
    local maximum is 0.9999 at x = 1.5705
    Optimization terminated successfully.
             Current function value: -0.999934
             Iterations: 13
             Function evaluations: 26
    local minimum is -0.9999 at x = 4.7126

    
    '''

    # Interpolate function of (x, y)
    spline_fit = interp1d(x, y, kind)

    # Interpolate function of (x, -y), which is used to find local maximum
    vspline_fit = interp1d(x, -1.0 * y, kind)

    # Generate data set of interpolate function 
    xx = np.linspace(x[0], x[-1], 50)
    yy = spline_fit(xx)

    # Arrays to store local minimum/maximum
    xmin = np.empty(0)
    ymin = np.empty(0)
    xmax = np.empty(0)
    ymax = np.empty(0)

    # count numbers
    maxcount = 0
    mincount = 0

    print "Interpolating: "
    for x0i, mi in zip(x0, max):
        # Local maximum
        if mi:
            xmax = np.append(xmax, fmin(vspline_fit, x0i))
            ymax = np.append(ymax, spline_fit(xmax[maxcount]))
            print "local maximum is %6.4f at x = %6.4f" % (ymax[maxcount], xmax[maxcount])
            maxcount += 1
        # Local minimum
        else:
            xmin = np.append(xmin, fmin(spline_fit, x0i))
            ymin = np.append(ymin, spline_fit(xmin[mincount]))
            print "local minimum is %6.4f at x = %6.4f" % (ymin[mincount], xmin[mincount])
            mincount += 1

    if maxcount == 0 and mincount == 0 :
        print " No local minimum/maximum found in [%6.4f, %6.4f]" % (x[0], x[-1])
        print " Choose another interval space."
        exit()


    if PLOT:
        # Plot original data and interpolate function
        plt.plot(x, y, 'gs-', xx, yy, 'b--', linewidth = 2)
        # Plot local extrema
        plt.plot(xmin[:], ymin[:], xmax[:], ymax[:], 'ro', ms = 8)
        plt.show()



def poly_fit(x, y, deg = 5):
    '''
    Return the local extreme of the polynomial function fitting to data series (x, y)

    >>> x = np.linspace(0, 2.0 * np.pi, 10)
    >>> y = np.sin(x)
    >>> find_extrm(x, y)
    Polynomial fitting: 
    local minimum is -1.0046 at x = 4.7241
    local maximum is 1.0046 at x = 1.5591

    '''

    # First fit the data series to a polynomial function
    coef1 = np.polyfit(x, y, deg) 

    # Get the first derivative of the fitted function
    pcoef1 = np.polyder(coef1, 1)

    # Get the second derivative
    pcoef2 = np.polyder(coef1, 2)

    # Find the roots for the equation of first derivative
    proots = np.roots(pcoef1)
    
    # Arrays to store local extrema
    xm = np.empty(np.size(proots))
    ym = np.empty(np.size(proots))

    # Number of found local extrema
    count = 0

    print "Polynomial fitting: "
    # Check we get the real root in the inverval
    for xr in proots:
        # Calculate the value of second derivative for each root
        der2 = np.polyval(pcoef2, xr)
        # Ensure it's real number and locates in the interval
        if xr.imag == 0 and xr >= x[0] and xr <= x[-1] :
            xm[count] = xr
            # jocal maximum if 2nd derivative is negative
            if der2 < 0:
                ym[count] = float(np.polyval(coef1, xm[count]))
                print "local maximum is %6.4f at x = %6.4f" % (ym[count], xm[count])
            # Local minimum if 2nd derivative is positive
            elif der2 > 0:
                ym[count] = float(np.polyval(coef1, xm[count]))
                print "local minimum is %6.4f at x = %6.4f" % (ym[count], xm[count])
            count += 1

    # Exit if no local extreme found
    if not count:
        print " No local minimum/maximum found in [%6.4f, %6.4f]" % (x[0], x[-1])
        print " Choose another interval space."
        exit()


    # Generate data set of fitted function 
    xx = np.linspace(x[0], x[-1], 50)
    yy = np.polyval(coef1, xx)
    
    if PLOT:
        # Plot original and fitted function data
        plt.plot(x, y, 'gs-', xx, yy, 'b--', linewidth = 2)
        # Plot the local maximum/minimum
        plt.plot(xm[:count], ym[:count], 'ro', ms = 6)
        plt.show()

def test(path):

    r, g = loadbin(path)
    # r, g = loadtxt(path)
    rinter, ginter = get_interval(r, g, interval = (3.5, 4.5))
    poly_fit(rinter, ginter)
    interplt_fit(rinter, ginter, [4], [True])

if __name__ == '__main__':

    import doctest
    # doctest.testmod()
    test("/home/lib/BGY3DM/version_BL_py/python/rdf.npz")
    # test("/home/lib/BGY3DM/version_BL_py/python/g2H")

