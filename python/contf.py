
import os
import bgy3d as bgy
import numpy as np
import pylab as plt

def m2dat(path):
    '''
    Convert vec*.m file to plain txt data file, defaulted new filename is vec*.dat.

    '''
    
    if os.path.isfile(path): 
        # Get the file name, by removing the directory prefix
        filename = path.split(os.sep)[-1]

        # Remove the suffix of filename
        fileprefix = filename.split(".")[-2]

        # Default newfilename
        newfile = fileprefix + ".dat"
        # Check whether file with the same name exists
        if os.path.isfile(newfile):
            print "file '" + newfile + "' already exists"
            print "overwrite/continue to use/write to new? (o/c/n)"
            try:
                opt = raw_input()
                inputting = True
                while inputting:
                    if opt == 'o':
                        inputting = False
                        # Remove the first and last line of the .m file
                        # redirect the content to the new file.
                        # cmdstring = "sed -e '1d; $d' " + path + " > " + newfile
                        # os.system(cmdstring)
                        continue
                    elif opt == 'n':
                        try:
                            newfile = raw_input("Enter a new file name:")
                            inputting = False
                        except KeyboardInterrupt:
                            print "\nUser terminated."
                        except EOFError:
                            print "\nunexpected EOF."
                        except:
                            print '\nSome error/exception occurred.'
                    elif opt == 'c':
                        inputting = False
                    elif opt == 'q':
                        inputting = False
                        exit()
                    else:
                        opt = raw_input("Unknown options, enter 'y' or 'n' ('q' to quit)")


            except KeyboardInterrupt:
                print "\nUser terminated."
                exit()
            except EOFError:
                print "\nunexpected EOF."
                exit()
            #except:
            #    print '\nSome error/exception occurred.'

        # Remove the first and last line of the .m file
        # redirect the content to the new file.
        cmdstring = "sed -e '1d; $d' " + path + " > " + newfile
        os.system(cmdstring)

        # Load the plain txt file
        f = file(newfile,'r')
        vec = np.loadtxt(f)
        f.close()

        # Reshape the array
        n = bgy.root3(np.size(vec))
        vec.shape = (n, n, n)

        return vec



    else:
        print "'" + path + "' not exist." 
        exit()

def get_plane(g, plane = 3):
    '''
    Get the plane in which the center point of grid box locates
    1 : plane x = 0
    2 : plane y = 0
    3 : plane z = 0 (default)
    '''

    # Get the zero and first moments 
    m = bgy.moments1(1.0 - g)

    # Round the center of the grid to the nearest integer
    center = np.empty(3,dtype='float')
    gcenter = np.empty(3, dtype='int')
    center = m[1:4] / m[0]
    for i in range(3):
        gcenter[i] = int(np.round(center[i]))

    if plane == 3:
        # z plane
        pd = g[:, :, gcenter[2]]
    elif plane == 2:
        # y plane
        pd = g[:, gcenter[1], :]
    elif plane == 1:
        # x plane
        pd = g[gcenter[0], :, :]
    else:
        print "unknown plane."
        exit()

    return pd

def contrf_plt(pd, interval = (-10, 10)):
    '''
    Draw the contourf plot of selected plane, defaulted interval is (-10, 10)
    '''

    NX = pd.shape[1]
    NY = pd.shape[0]

    x = np.linspace(interval[0], interval[1], NX)
    y = np.linspace(interval[0], interval[1], NY)

    cs = plt.contourf(x, y, pd) 
    cbar = plt.colorbar(cs)
    plt.show()

def test(path):
    import MA

    g = m2dat(path)

    pd = get_plane(g)

    contrf_plt(pd)




if __name__ == "__main__":

    test("/home/lib/BGY3DM/version_BL_py/python/vec.m")



