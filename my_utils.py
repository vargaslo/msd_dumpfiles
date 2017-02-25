import numpy as np
import sys

# # Useful functions

def find_files(loc, pattern):
    matches = []
    if (sys.version_info > (3, 5)):
        print (sys.version_info)
        import glob
        for filename in glob.glob('{}/**/{}'.format(loc, pattern), recursive=True):
            matches.append(filename)
    else:
        print (sys.version_info)
        import fnmatch
        import os
        for root, dirname, filenames in os.walk(loc):
            for filename in fnmatch.filter(filenames, pattern):
                matches.append(os.path.join(root, filename))
    return matches

def cosd(deg):
    return np.cos(deg/180. * np.pi)

def sind(deg):
    return np.sin(deg/180. * np.pi)

def xyz2fracM(crystal):
    a = crystal['a']
    b = crystal['b']
    c = crystal['c']
    alpha = crystal['alpha']
    beta  = crystal['beta']
    gamma = crystal['gamma']

    v = np.sqrt(1-cosd(alpha)**2-cosd(beta)**2-cosd(gamma)**2 + 2*cosd(alpha)*cosd(beta)*cosd(gamma))
    r1 = [1./a, -cosd(gamma)/a/sind(gamma), (cosd(alpha)*cosd(gamma)-cosd(beta)) / a/v/sind(gamma)]
    r2 = [0, 1./b/sind(gamma), (cosd(beta)*cosd(gamma)-cosd(alpha)) / b/v/sind(gamma)]
    r3 = [0, 0, sind(gamma)/c/v]
    M = np.array([r1, r2, r3])
    return M

def frac2xyzM(crystal):
    a = crystal['a']
    b = crystal['b']
    c = crystal['c']
    alpha = crystal['alpha']
    beta  = crystal['beta']
    gamma = crystal['gamma']

    v = np.sqrt(1-cosd(alpha)**2-cosd(beta)**2-cosd(gamma)**2 + 2*cosd(alpha)*cosd(beta)*cosd(gamma))
    r1 = [a, b*cosd(gamma), c*cosd(beta)]
    r2 = [0, b*sind(gamma), c*(cosd(alpha)-cosd(beta)*cosd(gamma))/sind(gamma)]
    r3 = [0, 0, c*v/sind(gamma)]
    M = np.array([r1, r2, r3])
    return M

def wrapcoords(alldata, crystal):

    Mfwd = xyz2fracM(crystal)
    Mrev = frac2xyzM(crystal)

    frac = np.inner(Mfwd, alldata)
    frac_mod = np.rollaxis(frac, 0,3)

    frac_wrap = frac_mod % 1

    xyz_wrap = np.inner(Mrev, frac_wrap)
    xyz_wrap_mod = np.rollaxis(xyz_wrap, 0, 3)

    return xyz_wrap_mod
