import numpy as np
import itertools
import my_utils
import read_lammps_dump
import sys
import os
import matplotlib.pyplot as plt

# python 1_format_dump.py <NF> <timestep> <dumpfile>
#dumpfile = sys.argv[1]
#timestep = 0.5   # timestep in fs
#NF = 200       # integer or 'all'


NF, timestep, dumpfile = sys.argv[1:4]

basedir = os.path.dirname(dumpfile)
timestep = float(timestep)
try:
    NF = int(NF)
except:
    pass


# crystal parameters for NU-1000
crystal = {}
crystal['a'] = 39.97
crystal['b'] = 40.00
crystal['c'] = 16.58 * 2
crystal['alpha'] = 90
crystal['beta'] = 90
crystal['gamma'] = 120


def coord_transform(unwrapped, crystal):

    Mfwd = my_utils.xyz2fracM(crystal)
    Mrev = my_utils.frac2xyzM(crystal)

    wrapped = my_utils.wrapcoords(unwrapped, crystal)
    zw = wrapped[:,:,2]

    # get xyz coordinates of the four corners
    P1 = np.dot(Mrev, np.array([0.0, 0.0, 0]))[np.newaxis, np.newaxis, :]
    P2 = np.dot(Mrev, np.array([1.0, 0.0, 0]))[np.newaxis, np.newaxis, :]
    P3 = np.dot(Mrev, np.array([1.0, 1.0, 0]))[np.newaxis, np.newaxis, :]
    P4 = np.dot(Mrev, np.array([0.0, 1.0, 0]))[np.newaxis, np.newaxis, :]

    # calculate deltax and deltay to each corner
    deltas1 = (wrapped - P1)[:,:,:2]
    deltas2 = (wrapped - P2)[:,:,:2]
    deltas3 = (wrapped - P3)[:,:,:2]
    deltas4 = (wrapped - P4)[:,:,:2]

    # calculate total displacement (s2=x2+y2) to each corner
    s1 = np.sqrt(np.sum(deltas1**2, axis=2))
    s2 = np.sqrt(np.sum(deltas2**2, axis=2))
    s3 = np.sqrt(np.sum(deltas3**2, axis=2))
    s4 = np.sqrt(np.sum(deltas4**2, axis=2))

    # find shortest distance to nearest corner
    s = np.minimum(s1, s2)
    s = np.minimum(s, s3)
    s = np.minimum(s, s4)

    szw = np.concatenate((s[:,:,np.newaxis], zw[:,:,np.newaxis]), axis=2)

    return s


def groupby_s(unwrapped, crystal):

    s = coord_transform(unwrapped, crystal)

    unwrapped_s = np.concatenate((unwrapped, s[:,:,np.newaxis]), axis=2)
    NF, Nmolec, dims = np.shape(unwrapped_s)

    meso_trj = []
    micro_trj = []
    composite_trj = []

    for i in range(Nmolec):
        trj = unwrapped_s[:, i, :]
        composite_trj.append(trj)

        gb = itertools.groupby(trj, lambda x: x[3]<=17)
        for in_meso, items in gb:
            if in_meso==True:
                meso_trj.append(list(items))
            else:
                micro_trj.append(list(items))

    outfile = os.path.join(basedir, '_x_micro_NF{}.dat'.format(NF))
    with open(outfile, 'w') as fout:
        fout.write('# ith xu yu zu s\n')
        for i,xyz in enumerate(micro_trj):
            for xyz_ in xyz:
                fout.write('{} {} {} {} {}\n'.format(i, *xyz_))

    outfile = os.path.join(basedir, '_x_meso_NF{}.dat'.format(NF))
    with open(outfile, 'w') as fout:
        fout.write('# ith xu yu zu s\n')
        for i,xyz in enumerate(meso_trj):
            for xyz_ in xyz:
                fout.write('{} {} {} {} {}\n'.format(i, *xyz_))

    outfile = os.path.join(basedir, '_x_composite_NF{}.dat'.format(NF))
    with open(outfile, 'w') as fout:
        fout.write('# ith xu yu zu s\n')
        for i,xyz in enumerate(composite_trj):
            for xyz_ in xyz:
                fout.write('{} {} {} {} {}\n'.format(i, *xyz_))

    return


t_ps, unwrapped = read_lammps_dump.read_dumpfile(dumpfile, timestep_fs=timestep, NF=NF)
groupby_s(unwrapped, crystal)

