import numpy as np
import itertools
import sys
import os

import my_utils
import savin_doyle

# INPUT LIST OF TRAJECTORIES (i.e. LIST OF (NFx3) NDARRAYS)
# syntax: python get_msd.py <path-to-infile>
# does this require timestep?
# what should debug do?
# what is skip?


timestep = 0.5
debug = True
skip = 1


def read_infile(infile):
    # convert contents of infile to list of (NFx3) ndarrays
    # input format is:
    # ith_trj xu yu zu s
    data = np.genfromtxt(infile)
    all_trajs = []
    for ith, trj in itertools.groupby(data, lambda x:x[0]):
        xu_yu_zu = np.array(list(trj))[:, 1:4]
        all_trajs.append(xu_yu_zu)
    return all_trajs


def process_dump(infile, timestep_fs):

    try:
        # read dump file and find regions
        list_of_trajs = read_infile(infile)
        print ('trjs in file: {}'.format(len(list_of_trajs)))
    except:
        print ('read_file failed!')
        raise

    # get msd and errors
    outfile = '{}.msdout'.format(os.path.abspath(infile))
    x = savin_doyle.estimators_to_file(list_of_trajs, outfile)

    # write msd and errors to file
#    with open(os.path.join(dirname, outfile), 'w') as fout:
#        pass
#    print ('Files written in: {}'.format(dirname))

    return


# get path from command line argument
if 0:
    import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    root.withdraw()
    dir_path = tkFileDialog.askdirectory()
else:
    dir_path = sys.argv[1]

# find the infile by pattern matching
matches = my_utils.find_files(dir_path, '_x_*dat')

# process the file
for infile in matches:
    try:
        process_dump(infile, timestep)
    except:
        print ("Error processing {}".format(infile))
