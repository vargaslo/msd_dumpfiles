import numpy as np
import os
import numbers

# # Function to read the dump file

def read_dumpfile(infile, timestep_fs=None, debug=False, NF=None, timesteps_per_frame=None):

    # require user to input the timestep
    if not isinstance(timestep_fs, numbers.Real):
        raise TypeError('Timestep must be entered as argument to read_file function')

    t = []
    xyz = []
    tmp_xyz = []
    frame = -1
    with open(infile, 'r') as fin:
        for line in fin:
            line_chunks = line.split()
            if line_chunks[0] != '#':
                if len(line_chunks)==4:
                    # store value of previous timestep
                    prev_timestep = timestep
                    ind, x, y, z = [float(i) for i in line_chunks]
                    tmp_xyz.append([x, y, z])
                if len(line_chunks)==2:

                    # read new timestep and number of molec
                    timestep, N = [int(i) for i in line_chunks]

                    # check if delta t is correct and if correct number of lines have been read
                    if frame>-1:
                        delta_t = timestep - prev_timestep
                        if len(tmp_xyz)==N and delta_t==timesteps_per_frame:
                            xyz.append(tmp_xyz)
                            t.append(timestep)
                            tmp_xyz = []
                            # check if desired number of frames have been written
                            if isinstance(NF, numbers.Real) and len(t)>=NF:
                                break
                        else:
                            print (prev_timestep, timestep, np.shape(tmp_xyz))
                            break

                    frame+=1


        # append final group
        if isinstance(NF, numbers.Real) and len(t)!=NF:
            if len(tmp_xyz)==N and delta_t==timesteps_per_frame:
                xyz.append(tmp_xyz)
                t.append(timestep)



    # Convert timestep to time
#    timesteps_per_frame = 1000
    fs_per_timestep = timestep_fs
    ps_per_fs = 0.001
    t_ps = (np.array(t) - t[0]) * fs_per_timestep * ps_per_fs

    # Get number of frames and reshape data array
    numframes, N, dims = np.shape(np.array(xyz))

    print ('{} frames read for {} molecules from {}\n'.format(numframes, N, os.path.abspath(infile)))

    return t_ps, xyz
