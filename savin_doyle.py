import numpy as np
import os

# INPUT LIST OF TRAJECTORIES (i.e. LIST OF (NFx3) NDARRAYS)


def get_weights(ni_list, n, theta0):
    # the only valid trajectories are when ni>=n
    valid = [i for i,ni in enumerate(ni_list) if ni>=n]
    ni = np.array(ni_list)[valid]
    sum_ni = np.sum(ni, dtype=float)
    th = sum_ni / theta0

    # get qi, si
    qi = [(nii+1) - n for nii in ni]
    si = [n if (q >= n) else q for q in qi]

    # weighting factors
    wi_list = np.array([1.*i/v/sum_ni for i,v in zip(ni, qi)])

    # for blocking method
    q_ = [i//v for i,v in zip(qi, si)]
    wi_list_ = np.array([1.*i/v/sum_ni for i,v in zip(ni, q_)])

    w = {}
    w['valid'] = valid
    w['ni'] = ni
    w['sum_ni'] = sum_ni
    w['theta'] = th
    w['qi'] = qi
    w['si'] = si
    w['wi'] = wi_list
    w['qi_'] = q_
    w['wi_'] = wi_list_

    return w


def get_estimator_v1(xyz, n, w):
    xyz = np.array(xyz)
    numparticles = np.shape(xyz)[0]

    disp = xyz[:, n:, :] - xyz[:, :-n, :]

    D2 = np.power(disp, 2)
    D4 = np.power(disp, 4)
    M1 = D2.mean(axis=(0,1))
    mean_D4 = D4.mean(axis=(0,1))
    M2 = mean_D4 / 3. - np.power(M1, 2)

    # this is the actual blocking method
    newD2 = D2[:, 0:w['qi_'][0]*w['si'][0], :].reshape((numparticles, w['qi_'][0], w['si'][0], 3), order='C')
    D2_ = np.mean(newD2, axis=2)  # average over the blocks

    mean_D2_1 = D2_.mean(axis=(0, 1))
    mean_D2_2 = (np.power(D2_, 2)).mean(axis=(0, 1))

    # variance when all observations are identically distributed (assumption of homogeneity)
    term2 =  (mean_D2_2 - mean_D2_1**2)  * (w['wi_'][0] / (1.-w['wi_'][0]))

    # total expression for variance of msd
    var_M1 = M2 / numparticles + term2

    return M1, var_M1


def get_estimator_v2(xyz, n, w):

    # This is the xyz data for the valid trajectories
    valid_xyz = np.array(xyz)[w['valid']]

    Nb = (len(w['valid']))  # number of particles in view

    # reset these for every value of lagtime
    d2 = 0
    d4 = 0
    d2_1 = 0
    d2_2 = 0

    # loop over each sub-trajectory i
    for i in range(Nb):

        valid_xyzi = np.array(valid_xyz[i])

        di = valid_xyzi[n:, :] - valid_xyzi[:-n, :]
        di2 = np.power(di, 2)
        di4 = np.power(di, 4)
        d2 += np.sum(w['wi'][i] * di2, axis=0)
        d4 += np.sum(w['wi'][i] * di4, axis=0)

        # this is the actual blocking method
        d2_block = np.reshape(di2[0:w['qi_'][i]*w['si'][i],:], (w['qi_'][i], w['si'][i], 3), order='C')
        d2_block = (np.mean(d2_block, axis=1))
        d2_1 += np.sum(w['wi_'][i] * d2_block, axis=0)
        d2_2 += np.sum(w['wi_'][i] * np.power(d2_block, 2), axis=0)

    # Equation 8
    term1 = (d4/3 - d2**2) / Nb
    term2 = (d2_2 - d2_1**2)
    term3 = np.sum([i*v**2 for i,v in zip(w['qi_'], w['wi_'])])
    term4 = np.sum([i*v*(1-v) for i,v in zip(w['qi_'], w['wi_'])])
    var_M1 = (term1 + term2*term3/term4)

    return d2, var_M1

def estimators_to_file(list_of_trajs, outfile):

    ni_list = [len(i)-1 for i in list_of_trajs]
    max_ni = max(ni_list)
    theta0 = np.sum(ni_list)

    print ("Trajectories: {}  MaxDuration: {}".format(len(list_of_trajs), max_ni+1))

    def makelist():
        full = []
        full.extend(range(1, 100, 1))
        full.extend(range(100, 1000, 10))
        full.extend(range(1000, 10000, 100))
        full.extend(range(10000, 100000, 1000))
        full = [i for i in full if i<=max_ni]
        return full

    with open(outfile, 'w', buffering=1) as fout:
        fout.write("lagtime(frames)  MSDx MSDy MSDz  varMSDx varMSDy varMSDz  theta \n")

        # loop over lagtimes "n"
        for n in makelist():

            # if trajectory lengths are all equal, we can average and forego the loops
            if len(set(ni_list))==1:

                w = get_weights(ni_list, n, theta0)

                d2, var_M1 = get_estimator_v1(list_of_trajs, n, w)

                # write to file
                fout.write("{:6} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:6.3f}\n".format(n, d2[0],d2[1],d2[2], var_M1[0],var_M1[1],var_M1[2], w['theta']))

            # if trajectory lengths differ, then analyze each trajectory in a loop
            else:

                w = get_weights(ni_list, n, theta0)

                # Stop the loop if less than two molecules are left or if theta < 0.5
                Nb = (len(w['valid']))  # number of particles in view
                if (Nb < 2) or (w['theta'] < 0.5):
                    break

                d2, var_M1 = get_estimator_v2(list_of_trajs, n, w)

                # write to file
                fout.write("{:6} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:6.3f}\n".format(n, d2[0],d2[1],d2[2], var_M1[0],var_M1[1],var_M1[2], w['theta']))


    #out = np.genfromtxt(outfile, skip_header=1, names=['lagtime', 'M1x', 'M1y', 'M1z', 'var_M1x', 'var_M1y', 'var_M1z', 'theta'])

    return

