import numpy as np
import matplotlib.pyplot as plt


# conversion from frames to ps
timesteps_per_frame = 1000
fs_per_timestep = 0.5
ps_per_fs = 0.001
ps_per_frame = ps_per_fs * fs_per_timestep * timesteps_per_frame

names = ['frame', 'msdx', 'msdy', 'msdz', 'varx', 'vary', 'varz', 'theta']
data_all = np.genfromtxt('_x_composite_NF29608.dat.msdout', names=names)
data_tri = np.genfromtxt('_x_micro_NF29608.dat.msdout', names=names)
data_hex = np.genfromtxt('_x_meso_NF29608.dat.msdout', names=names)

plt.plot(ps_per_frame*data_hex['frame'], data_hex['msdz'], label='meso')
plt.plot(ps_per_frame*data_all['frame'], data_all['msdz'], label='composite')
plt.plot(ps_per_frame*data_tri['frame'], data_tri['msdz'], label='micro')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_aspect('equal')
plt.gca().set_xlabel('Time (ps)')
plt.gca().set_ylabel('MSD z (A2)')

plt.legend(loc='best')
plt.grid(True)

plt.show()
