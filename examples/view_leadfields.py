"""Sample code to explore the leadfields computed with compute_leadfields.py
"""

import numpy as np
import scipy.io as io
from mayavi import mlab
from mesh import Mesh

# Load data
cortex = Mesh("model/cortex.vtk")
eeg_electrodes = np.loadtxt('model/eeg_channels_locations.txt')
eit_electrodes = np.loadtxt('model/eit_locations.txt')
ecog_electrodes = np.loadtxt('model/ecog_electrodes_locations.txt')
squids = np.loadtxt('model/meg_channels_locations.squids')

##############################################################################
# Load 4 leadfields
import openmeeg as om
G_meg = om.loadmat('leadfields/meg_leadfield.mat')
G_eeg = om.loadmat('leadfields/eeg_leadfield.mat')
G_eit = om.loadmat('leadfields/eit_leadfield.mat')
G_ip = om.loadmat('leadfields/ip_leadfield.mat')
G_ecog = om.loadmat('leadfields/ecog_leadfield.mat')

###############################################################################
# EEG leadfield
mlab.figure(1)
mlab.clf()

eeg_chan_idx = 28
cortex.plot(opacity=1, scalars=G_eeg[eeg_chan_idx, :])

# view EEG electrodes
mlab.points3d(eeg_electrodes[[eeg_chan_idx],0], eeg_electrodes[[eeg_chan_idx],1],
            eeg_electrodes[[eeg_chan_idx],2],
            opacity=0.5, scale_factor=12, color=(1,0,0))

###############################################################################
# ECoG leadfield
mlab.figure(2)
mlab.clf()
ecog_chan_idx = 0
cortex.plot(opacity=1, scalars=G_ecog[ecog_chan_idx, :])

# view EEG electrodes
mlab.points3d(ecog_electrodes[[ecog_chan_idx],0],
              ecog_electrodes[[ecog_chan_idx],1],
              ecog_electrodes[[ecog_chan_idx],2],
              opacity=0.5, scale_factor=8, color=(1,0,0))

###############################################################################
# MEG leadfield
mlab.figure(3)
mlab.clf()

meg_chan_idx = 30
cortex.plot(opacity=1, scalars=G_meg[meg_chan_idx,:])

# view MEG squids
mlab.quiver3d(squids[[meg_chan_idx], 0], squids[[meg_chan_idx], 1],
              squids[[meg_chan_idx], 2], -squids[[meg_chan_idx], 3],
              -squids[[meg_chan_idx], 4], -squids[[meg_chan_idx], 5],
              opacity=0.5, scale_factor=10, mode='cone')

###############################################################################
# EIT leadfield
mlab.figure(4)
mlab.clf()

# Generate sample current injection set up
n_electrodes = eit_electrodes.shape[0]
j_eit = np.zeros(n_electrodes)
idx_in = 10
idx_out = 13
j_eit[idx_in] = 1 # +1 current enters in idx_in electrode
j_eit[idx_out] = -1 # -1 current leaves in idx_out electrode

# Get corresponding potential
v_eit = np.dot(G_eit, j_eit)

# View results
electrodes_mesh = Mesh("model/eeg_channels_mesh.vtk")
electrodes_mesh.plot(scalars=v_eit)

mlab.points3d(eeg_electrodes[[idx_in, idx_out], 0],
              eeg_electrodes[[idx_in, idx_out], 1],
              eeg_electrodes[[idx_in, idx_out],2],
              opacity=0.5, scale_factor=12, color=(1,0,0))

###############################################################################
# Internal potential leadfield
int_elecs = np.loadtxt('model/internal_electrodes_locations.txt')
mlab.figure(5)
mlab.clf()

# Generate sample current generator configuration
n_dipoles = cortex.points.shape[0]
j_dipoles = np.zeros(n_dipoles)
dip_idx = 8000
dipoles = cortex.points
j_dipoles[dip_idx] = 1

# Compute the potential of the locations of the internal electrodes
v_int_elecs = np.dot(G_ip, j_dipoles)

# View results
mlab.points3d(int_elecs[:, 0], int_elecs[:, 1], int_elecs[:, 2], v_int_elecs)
mlab.points3d(dipoles[[dip_idx], 0], dipoles[[dip_idx], 1],
              dipoles[[dip_idx], 2], scale_factor=12, color=(1, 0, 0))
cortex.plot(color=(0.68, 0.68, 0.68), opacity=0.3)

mlab.show()
