""" View BEM head model
    - cortex + 3 BEM layers with EEG and MEG sensors locations
"""
print __doc__

from mesh import Mesh
import numpy as np
from mayavi import mlab

# define some nice colors for the surfaces
head_col = (0.95, 0.83, 0.83) # light pink
skull_col = (0.91, 0.89, 0.67)
brain_col = (0.67, 0.89, 0.91) # light blue
cortex_col = (0.68, 0.68, 0.68) # grey

head = Mesh("model/head.vtk")
brain = Mesh("model/skull.vtk")
skull = Mesh("model/brain.vtk")
cortex = Mesh("model/cortex.vtk")

mlab.clf()

head.plot(color=head_col, opacity=0.3)
skull.plot(color=skull_col, opacity=0.3)
brain.plot(color=brain_col, opacity=0.3)
cortex.plot(color=cortex_col, opacity=1)

# View EEG electrodes
electrodes = np.loadtxt('model/eeg_channels_locations.txt')
mlab.points3d(electrodes[:,0], electrodes[:,1], electrodes[:,2], opacity=0.5, scale_factor=6)

# View MEG squids
squids = np.loadtxt('model/meg_channels_locations.squids')
mlab.quiver3d(squids[:,0], squids[:,1], squids[:,2],
              -squids[:,3], -squids[:,4], -squids[:,5],
              opacity=0.5, scale_factor=10, mode='cone')

mlab.show()
