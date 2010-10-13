""" View BEM head model
    - cortex + 3 BEM layers with EEG and MEG sensors locations
"""
print __doc__

from mesh import Mesh
import numpy as np
from enthought.mayavi import mlab

# define some nice colors for the surfaces
head_col = (0.95, 0.83, 0.83) # light pink
skull_col = (0.91, 0.89, 0.67)
brain_col = (0.67, 0.89, 0.91) # light blue
cortex_col = (0.68, 0.68, 0.68) # grey

head = Mesh("head.tri")
brain = Mesh("skull.tri")
skull = Mesh("brain.tri")
cortex = Mesh("cortex.tri")

mlab.clf()

head.plot(color=head_col, opacity=0.3)
skull.plot(color=skull_col, opacity=0.3)
brain.plot(color=brain_col, opacity=0.3)
cortex.plot(opacity=1, scalars=cortex.normals[:,0])

# View EEG electrodes
electrodes = np.loadtxt('eeg_channels_locations.txt')
mlab.points3d(electrodes[:,0], electrodes[:,1], electrodes[:,2], opacity=0.5,
                scale_factor=6)

# View MEG squids
squids = np.loadtxt('meg_channels_locations.squids')
mlab.quiver3d(squids[:,0], squids[:,1], squids[:,2],
              -squids[:,3], -squids[:,4], -squids[:,5],
              opacity=0.5, scale_factor=10, mode='cone')

