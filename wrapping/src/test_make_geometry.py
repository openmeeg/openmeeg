#!/usr/bin/env python3
#
# WIP: copy this file in wrapping/src directory of the build tree
# then execute it (python3 ./wrapping/src/poc_geometry)

import openmeeg as om
from os import path as path
from optparse import OptionParser

data_path = path.dirname(path.abspath(__file__))
parser = OptionParser()
parser.add_option("-p","--path",dest="data_path",help="path to data folder",metavar="FILE",default=data_path)

options,args = parser.parse_args()
data_path    = options.data_path

# Load mesh data to mimic Head1.geom + Head1.cond

subject = "Head1"
dirpath = path.join(data_path,subject)

cortex = om.Mesh(path.join(dirpath,"cortex.1.tri"))
cortex.setName("cortex")
skull = om.Mesh(path.join(dirpath,"skull.1.tri"))
skull.setName("skull")
scalp = om.Mesh(path.join(dirpath,"scalp.1.tri"))
scalp.setName("scalp")

# It should be possible to have multiple oriented meshes per interface. e.g.
# interface1 = [(m1,om.OrientedMesh.Normal), (m2,om.OrientedMesh.Opposite), (m3,om.OrientedMesh.Normal)]
# It should also be possible to have a name added at the beginning of the tuple.

interfaces = [[(cortex,om.OrientedMesh.Normal)],
              [(skull,om.OrientedMesh.Normal)],
              [(scalp,om.OrientedMesh.Normal)]]

domains = {
    "Scalp" : ([(interfaces[1],om.SimpleDomain.Outside), (interfaces[2],om.SimpleDomain.Inside)],1.0) ,
    "Brain" : ([(interfaces[0],om.SimpleDomain.Inside)],1.0),
    "Air"   : ([(interfaces[2],om.SimpleDomain.Outside)],0.0),
    "Skull" : ([(interfaces[1],om.SimpleDomain.Inside), (interfaces[0],om.SimpleDomain.Outside)],0.0125)
}

g1 = om.make_geometry(domains);
g2 = om.Geometry(path.join(dirpath,subject+".geom"),path.join(dirpath,subject+".cond"))

assert g1.__class__ == g2.__class__
