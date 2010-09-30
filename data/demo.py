#!/usr/bin/env python

import openmeeg as om
from os import path
import sys
from optparse import OptionParser

data_path = path.dirname(path.abspath(__file__))
parser = OptionParser()
parser.add_option("-p", "--path", dest="data_path",
                  help="path to data folder", metavar="FILE", default=data_path)
(options, args) = parser.parse_args()
data_path = options.data_path

# =============
# = Load data =
# =============

subject = 'Head1'
cond_file = path.join(data_path, 'Models', subject, subject + '.cond')
geom_file = path.join(data_path, 'Models', subject, subject + '.geom')
source_mesh_file = path.join(data_path, 'Models', subject, subject + '.tri')
dipole_file = path.join(data_path, 'Models', subject, subject + '.dip')
squidsFile = path.join(data_path, 'Computations', subject, subject + '.squids')
patches_file = path.join(data_path, 'Computations', subject, subject + '.patches')

geom = om.Geometry()
geom.read(geom_file, cond_file)

mesh = om.Mesh()
mesh.load(source_mesh_file)

dipoles = om.Matrix()
dipoles.load(dipole_file)

sensors = om.Sensors()
sensors.load(squidsFile)

patches = om.Sensors()
patches.load(patches_file)

# =================================================
# = Compute forward problem (Build Gain Matrices) =
# =================================================

gauss_order = 3
use_adaptive_integration = True

hm                  = om.HeadMat(geom, gauss_order)
hminv               = hm.inverse()
ssm                 = om.SurfSourceMat(geom, mesh, gauss_order)
ss2mm               = om.SurfSource2MEGMat(mesh, sensors)
dsm                 = om.DipSourceMat(geom, dipoles, gauss_order, use_adaptive_integration)
ds2mm               = om.DipSource2MEGMat(dipoles, sensors)
h2mm                = om.Head2MEGMat(geom, sensors)
h2em                = om.Head2EEGMat(geom, patches)
gain_meg_surf       = om.GainMEG(hminv, ssm, h2mm, ss2mm)
gain_eeg_surf       = om.GainEEG(hminv, ssm, h2em)
gain_meg_dip        = om.GainMEG(hminv, dsm, h2mm, ds2mm)
gain_eeg_dip        = om.GainEEG(hminv, dsm, h2em)
gainadjoint_eeg_dip = om.GainEEGadjoint(hminv, dsm, h2em)

print "hm                  : %d x %d" % (hm.nlin(), hm.ncol())
print "hminv               : %d x %d" % (hminv.nlin(), hminv.ncol())
print "ssm                 : %d x %d" % (ssm.nlin(), ssm.ncol())
print "ss2mm               : %d x %d" % (ss2mm.nlin(), ss2mm.ncol())
print "dsm                 : %d x %d" % (ssm.nlin(), ssm.ncol())
print "ds2mm               : %d x %d" % (ss2mm.nlin(), ss2mm.ncol())
print "h2mm                : %d x %d" % (h2mm.nlin(), h2mm.ncol())
print "h2em                : %d x %d" % (h2mm.nlin(), h2mm.ncol())
print "gain_meg_surf       : %d x %d" % (gain_meg_surf.nlin(), gain_meg_surf.ncol())
print "gain_eeg_surf       : %d x %d" % (gain_eeg_surf.nlin(), gain_eeg_surf.ncol())
print "gain_meg_dip        : %d x %d" % (gain_meg_dip.nlin(), gain_meg_dip.ncol())
print "gain_eeg_dip        : %d x %d" % (gain_eeg_dip.nlin(), gain_eeg_dip.ncol())
print "gainadjoint_eeg_dip : %d x %d" % (gain_eeg_dip.nlin(), gain_eeg_dip.ncol())

# Leadfield MEG in one line :

gain_meg_surf_one_line = om.GainMEG(om.HeadMat(geom, gauss_order).inverse(),
                                    om.SurfSourceMat(geom,mesh, gauss_order),
                                    om.Head2MEGMat(geom, sensors),
                                    om.SurfSource2MEGMat(mesh, sensors))

print "gain_meg_surf_one_line : %d x %d" % (gain_meg_surf_one_line.nlin(),
                                            gain_meg_surf_one_line.ncol())

# ========================
# = Compute forward data =
# ========================

srcFile = path.join(data_path,'Computations/'+subject+'/'+subject+'.src')
sources = om.Matrix()
sources.load(srcFile)

noiseLevel = 0.0
est_meg = om.Forward(gain_meg_dip, sources,noiseLevel)
print "est_meg    : %d x %d"%(est_meg.nlin(), est_meg.ncol())

est_eeg = om.Forward(gain_eeg_dip, sources,noiseLevel)
print "est_eeg    : %d x %d"%(est_eeg.nlin(), est_eeg.ncol())

est_eegadjoint = om.Forward(gainadjoint_eeg_dip, sources,noiseLevel)
print "est_eegadjoint    : %d x %d" % (est_eegadjoint.nlin(), est_eegadjoint.ncol())

# ============================
# = Compute inverse problems =
# ============================

smooth_weight = 0.0001
maxIter = 300
tol = 0

smooth_matrix = mesh.gradient()
aiVector = mesh.areas()

meg_inverse_mn   = om.MN_inverse(est_meg, gain_meg_dip, smooth_weight)
meg_inverse_heat = om.HEAT_inverse(est_meg, gain_meg_dip, smooth_matrix, smooth_weight)
meg_inverse_tv   = om.TV_inverse(est_meg, gain_meg_dip, smooth_matrix,aiVector, smooth_weight, maxIter, tol)

# ==================================
# = Example of basic manipulations =
# ==================================

v1 = om.Vect3(1,0,0)
v2 = om.Vect3(0,1,0)
v3 = om.Vect3(0,0,1)

print v1.norm()
print (v1+v2).norm()

normale = om.Vect3(1,0,0)
t = om.Triangle(0,1,2,normale)

hmFile = subject+'.hm'
hm.save(hmFile)

ssmFile = subject+'.ssm'
ssm.save(ssmFile)

m1 = om.SymMatrix()
m1.load(hmFile)
print m1(0,0)
print m1.nlin()
print m1.ncol()

m2 = om.Matrix()
m2.load(ssmFile)
print m2(0,0)
print m2.nlin()
print m2.ncol()

sys.exit(0)

# ===================
# = Numpy interface =
# ===================

# For a Vector

vec = om.asarray(aiVector)
print vec.size
print vec[0:5]
print vec

# For a Matrix

mat = om.asarray(m2)

print mat.shape
print mat.sum()
mat[0:2,1:3] = 0
print mat[0:5,0:5]
