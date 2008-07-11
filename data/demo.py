#!/usr/bin/env python

import openmeeg as om

# =============
# = Load data =
# =============

subject='Head1'
condFile='Models/'+subject+'/'+subject+'.cond'
geomFile='Models/'+subject+'/'+subject+'.geom'
sourceMeshFile='Models/'+subject+'/'+subject+'.tri'
squidsFile='Computations/'+subject+'/'+subject+'.squids'
patchesFile='Computations/'+subject+'/'+subject+'.patches'

print condFile

geom = om.Geometry()
geom.read(geomFile,condFile)

mesh = om.Mesh()
mesh.load(sourceMeshFile)

sensors = om.Sensors()
sensors.load(squidsFile)

patches = om.Matrix()
patches.load(patchesFile)

# =================================================
# = Compute forward problem (Build Gain Matrices) =
# =================================================

gaussOrder = 3;

hm      = om.HeadMat(geom,gaussOrder);
hminv   = hm.inverse()
ssm      = om.SurfSourceMat(geom,mesh,gaussOrder);
ss2mm    = om.SurfSource2MEGMat(mesh,sensors)
h2mm    = om.Head2MEGMat(geom,sensors)
h2em    = om.Head2EEGMat(geom,patches)
gain_meg = om.GainMEG(hminv,ssm,h2mm,ss2mm)
gain_eeg = om.GainEEG(hminv,ssm,h2em)

print "hm         : %d x %d"%(hm.nlin(),hm.ncol())
print "hminv      : %d x %d"%(hminv.nlin(),hminv.ncol())
print "ssm        : %d x %d"%(ssm.nlin(),ssm.ncol())
print "ss2mm      : %d x %d"%(ss2mm.nlin(),ss2mm.ncol())
print "h2mm       : %d x %d"%(h2mm.nlin(),h2mm.ncol())
print "h2em       : %d x %d"%(h2mm.nlin(),h2mm.ncol())
print "gain_meg   : %d x %d"%(gain_meg.nlin(),gain_meg.ncol())
print "gain_eeg   : %d x %d"%(gain_eeg.nlin(),gain_eeg.ncol())

# Leadfield MEG in one line :

surf_gain_meg = om.GainMEG(om.HeadMat(geom,gaussOrder).inverse(),om.SurfSourceMat(geom,mesh,gaussOrder), \
om.Head2MEGMat(geom,sensors),om.SurfSource2MEGMat(mesh,sensors));

print "gain_meg (one line) : %d x %d"%(gain_meg.nlin(),gain_meg.ncol())

# ========================
# = Compute forward data =
# ========================

srcFile = 'Computations/'+subject+'/'+subject+'.src'
sources = om.Matrix()
sources.load(srcFile)

noiseLevel = 0.0
est_meg = om.Forward(gain_meg,sources,noiseLevel)
print "est_meg    : %d x %d"%(est_meg.nlin(),est_meg.ncol())

est_eeg = om.Forward(gain_eeg,sources,noiseLevel)
print "est_eeg    : %d x %d"%(est_eeg.nlin(),est_eeg.ncol())

# ============================
# = Compute inverse problems =
# ============================

smoothWeight = 0.0001
maxIter = 300
tol = 0

smoothMatrix = mesh.gradient()
aiVector = mesh.areas()

meg_inverse_mn   = om.MN_inverse(est_meg,gain_meg,smoothWeight)
meg_inverse_heat = om.HEAT_inverse(est_meg,gain_meg,smoothMatrix,smoothWeight)
meg_inverse_tv   = om.TV_inverse(est_meg,gain_meg,smoothMatrix,aiVector,smoothWeight,maxIter,tol)

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
hm.saveBin(hmFile)

ssmFile = subject+'.ssm'
ssm.saveBin(ssmFile)

m1 = om.SymMatrix()
m1.loadBin(hmFile)
print m1(0,0)
print m1.nlin()
print m1.ncol()

m2 = om.Matrix()
m2.loadBin(ssmFile)
print m2(0,0)
print m2.nlin()
print m2.ncol()
