#!/usr/bin/env python

import openmeeg as om

# =============
# = Load data =
# =============

subject='Head1'
condFile='HeadModels/'+subject+'/'+subject+'.cond'
geomFile='HeadModels/'+subject+'/'+subject+'.geom'
sourceMeshFile='HeadModels/'+subject+'/'+subject+'.tri'
squidsFile='Computations/'+subject+'/'+subject+'.squids'
patchesFile='Computations/'+subject+'/'+subject+'.patches'

print condFile

geom = om.Geometry()
geom.read(geomFile,condFile)

mesh = om.Mesh()
mesh.load(sourceMeshFile)

sensors = om.Sensors()
sensors.load(squidsFile)

patches = om.matrice()
patches.loadTxt(patchesFile)

# =======================
# = Build Gain Matrices =
# =======================

lhs      = om.LHS_matrice(geom,3);
lhsinv   = lhs.inverse()
rhs      = om.RHS_matrice(geom,mesh,3);
s2meg    = om.sToMEG_matrice(mesh,sensors)
v2meg    = om.vToMEG_matrice(geom,sensors)
v2eeg    = om.vToEEG_matrice(geom,patches)
gain_meg = om.HMEG_matrice(lhsinv,rhs,v2meg,s2meg)
gain_eeg = om.HEEG_matrice(lhsinv,rhs,v2eeg)

print "lhs        : %d x %d"%(lhs.nlin(),lhs.ncol())
print "lhsinv     : %d x %d"%(lhsinv.nlin(),lhsinv.ncol())
print "rhs        : %d x %d"%(rhs.nlin(),rhs.ncol())
print "s2meg      : %d x %d"%(s2meg.nlin(),s2meg.ncol())
print "v2meg      : %d x %d"%(v2meg.nlin(),v2meg.ncol())
print "v2eeg      : %d x %d"%(v2meg.nlin(),v2meg.ncol())
print "gain_meg   : %d x %d"%(gain_meg.nlin(),gain_meg.ncol())
print "gain_eeg   : %d x %d"%(gain_eeg.nlin(),gain_eeg.ncol())

# ===========================
# = Compute forward problem =
# ===========================

srcFile = 'Computations/'+subject+'/'+subject+'.src'
sources = om.matrice()
sources.loadTxt(srcFile)

noiseLevel = 0.0
est_meg = om.Forward_matrice(gain_meg,sources,noiseLevel)
print "est_meg    : %d x %d"%(est_meg.nlin(),est_meg.ncol())

est_eeg = om.Forward_matrice(gain_eeg,sources,noiseLevel)
print "est_eeg    : %d x %d"%(est_eeg.nlin(),est_eeg.ncol())

# ============================
# = Compute inverse problems =
# ============================

smoothWeight = 0.0001
maxIter = 300
tol = 0

smoothMatrix = mesh.gradient()
aiVector = mesh.areas()

meg_inverse_mn   = om.MN_inverse_matrice(est_meg,gain_meg,smoothWeight)
meg_inverse_heat = om.HEAT_inverse_matrice(est_meg,gain_meg,smoothMatrix,smoothWeight)
meg_inverse_tv   = om.TV_inverse_matrice(est_meg,gain_meg,smoothMatrix,aiVector,smoothWeight,maxIter,tol)

# ==================================
# = Example of basic manipulations =
# ==================================

v1 = om.Vect3(1,0,0)
v2 = om.Vect3(0,1,0)
v3 = om.Vect3(0,0,1)

print v1.norme()
print (v1+v2).norme()

normale = om.Vect3(1,0,0)
t = om.Triangle(0,1,2,normale)

lhsFile = subject+'.lhs'
lhs.saveBin(lhsFile)

rhsFile = subject+'.rhs'
rhs.saveBin(rhsFile)

m1 = om.symmatrice()
m1.loadBin(lhsFile)
print m1(0,0)
print m1.nlin()
print m1.ncol()

m2 = om.matrice()
m2.loadBin(rhsFile)
print m2(0,0)
print m2.nlin()
print m2.ncol()
