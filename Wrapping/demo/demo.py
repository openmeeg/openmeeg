#!/usr/bin/env python

import openmeeg

v1 = openmeeg.Vect3(1,0,0)
v2 = openmeeg.Vect3(0,1,0)
v3 = openmeeg.Vect3(0,0,1)

print v1.norme()
print (v1+v2).norme()

normale = openmeeg.Vect3(1,0,0)
t = openmeeg.Triangle(0,1,2,normale)

m1 = openmeeg.symmatrice()
m1.loadBin('Head1.lhs')
print m1(0,0)
print m1.nlin()
print m1.ncol()

m2 = openmeeg.matrice()
m2.loadBin('Head1.Hmeg')
print m2(0,0)
print m2.nlin()
print m2.ncol()

geom = openmeeg.Geometry()
geom.read('Head1.geom','Head1.cond')
mesh = openmeeg.Mesh()
mesh.load('Head1.tri')
sensors = openmeeg.Sensors()
sensors.load('Head1.squids')

patches = openmeeg.matrice()
patches.loadTxt('Head1.patches')

lhs     = openmeeg.LHS_matrice(geom,3);
lhsinv  = lhs.inverse()
rhs     = openmeeg.RHS_matrice(geom,mesh,3);
s2meg   = openmeeg.sToMEG_matrice(mesh,sensors)
v2meg   = openmeeg.vToMEG_matrice(geom,sensors)
v2eeg   = openmeeg.vToEEG_matrice(geom,patches)

print "lhs        : %d x %d"%(lhs.nlin(),lhs.ncol())
print "lhsinv     : %d x %d"%(lhsinv.nlin(),lhsinv.ncol())
print "rhs        : %d x %d"%(rhs.nlin(),rhs.ncol())
print "s2meg      : %d x %d"%(s2meg.nlin(),s2meg.ncol())
print "v2meg      : %d x %d"%(v2meg.nlin(),v2meg.ncol())
print "v2eeg      : %d x %d"%(v2meg.nlin(),v2meg.ncol())

gain_meg = openmeeg.HMEG_matrice(lhsinv,rhs,v2meg,s2meg)
print "gain_meg   : %d x %d"%(gain_meg.nlin(),gain_meg.ncol())

gain_eeg = openmeeg.HEEG_matrice(lhsinv,rhs,v2eeg)
print "gain_eeg   : %d x %d"%(gain_eeg.nlin(),gain_eeg.ncol())

sources = openmeeg.matrice('Head1.src')

noiseLevel = 0.0
est_meg = openmeeg.Forward_matrice(gain_meg,sources,noiseLevel)
print "est_meg    : %d x %d"%(est_meg.nlin(),est_meg.ncol())

est_eeg = openmeeg.Forward_matrice(gain_eeg,sources,noiseLevel)
print "est_eeg    : %d x %d"%(est_eeg.nlin(),est_eeg.ncol())

smoothWeight = 0.0001
maxIter = 300
tol = 0

smoothMatrix = mesh.gradient()
aiVector = mesh.areas()

meg_inverse_mn   = openmeeg.MN_inverse_matrice(est_meg,gain_meg,smoothWeight)
meg_inverse_heat = openmeeg.HEAT_inverse_matrice(est_meg,gain_meg,smoothMatrix,smoothWeight)
meg_inverse_tv   = openmeeg.TV_inverse_matrice(est_meg,gain_meg,smoothMatrix,aiVector,smoothWeight,maxIter,tol)
