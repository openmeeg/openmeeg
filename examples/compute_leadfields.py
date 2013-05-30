#!/usr/bin/env python

"""Compute the 5 types of leadfields supported by OpenMEEG
    - EEG
    - MEG
    - EIT
    - Internal Potential for dipolar source
    - Internal Potential for boundary injected current
"""
print __doc__

import openmeeg as om

###############################################################################
# Load data
geom_file                    = 'head_model.geom'
cond_file                    = 'head_model.cond'
dipoles_file                 = 'cortex_dipoles.txt'
squids_file                  = 'meg_channels_locations.squids'
eeg_electrodes_file          = 'eeg_channels_locations.txt'
ecog_electrodes_file         = 'ecog_electrodes_locations.txt'
internal_electrodes_file     = 'internal_electrodes_locations.txt'

geom = om.Geometry()
geom.read(geom_file, cond_file)

dipoles = om.Matrix()
dipoles.load(dipoles_file)

meg_sensors = om.Sensors()
meg_sensors.load(squids_file)

eeg_electrodes = om.Sensors()
eeg_electrodes.load(eeg_electrodes_file)

ecog_electrodes = om.Sensors()
ecog_electrodes.load(ecog_electrodes_file)

int_electrodes = om.Matrix()
int_electrodes.load(internal_electrodes_file)

###############################################################################
# Compute Leadfields
gauss_order = 3
use_adaptive_integration = True
dipole_in_cortex = True

hm = om.HeadMat(geom, gauss_order)
hm.invert()
hminv = hm
# hminv = hm.inverse() # to also test the adjoint method: comment the 2 previous lines, and uncomment this line, and the two others containing 'adjoint'
dsm = om.DipSourceMat(geom, dipoles, gauss_order,
                                 use_adaptive_integration, dipole_in_cortex)

# For EEG
h2em = om.Head2EEGMat(geom, eeg_electrodes)

# For ECoG
h2ecogm = om.Head2ECoGMat(geom, ecog_electrodes)

# For MEG
ds2mm = om.DipSource2MEGMat(dipoles, meg_sensors)
h2mm = om.Head2MEGMat(geom, meg_sensors)

# For EIT (using the same electrodes as EEG)
eitsm = om.EITSourceMat(geom, eeg_electrodes, gauss_order)

# For Internal Potential
iphm = om.Surf2VolMat(geom, int_electrodes)
ipsm = om.DipSource2InternalPotMat(geom, dipoles, int_electrodes)

eeg_leadfield  = om.GainEEG(hminv, dsm, h2em)
ecog_leadfield = om.GainEEG(hminv, dsm, h2ecogm)
meg_leadfield  = om.GainMEG(hminv, dsm, h2mm, ds2mm)
eit_leadfield  = om.GainEEG(hminv, eitsm, h2em)
ip_leadfield   = om.GainInternalPot(hminv, dsm, iphm, ipsm)
sip_leadfield  = om.GainStimInternalPot(hminv, eitsm,  iphm)

# eeg_leadfield_adjoint = om.GainEEGadjoint(geom,dipoles,hm, h2em)

print "hm             : %d x %d" % (hm.nlin(), hm.ncol())
print "hminv          : %d x %d" % (hminv.nlin(), hminv.ncol())
print "dsm            : %d x %d" % (dsm.nlin(), dsm.ncol())
print "h2em           : %d x %d" % (h2em.nlin(), h2em.ncol())
print "h2ecogm        : %d x %d" % (h2ecogm.nlin(), h2ecogm.ncol())
print "ds2mm          : %d x %d" % (ds2mm.nlin(), ds2mm.ncol())
print "h2mm           : %d x %d" % (h2mm.nlin(), h2mm.ncol())
print "eeg_leadfield  : %d x %d" % (eeg_leadfield.nlin(), eeg_leadfield.ncol())
print "ecog_leadfield : %d x %d" % (ecog_leadfield.nlin(), ecog_leadfield.ncol())
print "meg_leadfield  : %d x %d" % (meg_leadfield.nlin(), meg_leadfield.ncol())
print "eit_leadfield  : %d x %d" % (eit_leadfield.nlin(), eit_leadfield.ncol())
print "ip_leadfield   : %d x %d" % (ip_leadfield.nlin(), ip_leadfield.ncol())
print "sip_leadfield  : %d x %d" % (sip_leadfield.nlin(), sip_leadfield.ncol())

eeg_leadfield.save('eeg_leadfield.mat')
# eeg_leadfield_adjoint.save('eeg_leadfield_adjoint.mat')
ecog_leadfield.save('ecog_leadfield.mat')
meg_leadfield.save('meg_leadfield.mat')
eit_leadfield.save('eit_leadfield.mat')
ip_leadfield.save('ip_leadfield.mat')
sip_leadfield.save('sip_leadfield.mat')
