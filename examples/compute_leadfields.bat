set GEOMETRY=model\head_model.geom
set CONDUCTIVITIES=model\head_model.cond
set DIPOLES=model\cortex_dipoles_small.txt
set EEG_ELECTRODES=model\eeg_channels_locations.txt
set EIT_ELECTRODES=model\eit_locations.txt
set ECOG_ELECTRODES=model\ecog_electrodes_locations.txt
set SQUIDS=model\meg_channels_locations.squids
set INTERNAL_ELECTRODES=model\internal_electrodes_locations.txt

rem Leadfields
set EEG_LEADFIELD=leadfields\eeg_leadfield.mat
set EEG_LEADFIELD_ADJOINT=leadfields\eeg_leadfield_adjoint.mat
set ECOG_LEADFIELD=leadfields\ecog_leadfield.mat
set MEG_LEADFIELD=leadfields\meg_leadfield.mat
set MEG_LEADFIELD_ADJOINT=leadfields\meg_leadfield_adjoint.mat
set EIT_LEADFIELD=leadfields\eit_leadfield.mat
set IP_LEADFIELD=leadfields\ip_leadfield.mat

rem Name temporary matrices
rem For EEG and MEG

set HM=tmp\tmp.hm
set HMINV=tmp\tmp.hm_inv
set DSM=tmp\tmp.dsm

rem for EEG
set H2EM=tmp\tmp.h2em

rem for ECoG
set H2ECOGM=tmp\tmp.h2ecogm

rem for MEG
set H2MM=tmp\tmp.h2mm
set DS2MEG=tmp\tmp.ds2mm

rem for EIT
set EITSM=tmp\tmp.eitsm

rem for IP (Internal Potential)
set IPHM=tmp\tmp.iphm
set IPSM=tmp\tmp.ipsm

mkdir tmp
mkdir leadfields

rem Compute EEG gain matrix
om_assemble -HM %GEOMETRY% %CONDUCTIVITIES% %HM%
om_minverser %HM% %HMINV%
om_assemble -DSM %GEOMETRY% %CONDUCTIVITIES% %DIPOLES% %DSM%
om_assemble -H2EM %GEOMETRY% %CONDUCTIVITIES% %EEG_ELECTRODES% %H2EM%
om_gain -EEG %HMINV% %DSM% %H2EM% %EEG_LEADFIELD%
rem with adjoint
om_gain -EEGadjoint %GEOMETRY% %CONDUCTIVITIES% %DIPOLES% %HM% %H2EM% %EEG_LEADFIELD_ADJOINT%

rem Compute ECoG gain matrix
om_assemble -H2ECOGM %GEOMETRY% %CONDUCTIVITIES% %ECOG_ELECTRODES% "Cortex" %H2ECOGM%
om_gain -EEG %HMINV% %DSM% %H2ECOGM% %ECOG_LEADFIELD%

rem Compute MEG gain matrix
om_assemble -H2MM %GEOMETRY% %CONDUCTIVITIES% %SQUIDS% %H2MM%
om_assemble -DS2MM %DIPOLES% %SQUIDS% %DS2MEG%
om_gain -MEG %HMINV% %DSM% %H2MM% %DS2MEG% %MEG_LEADFIELD%
rem with adjoint
om_gain -MEGadjoint %GEOMETRY% %CONDUCTIVITIES% %DIPOLES% %HM% %H2MM% %DS2MEG% %MEG_LEADFIELD_ADJOINT%

rem Compute EIT gain matrix
om_assemble -EITSM %GEOMETRY% %CONDUCTIVITIES% %EIT_ELECTRODES% %EITSM%
om_gain -EEG %HMINV% %EITSM% %H2EM% %EIT_LEADFIELD%

rem Compute Internal Potential ...
om_assemble -H2IPM %GEOMETRY% %CONDUCTIVITIES% %INTERNAL_ELECTRODES% %IPHM%
rem ...for internal dipoles
om_assemble -DS2IPM %GEOMETRY% %CONDUCTIVITIES% %DIPOLES% %INTERNAL_ELECTRODES% %IPSM%
om_gain -IP %HMINV% %DSM% %IPHM% %IPSM% %IP_LEADFIELD%
rem ...for boundary-injected current
om_gain -SIP %HMINV%  %EITSM% %IPHM% %SIP_LEADFIELD%

