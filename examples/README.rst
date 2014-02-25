Demo scripts to compute leadfields with OpenMEEG

- Supports EEG, MEG, EIT and Internal potential leadfields

This folder contains a sample realistic dataset for EEG, MEG, EIT
and internal potential forward modeling.

The head model is a 3 layers model with 3 nested meshes:
brain.vtk, skull.vtk and head.vtk

To run the computation you can use the scripts:

On windows (bat script):
------------------------
   compute_leadfields.bat

On Linux or Mac OS X (bash script):
-----------------------------------
	./compute_leadfields.sh

Or using Python:
------------------------------------
	python compute_leadfields.py


The leadfields computed are stored in (matlab format):

    eeg_leadfield.mat (for EEG)

    meg_leadfield.mat (for MEG)

    eit_leadfield.mat (for EIT)

    ip_leadfield.mat (for Internal Potential)

The files used during the BEM computation are stored in the "tmp" folder.

See sample_output.txt to see what the scripts output should look like.

On a recent workstation the computation takes about 
   	    744.66 s    344.12 s  390.41 s om_assemble -HM
	 +  431.74 s	188.5 s   246.5 s  om_inverser
	 +  2689.71	    545.1 s	  719.13 s om_assemble -DSM
	 +  0.13 s	     0.02 s   0.04 s   om_assemble -H2EM
	 + 80.92 s	     3.19 s	  16.42 s  om_gain -EEG
	 + 3.39 s.	     0.88 s   1.23 s   om_assemble -H2MM
	 + 2.33 s	     0.21 s	  0.95 s   om_assemble -DS2MM
	 + 84.84 s	     4.17 s	  17.65 s  om_gain -MEG
                             173.49 s  om_assemble -EITSM
                               1.15 s  om_gain -EEG

If you meet some difficulties running this example please contact:

openmeeg-info@lists.gforge.inria.fr

The OpenMEEG developers.

