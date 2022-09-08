import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from mne.fixes import bincount
import openmeeg as om


def _col_corrs(a, b):
    """Compute correlation between paired columns, being careful about 0."""
    a = a - a.mean(0)
    b = b - b.mean(0)
    num = (a * b).mean(0)
    a_std = np.sqrt((a * a).mean(0))
    b_std = np.sqrt((b * b).mean(0))
    all_zero = (a_std == 0) & (b_std == 0)
    num[all_zero] = 1.
    a_std[all_zero] = 1.
    b_std[all_zero] = 1.
    return num / (a_std * b_std)


def _rdm(a, b):
    """Compute the ratio of norms, being careful about 0."""
    a_norm = np.linalg.norm(a, axis=0)
    b_norm = np.linalg.norm(b, axis=0)
    all_zero = (a_norm == 0) & (b_norm == 0)
    a_norm[all_zero] = 1.
    b_norm[all_zero] = 1.
    return a_norm / b_norm


###############################################################################
# Load data
def get_meg_leadfield_om(model, head_model):
    geom_file = f'data/{model}/{head_model}.geom'
    cond_file = f'data/{model}/{head_model}.cond'
    dipoles_file = f'data/{model}/{model}.dip'
    squids_file = f'data/{model}/{model}.squids'

    geom = om.Geometry(geom_file, cond_file)
    dipoles = om.Matrix(dipoles_file)
    meg_sensors = om.Sensors(squids_file)

    hm = om.HeadMat(geom)
    hm.invert()
    hminv = hm

    dsm = om.DipSourceMat(geom, dipoles, "Brain")

    # For MEG
    ds2mm = om.DipSource2MEGMat(dipoles, meg_sensors)
    h2mm = om.Head2MEGMat(geom, meg_sensors)

    meg_leadfield = om.GainMEG(hminv, dsm, h2mm, ds2mm)
    return meg_leadfield.array()


def _apply_weights(lf, ws, bins):
    B = np.array([bincount(bins, ws * x, bins[-1] + 1)
                  for x in lf.T], float)
    return B


def get_meg_leadfield_mne():
    import mne
    from mne.datasets import sample
    data_path = sample.data_path()

    # the raw file containing the channel location + types
    sample_dir = data_path / 'MEG' / 'sample'
    raw_fname = sample_dir / 'sample_audvis_raw.fif'
    # The paths to Freesurfer reconstructions
    subjects_dir = data_path / 'subjects'
    subject = 'sample'
    trans = sample_dir / 'sample_audvis_raw-trans.fif'
    src = mne.setup_source_space(subject, spacing='ico3', add_dist='patch',
                                 subjects_dir=subjects_dir)
    conductivity = (0.3,)  # for single layer
    model = mne.make_bem_model(subject='sample', ico=3,
                               conductivity=conductivity,
                               subjects_dir=subjects_dir)
    bem = mne.make_bem_solution(model)
    fwd = mne.make_forward_solution(raw_fname, trans=trans, src=src, bem=bem,
                                    meg=True, eeg=False, mindist=5.0,
                                    verbose=True)
    return fwd['sol']['data'].T


if __name__ == '__main__':
    leadfield_mne = get_meg_leadfield_mne()
    leadfield_om = get_meg_leadfield_om("mne_sample_ico3", "mne_sample_ico3_1layer")
    leadfield_ft = loadmat("data/mne_sample_ico3/meg_forward_nolte.mat")['lf_singleshell']
    leadfield_ft[np.isnan(leadfield_ft)] = 0

    weights = pd.read_csv("data/mne_sample_ico3/weights.csv")
    ws, bins = weights['ws'].values, weights['bins'].values

    leadfield_om = _apply_weights(leadfield_om, ws, bins)
    leadfield_ft = _apply_weights(leadfield_ft, ws, bins)

    corrs_om = _col_corrs(leadfield_om, leadfield_mne)
    corrs_om = corrs_om[~np.isnan(corrs_om)]

    corrs_ft = _col_corrs(leadfield_ft, leadfield_mne)
    corrs_ft = corrs_ft[~np.isnan(corrs_om)]

    plt.hist(corrs_om, bins=100, label="openmeeg vs mne")
    plt.hist(corrs_ft, bins=100, label="fieldtrip vs mne")
    plt.legend()
    plt.show()

    print(np.min(corrs_om), np.max(corrs_om), np.mean(corrs_om), np.median(corrs_om))
    print(np.min(corrs_ft), np.max(corrs_ft), np.mean(corrs_ft), np.median(corrs_ft))
