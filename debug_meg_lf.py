import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
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
def get_meg_leadfield(model, head_model):
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


if __name__ == '__main__':
    # leadfield = get_meg_leadfield("mne_sample_ico3", "mne_sample_ico3")
    leadfield = get_meg_leadfield("mne_sample_ico3", "mne_sample_ico3_1layer")

    leadfield_ft = loadmat("data/mne_sample_ico3/meg_forward_nolte.mat")['lf_singleshell']

    leadfield_ft[np.isnan(leadfield_ft)] = 0
    
    corrs = _col_corrs(leadfield, leadfield_ft)

    corrs = corrs[~np.isnan(corrs)]

    plt.hist(corrs, bins=100)

    print(np.min(corrs), np.max(corrs), np.mean(corrs), np.median(corrs))

    # leadfield_3l = get_meg_leadfield("Head1", "Head1")
    # leadfield_1l = get_meg_leadfield("Head1", "Head1_1_layer")

    # corrs = [np.corrcoef(x, y)[0, 1] for x, y in zip(leadfield_3l.T, leadfield_1l.T)]
    # print(np.array(corrs))
