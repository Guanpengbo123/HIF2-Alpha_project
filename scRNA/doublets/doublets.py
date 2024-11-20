
import numpy as np
import doubletdetection
import tarfile
import os
import re

!cd /data1/02.private/dengyj/analysis/Hif2a/matrix/doublets

dir = '/data1/02.private/dengyj/analysis/Hif2a/matrix/'
sub_dir = '/filtered_feature_bc_matrix/matrix.mtx.gz'

sample_dir = os.listdir(dir)
for i in range(len(sample_dir)):
  sample_dir[i] = dir + sample_dir[i] + sub_dir


for d in sample_dir:
    raw_counts = doubletdetection.load_mtx(d)
    zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
    raw_counts = raw_counts[:, ~zero_genes]
    clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
    doublets = clf.fit(raw_counts).predict(p_thresh=1e-16, voter_thresh=0.5)
    iter_obj = d.split('/')
    file_name = iter_obj[len(iter_obj) - 3] + '_doublets.txt'
    np.savetxt(file_name, doublets, fmt="%d", delimiter=",")


#########################2020.12.1 second sequencing (O1-10X5 I3-10X5 I4-10X5)

sample_dir = ['O1-10X5','I3-10X5','I4-10X5']
for i in range(len(sample_dir)):
  sample_dir[i] = dir + sample_dir[i] + sub_dir


for d in sample_dir:
    raw_counts = doubletdetection.load_mtx(d)
    zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
    raw_counts = raw_counts[:, ~zero_genes]
    clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
    doublets = clf.fit(raw_counts).predict(p_thresh=1e-16, voter_thresh=0.5)
    iter_obj = d.split('/')
    file_name = iter_obj[len(iter_obj) - 3] + '_doublets.txt'
    np.savetxt(file_name, doublets, fmt="%d", delimiter=",")
