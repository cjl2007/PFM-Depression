**Precision functional mapping of indidivudals with depression**
------ 

This repository contains code related to "Expansion of a frontostriatal salience network in individuals with depression” (in revision),
including routines for mapping functional brain networks in individual subjects. 

Usage
------ 

The "pfm_example_use.m" script implements the entire precision funcional mapping (PFM) procedure, including concatenation of the individual scans ("concatenate_scans.m"), generating distance matrix ("make_distance_matrix.m"), applying the infomap community detection algorithm ("pfm_wrapper.m"), and algorithmic assignment of network identities to infomap communities ("pfm_identify_networks.m"). This code assumes that the resting-state fMRI dataset is fully preprocessed, denoised, and surface-registered (see https://github.com/cjl2007/Liston-Laboratory-MultiEchofMRI-Pipeline). 

If you have issues, please email me (Chuck Lynch; cjl2007@med.cornell.edu), and I will try my best to help you.
