function smooth_cifti(Subdir,Cifti,SmoothedCifti,CorticalKernel,SubcorticalKernel)

% infer subject name;
str = strsplit(Subdir,'/');
Subject = str{end};

MIDTHICK_32k{1} = [Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii'];
MIDTHICK_32k{2} = [Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii'];

% smooth with geodesic (for surface data) and Euclidean (for volumetric data) Gaussian kernels;
system(['wb_command -cifti-smoothing ' Cifti ' ' num2str(CorticalKernel) ' ' num2str(SubcorticalKernel) ' COLUMN ' SmoothedCifti ' -left-surface ' MIDTHICK_32k{1} ' -right-surface ' MIDTHICK_32k{2} ' -merged-volume']);

end