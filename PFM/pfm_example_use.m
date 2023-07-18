%% An example precision functional mapping routine

% define some paths 
Paths{1} = '~/res0urces/'; % this folder contains ft_read / gifti functions for reading and writing cifti files (e.g., https://github.com/MidnightScanClub/MSCcodebase). 

% define subject 
% directory and surface files;
Subdir = '/path/to/subject/';
str = strsplit(Subdir,'/'); Subject = str{end};
MidthickSurfs{1} = [Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii'];
MidthickSurfs{2} = [Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii'];

% denoising QC;
grayplot_qa_func(Subdir);

% concatenate and smooth resting-state fMRI datasets;
nSessions = length(dir([Subdir '/func/rest/session_*']));
[C,ScanIdx,FD] = concatenate_scans(Subdir,'Rest_OCME+MEICA+MGTR',1:nSessions);
mkdir([Subdir '/func/rest/ConcatenatedCiftis']);
cd([Subdir '/func/rest/ConcatenatedCiftis']);

% make distance matrix and then regress
% adjacent cortical signals from subcortical voxels;
make_distance_matrix(C,MidthickSurfs,[Subdir '/anat/T1w/fsaverage_LR32k/'],8);
[C] = regress_cortical_signals(C,[Subdir '/anat/T1w/fsaverage_LR32k/DistanceMatrix.mat'],20);
ft_write_cifti_mod([Subdir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_Concatenated+SubcortRegression.dtseries.nii'],C);
save([Subdir '/func/rest/ConcatenatedCiftis/ScanIdx'],'ScanIdx');
save([Subdir '/func/rest/ConcatenatedCiftis/FD'],'FD');
clear C % clear intermediate file
    
% sweep a range of
% smoothing kernels;
for k = [0.85 1.7 2.55]
    smooth_cifti(Subdir,'Rest_OCME+MEICA+MGTR_Concatenated+SubcortRegression.dtseries.nii',['Rest_OCME+MEICA+MGTR_Concatenated+SubcortRegression+SpatialSmoothing' num2str(k) '.dtseries.nii'],k,k);
end

% load your concatenated resting-state dataset, here we selected the highest level of spatial smoothing
C = ft_read_cifti_mod([Subdir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_Concatenated+SubcortRegression+SpatialSmoothing2.55.dtseries.nii']);
C.data(:,FD < 0.3)=[]; % remove high motion volumes

% run precision mapping; (Note: there are hard set paths related to the resources folder and infomap in "pfm.m" user must first set)
pfm(C,[Subdir '/anat/T1w/fsaverage_LR32k/DistanceMatrix.mat'],[Subdir '/pfm/'],flip([0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05]),[1 5 10 50 50 50 50 50 50],10,[],{'CORTEX_LEFT','CEREBELLUM_LEFT','ACCUMBENS_LEFT','CAUDATE_LEFT','PALLIDUM_LEFT','PUTAMEN_LEFT','THALAMUS_LEFT','HIPPOCAMPUS_LEFT','AMYGDALA_LEFT','ACCUMBENS_LEFT','CORTEX_RIGHT','CEREBELLUM_RIGHT','ACCUMBENS_RIGHT','CAUDATE_RIGHT','PALLIDUM_RIGHT','PUTAMEN_RIGHT','THALAMUS_RIGHT','HIPPOCAMPUS_RIGHT','AMYGDALA_RIGHT','ACCUMBENS_RIGHT'},5);
spatial_filtering([Subdir '/pfm/Bipartite_PhysicalCommunities.dtseries.nii'],[Subdir '/pfm/'],['Bipartite_PhysicalCommunities+SpatialFiltering.dtseries.nii'],MidthickSurfs);

% remove some intermediate files;
system(['rm ' Subdir '/pfm/*.net']);
system(['rm ' Subdir '/pfm/*.clu']);
system(['rm ' Subdir '/pfm/*Log*']);

% multiplex 

% load the resting-state data;
C = ft_read_cifti_mod([Subdir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_Concatenated+SubcortRegression+SpatialSmoothing2.55.dtseries.nii']);
load([Subdir '/func/rest/ConcatenatedCiftis/ScanIdx.mat']);
load([Subdir '/func/rest/ConcatenatedCiftis/FD.mat']);
C.data = C.data(:,FD < 0.3); % apply motion censoring;
LayerIdx = ScanIdx(:,1); % study visit index
LayerIdx(FD > 0.3) = [];

% run multiplex precision mapping routine;
multiplex_pfm_wrapper(C,LayerIdx,[Subdir '/anat/T1w/fsaverage_LR32k/DistanceMatrix.mat'],[Subdir '/pfm/'],0.001,10,0.15,1,10,{'CORTEX_LEFT','CEREBELLUM_LEFT','ACCUMBENS_LEFT','CAUDATE_LEFT','PALLIDUM_LEFT','PUTAMEN_LEFT','THALAMUS_LEFT','HIPPOCAMPUS_LEFT','AMYGDALA_LEFT','ACCUMBENS_LEFT','CORTEX_RIGHT','CEREBELLUM_RIGHT','ACCUMBENS_RIGHT','CAUDATE_RIGHT','PALLIDUM_RIGHT','PUTAMEN_RIGHT','THALAMUS_RIGHT','HIPPOCAMPUS_RIGHT','AMYGDALA_RIGHT','ACCUMBENS_RIGHT'},[]);
spatial_filtering([Subdir '/pfm/MultiPlex_StateCommunities_Density0.001_RelaxRate0.15.dtseries.nii'],[Subdir '/pfm/'],['MultiPlex_StateCommunities_Density0.001_RelaxRate0.15+SpatialFiltering.dtseries.nii'],MidthickSurfs,50,50);
