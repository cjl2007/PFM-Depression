%% A tutorial covering precision functional mapping using an example dataset. 

%% Before you begin.

% add dependencies to Matlab search path
addpath(genpath([pwd '/PFM-Tutorial/Utilities']));

% define path to some software packages that will be needed
InfoMapBinary = '/home/charleslynch/miniconda3/bin/infomap'; % path to infomap binary; code tested on version 2.0.0 
WorkbenchBinary = '/usr/local/workbench/bin_linux64/wb_command'; % path to workbench binary; code tested on version 1.4.2

% number of 
% workers
nWorkers = 5;

%% Step 1: Temporal Concatenation of fMRI data from all sessions.

% define subject directory and name;
Subdir = [pwd '/WCM-ME/derivatives/sub-ME01/'];
Subject = 'ME01';

% define & create
% the pfm directory;
PfmDir = [Subdir '/pfm/'];
mkdir(PfmDir);

% count the number of imaging sessions;
nSessions = length(dir([Subdir '/processed_restingstate_timecourses/ses-func*']));

% preallocate;
ConcatenatedData = [];

% sweep through
% the sessions;
for i = 1:nSessions
    
    % count the number of runs in this session
    nRuns = length(dir([Subdir '/processed_restingstate_timecourses/ses-func' sprintf('%02d',i) '/*run-*.dtseries.nii']));
    
    % sweep 
    % through
    % the runs;
    for ii = 1:nRuns
        
        % load the denoised & fs_lr_32k surface-registered CIFTI file for run "ii" from session "i"...
        Cifti = ft_read_cifti_mod([Subdir '/processed_restingstate_timecourses/ses-func' sprintf('%02d',i) '/sub-' Subject '_ses-func' sprintf('%02d',i) '_task-rest_run-' sprintf('%02d',ii) '_bold_32k_fsLR.dtseries.nii']);
        Cifti.data = Cifti.data - mean(Cifti.data,2); % demean
        Tmask = load([Subdir '/processed_restingstate_timecourses/ses-func' sprintf('%02d',i) '/sub-' Subject '_ses-func' sprintf('%02d',i) '_task-rest_run-' sprintf('%02d',ii) '_bold_32k_fsLR_tmask.txt']);
        ConcatenatedData = [ConcatenatedData Cifti.data(:,Tmask==1)]; % 1 (Low motion timepoints) == FD < 0.3mm, 0 (High motion timepoints) == FD > 0.3mm

    end
    
end

% make a single CIFTI containing 
% time-series from all scans;
ConcatenatedCifti = Cifti;
ConcatenatedCifti.data = ConcatenatedData;

%% Step 2: Make a distance matrix.

% define fs_lr_32k midthickness surfaces;
MidthickSurfs{1} = [Subdir '/fs_LR/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii'];
MidthickSurfs{2} = [Subdir '/fs_LR/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii'];

% make the distance matrix;
pfm_make_dmat(ConcatenatedCifti,MidthickSurfs,PfmDir,nWorkers,WorkbenchBinary); %

% optional: regress adjacent cortical signal from subcortex to reduce artifactual coupling 
% (for example, between cerebellum and visual cortex, or between putamen and insular cortex)
[ConcatenatedCifti] = pfm_regress_adjacent_cortex(ConcatenatedCifti,[PfmDir '/DistanceMatrix.mat'],20);

% write out the CIFTI file;
ft_write_cifti_mod([Subdir '/pfm/sub-ME01_task-rest_concatenated_32k_fsLR.dtseries.nii'],ConcatenatedCifti);
 
%% Step 3: Apply spatial smoothing.

% define a range of gaussian 
% smoothing kernels (in sigma)
KernelSizes = [0.85 1.7 2.55];

% sweep a range of
% smoothing kernels;
for k = KernelSizes
    
    % smooth with geodesic (for surface data) and Euclidean (for volumetric data) Gaussian kernels;
    system([WorkbenchBinary ' -cifti-smoothing ' PfmDir '/sub-ME01_task-rest_concatenated_32k_fsLR.dtseries.nii '...
    num2str(k) ' ' num2str(k) ' COLUMN ' PfmDir '/sub-ME01_task-rest_concatenated_smoothed' num2str(k) '_32k_fsLR.dtseries.nii -left-surface ' MidthickSurfs{1} ' -right-surface ' MidthickSurfs{2} ' -merged-volume']);
    
end

%% Step 4: Run infomap.

% load your concatenated resting-state dataset, pick whatever level of spatial smoothing you want
ConcatenatedCifti = ft_read_cifti_mod([PfmDir '/sub-ME01_task-rest_concatenated_smoothed2.55_32k_fsLR.dtseries.nii']);

% define inputs;
DistanceMatrix = [Subdir '/pfm/DistanceMatrix.mat']; % can be path to file
DistanceCutoff = 10; % in mm; usually between 10 to 30 mm works well.
GraphDensities = flip([0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05]); % 
NumberReps = 50; % number of times infomap is run;
BadVertices = []; % optional, but you could include regions to ignore, if you know there is bad signal there.
Structures = {'CORTEX_LEFT','CEREBELLUM_LEFT','ACCUMBENS_LEFT','CAUDATE_LEFT','PALLIDUM_LEFT','PUTAMEN_LEFT','THALAMUS_LEFT','HIPPOCAMPUS_LEFT','AMYGDALA_LEFT','ACCUMBENS_LEFT','CORTEX_RIGHT','CEREBELLUM_RIGHT','ACCUMBENS_RIGHT','CAUDATE_RIGHT','PALLIDUM_RIGHT','PUTAMEN_RIGHT','THALAMUS_RIGHT','HIPPOCAMPUS_RIGHT','AMYGDALA_RIGHT','ACCUMBENS_RIGHT'};

% run infomap
pfm_infomap(ConcatenatedCifti,DistanceMatrix,PfmDir,GraphDensities,NumberReps,DistanceCutoff,BadVertices,Structures,nWorkers,InfoMapBinary);

% remove some intermediate files (optional)
system(['rm ' Subdir '/pfm/*.net']);
system(['rm ' Subdir '/pfm/*.clu']);
system(['rm ' Subdir '/pfm/*Log*']);

% define inputs;
Input = [PfmDir '/Bipartite_PhysicalCommunities.dtseries.nii'];
Output = 'Bipartite_PhysicalCommunities+SpatialFiltering.dtseries.nii';
MinSize = 50; % in mm^2

% perform spatial filtering
pfm_spatial_filtering(Input,PfmDir,Output,MidthickSurfs,MinSize,WorkbenchBinary);

%% Step 5: Algorithmic assignment of network identities to infomap communities.

% load the priors;
load('priors.mat');

% define inputs;
Ic = ft_read_cifti_mod([PfmDir '/Bipartite_PhysicalCommunities+SpatialFiltering.dtseries.nii']);
Output = 'Bipartite_PhysicalCommunities+AlgorithmicLabeling';
Column = 6; % column 6, representing graph density 0.01% in this example.

% run the network identification algorithm;
pfm_identify_networks(ConcatenatedCifti,Ic,MidthickSurfs,Column,Priors,Output,PfmDir,WorkbenchBinary);

%% Step 6: Review algorithmic network assignments, optionally adjust labels manually if needed.

% define inputs
XLS = [PfmDir '/Bipartite_PhysicalCommunities+AlgorithmicLabeling_NetworkLabels+ManualDecisions.xls']; 
Output = 'Bipartite_PhysicalCommunities+FinalLabeling';

% OPTIONAL: update network assignments according to manual decisions;
pfm_parse_manual_decisions(Ic,Column,MidthickSurfs,Priors,XLS,Output,PfmDir,WorkbenchBinary);

%% Step 7: Calculate size of each functional brain network

% define inputs
FunctionalNetworks = ft_read_cifti_mod([PfmDir '/Bipartite_PhysicalCommunities+FinalLabeling.dlabel.nii']);
VA = ft_read_cifti_mod([Subdir '/fs_LR/fsaverage_LR32k/' Subject '.midthickness_va.32k_fs_LR.dscalar.nii']);
Structures = {'CORTEX_LEFT','CORTEX_RIGHT'}; % in this case, cortex only.

% calculate the size of each functional brain network
NetworkSize = pfm_calculate_network_size(FunctionalNetworks,VA,Structures);

close all; % blank slate
H = figure; % prellocate parent figure
set(H,'position',[1 1 325 400]); hold;

% unique functional networks;
uCi = unique(nonzeros(FunctionalNetworks.data));

% sweep through
% the networks;
for i = 1:length(uCi)
    Tmp = nan(1,length(Priors.NetworkLabels));
    Tmp(i) = NetworkSize(i);
    barh(Tmp,'FaceColor',Priors.NetworkColors(i,:));
    text((NetworkSize(i)+0.1),i,[num2str(NetworkSize(i),3) '%']);
end

% make it pretty;
yticklabels(Priors.NetworkLabels); 
yticks(1:length(uCi)); ylim([0 21]);
xlim([0 20]); xticks(0:5:20);
set(gca,'fontname','arial','fontsize',10,'TickLength',[0 0],'TickLabelInterpreter','none');
xlabel('% of Cortical Surface');
print(gcf,[PfmDir '/FunctionalNetworkSizes'],'-dpdf');

