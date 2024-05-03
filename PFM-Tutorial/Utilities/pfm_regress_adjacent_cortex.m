function [Output] = pfm_regress_adjacent_cortex(Input,DistanceMatrix,Distance)

% count the number of cortical vertices (should be 59412);
nCorticalVertices = nnz(Input.brainstructure==1) + nnz(Input.brainstructure==2);

% load distance matrix;
D = smartload(DistanceMatrix);

% index of subcortical voxels
SubcortVoxels = (nCorticalVertices+1):size(D,1);

% trim to be subcortex x cortex;
d = D(SubcortVoxels,1:nCorticalVertices);

% find all voxels adjacent to cortex;
idx = find(min(d,[],2) <= Distance);
clear d % clear intermediate file;

% preallocate;
Output = Input;

% sweep all
% subcortical voxels
% nearby gray matter;
for i = 1:length(idx)
    
    % extract nearby gm signals;
    nb_gm_ts = Input.data(D(idx(i),:)<=Distance,:); %
    
    % average; if needed
    if size(nb_gm_ts,1) > 1
        nb_gm_ts = mean(nb_gm_ts);
    end
    
    % remove (possible) contamination of nearby cortical signals via linear regression
    [~,~,Output.data(SubcortVoxels(idx(i)),:),~,~] = regress(Input.data(SubcortVoxels(idx(i)),:)',[nb_gm_ts' ones(size(Input.data,2),1)]);
    
end

end

% subfunctions
function out = smartload(matfile)
out = load(matfile);
names = fieldnames(out);
out = eval(['out.' names{1}]);
end