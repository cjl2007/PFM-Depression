function pfm_infomap(C,D,OutDir,GraphDensities,NumberReps,MinDistance,BadVerts,Structures,NumberCores,InfoMapBinary)
% cjl; cjl2007@med.cornell.edu;
rng(44); % for reproducibility.

% define the infomap binary location (v2.0.0 used here);
% see https://www.mapequation.org/infomap/#Install
if isempty(InfoMapBinary)
    InfoMapBinary = '/home/charleslynch/miniconda3/bin/infomap';
end

% define a list of regions 
% to be considered in the
% community detection routine;
if isempty(Structures)
    Structures = unique(C.brainstructurelabel);  
end

% number of algorithm repetitions 
% performed at each graph density 
if length(NumberReps) ==1  
NumberReps = repmat(NumberReps,length(GraphDensities));      
end

% make output
% directory ;
mkdir(OutDir);

% change
% directory;
cd(OutDir);

% extract brain structure;
BrainStructure = C.brainstructure;
BrainStructure(BrainStructure < 0) = [];
BrainStructureLabels = C.brainstructurelabel;

% count the number of cortical vertices;
nCorticalVertices = nnz(C.brainstructure==1) + nnz(C.brainstructure==2);

% index of relevant vertices and voxels;
GoodVerts = find(ismember(BrainStructure,find(ismember(BrainStructureLabels,Structures))));
GoodVerts(ismember(GoodVerts,BadVerts)) = [];

% load
% distance
% matrix;
if ischar(D)
    D = smartload(D); 
end

% trim distance matrix
D((nCorticalVertices + 1):end,(nCorticalVertices + 1):end) = 0; % note: this will ensure that subcortical-subcortical edges are set to zero;
D = D(GoodVerts,GoodVerts);

% fc matrix with 
% local edges removed;
m = corr(C.data(GoodVerts,:)');
m(eye(size(m,1))==1) = 0; % remove the diagonal;
m(D<=MinDistance) = 0; % remove local edges and set subcortical-subcortical FC to zero (see note above);
m(isnan(m)) = 0; % remove nans

% calculate the
% number of nodes;
nNodes = size(m,1);

% sweep through the densities
for d = 1:length(GraphDensities)
    
    % preallocate;
    ConnectionMatrix = false(nNodes,nNodes); %
    
    % sweep all
    % of the nodes;
    for i = 1:nNodes
        
        % if nonzero
        % values exist;
        if any(m(:,i))
            
            % sort node fc from
            % strongest to weakest;
            [~,idx] = sort(m(:,i),'descend');
            ConnectionMatrix(idx(1:ceil(nNodes .* GraphDensities(d))),i) = true;
            ConnectionMatrix(i,idx(1:ceil(nNodes .* GraphDensities(d)))) = true;
            
        end
        
    end
    
    % index of relevant FC features;
    ind = find(triu(ConnectionMatrix,1)); 
    
    % determine the edges and values
    [x,y] = ind2sub(size(m),ind);
    z = m(ind); 
    towrite = [x y z];
    
    % make the bipartite pajek file;
    fid = fopen([OutDir '/Bipartite_Density' num2str(GraphDensities(d)) '.net'],'W');
    fprintf(fid,'*Vertices %d\n',size(m,1));
    fprintf(fid,'%d "%d"\n',[1:size(m,1); 1:size(m,1)]);
    fprintf(fid,'*Edges %d\n',length(x));
    fprintf(fid,'%d %d %f\n',towrite');
    fclose(fid);
    
end

% clear
% large
% variable;
clear m ;

% only if multiple 
% workers are requested
if NumberCores > 1
    
    % start a parpool;
    pool = parpool('local',...
    NumberCores);
    
    % sweep through the densities;
    parfor d = 1:length(GraphDensities)
        
        % run InfoMap;
        system([InfoMapBinary ' ' OutDir '/Bipartite_Density' num2str(GraphDensities(d)) '.net ' OutDir '/ --clu -2 -s 42 -N ' num2str(NumberReps(d)) ' --no-self-links >> ' OutDir '/Bipartite_Density' num2str(GraphDensities(d)) '_LogFile_' datestr(datetime) '.txt']);
        
    end
    
    % delete
    % parpool
    delete(pool);
    
else
    
    % sweep through the densities;
    for d = 1:length(GraphDensities)
        
        % run InfoMap;
        system([InfoMapBinary ' ' OutDir '/Bipartite_Density' num2str(GraphDensities(d)) '.net ' OutDir '/ --clu -2 -s 42 -N ' num2str(NumberReps(d)) ' --no-self-links >> ' OutDir '/Bipartite_Density' num2str(GraphDensities(d)) '_LogFile_' datestr(datetime) '.txt']);
        
    end
    
end

InfoMap_Ci = C; % preallocate output Cifti file;
InfoMap_Ci.data = zeros(size(C.data,1),length(GraphDensities)); % vertices x graph densities;

% index of relevant vertices and voxels
idx = find(ismember(BrainStructure,find(ismember(BrainStructureLabels,Structures))));
idx(ismember(idx,BadVerts)) = [];

% sweep through densities;
for d = 1:length(GraphDensities)
    
    % load community partition & reorder;
    output = readmatrix([OutDir '/Bipartite_Density' num2str(GraphDensities(d)) '.clu'],'Delimiter',' ','NumHeaderLines',9,'FileType','text');
    
    % reorder by node;
    [~,reorder] = sort(output(:,1));
    Tmp = output(reorder,2);
    
    % log InfoMap community affiliations;
    InfoMap_Ci.data(idx,d) = Tmp;
   
end

% sweep through densities;
for d = 1:length(GraphDensities)
    
    % define unique communities
    uCi = unique(nonzeros(InfoMap_Ci.data(:,d)));
    
    % preallocate some variables;
    rm_idx = zeros(length(uCi),1); % "remove index"
    uCi_idx = cell(length(uCi),1); % community indices
    
    % sweep through
    % unique communities;
    for i = 1:length(uCi)
        
        % save an index of vertices
        % affiliated with this community;
        uCi_idx{i} = find(InfoMap_Ci.data(:,d)==uCi(i));
        
        % mark this community for removal if < 10 vertices;
        if length( find(InfoMap_Ci.data(:,d)==uCi(i)) ) < 10
            rm_idx(i) = 1;
        end
        
    end
    
    % rm. small communities;
    for i = 1:length(rm_idx)
        if rm_idx(i)==1
            InfoMap_Ci.data(uCi_idx{i},d) = 0;
        end
    end
    
end

% write out the physical bipartite communities;
ft_write_cifti_mod([OutDir '/Bipartite_PhysicalCommunities.dtseries.nii'],InfoMap_Ci);

end

% subfunctions
function out = smartload(matfile)
out = load(matfile);
names = fieldnames(out);
out = eval(['out.' names{1}]);
end