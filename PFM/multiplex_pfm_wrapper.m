function multiplex_pfm_wrapper(C,LayerIdx,DistanceMatrix,OutDir,Densities,NumberReps,RelaxRate,RelaxLimit,MinDistance,BrainStructures,BadVerts)
% cjl; cjl2007@med.cornell.edu;
rng(44); % for reproducibility.

% define "resource" directories;
ResourceDir = '/home/charleslynch/MultiEchofMRI-Pipeline/res0urces';
addpath(genpath(ResourceDir));

% define the infomap binary location;
InfoMap = '/home/charleslynch/miniconda3/bin/infomap';

% define a list of regions 
% to be considered in the
% community detection routine;
if isempty(BrainStructures)
    BrainStructures = unique(C.brainstructurelabel);  
end

% infer the number of sessions;
nSessions = length(unique(LayerIdx));

% make output
% directory;
mkdir(OutDir);

% change
% directory;
cd(OutDir);

% extract brain structure;
brain_structure = C.brainstructure;
brain_structure(brain_structure < 0) = [];
brain_structure_labels = C.brainstructurelabel;

% count the number of cortical vertices;
nCorticalVertices = nnz(C.brainstructure==1) + nnz(C.brainstructure==2);

% index of relevant vertices and voxels;
GoodVerts = find(ismember(brain_structure,find(ismember(brain_structure_labels,BrainStructures))));
GoodVerts(ismember(GoodVerts,BadVerts)) = [];

% load and trim
% the distance matrix;
D = smartload(DistanceMatrix);
D((nCorticalVertices + 1):end,(nCorticalVertices + 1):end) = 0; % note: this will ensure that subcortical-subcortical edges are set to zero;
D = D(GoodVerts,GoodVerts);

% sweep through
% the graph densities;
for d = 1:length(Densities)
    
    % preallocate;
    intra = []; % intra module links;
    
    % sweep through
    % the sessions;
    for s = 1:nSessions
        
        % fc matrix with local edges removed;
        m = paircorr_mod(C.data(GoodVerts,LayerIdx==s)');
        m(eye(size(m,1))==1) = 0; % remove the diagonal;
        m(D<=MinDistance) = 0; % remove local edges and set subcortical-subcortical FC to zero;
        m(isnan(m)) = 0; % remove nans
        
        % calculate the
        % number of nodes;
        num_nodes = size(m,1);
        
        % preallocate;
        connection_mat = false(num_nodes,num_nodes); %
        
        % sweep all
        % of the nodes;
        for i = 1:num_nodes
            
            % if nonzero
            % values exist;
            if any(m(:,i))
                
                % sort node fc from
                % strongest to weakest;
                [~,idx] = sort(m(:,i),'descend');
                connection_mat(idx(1:ceil(num_nodes .* Densities(d))),i) = true;
                connection_mat(i,idx(1:ceil(num_nodes .* Densities(d)))) = true;
                
            end
            
        end
        
        m = triu(m,1); % only the upper triangle;
        ind = find(triu(connection_mat,1)); % index of relevant FC features;
        
        % get vertex pair numbers;
        [x,y] = ind2sub(size(m),ind);
        
        % log the intra layer edges;
        intra = [intra ; ones(length(x),1)*s x ones(length(x),1)*s y m(ind)];

    end
    
    O = C; % preallocate an output variable;
    O.data = zeros(size(O.data,1),1);
    
    % make the
    % multiplex
    % pajek file;
    nodes = size(m,1);
    nodenum = 1:nodes;
    fid = fopen([OutDir '/MultiPlex_Density' num2str(Densities(d)) '.net'],'W');
    fprintf(fid,'*Vertices %d\n',nodes);
    fprintf(fid,'%d "%d"\n',[nodenum; nodenum]);
    fprintf(fid,'*Multilayer\n');
    fprintf(fid,'# intra\n');
    fprintf(fid,'# layer_id node_id layer_id node_id weight\n');
    fprintf(fid,'%d %d %d %d %f\n',intra');
    fclose(fid);
    
end

O = C;
O.data = zeros(size(C.data,1),nSessions); % vertices/voxels x nSessions ;
clear C; % clear large intermediate file;

% sweep through
% the graph densities;
for d = 1:length(Densities)
    
    % sweep through
    % the relax rates;
    for r = 1:length(RelaxRate)
        
        % run the map equation (this may take awhile);
        % note: multilayer relax up/down set to a single layer to respect temporal nature of the network;
        system([InfoMap ' ' OutDir '/MultiPlex_Density' num2str(Densities(d)) '.net ' OutDir ' --clu -2 -s 42 -N ' num2str(NumberReps) ' --multilayer-relax-limit ' num2str(RelaxLimit) ' --multilayer-relax-rate ' num2str(RelaxRate(r)) ' --no-self-links >> ' OutDir '/MultiPlex_Density' num2str(Densities(d)) '_RelaxRate' num2str(RelaxRate(r)) '_LogFile_' datestr(datetime) '.txt']); %
        
        % rename the state communities output file (purely for convenience) and read in the results;
        system(['mv ' OutDir '/MultiPlex_Density' num2str(Densities(d)) '.clu ' OutDir '/MultiPlex_PhysicalCommunities_Density' num2str(Densities(d)) '_RelaxRate' num2str(RelaxRate(r)) '.clu']); % ocd, but whatever.
        system(['mv ' OutDir '/MultiPlex_Density' num2str(Densities(d)) '_states.clu ' OutDir '/MultiPlex_StateCommunities_Density' num2str(Densities(d)) '_RelaxRate' num2str(RelaxRate(r)) '.clu']); %
        Output = readmatrix([OutDir '/MultiPlex_StateCommunities_Density' num2str(Densities(d)) '_RelaxRate' num2str(RelaxRate(r)) '.clu'],'Delimiter',' ','NumHeaderLines',9,'FileType','text');
        [~,reorder] = sort(Output(:,4)); % reorder by node (vertex / voxel);
        Output = Output(reorder,:);
        
        InfoMap_Ci = O; % preallocate output file;
        
        % sweep all the nodes;
        for i = 1:length(GoodVerts)
            
            % index of vertex / voxel "i"
            idx = find(Output(:,4)==i); % note: index length should == nSessions;
            
            % sweep the layers
            % (the study visits);
            for ii = 1:length(idx)
                InfoMap_Ci.data(GoodVerts(i),Output(idx(ii),5)) = Output(idx(ii),2);
            end
            
        end
        
        % remove small network pieces
        % (those with < 10 vertices)
        
        % sweep the layers
        for i = 1:nSessions
            
            % define unique communities
            uCi = unique(nonzeros(InfoMap_Ci.data(:,i)));
            
            % preallocate some variables;
            rm_idx = zeros(length(uCi),1); % "remove index"
            uCi_idx = cell(length(uCi),1); % community indices
            
            % sweep through
            % unique communities;
            for ii = 1:length(uCi)
                
                % save an index of vertices
                % affiliated with this community;
                uCi_idx{ii} = find(InfoMap_Ci.data(:,i)==uCi(ii));
                
                % mark this community for removal if < 10 vertices;
                if length( find(InfoMap_Ci.data(:,i)==uCi(ii)) ) < 10
                    rm_idx(ii) = 1;
                end
                
            end
            
            % rm. small communities;
            for ii = 1:length(rm_idx)
                if rm_idx(ii)==1
                    InfoMap_Ci.data(uCi_idx{ii},i) = 0;
                end
            end
            
        end
        
        % write out the multiplex state communities;
        ft_write_cifti_mod([OutDir '/MultiPlex_StateCommunities_Density' num2str(Densities(d))...
        '_RelaxRate' num2str(RelaxRate(r)) '.dtseries.nii'],InfoMap_Ci);
        
    end
    
end

end