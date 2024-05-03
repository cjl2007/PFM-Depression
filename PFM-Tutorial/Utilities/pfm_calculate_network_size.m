function [NetworkSize] = pfm_calculate_network_size(FunctionalNetworks,VA,Structures)

% read in the
% network CIFTI;
if ischar(FunctionalNetworks)
    FunctionalNetworks = ft_read_cifti_mod(FunctionalNetworks);
end

% read in the
% surface areas;
if ischar(VA)
    VA = ft_read_cifti_mod(VA);
end

VA = VA.data; % extract data; tack on uniform values for subcortex
VA = [VA ; ones((size(FunctionalNetworks.data,1)-size(VA,1)),1)];

% extract brain structure;
BrainStructure = FunctionalNetworks.brainstructure;
BrainStructure(BrainStructure < 0) = [];
BrainStructureLabels = FunctionalNetworks.brainstructurelabel;

% index of relevant vertices and/or voxels;
Idx = find(ismember(BrainStructure,find(ismember(BrainStructureLabels,Structures))));
FunctionalNetworks = FunctionalNetworks.data(Idx); % networks
VA = VA(Idx); % surface areas

% unique functional networks
uCi = unique(nonzeros(FunctionalNetworks));

% preallocate;
NetworkSize = zeros(1,length(uCi));

% sweep through
% the networks;
for i = 1:length(uCi)
    NetworkSize(i) = sum(VA(FunctionalNetworks==uCi(i))) / sum(VA) * 100;
end

end



