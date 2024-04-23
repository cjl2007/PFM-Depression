function make_distance_matrix(RefCifti,MidthickSurfs,OutDir,nslots)
% cjl2007@med.cornell.edu; 

% start parpool;
pool = parpool('local',nslots);

try % make hidden directory
    mkdir([OutDir '/tmp/']);
catch
end

% load
% reference CIFTI
if ischar(RefCifti)
    RefCifti = ft_read_cifti_mod(RefCifti);
end

RefCifti.data=[]; % remove data, not needed

% load midthickness surfaces 
LH = gifti(MidthickSurfs{1});
RH = gifti(MidthickSurfs{2});

% find cortical vertices on surface cortex (not medial wall)
lh_idx = RefCifti.brainstructure(1:length(LH.vertices))~=-1;
rh_idx = RefCifti.brainstructure((length(LH.vertices)+1):(length(LH.vertices)+length(RH.vertices)))~=-1;

% preallocate "reference verts"
LH_verts=1:length(LH.vertices);
RH_verts=1:length(RH.vertices);

% cortical vertices only
LH_verts=LH_verts(lh_idx);
RH_verts=RH_verts(rh_idx);

% sweep through vertices
parfor i = 1:length(LH_verts)
    
    % calculate geodesic distances from vertex i
    system(['wb_command -surface-geodesic-distance ' MidthickSurfs{1} ' ' num2str(LH_verts(i)-1) ' ' OutDir '/tmp/temp_' num2str(i) '.shape.gii']);
    temp = gifti([OutDir '/tmp/temp_' num2str(i) '.shape.gii']);
    system(['rm ' OutDir '/tmp/temp_' num2str(i) '.shape.gii']);
    lh(:,i) = temp.cdata(lh_idx); % log distances
        
end

% convert to uint8
lh = uint8(lh);

% sweep through vertices
parfor i = 1:length(RH_verts)
    
    % calculate geodesic distances from vertex i
    system(['wb_command -surface-geodesic-distance ' MidthickSurfs{2} ' ' num2str(RH_verts(i)-1) ' ' OutDir '/tmp/temp_' num2str(i) '.shape.gii']);
    temp = gifti([OutDir '/tmp/temp_' num2str(i) '.shape.gii']);
    system(['rm ' OutDir '/tmp/temp_' num2str(i) '.shape.gii']);
    rh(:,i) = temp.cdata(rh_idx); % log distances
    
end

% delete 
% parpool
delete(pool);

% remove temp dir.;
[~,~]=system(['rm -rf ' OutDir '/tmp/']);

% convert to uint8
rh = uint8(rh);

% piece together results (999 = inter-hemispheric)
top = [lh ones(length(lh),length(rh))*999]; % lh & dummy rh
bottom = [ones(length(rh),length(lh))*999 rh]; % dummy lh & rh
D = uint8([top;bottom]); % combine hemispheres; cortical surface only so far 

% save distance matrix;
save([OutDir '/DistanceMatrixCortexOnly'],'D','-v7.3');

% extract coordinates for all cortical vertices 
coords_surf=[LH.vertices; RH.vertices]; % combine hemipsheres 
surf_indices_incifti = RefCifti.brainstructure > 0 & RefCifti.brainstructure < 3;
surf_indices_incifti = surf_indices_incifti(1:size(coords_surf,1));
coords_surf = coords_surf(surf_indices_incifti,:);
coords_subcort = RefCifti.pos(RefCifti.brainstructure>2,:);
coords = [coords_surf;coords_subcort]; % combine 

% compute euclidean distance 
% between all vertices & voxels 
D2 = uint8(pdist2(coords,coords));

% combine distance matrices; geodesic & euclidean  
D = [D ; D2(size(D,1)+1:end,1:size(D,2))]; % vertcat
D = [D  D2(1:size(D,1),size(D,2)+1:end)]; % horzcat 
clear D2;

% save distance matrix;
save([OutDir '/DistanceMatrix'],'D','-v7.3');

% clear 
% distances
clear D;

end