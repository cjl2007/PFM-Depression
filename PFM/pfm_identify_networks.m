function pfm_identify_networks(C,Ic,MidthickSurfs,Col,Priors,OutFile,OutDir)
% cjl; cjl2007@med.cornell.edu;
rng(44); % for reproducibility.

% make output directory;
if ~exist(OutDir,'dir')
    mkdir(OutDir);
end

% read in the 
% resting-state data
if ischar(C)
    C = ft_read_cifti_mod(C);
end

% count the number of cortical vertices; note: this should be 59412;
nCorticalVertices = nnz(C.brainstructure==1) + nnz(C.brainstructure==2);

% calculate the functional connectivity matrix;
m = paircorr_mod(C.data(1:nCorticalVertices,:)');
m(eye(size(m,1))==1) = 0; % remove the diagonal;
m(isnan(m)) = 0; % remove nans

% read in 
% infomap
% communities 
if ischar(Ic)
    Ic = ft_read_cifti_mod(Ic); % Ic == infomap communities
end

% extract
% graph density
% of interest ;
if ~isempty(Col)
    Ic.data = Ic.data(:,Col);
end

O = Ic; % preallocate output
O.data = zeros(size(Ic.data));

% unique infomap communities;
uCi = unique(nonzeros(Ic.data));

% preallocate functional connectivity of each community
uCi_FC = zeros(nCorticalVertices,length(uCi)); %

% sweep each of the
% unique communities
for i = 1:length(uCi)
    
    % calculate functional connectivity profile of community "i";
    uCi_FC(:,i) = nanmean(m(:,Ic.data(1:nCorticalVertices)==uCi(i)),2);
    
end

% calculate the spatial
% similarity of FC with priors;
uCi_rho = corr(uCi_FC,Priors.FC);

% preallocate variable representing distance of 
% community spatial locations relative to spatial priors;
uCi_Spatial = zeros(length(uCi),size(Priors.Spatial,2)); %
 
% sweep each of the
% unique communities
for i = 1:length(uCi)
    
    % loop through the networks;
    for ii = 1:size(Priors.Spatial,2)
      
        % average probability of community "i" belonging to functional network "ii";
        uCi_Spatial(i,ii) = mean(Priors.Spatial(Ic.data(1:nCorticalVertices)==uCi(i),ii));
        
    end
    
end

% Loop through the communities and assign network labels based on combination 
% of functional connectivity and spatial locations relative to the priors.  

% preallocate lots of variables;
S.Community = zeros(length(uCi),1);
S.Network = cell(length(uCi),1);
S.Network_ManualDecision = cell(length(uCi),1);
S.R = zeros(length(uCi),1);
S.G = zeros(length(uCi),1);
S.B = zeros(length(uCi),1);
S.FC_Similarity = zeros(length(uCi),1);
S.Spatial_Score = zeros(length(uCi),1);
S.Confidence = zeros(length(uCi),1);

% sweep through the networks;
for i = 1:size(Priors.Spatial,2)-1
    S.(['Alt_' num2str(i) '_Network']) = cell(length(uCi),1);
    S.(['Alt_' num2str(i) '_FC_Similarity']) = zeros(length(uCi),1);
    S.(['Alt_' num2str(i) '_Spatial_Score']) = zeros(length(uCi),1);
end

% sweep through
% the communities;
for i = 1:length(uCi)
    
    % sort from best match to worst match;
    [x,y] = sort(uCi_rho(i,:) .* uCi_Spatial(i,:),'Descend'); 
    % notes: x = similarity values, y = index of networks

    % log all
    % the info;
    S.Community(i) = i;
    S.R(i) = Priors.NetworkColors(y(1),1);
    S.G(i) = Priors.NetworkColors(y(1),2);
    S.B(i) = Priors.NetworkColors(y(1),3);
    S.Network{i} = Priors.NetworkLabels{y(1)};
    S.FC_Similarity(i) = uCi_rho(i,y(1));
    S.Spatial_Score(i) = uCi_Spatial(i,y(1));
    S.Confidence(i) = (x(1) - x(2)) / x(2);
    
    % sweep through the networks;
    for ii = 1:size(Priors.Spatial,2)-1
        S.(['Alt_' num2str(ii) '_Network']){i,1} = Priors.NetworkLabels{y(ii+1)};
        S.(['Alt_' num2str(ii) '_FC_Similarity'])(i) = uCi_rho(i,y(ii+1));
        S.(['Alt_' num2str(ii) '_Spatial_Score'])(i) = uCi_Spatial(i,y(ii+1));
    end
    
end

O = Ic; % preallocate output
O.data = zeros(size(Ic.data));

% sweep through
% the communities;
for i = 1:length(uCi)
    
    % log the wta network identity;
    O.data(Ic.data(:,1)==uCi(i),1) = find(strcmp(Priors.NetworkLabels,S.Network{i}));
    
end

T = struct2table(S); % write out .xls sheet;
writetable(T,[OutDir '/' OutFile '_NetworkLabels.xls']);

% write out the temporary cifti file;
ft_write_cifti_mod([OutDir '/Tmp.dtseries.nii'],O);

% write out the first network;
system(['echo ' char(Priors.NetworkLabels{1}) ' > ' OutDir '/LabelListFile.txt ']);
system(['echo 1 ' num2str(round(Priors.NetworkColors(1,1)*255)) ' ' num2str(round(Priors.NetworkColors(1,2)*255)) ' ' num2str(round(Priors.NetworkColors(1,3)*255)) ' 255 >> ' OutDir '/LabelListFile.txt ']);

% sweep through the networks;
for i = 2:length(Priors.NetworkLabels)
    
    system(['echo ' char(Priors.NetworkLabels{i}) ' >> ' OutDir '/LabelListFile.txt ']);
    system(['echo ' num2str(i) ' ' num2str(round(Priors.NetworkColors(i,1)*255)) ' ' num2str(round(Priors.NetworkColors(i,2)*255)) ' ' num2str(round(Priors.NetworkColors(i,3)*255)) ' 255 >> ' OutDir '/LabelListFile.txt ']);
    
end

% make dense label file + network borders;
system(['wb_command -cifti-label-import ' OutDir '/Tmp.dtseries.nii ' OutDir '/LabelListFile.txt ' OutDir '/' OutFile '.dlabel.nii -discard-others']);
system(['wb_command -cifti-label-to-border ' OutDir '/' OutFile '.dlabel.nii -border ' MidthickSurfs{1} ' ' OutDir '/' OutFile '.L.border']); % LH
system(['wb_command -cifti-label-to-border ' OutDir '/' OutFile '.dlabel.nii -border ' MidthickSurfs{2} ' ' OutDir '/' OutFile '.R.border']); % RH

% remove some intermediate files;
system(['rm ' OutDir '/Tmp.dtseries.nii']);
system(['rm ' OutDir '/LabelListFile.txt']);

O = Ic; % preallocate output;
O.data = zeros(size(O.data,1),size(uCi_FC,2)); % blank slate
O.data(1:nCorticalVertices,:) = uCi_FC; % log community FC profiles
ft_write_cifti_mod([OutDir '/' OutFile '_FC_WholeBrain.dtseries.nii'],O);

O = Ic; % preallocate output;
O.data = zeros(size(O.data,1),length(uCi)); 

% sweep through
% the communities;
for i = 1:length(uCi)
    O.data(Ic.data==uCi(i),i) = 1;
end

% write out the communities
ft_write_cifti_mod([OutDir '/'...
OutFile '_InfoMapCommunities.dtseries.nii'],O);

O = Ic; % preallocate output;
O.data = zeros(size(O.data,1),length(uCi)); 

% preallocate;
FC = zeros(length(uCi));
Ci = zeros(length(uCi),1);

% sweep infomap
% communities
for i = 1:length(uCi)
    
    % sweep infomap
    % communities
    for ii = 1:length(uCi)

        % calculate functional connectivity strength between communities "i" and "ii"
        tmp = m(Ic.data(1:nCorticalVertices)==uCi(i),Ic.data(1:nCorticalVertices)==uCi(ii));
        FC(i,ii) = nanmean(tmp(:));
        O.data(Ic.data==uCi(ii),i) = FC(i,ii);
        
    end
    
    % log the network assignment for community "i";
    Ci(i) = find(strcmp(Priors.NetworkLabels,S.Network{i}));
    
end

% write out cifti describing avgerage functional connectivity  between communities;
ft_write_cifti_mod([OutDir '/' OutFile '_FC_btwn_InfoMapCommunities.dtseries.nii'],O);

% generate some visualizations showing how well the initial labeling has
% separated communities into clusters with similar functional connectivity

H = figure; % prellocate parent figure (these are hard set parameters);
subaxis(10,10,1:10:81,'MB',0.005,'MT',0.095,'ML',0.05,'MR',0.05); % left subplot

% preallocate variable that controls 
% size of network representation in stacked bar;
PercentageOfNodes = zeros(1,size(Priors.Spatial,2));

% calculate percentage of cortical
% vertices belong to network "i"
for i = 1:size(Priors.Spatial,2)
    PercentageOfNodes(i) = sum(Ci==i)/(length(Ci));
end

% create a stacked bar graph where colors indicate 
% network membership of nodes in the functional connectivity matrix
Tmp = bar([flip(PercentageOfNodes) ; nan(size(PercentageOfNodes))],"stacked"); % NaN work around not needed for matlab 2019b and later;
xlim([0.9 1.1]); % these values work okay on scully;
ylim([0 sum(PercentageOfNodes)]);
axis('off');

% sweep through
% the networks 
for i = 1:size(Priors.Spatial,2)
    Idx = (size(Priors.Spatial,2)+1)-i;
    Tmp(i).FaceColor = Priors.NetworkColors(Idx,:);
end

% bottom subplot;
subaxis(10,10,92:1:100,'MB',0.05,'MT',0.05,'ML',0,'MR',0.05);

% create another stacked bar graph where colors indicate 
% network membership of nodes in the functional connectivity matrix
Tmp = barh([PercentageOfNodes ; nan(size(PercentageOfNodes))],"stacked");
ylim([0.9 1.1]);
xlim([0 sum(PercentageOfNodes)]);
axis('off');

% sweep throughthe networks; 
for i = 1:size(Priors.Spatial,2)
    Tmp(i).FaceColor = Priors.NetworkColors(i,:);
end

% plot the functional connectivity 
% matrix (sorted by network);
Tmp = zeros(10); Tmp(:,1) = 1;
subaxis(10,10,find(Tmp==0),'MB',0.1,'MT',0,'ML',0.1,'MR',0.05);
[~,sort_idx] = sort(Ci); % sort by network membership
imagesc(FC(sort_idx,sort_idx)); hold;

% add
% grid
count = 0;
for i = 1:19
  count = count + length(find(Ci==i));  
  hline(count+0.5,'k');
  vline(count+0.5,'k');
end

% make the plot pretty;
set(gca,'TickLength',[0 0])
xticklabels('');
yticklabels('');
colormap(redbluecmap);
set(H,'position',[1 1 315 315]);
caxis([-1 1]);
print(gcf,[OutDir '/' OutFile...
'_FC_btwn_InfoMapCommunities'],'-dpdf');
close;

end

