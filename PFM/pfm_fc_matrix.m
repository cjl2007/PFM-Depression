function pfm_fc_matrix(Ci,FC,Priors,OutFile)

% Ci == vector of network memberships (1  x N)
% FC == Functional connectivity matrix (N x N); not yet sorted by network memberships
% Priors == .mat file used for PFM-ing; contains network colors (RGB codes)
% OutFile == file name 

H = figure; % prellocate parent figure (these are hard set parameters);
subaxis(10,10,1:10:81,'MB',0.005,'MT',0.095,'ML',0.05,'MR',0.05); % left subplot

% unique networks;
uCi = unique(Ci);

% preallocate variable that controls 
% size of network representation in stacked bar;
uCi_PercentageOfNodes = zeros(1,length(uCi));

% calculate percentage of cortical
% vertices belong to network "i"
for i = 1:length(uCi)
    uCi_PercentageOfNodes(i) = sum(Ci==i)/(length(Ci));
end

% create a stacked bar graph where colors indicate 
% network membership of nodes in the functional connectivity matrix
Tmp = bar([flip(uCi_PercentageOfNodes) ; nan(size(uCi_PercentageOfNodes))],"stacked"); % NaN work around not needed for matlab 2019b and later;
xlim([0.9 1.1]); % these values work okay on scully;
ylim([0 sum(uCi_PercentageOfNodes)]);
axis('off');

% sweep through
% the networks 
for i = 1:length(uCi)
    Idx = (length(uCi)+1)-i;
    Tmp(i).FaceColor = Priors.NetworkColors(Idx,:);
end

% bottom subplot;
subaxis(10,10,92:1:100,'MB',0.05,'MT',0.05,'ML',0,'MR',0.05);

% create another stacked bar graph where colors indicate 
% network membership of nodes in the functional connectivity matrix
Tmp = barh([uCi_PercentageOfNodes ; nan(size(uCi_PercentageOfNodes))],"stacked");
ylim([0.9 1.1]);
xlim([0 sum(uCi_PercentageOfNodes)]);
axis('off');

% sweep through
% the networks; 
for i = 1:length(uCi)
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
  hline(count,'k');
  vline(count,'k');
end

axis('off');
colormap(redbluecmap);
caxis([-1 1]);
set(H,'position',[1 1 315 315]);
print(gcf,OutFile,'-dpdf');
close;

end