function pfm_parse_manual_decisions(Ic,Col,MidthickSurfs,Priors,XLS,OutFile,OutDir,WorkbenchBinary)
% cjl; cjl2007@med.cornell.edu;

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

% unique infomap communities;
uCi = unique(nonzeros(Ic.data));

% read in the labels;
Labels = readtable(XLS);

O = Ic; % preallocate output
O.data = zeros(size(Ic.data));

% sweep through
% the communities;
for i = 1:length(uCi)
    
    % check to see if a manual decision was made;
    if isempty(Labels.Network_ManualDecision{i}) 
        
        % log the network identity;
        O.data(Ic.data(:,1)==uCi(i),1) = find(strcmp(Priors.NetworkLabels,Labels.Network{i}));
        
    else
        
        % only if the specified name matches a prior label;
        if ismember(Labels.Network_ManualDecision{i},Priors.NetworkLabels)
            
            % log the network identity;
            O.data(Ic.data(:,1)==uCi(i),1) = find(strcmp(Priors.NetworkLabels,Labels.Network_ManualDecision{i}));
            
        end
        
    end
    
end

% write out the temporary cifti file;
%O.data(O.data==find(strcmp(Priors.NetworkLabels,'Noise')))=0;
ft_write_cifti_mod([OutDir '/Tmp.dtseries.nii'],O);
system([WorkbenchBinary ' -cifti-dilate ' OutDir '/Tmp.dtseries.nii COLUMN 10 10 ' OutDir '/Tmp.dtseries.nii -left-surface ' MidthickSurfs{1} ' -right-surface ' MidthickSurfs{1} ' -nearest']);

% write out the first network;
system(['echo ' char(Priors.NetworkLabels{1}) ' > ' OutDir '/LabelListFile.txt ']);
system(['echo 1 ' num2str(round(Priors.NetworkColors(1,1)*255)) ' ' num2str(round(Priors.NetworkColors(1,2)*255)) ' ' num2str(round(Priors.NetworkColors(1,3)*255)) ' 255 >> ' OutDir '/LabelListFile.txt ']);

% sweep through the networks;
for i = 2:length(Priors.NetworkLabels)
    
    system(['echo ' char(Priors.NetworkLabels{i}) ' >> ' OutDir '/LabelListFile.txt ']);
    system(['echo ' num2str(i) ' ' num2str(round(Priors.NetworkColors(i,1)*255)) ' ' num2str(round(Priors.NetworkColors(i,2)*255)) ' ' num2str(round(Priors.NetworkColors(i,3)*255)) ' 255 >> ' OutDir '/LabelListFile.txt ']);
    
end

% make dense label file + network borders;
system([WorkbenchBinary ' -cifti-label-import ' OutDir '/Tmp.dtseries.nii ' OutDir '/LabelListFile.txt ' OutDir '/' OutFile '.dlabel.nii -discard-others']);
system([WorkbenchBinary ' -cifti-label-to-border ' OutDir '/' OutFile '.dlabel.nii -border ' MidthickSurfs{1} ' ' OutDir '/' OutFile '.L.border']); % LH
system([WorkbenchBinary ' -cifti-label-to-border ' OutDir '/' OutFile '.dlabel.nii -border ' MidthickSurfs{2} ' ' OutDir '/' OutFile '.R.border']); % RH

% remove some intermediate files;
system(['rm ' OutDir '/Tmp.dtseries.nii']);
system(['rm ' OutDir '/LabelListFile.txt']);

end



