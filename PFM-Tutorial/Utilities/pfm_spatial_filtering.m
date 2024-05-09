function pfm_spatial_filtering(Input,OutDir,Output,MidthickSurfs,MinSize,WorkbenchBinary)

% read in input and 
% preallocate output cifti;
Input = ft_read_cifti_mod(Input);
O = Input; O.data = zeros(size(Input.data)); % blank slate;

% sweep through the graph densities
for s = 1:size(Input.data,2)
    
    A = Input; % preallocate
    A.data = Input.data(:,s);
    
    % define the unique communities;
    uCi = unique(nonzeros(A.data(:,1)));
    uCi(isnan(uCi))=[];
    
    % sweep communities;
    for i = 1:length(uCi)
        
        B = A; % preallocate;
        
        % set all other
        % communities to zero;
        B.data(B.data~=uCi(i)) = 0;
        B.data(B.data~=0) = 1;
        
        % write out temporary CIFTI;
        ft_write_cifti_mod([OutDir '/temp_' num2str(i)],B);
        
        % remove small clusters;
        system([WorkbenchBinary ' -cifti-find-clusters ' OutDir '/temp_' num2str(i) '.dtseries.nii 0 ' num2str(MinSize) ' 0 ' num2str(MinSize) ' COLUMN ' OutDir '/temp_' num2str(i) '.dtseries.nii -left-surface ' MidthickSurfs{1} ' -right-surface ' MidthickSurfs{2} ' -merged-volume']);
        C = ft_read_cifti_mod([OutDir '/temp_' num2str(i) '.dtseries.nii']); % read in CIFTI post-removal of small networks;
        O.data(C.data~=0,s) = uCi(i); % log "large" clusters;
        
    end
    
    % remove some
    % intermediate files;
    system(['rm ' OutDir...
    '/temp_*.dtseries.nii']);
    
end

% write out a cifti where we have removed small network pieces and any blank vertices are assigned to the closest community;
ft_write_cifti_mod([OutDir '/' Output],O);
system([WorkbenchBinary ' -cifti-dilate ' OutDir '/' Output ' COLUMN 50 50 -left-surface ' MidthickSurfs{1} ' -right-surface ' MidthickSurfs{2} ' ' OutDir '/' Output ' -nearest']);

end