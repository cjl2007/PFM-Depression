function grayplot_qa_func(Subdir)

% make grayplot dir. 
system(['mkdir ' Subdir '/func/qa/GrayPlotQA/']);

% define number of runs;
sessions = dir([Subdir '/func/rest/session_*']);

% sweep through sessions;
for s = 1:length(sessions)
    
    % define number of runs;
    runs = dir([Subdir '/func/rest/session_' num2str(s) '/run_*']);
    
    % sweep through runs;
    for r = 1:length(runs)
        
        % try
        
        H = figure; % prellocate parent figure
        set(H,'position',[1 1 1200 800]);
        [~] = evalc('hold');
        title(['Session:' num2str(s) ', Run:' num2str(r)]);
        
        % load various CIFTIs;
        a = ft_read_cifti_mod([Subdir '/func/rest/session_' num2str(s) '/run_' num2str(r) '/Rest_OCME.dtseries.nii']);
        b = ft_read_cifti_mod([Subdir '/func/rest/session_' num2str(s) '/run_' num2str(r) '/Rest_OCME+MEICA.dtseries.nii']);
        c = ft_read_cifti_mod([Subdir '/func/rest/session_' num2str(s) '/run_' num2str(r) '/Rest_OCME+MEICA+MGTR.dtseries.nii']);

        % remove mean and make std=1
        A_tSNR = nanmean(mean(a.data,2) ./ std(a.data,[],2));
        A = detrend(a.data','constant');
        A = A./repmat(std(A),size(A,1),1);
        A = A';
        
        % remove mean and make std=1
        B_tSNR = nanmean(mean(b.data,2) ./ std(b.data,[],2));
        B = detrend(b.data','constant');
        B = B./repmat(std(B),size(B,1),1);
        B = B';
        
        % remove mean and make std=1
        C_tSNR = nanmean(mean(c.data,2) ./ std(c.data,[],2));
        C = detrend(c.data','constant');
        C = C./repmat(std(C),size(C,1),1);
        C = C';
        
        % load motion estimates;  (filtered FD; 2 TRs);
        rp = load([Subdir '/func/rest/session_' num2str(s) '/run_' num2str(r) '/MCF.par']);
        tr = load([Subdir '/func/rest/session_' num2str(s) '/run_' num2str(r) '/TR.txt']);
        [fd,rp]=calc_fd(rp,tr);

        % plot head motion;
        subaxis(12,1,[1 2],'MB',0.05,'MT',0.05,'ML',0.05,'MR',0.05)
        plot(fd,'r'); hold;
        plot(normalize(rp-mean(rp)),'Color',[.5 .5 .5 .5]);
        xlim([0 length(fd)]);
        ylabel('mm');
        ylim([0 1]);
        yticks(0:2)
        xticks('');
        set(gca,'FontName','Arial','FontSize',10,'TickLength',[0 0]);
        title(['Head Motion (Median FD: ' num2str(median(fd),'%.2f') ', Max FD: ' num2str(max(fd),'%.2f') ', ' num2str(100 * length(find(fd < 0.3)) / length(fd),'%.2f') '% < 0.3mm)']);

        try % extract physio information;
        [resp] = extract_resp([Subdir '/physio/unprocessed/rest/session_' num2str(s) '/run_' num2str(r) '/'],50);
        catch
        end
        
        % plot respiration;
        subaxis(12,1,[3 4],'MB',0.05,'MT',0.05,'ML',0.05,'MR',0.05)
        if exist('resp','var')
        plot(smoothdata(zscore(resp)),'b');
        xlim([0 length(resp)]);
        end
        ylabel('z-score');
        yticklabels('');
        xlabel('');
        xticks('');
        xticklabels('');
        set(gca,'FontName','Arial','FontSize',10,'TickLength',[0 0]);
        title('Respiration Belt');
        ylim([-3 3]);
        clear resp;
        
        % plot optimally-combined time-series
        subaxis(12,1,[5 6],'MB',0.05,'MT',0.05,'ML',0.05,'MR',0.05)
        imagesc(A);
        hline(59412,'g');
        text(0.005,0.4,'Cortex','Units','normalized','Color','white','FontSize',10,'FontName','Arial');
        text(0.005,0.15,'Subcortex + Cerebellum','Units','normalized','Color','white','FontSize',10,'FontName','Arial');
        colormap(gray);
        caxis([-3 3 ]);
        yticks('');        
        xticks('');
        set(gca,'FontName','Arial','FontSize',10,'TickLength',[0 0]);
        title(['Optimally-Combined (OC-ME) (tSNR:' num2str(A_tSNR) ')']);
        
        % plot variance discarded after ME-ICA
        subaxis(12,1,[7:8],'MB',0.05,'MT',0.05,'ML',0.05,'MR',0.05)
        imagesc(A-B);
        hline(59412,'g');
        text(0.005,0.4,'Cortex','Units','normalized','Color','white','FontSize',10,'FontName','Arial');
        text(0.005,0.2,'Subcortex + Cerebellum','Units','normalized','Color','white','FontSize',10,'FontName','Arial');
        colormap(gray);
        caxis([-3 3 ]);
        yticks('');
        xticks('');
        set(gca,'FontName','Arial','FontSize',10,'TickLength',[0 0]);
        title('Discarded by ME-ICA');
        
        % plot variance retained by ME-ICA;
        subaxis(12,1,[9:10],'MB',0.05,'MT',0.05,'ML',0.05,'MR',0.05)
        imagesc(B-C);
        hline(59412,'g');
        text(0.005,0.4,'Cortex','Units','normalized','Color','white','FontSize',10,'FontName','Arial');
        text(0.005,0.2,'Subcortex + Cerebellum','Units','normalized','Color','white','FontSize',10,'FontName','Arial');
        colormap(gray);
        caxis([-3 3 ]);
        yticks('');
        xticks('');
        set(gca,'FontName','Arial','FontSize',10,'TickLength',[0 0]);
        title('Discarded by MGTR');
        
        % fully denoised;
        subaxis(12,1,[11:12],'MB',0.05,'MT',0.05,'ML',0.05,'MR',0.05)
        imagesc(C);
        hline(59412,'g');
        text(0.005,0.4,'Cortex','Units','normalized','Color','white','FontSize',10,'FontName','Arial');
        text(0.005,0.2,'Subcortex + Cerebellum','Units','normalized','Color','white','FontSize',10,'FontName','Arial');
        colormap(gray);
        caxis([-3 3 ]);
        yticks('');
        set(gca,'FontName','Arial','FontSize',10,'TickLength',[0 0]);
        title(['Retained by ME-ICA + MGTR (tSNR:' num2str(C_tSNR) ')']);
        xlabel('Time (TRs)');
        
        % save grayplot;
        saveas(gcf,[Subdir '/func/qa/GrayPlotQA/GrayPlotQA_Session' num2str(s) '_Run' num2str(r) '.jpg']); % save image
        close all;
        
        %catch
        %end
        
    end
    
end

end

function [fd,rp]=calc_fd(rp,tr)

% nyquist freq.
nyq = (1/tr)/2;

% create a tailored
% stop band filter;
stopband = [0.2 (nyq-0.02)];
[B,A] = butter(10,stopband/nyq,'stop');

% apply stopband filter 
for i = 1:size(rp,2)
    rp(:,i) = filtfilt(B,A,rp(:,i));
end

% calc. backward difference;
n_trs = round(2.5 / tr);

fd = rp; % preallocate
fd(1:n_trs,:) = 0; % by convention

% sweep the columns;
for i = 1:size(rp,2)
    for ii = (n_trs+1):size(fd,1)
        fd(ii,i) = abs(rp(ii,i)-rp(ii-n_trs,i));
    end
end

fd_ang = fd(:,1:3); % convert rotation columns into angular displacement...
fd_ang = fd_ang / (2 * pi); % fraction of circle
fd_ang = fd_ang * 100 * pi; % multiplied by circumference

fd(:,1:3) = []; % delete rotation columns,
fd = [fd fd_ang]; % add back in as angular displacement
fd = sum(fd,2); % sum

rp(:,1:3) = [];
rp = [rp fd_ang];

end

function [resp] = extract_resp(pdir,freq)
% cjl;

% define various files 
info = dir([pdir '/*Info.log']);
resp = dir([pdir '/*RESP.log']);

% log information.log file 
info = importdata([pdir '/' info.name],' ',10);
acq_start = info.data(1,3); % start of acquisition  
acq_end = info.data(end,4); % end of acquisition

% notes: the start and stop time is determined by the info file;

% load respiration data 
resp = importdata([pdir '/' resp.name],' ',8);

% reformat respiration data;
resp_time = resp.textdata([7:size(resp.textdata,1)],1); % ignore the first 6 rows; which contain text
resp.data(strcmp(resp_time,'PULS_TRIGGER'))=[]; % remove time points contaminated by the trigger;
resp_time(strcmp(resp_time,'PULS_TRIGGER'))=[]; % remove time points contaminated by the trigger;
resp_time = str2num(cell2mat(resp_time)); % convert data type (cell --> double); 
resp_full = horzcat(resp_time,resp.data); % combine respiratory data with time stamps

% at this point; resp_full is too long relative to the info file. 
% if this is chuck's 14.43 min. 5 echo resting-state scan, it is probably off by about 30 seconds; 
% I suspect this is related to the single-band reference images that
% collected prior to the main functional data. Specifically, it seems that
% during the ~30 second period of time that the SBrefs are collected, data
% is being written to _Resp.log but not _Info.log. Note also that the first
% trigger occurs after all the SBrefs (1 per TE) are collected.

% the takeaway from above is that we need to find the point in
% resp_full(:,1) that is closest in time to the first trigger (acq_start).
% Everything in before acq_start is generally not of interest, so we can discard it.  

% find the row number in resp_full 
% containing a timestamp closest in
% its absolute value to acq_start.
[~,idx] = min(abs(acq_start - resp_full(:,1)));
resp_full([1:(idx-1)],:) = [];  % 

% do the same thing for the end; in practice this probably isnt needed but good to double check;
[~,idx] = min(abs(acq_end- resp_full(:,1)));
resp_full([(idx+1):size(resp_full,1)],:) = [];

% fill outliers and lightly smooth data 
resp = filloutliers(resp_full(:,2),'linear','movmedian',100);
resp = smoothdata(resp,'sgolay',freq);

end


