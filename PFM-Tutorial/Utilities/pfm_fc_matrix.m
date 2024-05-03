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

% subfunctions
function h=subaxis(varargin)
%SUBAXIS Create axes in tiled positions. (just like subplot)
%   Usage:
%      h=subaxis(rows,cols,cellno[,settings])
%      h=subaxis(rows,cols,cellx,celly[,settings])
%      h=subaxis(rows,cols,cellx,celly,spanx,spany[,settings])
%
% SETTINGS: Spacing,SpacingHoriz,SpacingVert
%           Padding,PaddingRight,PaddingLeft,PaddingTop,PaddingBottom
%           Margin,MarginRight,MarginLeft,MarginTop,MarginBottom
%           Holdaxis
%
%           all units are relative (i.e. from 0 to 1)
%
%           Abbreviations of parameters can be used.. (Eg MR instead of MarginRight)
%           (holdaxis means that it wont delete any axes below.)
%
%
% Example:
%
%   >> subaxis(2,1,1,'SpacingVert',0,'MR',0); 
%   >> imagesc(magic(3))
%   >> subaxis(2,'p',.02);
%   >> imagesc(magic(4))
%
% 2001-2014 / Aslak Grinsted  (Feel free to modify this code.)

f=gcf;

UserDataArgsOK=0;
Args=get(f,'UserData');
if isstruct(Args) 
    UserDataArgsOK=isfield(Args,'SpacingHorizontal')&isfield(Args,'Holdaxis')&isfield(Args,'rows')&isfield(Args,'cols');
end
OKToStoreArgs=isempty(Args)|UserDataArgsOK;

if isempty(Args)&&(~UserDataArgsOK)
    Args=struct('Holdaxis',0, ...
        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
        'rows',[],'cols',[]); 
end
Args=parseArgs(varargin,Args,{'Holdaxis'},{'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});

if (length(Args.NumericArguments)>2)
    Args.rows=Args.NumericArguments{1};
    Args.cols=Args.NumericArguments{2};
%remove these 2 numerical arguments
    Args.NumericArguments={Args.NumericArguments{3:end}};
end

if OKToStoreArgs
    set(f,'UserData',Args);
end


switch length(Args.NumericArguments)
   case 0
       return % no arguments but rows/cols.... 
   case 1
       if numel(Args.NumericArguments{1}) > 1 % restore subplot(m,n,[x y]) behaviour
           [x1 y1] = ind2sub([Args.cols Args.rows],Args.NumericArguments{1}(1)); % subplot and ind2sub count differently (column instead of row first) --> switch cols/rows
           [x2 y2] = ind2sub([Args.cols Args.rows],Args.NumericArguments{1}(end));
       else
           x1=mod((Args.NumericArguments{1}-1),Args.cols)+1; x2=x1;
           y1=floor((Args.NumericArguments{1}-1)/Args.cols)+1; y2=y1;
       end
%       x1=mod((Args.NumericArguments{1}-1),Args.cols)+1; x2=x1;
%       y1=floor((Args.NumericArguments{1}-1)/Args.cols)+1; y2=y1;
   case 2
      x1=Args.NumericArguments{1};x2=x1;
      y1=Args.NumericArguments{2};y2=y1;
   case 4
      x1=Args.NumericArguments{1};x2=x1+Args.NumericArguments{3}-1;
      y1=Args.NumericArguments{2};y2=y1+Args.NumericArguments{4}-1;
   otherwise
      error('subaxis argument error')
end
    

cellwidth=((1-Args.MarginLeft-Args.MarginRight)-(Args.cols-1)*Args.SpacingHorizontal)/Args.cols;
cellheight=((1-Args.MarginTop-Args.MarginBottom)-(Args.rows-1)*Args.SpacingVertical)/Args.rows;
xpos1=Args.MarginLeft+Args.PaddingLeft+cellwidth*(x1-1)+Args.SpacingHorizontal*(x1-1);
xpos2=Args.MarginLeft-Args.PaddingRight+cellwidth*x2+Args.SpacingHorizontal*(x2-1);
ypos1=Args.MarginTop+Args.PaddingTop+cellheight*(y1-1)+Args.SpacingVertical*(y1-1);
ypos2=Args.MarginTop-Args.PaddingBottom+cellheight*y2+Args.SpacingVertical*(y2-1);

if Args.Holdaxis
    h=axes('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
else
    h=subplot('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
end


set(h,'box','on');
%h=axes('position',[x1 1-y2 x2-x1 y2-y1]);
set(h,'units',get(gcf,'defaultaxesunits'));
set(h,'tag','subaxis');



if (nargout==0), clear h; end;
end
function ArgStruct=parseArgs(args,ArgStruct,varargin)
% Helper function for parsing varargin. 
%
%
% ArgStruct=parseArgs(varargin,ArgStruct[,FlagtypeParams[,Aliases]])
%
% * ArgStruct is the structure full of named arguments with default values.
% * Flagtype params is params that don't require a value. (the value will be set to 1 if it is present)
% * Aliases can be used to map one argument-name to several argstruct fields
%
%
% example usage: 
% --------------
% function parseargtest(varargin)
%
% %define the acceptable named arguments and assign default values
% Args=struct('Holdaxis',0, ...
%        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
%        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
%        'rows',[],'cols',[]); 
%
% %The capital letters define abrreviations.  
% %  Eg. parseargtest('spacingvertical',0) is equivalent to  parseargtest('sv',0) 
%
% Args=parseArgs(varargin,Args, ... % fill the arg-struct with values entered by the user
%           {'Holdaxis'}, ... %this argument has no value (flag-type)
%           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
%
% disp(Args)
%
%
%
%
% Aslak Grinsted 2004

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

persistent matlabver

if isempty(matlabver)
    matlabver=ver('MATLAB');
    matlabver=str2double(matlabver.Version);
end

Aliases={};
FlagTypeParams='';

if (length(varargin)>0) 
    FlagTypeParams=lower(strvcat(varargin{1}));  %#ok
    if length(varargin)>1
        Aliases=varargin{2};
    end
end
 

%---------------Get "numeric" arguments
NumArgCount=1;
while (NumArgCount<=size(args,2))&&(~ischar(args{NumArgCount}))
    NumArgCount=NumArgCount+1;
end
NumArgCount=NumArgCount-1;
if (NumArgCount>0)
    ArgStruct.NumericArguments={args{1:NumArgCount}};
else
    ArgStruct.NumericArguments={};
end 


%--------------Make an accepted fieldname matrix (case insensitive)
Fnames=fieldnames(ArgStruct);
for i=1:length(Fnames)
    name=lower(Fnames{i,1});
    Fnames{i,2}=name; %col2=lower
    Fnames{i,3}=[name(Fnames{i,1}~=name) ' ']; %col3=abreviation letters (those that are uppercase in the ArgStruct) e.g. SpacingHoriz->sh
    %the space prevents strvcat from removing empty lines
    Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams)); %Does this parameter have a value?
end
FnamesFull=strvcat(Fnames{:,2}); %#ok
FnamesAbbr=strvcat(Fnames{:,3}); %#ok

if length(Aliases)>0  
    for i=1:length(Aliases)
        name=lower(Aliases{i,1});
        FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(name,FnamesFull); %&??????? exact or not? 
        end
        Aliases{i,2}=FieldIdx;
        Aliases{i,3}=[name(Aliases{i,1}~=name) ' ']; %the space prevents strvcat from removing empty lines
        Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
    end
    %Append aliases to the end of FnamesFull and FnamesAbbr
    FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); %#ok
    FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3})); %#ok
end

%--------------get parameters--------------------
l=NumArgCount+1; 
while (l<=length(args))
    a=args{l};
    if ischar(a)
        paramHasValue=1; % assume that the parameter has is of type 'param',value
        a=lower(a);
        FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(a,FnamesFull); 
        end
        if (length(FieldIdx)>1) %shortest fieldname should win 
            [mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));%#ok
            FieldIdx=FieldIdx(mxi);
        end
        if FieldIdx>length(Fnames) %then it's an alias type.
            FieldIdx=Aliases{FieldIdx-length(Fnames),2}; 
        end
        
        if isempty(FieldIdx) 
            error(['Unknown named parameter: ' a])
        end
        for curField=FieldIdx' %if it is an alias it could be more than one.
            if (Fnames{curField,4})
                if (l+1>length(args))
                    error(['Expected a value for parameter: ' Fnames{curField,1}])
                end
                val=args{l+1};
            else %FLAG PARAMETER
                if (l<length(args)) %there might be a explicitly specified value for the flag
                    val=args{l+1};
                    if isnumeric(val)
                        if (numel(val)==1)
                            val=logical(val);
                        else
                            error(['Invalid value for flag-parameter: ' Fnames{curField,1}])
                        end
                    else
                        val=true;
                        paramHasValue=0; 
                    end
                else
                    val=true;
                    paramHasValue=0; 
                end
            end
            if matlabver>=6
                ArgStruct.(Fnames{curField,1})=val; %try the line below if you get an error here
            else
                ArgStruct=setfield(ArgStruct,Fnames{curField,1},val); %#ok <-works in old matlab versions
            end
        end
        l=l+1+paramHasValue; %if a wildcard matches more than one
    else
        error(['Expected a named parameter: ' num2str(a)])
    end
end
end
function vline(x, varargin)

% VLINE plot a vertical line in the current graph

abc = axis;
x = [x x];
y = abc([3 4]);
if length(varargin)==1
  varargin = {'color', varargin{1}};
end
h = line(x, y);
set(h, varargin{:});
end
function hhh=hline(y,in1,in2)
% function h=hline(y, linetype, label)
%
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% hline also supports vector inputs to draw multiple lines at once.  For example,
%
% hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
%
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(y)>1  % vector input
    for I=1:length(y)
        switch nargin
            case 1
                linetype='r:';
                label='';
            case 2
                if ~iscell(in1)
                    in1={in1};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                label='';
            case 3
                if ~iscell(in1)
                    in1={in1};
                end
                if ~iscell(in2)
                    in2={in2};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                if I>length(in2)
                    label=in2{end};
                else
                    label=in2{I};
                end
        end
        h(I)=hline(y(I),linetype,label);
    end
else
    switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            linetype=in1;
            label='';
        case 3
            linetype=in1;
            label=in2;
    end
    
    g=ishold(gca);
    hold on
    
    x=get(gca,'xlim');
    h=plot(x,[y y],linetype);
    if ~isempty(label)
        yy=get(gca,'ylim');
        yrange=yy(2)-yy(1);
        yunit=(y-yy(1))/yrange;
        if yunit<0.2
            text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
        else
            text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
        end
    end
    
    if g==0
        hold off
    end
    set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
end % else

if nargout
    hhh=h;
end

end