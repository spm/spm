function [po,freq,indices] = spm_opm_psd(S)
% Compute PSD for OPM data (for checking noise floor)
% FORMAT [po,freq,sel] = spm_opm_psd(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                       - Default: no Default
%   S.triallength   - window size (ms)                      - Default: 1000
%   S.bc            - boolean to dc correct                 - Default: 0
%   S.channels      - channel  to estimate PSD from         - Default: 'ALL'
%   S.plot          - boolean to plot or not                - Default: 0
%   S.units         - units of measurement                  - Default: 'fT'
%   S.constant      - constant line to draw as reference    - Default: 15
%   S.wind          - function handle for window            - Default: @hanning
%	S.selectbad		- highlights and enables selection of 
%                     bad channels			                - Default: 0
%
%
% Output:
%   psd             - power spectral density
%   f               - frequencies psd is sampled at
%	indices			- selected channel index
%						To get labels use:
%							plotted_lab = chanlabels(S.D,S.channels);
%							sel_lab = plotted_lab(sel);
%__________________________________________________________________________

% Tim Tierney
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


%-ArgCheck
%--------------------------------------------------------------------------

if ~isfield(S, 'units'),         S.units = 'fT'; end
if ~isfield(S, 'triallength'),   S.triallength = 1000; end
if ~isfield(S, 'constant'),      S.constant = 15; end
if ~isfield(S, 'bc'),            S.bc = 0; end
if ~isfield(S, 'channels'),      S.channels = 'ALL'; end
if ~isfield(S, 'plot'),          S.plot = 0; end
if ~isfield(S, 'D'),             error('D is required'); end
if ~isfield(S, 'trials'),        S.trials=0; end
if ~isfield(S, 'wind'),          S.wind=@hanning; end
if ~isfield(S, 'selectbad'),     S.selectbad = 0; end

%-channel Selection
%--------------------------------------------------------------------------
indices = [];
labs = S.channels;
regex = cell(1,length(labs));
for i = 1:length(labs)
    if isa(labs,'cell')
        regex{i} = ['regexp_(',labs{i},')'];
    else
        regex{i} = ['regexp_(',labs,')'];
    end
end
chans = [S.D.selectchannels(regex), indchantype(S.D,labs)];
labs = chanlabels(S.D,chans);
%- set window
%--------------------------------------------------------------------------
fs = S.D.fsample();

if (size(S.D,3)>1)
    N=size(S.D,2);
    nepochs = size(S.D,3);
else
    N = round(S.triallength/1000*S.D.fsample);
    begs = 1:N:size(S.D,2);
    ends= begs+N-1;
    
    if(ends(end)>size(S.D,2))
        ends(end)=[];
        begs(end)=[];
    end
    nepochs=length(begs);
end

try
  wind  = window(S.wind,N);
catch
  % catch for environments without signal toolbox
  wind = ones(N,1);
end

% Nf= ceil((N+1)/2);
coFac= max(wind)/mean(wind);
wind = repmat(wind,1,length(chans));

%- create PSD
%--------------------------------------------------------------------------
freq = 0:fs/N:fs/2;
odd=mod(N,2)==1;
psdx= zeros(length(freq),length(chans));

for j = 1:nepochs
    
    % if data is already epoched then extract epochs
    if (size(S.D,3)>1)
        Btemp=S.D(chans,:,j)';
    else % otherwise use triallength
        inds = begs(j):ends(j);
        Btemp=S.D(chans,inds,1)';
    end
    
    % window data and correct for amplitude loss;
    Btemp = Btemp.*wind*coFac;
    
    % baseline correct if required 
    if(S.bc)
        mu=median(Btemp);
        zf = bsxfun(@minus,Btemp,mu);
        fzf = zf;
    else
        fzf=Btemp;
    end
    
    % fourier transform data and get RMS
    xdft = fft(fzf);
    xdft = xdft(1:floor(N/2+1),:);
    tmppsd = abs(xdft)./sqrt(N*fs);
    
    
    if(odd)
        tmppsd(2:end) = tmppsd(2:end);
    else
        tmppsd(2:end-1) = tmppsd(2:end-1);
    end
    
    % accumulate avearge PSD to limit memory usage
    psdx=psdx + tmppsd/nepochs;
end

%- plot
%--------------------------------------------------------------------------

po = psdx;
% Identify outliers
med = repmat(median(po,2),1,size(po,2));
ma = repmat(median(abs(po-med),2)*1.48,1,size(po,2));
Z = (po-med)./ma;
bin = Z>abs(spm_invNcdf(.025/(size(po,2)*100),0,1));

if(S.plot)
    figure();

	h = plot(freq,po,'LineWidth',2);
	hold on
	for i = 1:numel(h)
		tag = [labs{i}, ', Index: ' num2str(indchannel(S.D,labs{i}))];
		set(h(i), 'UserData', false, 'tag',tag);
	end

    set(gca,'yscale','log')
    
    xp2 =0:round(freq(end));
    yp2=ones(1,round(freq(end))+1)*S.constant;
    p2 =plot(xp2,yp2,'--k');
    p2.LineWidth=2;
    p3=semilogy(freq,median(po,2),'LineWidth',2,'tag','Median');
    p3.Color='k';
    xlabel('Frequency (Hz)')
    labY = ['$$PSD (' S.units ' \sqrt[-1]{Hz}$$)'];
    ylabel(labY,'interpreter','latex')
    
    grid on
    ax = gca; % current axes
    ax.FontSize = 13;
    ax.TickLength = [0.02 0.02];
    fig= gcf;
    fig.Color=[1,1,1];
    xlim([0,100]);

	% Highlight bad channel/frequency combinations
	if(S.selectbad)
    indices =[];
		g = gobjects(1, numel(labs));
		hold on
		for i = 1:numel(labs)
    		x = freq(bin(:,i));
    		y = po(bin(:,i),i);
		
    		if ~isempty(x) && ~isempty(y)
        		tag = [labs{i}, ', Index: ' num2str(indchannel(S.D, labs{i}))];
        		g(i) = plot(x, y, 'k*', 'Tag', tag);
        		set(g(i), 'UserData', false);
    		end
		end
		sel = false(1,length(labs));
		numChannels = size(po, 2);
		anno = [];
		set(fig, 'WindowButtonDownFcn', @startSelection);
		
		% Wait for user to interact and close the figure
		waitfor(fig);
	else
		datacursormode on
		dcm = datacursormode(gcf);
 		set(dcm,'UpdateFcn',@getLabel)
	end
end

function startSelection(~, ~)
point1 = get(ax, 'CurrentPoint');
rbbox; 
point2 = get(ax, 'CurrentPoint');

% Determine if user clicked or clicked and dragged
if norm(point1(1,1:2) - point2(1,1:2)) < 2
	xClick = point1(1,1);
    yClick = point1(1,2);
    clickedObj = hittest(fig);
	
	% Make a text box (similar to datacursor)
	if isgraphics(clickedObj, 'line')
		if ~isempty(anno)
			set(anno,'Visible','off')
		end
		
		tagText = get(clickedObj, 'Tag');
        txt = {tagText;
            ['Frequency: ', num2str(xClick)];
            ['RMS: ', num2str(yClick)];
        };
 
		% Position on ylog axis
        axPos = get(ax, 'Position');
        xlims = get(ax, 'XLim');
        ylims = get(ax, 'YLim');
        yNorm = (log10(yClick) - log10(ylims(1))) / (log10(ylims(2)) - log10(ylims(1)));
		yNorm = axPos(2) + yNorm * axPos(4);
        xNorm = axPos(1) + (xClick - xlims(1)) / diff(xlims) * axPos(3);

        anno = annotation(fig, 'textbox', [xNorm, yNorm, 0.1, 0], ...
                    'String', txt, ...
                    'FitBoxToText', 'on', ...
                    'BackgroundColor', 'w', ...
                    'EdgeColor', [0.5 0.5 0.5], ...
                    'Margin', 2);
	end
else
	if ~isempty(anno)
		set(anno,'Visible','off')
	end
	% Get the x and y limits of the selection box
	xLimits = [min(point1(1,1), point2(1,1)), max(point1(1,1), point2(1,1))];
	yLimits = [min(point1(1,2), point2(1,2)), max(point1(1,2), point2(1,2))];
	
	% Check which lines are inside the box
	for chanIdx = 1:numChannels
		xData = get(h(chanIdx), 'XData');
		yData = get(h(chanIdx), 'YData');
		
		% Find any points within the x and y limits of the selection box
		inX = (xData >= xLimits(1)) & (xData <= xLimits(2));
		inY = (yData >= yLimits(1)) & (yData <= yLimits(2));
	
		% If any point of the line is within the selection box, toggle selection
		if any(inX & inY)
			if get(h(chanIdx), 'UserData') == false
				set(h(chanIdx), 'UserData', true); 
				set(h(chanIdx), 'Color', [h(chanIdx).Color(1:3), 0]);
				if isgraphics(g(chanIdx))
					set(g(chanIdx), 'Visible', 'off');
				end
				sel(chanIdx) = true;
				ta = get (h(chanIdx),'tag');
        ind = str2double(ta((strfind(ta,':')+1):end));
      	indices = [indices, ind];
			end
		end
	end
	
	% Exclude selected channels from the median calculation
	set(p3, 'YData', median(po(:, ~sel), 2));
	maxSelPo = max(max(po(:,~sel)));
	minSelPo = min(min(po(:,~sel)));
	try
		set(gca,'YLim',[minSelPo,maxSelPo])
	catch
	end
end
end

function txt = getLabel(~,event)
pos = get(event,'Position');
dts = get(event.Target,'Tag');
txt = {dts,...
		['Frequency: ',num2str(pos(1))],...
		['RMS: ',num2str(pos(2))]};
end
end


