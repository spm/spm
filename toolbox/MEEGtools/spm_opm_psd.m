function [po,freq,sel] = spm_opm_psd(S)
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
%   S.plotbad       - place asterisk over unusual channels  - Default: 0
%	S.interact		- allow inspection of channels			- Default: 0
%	S.select		- enable selection of channels			- Default: 0

%
% Output:
%   psd             - power spectral density
%   f               - frequencies psd is sampled at
%	sel				- selected channel index (see S.interactive)
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
if ~isfield(S, 'interact'),      S.interact =1; end
if ~isfield(S, 'plotbad'),       S.plotbad = 0; end
if ~isfield(S, 'select')		
	S.select = 0; 
	sel = [];
else
	S.plot = 1;
	S.interact = 0;
end

%-channel Selection
%--------------------------------------------------------------------------

labs = S.channels;
regex = {};
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

Nf= ceil((N+1)/2);
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
med = repmat(median(po,2),1,size(po,2));
ma = repmat(median(abs(po-med),2)*1.48,1,size(po,2));
Z = (po-med)./ma;

bin = Z>abs(spm_invNcdf(.025/(size(po,2)*100),0,1));

if(S.plot)
    f= figure();
	if S.select
		h = plot(freq,po,'LineWidth',2);
		hold on
		for i = 1:numel(h)
			set(h(i), 'UserData', false);
		end
	else
		hold on
		for i = 1:size(po,2)
			tag = [labs{i}, ', Index: ' num2str(indchannel(S.D,labs{i}))];
			plot(freq,po(:,i)','LineWidth',2,'tag',tag);
		end
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
    if(S.plotbad)
      for i = 1:size(po,2)
        tag = [labs{i}, ', Index: ' num2str(indchannel(S.D,labs{i}))];
        plot(freq(bin(:,i)), po(bin(:,i),i),'k*','tag',tag);
      end
    end
    
    
	if(S.interact)
		datacursormode on
		dcm = datacursormode(gcf);
		set(dcm,'UpdateFcn',@getLabel)
	end
	if(S.select)
		sel = false(1,length(labs));
		numChannels = size(po, 2);
		
		% Set up selection mechanism
		set(fig, 'WindowButtonDownFcn', @startSelection);
		
		% Wait for user to interact and close the figure
		waitfor(fig);
	end
end

function startSelection(~, ~)
% Draw selection box
point1 = get(ax, 'CurrentPoint');
rbbox; 
point2 = get(ax, 'CurrentPoint');

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
			sel(chanIdx) = true;
			disp(['Selected ',labs{chanIdx}])
		end
	end
end
% Exclude selected channels from the mean calculation
set(p3, 'YData', median(po(:, ~sel), 2));
maxSelPo = max(max(po(:,~sel)));
minSelPo = min(min(po(:,~sel)));
set(gca,'YLim',[minSelPo,maxSelPo])
end

end

function txt = getLabel(trash,event)
pos = get(event,'Position');
dts = get(event.Target,'Tag');
txt = {dts,...
       ['Frequency: ',num2str(pos(1))],...
     ['RMS: ',num2str(pos(2))]};
end