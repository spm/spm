function [po,freq] = spm_opm_psd(S)
% Compute PSD for OPM data(for checking noise floor)
% FORMAT D = spm_opm_psd(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                       - Default: no Default
%   S.triallength   - window size (ms)                      - Default: 1000
%   S.bc            - boolean to dc correct                 - Default: 0
%   S.channels      - channels to estimate PSD from         - Default: 'ALL'
%   S.plot          - boolean to plot or not                - Default: 0
%   S.units         - units of measurement                  - Default: 'fT'
%   S.constant      - constant line to draw as reference    - Default: 15
%   S.wind          - function handle for window            - Default: @hanning
%
% Output:
%   psd             - power spectral density
%   f               - frequencies psd is sampled at
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
if ~isfield(S, 'interact'),      S.interact=1; end


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

if(S.plot)
    f= figure();
    hold on
    for i = 1:size(po,2)
        tag = [labs{i}, ', Index: ' num2str(indchannel(S.D,labs{i}))];
        plot(freq,po(:,i)','LineWidth',2,'tag',tag);
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
    if(S.interact)
        datacursormode on
        dcm = datacursormode(gcf);
        set(dcm,'UpdateFcn',@getLabel)
    end
end

end

function txt = getLabel(trash,event)
pos = get(event,'Position');
dts = get(event.Target,'Tag');
txt = {dts,...
       ['Frequency: ',num2str(pos(1))],...
     ['RMS: ',num2str(pos(2))]};
end