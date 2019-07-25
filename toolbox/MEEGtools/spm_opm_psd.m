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
% Output:
%   psd             - power spectral density
%   f               - frequencies psd is sampled at
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_psd.m 7647 2019-07-25 13:58:53Z tim $

%-ArgCheck
%--------------------------------------------------------------------------

if ~isfield(S, 'units'),         S.units = 'fT'; end
if ~isfield(S, 'triallength'),   S.triallength = 1000; end
if ~isfield(S, 'constant'),      S.constant = 15; end
if ~isfield(S, 'bc'),            S.bc = 0; end
if ~isfield(S, 'channels'),      S.channels = 'ALL'; end
if ~isfield(S, 'plot'),          S.plot = 0; end
if ~isfield(S, 'D'),             error('D is required'); end


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
%-epoch dataset
%--------------------------------------------------------------------------
if size(S.D,3)>1
    eD=S.D;
else
    nsamps = S.triallength/1000*S.D.fsample;
    beg = 1:nsamps:size(S.D,2);
    endsamp =  beg+(nsamps-1);
    inRange = ~(beg>size(S.D,2)|endsamp>size(S.D,2));
    eD = zeros(length(chans),nsamps,sum(inRange));
    for i =1:length(inRange)
        if(inRange(i))
            eD(:,:,i)=S.D(chans,beg(i):endsamp(i),:);
        end
    end
    
end
%- set window
%--------------------------------------------------------------------------
fs = S.D.fsample();
N =size(eD,2);
Nf= ceil((N+1)/2);
nepochs=size(eD,3);
pow = zeros(Nf,size(eD,1),nepochs);
wind  = window(@flattopwin,size(eD,2));
wind = repmat(wind,1,size(eD,1));

%- create PSD
%--------------------------------------------------------------------------

for j = 1:nepochs
    Btemp=eD(:,:,j)';
    Btemp = Btemp.*wind;
    mu=mean(Btemp);
    zf = bsxfun(@minus,Btemp,mu);
    if(S.bc)
        fzf = zf;
    else
        fzf=Btemp;
    end
    
    N= length(fzf);
    xdft = fft(fzf);
    xdft=xdft(1:floor(N/2+1),:);
    psdx = abs(xdft)./sqrt(N*fs);
    freq = 0:fs/size(fzf,1):fs/2;
    odd=mod(size(fzf,1),2)==1;
    if(odd)
        psdx(2:end) = sqrt(2)*psdx(2:end);
    else
        psdx(2:end-1) = sqrt(2)*psdx(2:end-1);
    end
    pow(:,:,j) =psdx;
end

%- plot
%--------------------------------------------------------------------------
po = median(pow(:,:,:),3);
if(S.plot)
    figure()
    semilogy(freq,po,'LineWidth',2);
    hold on
    xp2 =0:round(freq(end));
    yp2=ones(1,round(freq(end))+1)*S.constant;
    p2 =plot(xp2,yp2,'--k');
    p2.LineWidth=2;
    xlabel('Frequency (Hz)')
    labY = ['$$PSD (' S.units ' \sqrt[-1]{Hz}$$)'];
    ylabel(labY,'interpreter','latex')
    grid on
    ax = gca; % current axes
    ax.FontSize = 13;
    ax.TickLength = [0.02 0.02];
    fig= gcf;
    fig.Color=[1,1,1];

end

end
