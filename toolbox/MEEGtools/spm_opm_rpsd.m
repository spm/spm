function [shield,f] = spm_opm_rpsd(S)
% Compute relative PSD of two OPM datasets (for checking shielding factors)
% FORMAT D = spm_opm_rpsd(S)
%   S               - input structure
%  fields of S:
%   S.D1            - SPM MEEG object                      - Default: no Default
%   S.D2            - SPM MEEG object                      - Default: no Default
%   S.triallength   - window size (ms)                      - Default: 1000
%   S.bc            - boolean to dc correct                 - Default: 0
%   S.channels      - channels to estimate PSD from         - Default: 'ALL'
%   S.dB            - boolean to return decibels            - Default: 0
%   S.plot          - boolean to plot or not                - Default: 0
% Output:
%   sf              - Shielding factor ( in data units or decibels)
%   f               - frequencies psd is sampled at
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_rpsd.m 7782 2020-02-11 11:33:49Z tim $

%-ArgCheck
%--------------------------------------------------------------------------

if ~isfield(S, 'units'),         S.units = 'fT'; end
if ~isfield(S, 'triallength'),   S.triallength = 1000; end
if ~isfield(S, 'constant'),      S.constant = 15; end
if ~isfield(S, 'bc'),            S.bc = 0; end
if ~isfield(S, 'channels'),      S.channels = 'ALL'; end
if ~isfield(S, 'plot'),          S.plot = 0; end
if ~isfield(S, 'dB'),            S.dB = 0; end
if ~isfield(S, 'D1'),            error('D1 is required'); end
if ~isfield(S, 'D2'),            error('D2 is required'); end

%- Checks
%--------------------------------------------------------------------------
fsEqual = S.D1.fsample==S.D2.fsample;
if (~fsEqual),error('Sample rates should be identical');end

%-First PSD
%--------------------------------------------------------------------------
args=[];
args.D = S.D1;
args.triallength=S.triallength;
args.bc=S.bc;
args.channels=S.channels;
args.plot=S.plot;
args.trials=1;
[p1,~]=spm_opm_psd(args);

%- Second PSD
%--------------------------------------------------------------------------
args=[];
args.D = S.D2;
args.triallength=S.triallength;
args.bc=S.bc;
args.channels=S.channels;
args.plot=S.plot;
args.trials=1;
[p2,f]=spm_opm_psd(args);

%- Get same number of trials for each PSD
%--------------------------------------------------------------------------
ntrials = min([size(p1,3),size(p2,3)]);
p1 = p1(:,:,1:ntrials);
p2 = p2(:,:,1:ntrials);
p1 = median(p1,3);
p2 = median(p2,3);

%- Ratio (or difference) of corresponding sensor in each dataset
%--------------------------------------------------------------------------
chan1 = chanlabels(S.D1);
chan2 = chanlabels(S.D2);
inCommon=zeros(size(p1));
keep = zeros(size(p1,2),1);

for i =1:length(chan1)
    index= strcmp(chan1{i},chan2);
    if (any(index))
        keep(i)=1;
        if(S.dB)
            inCommon(:,i)=20*log10(p1(:,i)./p2(:,index));
        else
            inCommon(:,i)=p1(:,i)-p2(:,index);
        end
    end
end
shield = inCommon(:,boolean(keep));
%- Plot
%--------------------------------------------------------------------------
if(S.plot)
    if(S.dB)
        figure()
        plot(f,shield,'LineWidth',2);
        hold on
        xp2 =0:round(f(end));
        yp2=ones(1,round(f(end))+1)*0;
        p2 =plot(xp2,yp2,'--k');
        p2.LineWidth=2;
        xlabel('Frequency (Hz)')
        labY = ['Shielding Factor (dB)'];
        ylabel(labY)
        grid on
        ax = gca; % current axes
        ax.FontSize = 13;
        ax.TickLength = [0.02 0.02];
        fig= gcf;
        fig.Color=[1,1,1];
        legend(chan1{boolean(keep)});
    else
        figure()
        plot(f,shield,'LineWidth',2);
        hold on
        xp2 =0:round(f(end));
        yp2=ones(1,round(f(end))+1)*0;
        p2 =plot(xp2,yp2,'--k');
        p2.LineWidth=2;
        xlabel('Frequency (Hz)')
        labY = ['$$PSD (fT' ' \sqrt[-1]{Hz}$$)'];
        ylabel(labY,'interpreter','latex')
        grid on
        ax = gca; % current axes
        ax.FontSize = 13;
        ax.TickLength = [0.02 0.02];
        fig= gcf;
        fig.Color=[1,1,1];
        legend(chan1{boolean(keep)});
    end
end

end
