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
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney

%-ArgCheck
%--------------------------------------------------------------------------
if ~isfield(S, 'triallength'),   S.triallength = 1000;                  end
if ~isfield(S, 'bc'),            S.bc = 0;                              end
if ~isfield(S, 'channels'),      S.channels = 'ALL';                    end
if ~isfield(S, 'plot'),          S.plot = 0;                            end
if ~isfield(S, 'dB'),            S.dB = 1;                              end
if ~isfield(S, 'D1'),            error('D1 is required');               end
if ~isfield(S, 'D2'),            error('D2 is required');               end
if ~isfield(S, 'interact'),      S.interact = 1;                        end

%- Checks
%--------------------------------------------------------------------------
fsEqual = S.D1.fsample==S.D2.fsample;
if (~fsEqual),error('Sample rates should be identical');end

%-First PSD
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Processing First PSD ',spm('time'));

args=[];
args.D = S.D1;
args.triallength=S.triallength;
args.bc=S.bc;
args.channels=S.channels;
[p1,~]=spm_opm_psd(args);


fprintf('%-40s: %30s\n','Processing Second PSD ',spm('time'));

%- Second PSD
%--------------------------------------------------------------------------
args=[];
args.D = S.D2;
args.triallength=S.triallength;
args.bc=S.bc;
args.channels=S.channels;
[p2,f]=spm_opm_psd(args);

%- Ratio (or difference) of corresponding sensor in each dataset
%--------------------------------------------------------------------------

if(S.dB)
    shield=20*log10(p1./p2);
else
    shield=p1-p2;
end

%- Plot
%--------------------------------------------------------------------------
if(S.plot)
    % get labels (in case S.channels = 'ALL')
    cid = S.D1.selectchannels(S.channels);
    labs = S.D1.chanlabels(cid);
    figure()
    hold on
    for i = 1:size(shield,2)
        tag = [labs{i}, ', Index: ' num2str(cid(i))];
        plot(f,shield(:,i)','LineWidth',2,'tag',tag);
    end
    % plot(f,shield,'LineWidth',2);
    xp2 =0:round(f(end));
    yp2=ones(1,round(f(end))+1)*0;
    p2 =plot(xp2,yp2,'--k');
    p2.LineWidth=2;
    if(S.dB)
        xlabel('Frequency (Hz)')
        labY = ['Shielding Factor (dB)'];
        ylabel(labY)
        % legend(chan1{boolean(keep)});
    else
        xlabel('Frequency (Hz)')
        labY = ['$$PSD (fT' ' \sqrt[-1]{Hz}$$)'];
        ylabel(labY,'interpreter','latex')
        %legend(chan1{boolean(keep)});
    end
    grid on
    ax = gca; % current axes
    ax.FontSize = 13;
    ax.TickLength = [0.02 0.02];
    fig= gcf;
    fig.Color=[1,1,1];

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
