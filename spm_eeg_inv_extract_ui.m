function [D] = spm_eeg_inv_extract_ui(varargin)
% GUI for contrast of evoked responses and power for an MEG-EEG model
% FORMAT [D] = spm_eeg_inv_extract_ui(D, val, XYZ)
% Sets:
%
%     D.contrast.woi   - time (ms) window of interest
%     D.contrast.fboi  - freq (Hz) window of interest
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_extract_ui.m 5640 2013-09-18 12:02:29Z vladimir $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{1:2});

D.inv{val}.source = [];

if nargin<3
    warndlg('Position in MNI coordinates should be provided');
    return
elseif ~isequal(size(varargin{3}), [1, 3])
    warndlg('A single point in MNI coordinates should be provided');
    return
end

D.inv{val}.source.XYZ = varargin{3};


% get source label
%--------------------------------------------------------------------------
D.inv{val}.source.label{1} = spm_input('Source label','+1','s');

% get VOI radius
%--------------------------------------------------------------------------
D.inv{val}.source.rad = spm_input('VOI radius','+1','r', 5, [1 1]);

% trials or or evoked?
%--------------------------------------------------------------------------
if strcmp(D.type, 'single')
    str  = {'evoked','trials'};
    type = spm_input('What to save?','+1','b',str,[],1);
else
    type = 'evoked';
end
D.inv{val}.source.type = type;

% output file name
%--------------------------------------------------------------------------
fname = spm_input('File name','+1','s', ['i' spm_file(D.fname,'basename')]);
D.inv{val}.source.fname = fullfile(D.path, [fname '.mat']);

% extract
%--------------------------------------------------------------------------
[Ds, D] = spm_eeg_inv_extract(D);

% display
%==========================================================================
if ~isempty(Ds)
    
    try, con = D.con; catch, con = 1; end
    
    if con == 0
        con = 1;
    end
    
    Fgraph = spm_figure('GetWin','Graphics');
    clf(Fgraph)
    figure(Fgraph)
    
    cl = Ds.condlist;
    
    for i = 1:numel(cl)
        if i == con
            Color = [1 0 0];
        else
            Color = [1 1 1]*(1 - 1/4);
        end
        
        plot(time(Ds, [], 'ms'), squeeze(Ds(1, :, Ds.indtrial(cl{i}, 'GOOD'))), 'Color', Color);                
        
        hold on
    end
    
    axis square
    
    title({sprintf('Extracted data from source %s', Ds.chanlabels{1}),...
        sprintf('at %d %d %d mm', round(D.inv{val}.source.XYZ))});
    
    xlabel('PST {ms}')
    ylabel('source amplitude')
    drawnow
end