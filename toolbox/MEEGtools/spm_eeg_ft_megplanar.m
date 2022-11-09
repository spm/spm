function D = spm_eeg_ft_megplanar(S)
% Function for transforming MEG data to planar gradient
%
% FORMAT  D = spm_eeg_ft_megplanar(S)
%
% S           - input structure (optional)
% (optional) fields of S:
%   S.D       - filename, or M/EEG object
%   S.prefix  - prefix (default L)
%
% Output
%   D - dataset converted to planar gradient
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename);
spm('FigName','Fieldtrip megplanar'); 

%-Load SPM EEG file and verify consistency
%--------------------------------------------------------------------------
if nargin == 0
    S = [];
end

if ~isfield(S, 'prefix'),       S.prefix   = 'L';           end

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
    S.D = D;
end

if ischar(D)
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

%-Configure the spectral analysis
%--------------------------------------------------------------------------

cfg = [];
cfg.neighbours = ft_prepare_neighbours(struct('method', 'template', 'template', 'CTF275_neighb.mat'), D.sensors('MEG'));
cfg.channel = 'all';

pl = ft_megplanar(cfg,  D.ftraw);

pl = rmfield(pl, 'hdr');
pl.grad.balance = [];

Dout = spm_eeg_ft2spm(pl, spm_file(fullfile(D), 'prefix', S.prefix));

%-Copy some additional information from the original file
%--------------------------------------------------------------------------
Dout  = conditions (Dout, ':', D.conditions);

Dout  = chantype(Dout, indchannel(Dout, setdiff(Dout.chanlabels, D.chanlabels)), 'MEGPLANAR');

[sel1, sel2] = spm_match_str(Dout.chanlabels, D.chanlabels);

Dout = chantype(Dout, sel1, chantype(D, sel2));
Dout = badchannels(Dout, sel1, badchannels(D, sel2));
Dout = coor2D(Dout, sel1, coor2D(D, sel2));

Dout = badtrials(Dout, badtrials(D), 1);
Dout = events(Dout, ':', D.events(':'));
Dout = history(Dout, history(D));


%-Update history
%--------------------------------------------------------------------------
Dout = history(Dout, mfilename, S);

%-Save
%--------------------------------------------------------------------------
save(Dout);

D = Dout;
