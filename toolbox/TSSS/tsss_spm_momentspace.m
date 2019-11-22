function D = tsss_spm_momentspace(S)
% Switch a dataset to SSS space using virtual montage
% FORMAT D = spm_eeg_crop(S)
%
% S        - input struct
%  fields of S:
%   D          - MEEG object or filename of M/EEG mat-file with data after
%                TSSS tool
%   condthresh - threshold on condition number for regularisation
%
% Output:
% D        - MEEG object (also written on disk)
%
% Reference: Vrba J, Taulu S, Nenonen J, Ahonen A. Signal space separation
% beamformer. Brain Topogr. 2010 Jun;23(2):128-33.
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: tsss_spm_momentspace.m 7703 2019-11-22 12:06:29Z guillaume $

SVNrev = '$Rev: 7703 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','TSSS momenspace transformation'); spm('Pointer','Watch');

if ~isfield(S, 'condthresh'),       S.condthresh   = 80;           end

D = spm_eeg_load(S.D);

if ~isfield(D, 'SSS')
    error('Run the TSSS tool first');
end

D = montage(D, 'clear');

SSS = D.SSS;

[SN_new, sss_indices , nmodes] = basis_condition_adjustment(SSS(1).SN, size(SSS.SNin, 2), S.condthresh);
pSN = pinv(SN_new);
pSN = pSN(1:nmodes, :);

labelnew = {};
for i = 1:nmodes
    labelnew{i} = ['moment' num2str(i)];
end

mont = [];
mont.labelorg = D.chanlabels(D.indchantype('MEGANY'))';
mont.labelnew = labelnew;
mont.tra = pSN;

if ~isempty(S.addchannels)
    chantypeorg = D.chantype(D.indchannel(S.addchannels));
    chanunitorg = D.units(D.indchannel(S.addchannels));
    for c = 1:numel(S.addchannels)
        mont.labelorg(end+1) = S.addchannels(c);
        mont.labelnew(end+1) = S.addchannels(c);
        mont.tra(end+1, end+1) = 1;
    end
end

D = montage(D, 'add', mont);
D = chantype(D, strmatch('moment', D.chanlabels), 'MEG');
D = units(D, strmatch('moment', D.chanlabels), 'fT');

if ~isempty(S.addchannels)
    D = chantype(D, D.indchannel(S.addchannels), chantypeorg);
    D = units(D, D.indchannel(S.addchannels), chanunitorg);
end

D = badchannels(D, D.indchantype('MEGANY'), 0);

save(D);

spm('FigName','TSS momenspace transformation: done'); spm('Pointer','Arrow');