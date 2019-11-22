function res = bf_inverse_minimumnorm(BF, S)
% Computes Minimum Norm projectors
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, using the code from Matti Stenroos and Olaf Hauk
% http://imaging.mrc-cbu.cam.ac.uk/meg/AnalyzingData/DeFleCT_SpatialFiltering_Tools
%
% Please cite:
% Hauk O, Stenroos M.
% A framework for the design of flexible cross-talk functions for spatial filtering of EEG/MEG data: DeFleCT.
% Human Brain Mapping 2013
% $Id: bf_inverse_minimumnorm.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    snr = cfg_entry;
    snr.tag = 'snr';
    snr.name = 'SNR';
    snr.strtype = 'r';
    snr.num = [1 1];
    snr.val = {5};
    snr.help = {'The assumed ratio of variances of signal and noise,',...
        'used for setting the regularisation parameter.'};
    
    trunc = cfg_entry;
    trunc.tag = 'trunc';
    trunc.name = 'Truncation parameter';
    trunc.strtype = 'w';
    trunc.num = [1 1];
    trunc.val = {0};
    trunc.help = {'The number of (smallest) singular values of the covariance matrix that are set to ',...
        'zero before making the whitener. For example, if the data has been SSP-projected, it needs to be at least the number of ',...
        'components projected away.'};
    
    minimumnorm      = cfg_branch;
    minimumnorm.tag  = 'minimumnorm';
    minimumnorm.name = 'Minimum norm';
    minimumnorm.val  = {snr, trunc};
    minimumnorm.help = {'Minimum norm implementation by Matti Stenroos and Olaf Hauk',...
        'http://imaging.mrc-cbu.cam.ac.uk/meg/AnalyzingData/DeFleCT_SpatialFiltering_Tools',...
        'Please cite:'...
        'Hauk O, Stenroos M.',...
        'A framework for the design of flexible cross-talk functions for spatial filtering of EEG/MEG data: DeFleCT.',...
        'Human Brain Mapping 2013'};
    res = minimumnorm;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end


C     =  BF.features.(S.modality).C;
U     =  BF.features.(S.modality).U;

L = S.L;
W = cell(size(L));

nvert = numel(W);

nori  = size(L{1}, 2);

LL = [];

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, ['Preparing ' S.modality ' leadfields']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    lf = U'*L{i};     
    
    LL = cat(2, LL, lf);
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

WW = MNestimator(LL,C, S.snr,S.trunc);

spm_progress_bar('Init', nvert, ['Preparing ' S.modality ' filters']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    W{i} = WW(((i-1)*nori+1):(i*nori), :); 
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


spm_progress_bar('Clear');

res.W = W;
