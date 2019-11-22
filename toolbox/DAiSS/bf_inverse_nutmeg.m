function res = bf_inverse_nutmeg(BF, S)
% Interface to NUTMEG inverse methods 
% http://www.nitrc.org/plugins/mwiki/index.php/nutmeg:MainPage
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging


% $Id: bf_inverse_nutmeg.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    method = cfg_menu;
    method.tag = 'method';
    method.name = 'Method';
    method.labels = {
        'sLORETA'
        'swLORETA'
        'dSPM'};       
    method.values = method.labels;
    method.help = {'Select one of NUTMEG methods'};
    
    snr = cfg_entry;
    snr.tag = 'snr';
    snr.name = 'SNR';
    snr.strtype = 'e';
    snr.num = [0 0];
    snr.val = {[]};
    snr.help = {'The assumed ratio of variances of signal and noise,',...
        'used for setting the regularisation parameter.'};
    
    regularisation = cfg_entry;
    regularisation.tag = 'regularisation';
    regularisation.name = 'Regularisation parameter';
    regularisation.strtype = 'e';
    regularisation.num = [0 0];
    regularisation.val = {[]};
    regularisation.help = {'Optional regularization parameter.'};
    
    
    nutmeg      = cfg_branch;
    nutmeg.tag  = 'nutmeg';
    nutmeg.name = 'NUTMEG methods';
    nutmeg.val  = {method, snr, regularisation};
    nutmeg.help = {'Methods ported from NUTMEG toolbox',...
        'See http://www.nitrc.org/plugins/mwiki/index.php/nutmeg:MainPage'};
    
    res = nutmeg;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end


C     =  BF.features.(S.modality).C;
U     =  BF.features.(S.modality).U;

L = S.L;
W = cell(size(L));
nvert = numel(W);

data.Ryy = C;
flags    = [];
if ~isempty(S.snr)
    flags.snr = S.snr;
end
if ~isempty(S.regularisation)
    flags.gamma = S.regularisation;
end

LL = [];
spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, ['Preparing ' S.modality ' leadfields']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    lf = U'*L{i};       
    
    LL = cat(3, LL, lf);
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


switch S.method
    case 'sLORETA'
        w = nut_sLORETA(LL,data,flags);
    case 'swLORETA'
        w = nut_swLORETA(LL,data,flags);
    case  'dSPM'
        w = nut_dSPM(LL,data,flags);
end


spm_progress_bar('Init', nvert, ['Preparing ' S.modality ' filters']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    W{i} = spm_squeeze(w(:, :, i), 3)';
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


spm_progress_bar('Clear');

res.W = W;