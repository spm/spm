function res = bf_inverse_eloreta(BF, S)
% Computes eLORETA projectors
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, using the code from Guido Nolte
%
% please cite
% R.D. Pascual-Marqui: Discrete, 3D distributed, linear imaging methods of electric neuronal activity. Part 1: exact, zero
% error localization. arXiv:0710.3341 [math-ph], 2007-October-17, http://arxiv.org/pdf/0710.3341

% $Id: bf_inverse_eloreta.m 7706 2019-11-22 16:30:29Z spm $

%--------------------------------------------------------------------------
if nargin == 0
    regularisation = cfg_entry;
    regularisation.tag = 'regularisation';
    regularisation.name = 'Regularisation parameter';
    regularisation.strtype = 'r';
    regularisation.num = [1 1];
    regularisation.val = {0.05};
    regularisation.help = {'Optional regularization parameter (default is .05 corresponding ',...
        'to 5% of the average of the eigenvalues of some matrix to be inverted.)'};
    
    
    eloreta      = cfg_branch;
    eloreta.tag  = 'eloreta';
    eloreta.name = 'eLORETA';
    eloreta.val  = {regularisation};
    eloreta.help = {'eLORETA implementation by Guido Nolte',...
        'Please cite',...
        'R.D. Pascual-Marqui: Discrete, 3D distributed, linear imaging methods of electric neuronal activity. Part 1: exact, zero',...
        'error localization. arXiv:0710.3341 [math-ph], 2007-October-17, http://arxiv.org/pdf/0710.3341'};
    
    res = eloreta;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

U     =  BF.features.(S.modality).U;

L = S.L;
W = cell(size(L));
nvert = numel(W);

LL = [];

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, ['Preparing ' S.modality ' leadfields']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    lf = U'*L{i};
    
    lf = reshape(lf, [size(lf, 1), 1, size(lf, 2)]);
    
    LL = cat(2, LL, lf);
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


% construct the spatial filter
w = mkfilt_eloreta_v2(LL,S.regularisation);

spm_progress_bar('Init', nvert, ['Preparing ' S.modality ' filters']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    W{i} = spm_squeeze(w(:, i, :), 2)';
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


spm_progress_bar('Clear');

res.W = W;
