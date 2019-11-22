function res = bf_inverse_champagne(BF, S)
% Computes Champagne filters
% See Owen et al. Performance evaluation of the Champagne source 
% reconstruction algorithm on simulated and real M/EEG data. Neuroimage. 2012 Mar;60(1):305-23
% Code contributed by Sri Nagarajan
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_inverse_champagne.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    nem = cfg_entry;
    nem.tag = 'nem';
    nem.name = 'Number of EM iterations';
    nem.strtype = 'n';
    nem.num = [1 1];
    nem.val = {100};
       
    vcs         = cfg_menu;
    vcs.tag     = 'vcs';
    vcs.name    = 'Voxel covariance structure';    
    vcs.labels  = {
        'scalar'
        'diagonal'
        'general'
        }';
    vcs.values  = {0, 1, 2};
    vcs.val = {2};
    
    nupd         = cfg_menu;
    nupd.tag     = 'nupd';
    nupd.name    = 'Noise covariance';
    nupd.labels  = {
        'use provided'
        'learn scalar'
        'learn diagonal'
        }';
    nupd.values  = {0, 1, 2};
    nupd.val = {0};
    
    champagne      = cfg_branch;
    champagne.tag  = 'champagne';
    champagne.name = 'Champagne';
    champagne.val  = {nem, vcs, nupd};
    res = champagne;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end


C     =  BF.features.(S.modality).C;
Y     =  BF.features.(S.modality).Y;

L = S.L;

W = cell(size(L));
nvert = numel(W);
nd = size(L{1}, 2);

LL = zeros(size(L{1}, 1), nvert*nd);
ind = 1;

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, ['Preparing ' S.modality ' leadfields']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert       
    
    for j = 1:nd
        lf = L{i}(:, j);
        LL(:, ind) = lf./norm(lf);
        ind = ind+1;
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

Fgraph  = spm_figure('GetWin', 'Champagne'); figure(Fgraph); clf

%****************************************************************************************
[gamma,x,w,sigu,like]=champagne_aug2015(Y, LL, C, S.nem, nd, S.vcs, S.nupd, [], 0, Fgraph);
%****************************************************************************************

spm_progress_bar('Init', nvert, ['Preparing ' S.modality ' filters']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    W{i} = w((i-1)*nd+(1:3), :);
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


spm_progress_bar('Clear');

spm('Pointer', 'Arrow');

res.W = W;