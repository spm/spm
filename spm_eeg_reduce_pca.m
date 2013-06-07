function res = spm_eeg_reduce_pca(S)
% Plugin for data reduction using PCA
% FORMAT res = spm_eeg_reduce_pca(S)
%
% S                     - input structure
% fields of S:
%    S.ncomp            - number of PCA components
%                            
% Output:
%  res - 
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided:
%      montage struct implementing projection to PCA subspace
%______________________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak 
% $Id: spm_eeg_reduce_pca.m 5528 2013-06-07 11:47:27Z vladimir $


if nargin == 0 
    ncomp = cfg_entry;
    ncomp.tag = 'ncomp';
    ncomp.name = 'Number of components';
    ncomp.strtype = 'n';
    ncomp.num = [1 1];
    ncomp.val = {1};
    ncomp.help = {'Number of PCA components'};        
    
    pca = cfg_branch;
    pca.tag = 'pca';
    pca.name = 'PCA data reduction';
    pca.val = {ncomp};
    
    res = pca;
    
    return
end

D = S.D;

YY = 0;
ns = 0;

ntrials = D.ntrials;

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

for i = 1:ntrials
    Y  = squeeze(D(S.chanind, :, i));
    Y  = detrend(Y', 'constant');
    YY = YY+(Y'*Y);
    ns = ns + D.nsamples-1;
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

C = YY/ns;

[U,lambda,V] = svd(C);

% Assuming projecting to columns
montage = [];
montage.labelorg = D.chanlabels(S.chanind);
montage.tra = U(:, 1:S.ncomp)';
for i = 1:S.ncomp
    montage.labelnew{i, 1} = ['comp' num2str(i)];
    montage.chantypenew{i}='MEGPCACOMP';
end

res = montage;
