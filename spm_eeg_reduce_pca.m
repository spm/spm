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
% $Id: spm_eeg_reduce_pca.m 6813 2016-06-18 18:41:31Z vladimir $


if nargin == 0
    ncomp = cfg_entry;
    ncomp.tag = 'ncomp';
    ncomp.name = 'Number of components';
    ncomp.strtype = 'n';
    ncomp.num = [1 1];
    ncomp.val = {1};
    ncomp.help = {'Number of PCA components'};
    
    threshold = cfg_entry;
    threshold.tag = 'threshold';
    threshold.name = 'Variance threshold';
    threshold.strtype = 'r';
    threshold.num = [1 1];
    threshold.val = {0};
    threshold.help = {'Threshold (1 > u > 0) for normalized eigenvalues'};
    
    pca = cfg_branch;
    pca.tag = 'pca';
    pca.name = 'PCA data reduction';
    pca.val = {ncomp, threshold};
    
    res = pca;
    
    return
end

YY = 0;
ns = 0;

for f = 1:numel(S.D)
    
    if iscell(S.D)
        D = S.D{f};
    else
        D = S.D;
    end
    
    isTF = strncmp(D.transformtype, 'TF', 2);
    
    
    iscont = isequal(D.type, 'continuous');
    
    if iscont
        np      = round(D.fsample);
        ntrials = ceil(D.nsamples/np);
    else
        ntrials = D.ntrials;
    end
    
    spm('Pointer', 'Watch');drawnow;
    spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
    if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
    else Ibar = 1:ntrials; end
    
    for i = 1:ntrials
        if iscont
            if i<ntrials
                if isTF
                    Y  = squeeze(D(S.chanind, :, ((i-1)*np+1):i*np, 1));
                    Y  = reshape(Y, size(Y, 1), []);
                else
                    Y  = squeeze(D(S.chanind, ((i-1)*np+1):i*np, 1));
                end
            else
                if isTF
                    Y  = squeeze(D(S.chanind, :, ((i-1)*np+1):end, 1));
                    Y  = reshape(Y, size(Y, 1), []);
                else
                    Y  = squeeze(D(S.chanind, ((i-1)*np+1):end, 1));
                end
            end
        else
            if isTF
                Y  = squeeze(D(S.chanind, :, :, i));
                Y  = reshape(Y, size(Y, 1), []);
            else
                Y  = squeeze(D(S.chanind, :, i));
            end
        end
        Y  = detrend(Y', 'constant');
        YY = YY+(Y'*Y);
        ns = ns + size(Y, 2);
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end
    
    spm_progress_bar('Clear');    
end

C = YY/ns;

[U,lambda,V] = spm_svd(C, S.threshold);

ncomp = min(S.ncomp, size(U, 2));

% Assuming projecting to columns
montage = [];
montage.labelorg = D.chanlabels(S.chanind);
montage.tra = U(:, 1:ncomp)';
for i = 1:ncomp
    montage.labelnew{i, 1} = ['comp' num2str(i)];
    montage.chantypenew{i}= 'PHYS';
end

res = montage;
