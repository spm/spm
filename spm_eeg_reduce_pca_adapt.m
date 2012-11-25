function res = spm_eeg_reduce_pca_adapt(S)
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
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Mark Woolrich
% $Id: spm_eeg_reduce_pca_adapt.m 5079 2012-11-25 18:38:18Z vladimir $


if nargin == 0 
    ncomp = cfg_entry;
    ncomp.tag = 'ncomp';
    ncomp.name = 'Number of components';
    ncomp.strtype = 'n';
    ncomp.num = [1 1];
    ncomp.val = {1};
    ncomp.help = {'Number of PCA components'};        
    
    pca_adapt = cfg_branch;
    pca_adapt.tag = 'pca_adapt';
    pca_adapt.name = 'Adaptive PCA';
    pca_adapt.val = {ncomp};
    
    res = pca_adapt;
    
    return
end

D = S.D;

YY = 0;
ns = 0;

ntrials=length(S.trials); %MWW

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

for i = 1:ntrials %MWW
    for j = 1:numel(S.samples) %MWW
        Y  = D(S.chanind, S.samples{j}, S.trials(i)); %MWW
        Y  = detrend(Y', 'constant');
        YY = YY+(Y'*Y);
        ns = ns + D.nsamples-1;
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end;

spm_progress_bar('Clear');

C = YY/ns;

[U,dum] = svd(C);

%%%%%%%
% MWW 
ncomp=S.ncomp;
ncomp_adapt=spm_pca_order(C);
if(ncomp==-1 || ncomp>ncomp_adapt)
   ncomp = ncomp_adapt; 
end
%%%%%%%

% Assuming projecting to columns
montage = [];
montage.labelorg = D.chanlabels(S.chanind);
montage.tra = U(:, 1:ncomp)'; 

%%%%%%%
% MWW 
%Y2=D(S.chanind, S.samples{1}, S.trials);
%Y2=reshape(permute(Y2,[2 3 1]),size(Y2,2)*size(Y2,3),size(Y2,1));
%[allsvd,Apca]=pca(Y2,ncomp);
%pinvApca=pinv(Apca);
%montage.tra = pinvApca;
%%%%%%%

for i = 1:S.ncomp
    montage.labelnew{i, 1} = ['comp' num2str(i)];

    % MWW added:
    montage.chantypenew{i}='MEGPCACOMP';
end

res = montage;
