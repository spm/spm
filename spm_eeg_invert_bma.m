function [Jbma,qCbma,PostMax]=spm_eeg_invert_bma(manyinverse,F)
% Compute Bayesian Model Average given a set of current distributions and model evidences
% FORMAT [Jbma,qCbma,PostMax]=spm_eeg_invert_bma(manyinverse,F)
%
% At the moment makes an estimate of posterior covariance based on relative
% weighting of the input posteriors but this could be changed in future.
% At the moment adds a random DC offset (and not a random time series) to
% the estimated current distribution at each vertex.
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% Jose David Lopez and Gareth Barnes
% $Id: spm_eeg_invert_bma.m 6125 2014-07-30 14:07:04Z guillaume $


%-1. Load F with the free energies and take the probabilities
%--------------------------------------------------------------------------
N = numel(manyinverse);
for j=1:N
    F(j) = manyinverse{j}.F;
end

Fx = F - mean(F);
Fx = exp(Fx);
sc = sum(Fx);
posmy = Fx/sc;

%-2. Perform the Occam's window
%--------------------------------------------------------------------------
%Fxm = max(Fx);
%index = find(posmy>Fxm/20);

%-3. With the index of the accepted values perform the BMA
%--------------------------------------------------------------------------
iter    = 2000;
Jbma    = manyinverse{1}.J{1}.*0;
PostMax = zeros(length(Jbma),1);
qCbma   = manyinverse{1}.qC.*0;

T = manyinverse{1}.T;

tmp = zeros(size(Jbma,1),size(T,1));

for j=1:N
    qC = manyinverse{j}.qC*diag(manyinverse{j}.qV)'; % variance over time
    sdmean(j,:)=mean(sqrt(qC')); % get mean sd per source (too time consuming to account for changes over time
end

disp('running BMA');
for i = 1:iter
    b = spm_multrnd(posmy,1);            % Take a sample from posmy
    
    J = cell2mat(manyinverse{b}.J);
    T = manyinverse{b}.T;
    
    x = spm_normrnd(zeros(size(sdmean,2),1),diag(sdmean(b,:).^2),1); % Generate a Gaussian sample
    
    tmp = x*ones(1,size(T,1));
    
    J=J+tmp*T;
    [dum,maxind]=max(var((J*T')'));      % Get current estimate with largest variance
    PostMax(maxind)=PostMax(maxind)+1/iter; % Add to posterior distribution of maximum
    
    Jbma=Jbma+J./iter;
    qCbma=qCbma+manyinverse{b}.qC./iter; % Make a weighted covariance
    if i/100==round(i/100)
        disp(sprintf('BMA %d percent done',round(100*i/iter)));
    end
end
