function xp = spm_beta_compare(alpha1,alpha2,Nsamp)
% Compute probability that r1 > r2
% FORMAT xp = spm_beta_compare(alpha1,alpha2,Nsamp)
% 
% Input:
% alpha1    - Beta parameters for first density
% alpha2    - Beta parameters for second density
% Nsamp     - number of samples used to compute xp [default = 1e4]
% 
% Output:
% xp        - exceedance probability
%
% Compute probability that r1 > r2 where p(r1)=Beta(r1|alpha1), 
% p(r2)=Beta(r2|alpha2). Uses sampling. 
% Useful for comparing groups in RFX model inference
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_beta_compare.m 3458 2009-10-13 10:05:35Z maria $

if nargin < 3
    Nsamp = 1e4;
end

Nk = length(alpha1);

% Perform sampling in blocks
%--------------------------------------------------------------------------
blk = ceil(Nsamp*Nk*8 / 2^28);
blk = floor(Nsamp/blk * ones(1,blk));
blk(end) = Nsamp - sum(blk(1:end-1));

xp = 0;
for i=1:length(blk)
    
    % Sample from univariate gamma densities then normalise
    % (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
    % 209-230)
    %----------------------------------------------------------------------
    r1 = zeros(blk(i),Nk);
    for k = 1:Nk
        r1(:,k) = spm_gamrnd(alpha1(k),1,blk(i),1);
    end
    sr1 = sum(r1,2);
    for k = 1:Nk
        r1(:,k) = r1(:,k)./sr1;
    end
    
    r2 = zeros(blk(i),Nk);
    for k = 1:Nk
        r2(:,k) = spm_gamrnd(alpha2(k),1,blk(i),1);
    end
    sr2 = sum(r2,2);
    for k = 1:Nk
        r2(:,k) = r2(:,k)./sr2;
    end
    
    % Exceedance probabilities:
    
    xp = xp + length(find (r1(:,1)>r2(:,1)));
    
end
xp = xp / Nsamp;
