function res = bf_features_cov_bysamples(BF, S)
% Simple covariance computation to handle variable width WOIs, 
% Requires S.samples as a [1 x samples x ntrials] matrix of logical indices
% indicating which data points should be used in the cov estimation
%
% Mark Woolrich 2014

%--------------------------------------------------------------------------
if nargin == 0,
    cov_bysamples      = cfg_const;
    cov_bysamples.tag  = 'cov_bysamples';
    cov_bysamples.name = 'Cov estimation with variable width WOIs';
    cov_bysamples.val  = {0};
    
    res = cov_bysamples;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

YY = 0;
ns = 0;

ntrials = length(S.trials);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

num_of_invalid_covs=0;

for i = 1:ntrials,
    nsamps=sum(S.samples(1,:,S.trials(i))>0);
    
    if(nsamps>1)
        Y  = D(S.channels, find(S.samples(1,:,S.trials(i))), S.trials(i));
        Y  = detrend(Y', 'constant');
        YY = YY+(Y'*Y);
        ns = ns + nsamps - 1;      
    else
        num_of_invalid_covs=num_of_invalid_covs+1;
    end;    
    
    if ismember(i, Ibar),        
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

if(ntrials-num_of_invalid_covs < 10),
    warning(['Only ' num2str(ntrials-num_of_invalid_covs) ' valid trial covariances']);
end;

C = YY/ns;

features.C = C;
features.N = ns;

res = features;