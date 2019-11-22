function res = bf_features_regmulticov(BF, S)
% Simple covariance computation with regularization
%
% Mark Woolrich

%--------------------------------------------------------------------------
if nargin == 0
    regmulticov      = cfg_const;
    regmulticov.tag  = 'regmulticov';
    regmulticov.name = 'Regularized multiple covariance';
    regmulticov.val  = {0};
    
    res = regmulticov;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

ntrials = length(S.trials);

classchanind = find(strcmp(D.chanlabels,'Class')); %MWW

if ~isempty(classchanind)
    NK = max(squash(D(classchanind,:,:)));
else
    NK = 1;
end

for k = 1:NK
    
    YY = 0;
    ns = 0;
    
    spm('Pointer', 'Watch');drawnow;
    spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
    if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
    else Ibar = 1:ntrials; end
    
    num_of_invalid_covs=0;
    
    for i = 1:ntrials,
        if ~isempty(classchanind)
            sampleind = S.samples(D(classchanind, S.samples, S.trials(i)) == k);
        else
            sampleind = S.samples;
        end
        
        nsamps = length(sampleind);
        
        if(nsamps>1)
            Y  = D(S.channels, sampleind, S.trials(i));
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
        warning(['Only ' num2str(ntrials-num_of_invalid_covs) ' valid trial covariances for class ' num2str(k)]);
    end;
    
    C = YY/ns;
    
    features.C{k} = C;
    features.N{k} = ns;
    
end

res = features;