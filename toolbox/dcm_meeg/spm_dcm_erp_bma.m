function [r,xp] = spm_dcm_erp_bma (BMS_file,stats,params)
% Compute posterior over connections, or modulatory gains in them, from BMA
% FORMAT [r,xp] = spm_dcm_erp_bma (BMS_file,stats,params)
%
% BMS_file      Name of Bayesian Model Selection .mat file
% stats         'ffx' or 'rfx' 
% params        Parameter data structure
%               .type          'A' (connection) or 'B' (gain in connection)
%               .hier          'forward', 'backward' or 'lateral' (for conn_type='A')
%               .ip            eg 1, 2, 3 indexes modulatory input (for conn_type='B')
%               .to            to region eg 3
%               .from          from region eg 1
%               .xt            exceedance threshold (typically set to 1)
%               .C             [nr x nr] contrast matrix where nr is the number of regions 
%
%
% r             posterior samples
% xp            exceedance probability
%               This is the posterior probability that the connection is
%               larger than params.xt. Alternatively, if you are looking at
%               a contrast of connections, its the posterior probability
%               that the contrast is greater than zero.
%
% The parameters returned by Bayesian Model Averaging (BMA) are the 'latent'
% variables A and B which are Gaussian (and consequently can be positive or
% negative).  
%
% The corresponding connection strengths (rA) or gains in connection
% strength (rB) are an exponential function of these latent variables. 
% These are the values we are interested in and want to make an inference
% about.
%
% This routine computes the posterior distribution over rA or rB by
% generating samples from the latent variables, and exponentiating each
% sample. 
%
% The probability that the rA or RB values are greater than some threshold
% xt (such as unity) is then just the proportion of posterior samples that
% are greater than xt.
%
% If a contrast matrix (C) is not specifed this function looks at a single
% connection or gain. To look at relative sizes of connection/gain values 
% enter a C matrix. eg. to test, in a 3 region DCM, is connection from 3
% to 2 bigger than 2 to 3 ? set C=[0 0 0; 0 0 1; 0 -1 0].
%
%--------------------------------------------------------------------------
% 
% Example usage: 
%
% 1. Look at a single connection value:
%
% params.type='A'; params.hier='forward'; 
% params.to=3; params.from=1; params.xt=1;
% spm_dcm_erp_bma([],'ffx',params);
%
% 2. Look at a single gain value:
%
% params.type='B'; params.ip=1; 
% params.to=1; params.from=1; params.xt=1;
% spm_dcm_erp_bma([],'ffx',params);
%
% 3. Look at a contrast of connection values:
%
% params.type='B'; params.ip=1;
% params.C=[0 0 0; 0 0 1; 0 -1 0];
% spm_dcm_erp_bma([],'ffx',params);
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_erp_bma.m 4529 2011-10-18 14:58:20Z guillaume $


%-Get parameters
%--------------------------------------------------------------------------
try
    load (BMS_file);
catch
    [BMS_file,sts] = spm_select(1,'^BMS\.mat$','Select BMS.mat file');
    if ~sts, return; end
    load (BMS_file);
end

if isfield(params,'C')
    con = 1;
else
    con = 0;
end
    
switch lower(stats)
    case 'ffx'
        bma = BMS.DCM.ffx.bma;
    case 'rfx'
        bma = BMS.DCM.rfx.bma;
    otherwise
        error('The stats field shout be set to ''ffx'' or ''rfx''');
end

%-Get mean and SD of effects
%--------------------------------------------------------------------------
switch params.type
    case 'A'
        hier_type = {'forward','backward','lateral'};
        h = find(ismember(hier_type,params.hier));
        M = bma.mEp.A{h};
        S = bma.sEp.A{h};
        if ~con
            str=sprintf('%s connection (rA) from region %d to %d:',hier_type{h},params.from,params.to);
        end
    case 'B'
        M = bma.mEp.B{params.ip};
        S = bma.sEp.B{params.ip};
        if ~con
            str=sprintf('Effect (rB) of modulatory input %d on connection from region %d to %d:',params.ip,params.from,params.to);
        end
    otherwise
        error('conn_type must be ''A'' or ''B''');
end
        

%-Generate N samples in latent space, then exponentiate to get estimate
% of effects
%--------------------------------------------------------------------------
N=10000;

switch con
    case 1
        % Look at contrast of connections
        nr = size(M,1);
        m  = M(:);
        s  = S(:);
        c  = params.C(:);
        
        % Full posterior covariances are not stored 
        % only univariate SD's so use these
        x  = (s*ones(1,N)).*randn(nr*nr,N)+m*ones(1,N);
        r  = c'*exp(x);
        rm = mean(r);
        xp = length(find(r>0))/length(r);
        switch params.type
            case 'A'
                str=sprintf('Contrast of %s connections (rA)',hier_type{h});
            case 'B'
                str=sprintf('Contrast of effects (rB) of modulatory input %d ',params.ip);
        end
    case 0
        % Just look at a single connection
        m  = M(params.to,params.from);
        s  = S(params.to,params.from);
        x  = s*randn(N,1)+m;
        r  = exp(x);
        rm = mean(r);
        xp = length(find(r>params.xt))/length(r);
end

%-Display
%--------------------------------------------------------------------------
figure
hist(r,20);
set(gca,'FontSize',18);
set(gca,'YTickLabel',[]);
ylabel(sprintf('p(r%s|Y)',params.type));
xlabel(sprintf('r%s',params.type));

%-Report
%--------------------------------------------------------------------------
fprintf('%s\n',str);
fprintf('Posterior mean = %1.2f\n',rm);
fprintf('Exceedance probability = %1.6f\n', xp);
