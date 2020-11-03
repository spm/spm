function Tab = spm_COVID_table(Ep,Cp,M)
% FORMAT Tab = spm_COVID_table(Ep,Cp,M)
% Ep  - conditional expectations
% Cp  - conditional covariance as
%
% Tab - table of conditional estimators
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_COVID_table.m 8001 2020-11-03 19:05:40Z karl $

% Get data for the United Kingdom (including total tests R)
%==========================================================================

% get parameter strings
%--------------------------------------------------------------------------
[pE,pC,str] = spm_SARS_priors;

% use current priors if supplied
%--------------------------------------------------------------------------
try
    pE = M.pE;
    pC = M.pC;
end

% create a table of various risks (please see below)
%--------------------------------------------------------------------------
s     = spm_invNcdf(1 - 0.05);
for i = 1:numel(str.field)
    
    % scale parameters
    %----------------------------------------------------------------------
    ep = getfield(Ep,str.field{i});    % posterior expectation
    pe = getfield(pE,str.field{i});    % prior expectation
    pc = getfield(pC,str.field{i});    % prior covariance
    cp = Cp(i,i);                      % posterior covariance
    
    % priors
    %----------------------------------------------------------------------
    Pe = exp(pe);                      % exponentiated
    lB = exp(pe - sqrt(pc)*s);         % lower confidence bound
    uB = exp(pe + sqrt(pc)*s);         % upper confidence bound  
    
    % posteriors
    %----------------------------------------------------------------------
    eP = exp(ep);                      % exponentiated
    lb = exp(ep - sqrt(cp)*s);         % lower confidence bound
    ub = exp(ep + sqrt(cp)*s);         % upper confidence bound

    Tab{i,1} = i;
    Tab{i,2} = str.field{i};
    Tab{i,3} = str.names{i};
    
    Tab{i,4} = Pe;
    Tab{i,5} = lB;
    Tab{i,6} = uB;
    
    Tab{i,7} = eP;
    Tab{i,8} = lb;
    Tab{i,9} = ub;
        
end
VariableNames{1} = 'number';
VariableNames{2} = 'name';
VariableNames{3} = 'description';
VariableNames{4} = 'prior';
VariableNames{5} = 'lower';
VariableNames{6} = 'upper';
VariableNames{7} = 'posterior';
VariableNames{8} = 'from';
VariableNames{9} = 'to';

Tab = cell2table(Tab);
Tab.Properties.Description  = 'parameter estimates';
Tab.Properties.VariableNames = VariableNames;








