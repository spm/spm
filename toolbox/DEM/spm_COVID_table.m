function Tab = spm_COVID_table(Ep,Cp,M)
% FORMAT Tab = spm_COVID_table(Ep,Cp,M)
% FORMAT Tab = spm_COVID_table(DCM)
% Ep  - conditional expectations
% Cp  - conditional covariances
% M   - model
%
% Tab - table of conditional estimators
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Get data for the United Kingdom (including total tests R)
%==========================================================================

% setup
%--------------------------------------------------------------------------
if nargin < 2
    Cp = Ep.Cp;
    M  = Ep.M;
    Ep = Ep.Ep;
end

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
Cp    = spm_unvec(diag(Cp),pC);
s     = spm_invNcdf(1 - 0.05);
for i = 1:numel(str.field)
    
    % scale parameters
    %----------------------------------------------------------------------
    ep = Ep.(str.field{i});            % posterior expectation
    pe = pE.(str.field{i});            % prior expectation
    pc = pC.(str.field{i});            % prior covariance
    cp = Cp.(str.field{i});            % posterior covariance
    
    ep = ep(1);            % posterior expectation
    pe = pe(1);            % prior expectation
    pc = pc(1);            % prior covariance
    cp = cp(1);            % posterior covariance
    
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

    Tab{i,1}  = sprintf('%g',i);
    Tab{i,2}  = sprintf('%s',lower(str.field{i}));
    Tab{i,3}  = sprintf('%s',str.names{i});
    
    Tab{i,4}  = sprintf('%g',eP);
    Tab{i,5}  = sprintf('%g',lb);
    Tab{i,6} = sprintf('%g',ub);
    
    Tab{i,7}  = sprintf('%g',Pe);
    Tab{i,8}  = sprintf('%g',1/pc);

    Tab{i,9}  = sprintf('%g',lB);
    Tab{i,10}  = sprintf('%g',uB);
   
end

VariableNames{1}  = 'number';
VariableNames{2}  = 'name';
VariableNames{3}  = 'description';

VariableNames{4}  = 'posterior';
VariableNames{5}  = 'from';
VariableNames{6} = 'to';

VariableNames{7}  = 'prior';
VariableNames{8}  = 'precision';
VariableNames{9}  = 'lower';
VariableNames{10}  = 'upper';

Tab = cell2table(Tab);
Tab.Properties.Description  = 'parameter estimates';
Tab.Properties.VariableNames = VariableNames;

writetable(Tab,'Parameters.csv','FileType','text','Delimiter',',')

% age specific parameters
%==========================================================================
clear tab
i     = 0;
for j = 1:numel(str.field)
    
    % scale parameters
    %----------------------------------------------------------------------
    ep = Ep.(str.field{j});                % posterior expectation  
    ep = ep(:,1);
    if size(ep,1) > 1
        
        i = i + 1;
        
        % posteriors
        %----------------------------------------------------------------------
        eP = exp(ep);                      % exponentiated
        
        tab{i,1}  = sprintf('%g',j);
        tab{i,2}  = sprintf('%s',lower(str.field{j}));
        tab{i,3}  = sprintf('%s',str.names{j});
        
        for k = 1:size(ep,1)
            tab{i,3 + k}  = sprintf('%g',eP(k,1));
        end
        
    end
    
end

tab = cell2table(tab);
tab.Properties.Description  = 'parameter estimates'

return







