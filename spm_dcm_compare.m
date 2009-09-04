function spm_dcm_compare(P)
% Compare two or more estimated models
% FORMAT spm_dcm_compare(P)
%
% P  - a char or cell array of DCM filenames
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Klaas Enno Stephan
% $Id: spm_dcm_compare.m 3363 2009-09-04 15:11:19Z christophe $


% Get DCM filenames
%--------------------------------------------------------------------------
if ~nargin
    P = spm_select([2 inf],'^DCM.*\.mat$','select DCM*.mat files');
end
if ~iscell(P), P = cellstr(P); end
num_models = numel(P);

% Load anc check all models and compute their evidence
%--------------------------------------------------------------------------
name = {};
for model_index=1:num_models
    
    try
        load(P{model_index},'DCM','-mat');
    catch
        error('Cannot load DCM file "%s".',P{model_index});
    end

    % Check that none of the models is an averaged DCM
    %----------------------------------------------------------------------
    if isfield(DCM,'averaged') && DCM.averaged
        str = {...
            ['Model ' P{model_index} ' is an averaged DCM. '],...
             'Please note that model comparison is not valid for averaged DCMs. ',...
             'Procedure aborted.'};
        spm('alert*',str,'DCM <Compare>',spm('Cmdline'));
        return
    end
    
    % Check that all models refer to the same set of VOIs
    %----------------------------------------------------------------------
    if model_index == 1
        VOIs = {DCM.xY.name};
    else
        if ~isequal(VOIs,{DCM.xY.name})
            str = {...
                'Selected models contain different sets of VOIs!',...
                'Please note that model comparison is only valid for models with identical time series (VOIs).',...
                'Procedure aborted.'};
            spm('alert*',str,'DCM <Compare>',spm('Cmdline'));
            return
        end
    end
    
    % Compute Model Evidence
    %----------------------------------------------------------------------
    [p,filename]           = fileparts(P{model_index});
    name{end + 1}          = filename;
    evidence(model_index)  = DCM.F;
    
end


% compute conditional probability of DCMs under flat priors.
%--------------------------------------------------------------------------
F    = evidence - min(evidence);
i    = F < (max(F) - 32);
P    = F;
P(i) = max(F) - 32;
P    = P - min(P);
P    = exp(P);
P    = P/sum(P);

   
% display results
%--------------------------------------------------------------------------   
Fgraph = spm_figure('GetWin','Graphics'); figure(Fgraph); clf

subplot(2,1,1)
n      = length(name);
barh(1:n,F), hold on
plot([1 1]*(max(F) - 3),[0 n],'LineWidth',4,'Color','r'), hold off
set(gca,'YTick',1:n)
set(gca,'YTickLabel',name)
xlabel('log-evidence (relative)')
title('Bayesian model comparison','FontSize',16)
axis square
grid on

subplot(2,1,2)
barh(1:n,P)
set(gca,'YTick',1:n)
set(gca,'YTickLabel',name)
title({'conditional model probability';'under uniform model priors'})
xlabel('posterior probability')
axis square
grid on


% Output model comparison details to MATLAB command window
%--------------------------------------------------------------------------
assignin('base','DCM_cond_prob',P)
assignin('base','DCM_log_ev',F)

fprintf('\n\n DCM Model Comparison \n\n Relative Log Model Evidences:')
for i = 1:length(F)
    fprintf( '\n %s:  %-6.2f ',name{i},F(i));
end
fprintf('\n\n Conditional Model Probabilities under flat priors:')
for i = 1:length(F)
    fprintf( '\n %s:  %-6.2f',name{i},P(i));
end
fprintf('\n\n These have been stored in the workspace as ''DCM_log_ev '' and ''DCM_cond_prob''\n')
