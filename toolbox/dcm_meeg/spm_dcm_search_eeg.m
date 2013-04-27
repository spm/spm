function spm_dcm_search_eeg(P)
% Post hoc optimisation of DCMs (under Laplace approximation)
% FORMAT spm_dcm_search_eeg([P])
%
% P - character/cell array of DCM filenames or cell array of DCM structures
%
%--------------------------------------------------------------------------
% spm_dcm_search_eeg operates on different DCMs of the same data to find
% the best model. It will invert the full model whose free-parameters are
% the union (superset) of all free parameters in each model specified. The
% routine then uses a post hoc selection procedure to evaluate the log-
% evidence and conditional density over free-parameters of each model
% specified.
%
% The outputs of this routine are graphics reporting the model space search
% (optimisation) and a DCM_optimum (in the first DCMs directory) for the
% best DCM. The structural and function (spectral embedding) graphs are
% based on this DCM.
%
% Conditional esimates (Ep, Cp and F values) in DCM_??? (specifed by P) are
% replaced by their reduced estimates (but only these esimates) in rDCM_???
%
% DCM_optimum  contains the fields:
%
%        DCM.P   - character/cell array of DCM filenames
%        DCM.PF  - their associated free energies
%        DCM.PP  - and posterior (model) probabilities
%
% In addition, the free energies and posterior estimates of each DCM in P
% are saved for subsequent searches over different partitions of model
% space.
%
% See alos: spm_dcm_post_hoc.m
%
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_search_eeg.m 5454 2013-04-27 10:46:41Z karl $

% get filenames
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end

if ischar(P), P = cellstr(P); end
N = numel(P);

%-Check models are compatible in terms of their data
%==========================================================================
for j = 1:N
    
    % get prior covariances
    %----------------------------------------------------------------------
    try, load(P{j}); catch, DCM = P{j}; end
    
    % and compare it with the first model
    %----------------------------------------------------------------------
    if j == 1
        Y = DCM.xY.y;
    else
        try
            if any(spm_vec(Y) - spm_vec(DCM.xY.y))
                fprintf('\nPlease check model %i for data compatibility\n',j)
                return
            end
        catch
            fprintf('\nPlease check model %i for compatibility\n',j)
            return
        end
    end
    
    % number of free parameters
    %----------------------------------------------------------------------
    par(:,j) = logical(spm_vec(DCM.A,DCM.B,DCM.C));

end

%-load full model
%==========================================================================
[i j] = max(sum(par));
try, load(P{j}); catch, DCM = P{j}; end

% check this is a full model and that is has been inverted
% -------------------------------------------------------------------------
if any(any(par,2) > par(:,j))
    fprintf('\nPlease ensure your models are nested\n')
    return
end
if ~all(isfield(DCM,{'Ep','Cp'}))
    fprintf('\nPlease invert model %i\n',j)
    return
end

% directory to save reduced [r]DCMs
% -------------------------------------------------------------------------
pathname = spm_file(DCM.name,'path');
if ~exist(pathname,'dir');
    pathname = pwd;
end


% Get full priors and posteriors
% -------------------------------------------------------------------------
qE    = DCM.Ep;
qC    = DCM.Cp;
pE    = DCM.M.pE;
pC    = DCM.M.pC;

%-Loop through models and get log-evidences
%==========================================================================
for j = 1:N
    
    % Get reduced model specification
    % ---------------------------------------------------------------------
    try, load(P{j}); catch, DCM = P{j}; end
    
    switch DCM.options.analysis
        
        % Get model (priors)
        % -----------------------------------------------------------------
        case {'ERP'}
            
            [rE,rC] = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
            
        otherwise
            
            fprintf('\nDCM analysis (%s) is not currently supported\n',DCM.options.analysis);
            
    end   
    
    % evaluate (reduced) free-energy and posteriors
    % ---------------------------------------------------------------------
    [F,Ep,Cp] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
    
    % Put reduced conditional estimates in DCM
    % =====================================================================
    
    % Store parameter estimates
    %----------------------------------------------------------------------
    DCM.M.pC = rC;
    DCM.Ep   = Ep;
    DCM.Cp   = Cp;
    DCM.F    = F;
    
    filename = spm_file(DCM.name,'basename');
    name{j}  = filename(5:end);
    filename = spm_file(DCM.name,'filename');
    filename = spm_file(filename,'prefix','r');
    filename = fullfile(pathname,filename);
    
    
    % Save DCM
    %======================================================================
    save(filename,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));
    P{j} = DCM;
    
    % Record free-energy and MAP estimates
    %----------------------------------------------------------------------
    P{j} = DCM;
    G(j) = F;
    
end

% Model evidences and best model
% =========================================================================
G     = G - min(G);
p     = exp(G - max(G));
p     = p/sum(p);

% Get selected model
%--------------------------------------------------------------------------
[q,j] = max(p);
try, load(P{j}); catch, DCM = P{j}; end

i   = spm_fieldindices(DCM.Ep,'A','B','C');
qE  = spm_vec(DCM.Ep);
Ep  = spm_vec(DCM.Ep);
qC  = DCM.Cp;
Cp  = DCM.Cp;
F   = DCM.F;

% Show results
% -------------------------------------------------------------------------
spm_figure('Getwin','Graphics'); clf

subplot(2,2,1)
if length(P) > 32, plot(G,'k'), else bar(diag(G),N), end
title('log-posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('log-probability','FontSize',12)
legend(name)
axis square

subplot(2,2,2)
if length(P) > 32, plot(p,'k'), else, bar(diag(p),N), end
title('model posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('probability','FontSize',12)
axis square

% Show full and reduced conditional estimates (for optimum DCM)
%--------------------------------------------------------------------------
subplot(2,2,3)
spm_plot_ci(qE(i),qC(i,i))
title('MAP connections (full)','FontSize',16)
axis square
a   = axis;

subplot(2,2,4)
spm_plot_ci(Ep(i),Cp(i,i))
title('MAP connections (optimum)','FontSize',16)
axis square
axis(a)

% Show structural and functional graphs
%--------------------------------------------------------------------------
spm_figure('Getwin','Graph'); clf

spm_dcm_graph(DCM,DCM.Ep.A)


%-Save optimum and full DCM
%==========================================================================
DCM.P  = P;
DCM.PF = G;
DCM.PP = p;

% Reduced model (optimum)
%--------------------------------------------------------------------------
filename = fullfile(pathname,'DCM_optimum.mat');
save(filename,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));


