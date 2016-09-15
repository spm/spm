function [GCM,PEB] = spm_dcm_peb_to_gcm(PEB, GCM_template, options)
% Generate an array of DCMs based on parameters from a PEB model.
% 
% Any parameters not included in the PEB model will be fixed at their 
% prior means (alternative fixed values can be selected).
%
% -------------------------------------------------------------------------
% Inputs:
%
% PEB - PEB model containing at least the following fields:
%
%                PEB.Ep   - group level parameter means
% 
%                PEB.Pind - indices within the DCM of parameters included
%                           in the PEB.
%
%                PEB.beta - between-subjects variance for each parameter is
%                           set to a fraction of the within-subject DCM
%                           priors: GCM_template{x}.M.pC / beta
%
%                PEB.Ce   - alternatively, a between-subjects covariance
%                           matrix can be provided, in which case beta is 
%                           ignored
%
% GCM_template - cell array of dimension [1 x models], where each element
%                is a DCM. These DCMs provide the structure of the models
%                that will be simulated (so don't need to be estimated).
%
%                If any parameters are not included in the PEB, they will
%                be fixed at values in GCM_template{m}.Ep. If this is not
%                present, they will be fixed at the priors in
%                GCM_template.M.pE{1}.
%
%                Alternatively, a matrix of size [subjects x models] can be
%                given, allowing subject-specific values for parameters not
%                in included in the PEB.
%
% options - settings structure for the simulation with fields:
%
%           options.nsubjects - number of subjects to generate
%
%           options.ratio - the ratio of posterior:prior covariance for
%                           the PEB model. [default 2]
%
% Returns:
%
% GCM  - DCM array populated with simulated data
% PEB  - PEB structure updated with PEB.Ce if not already present
%
% -------------------------------------------------------------------------
% Example:
% GCM = spm_dcm_peb_to_gcm(PEB, DCM);
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id: spm_dcm_peb_to_gcm.m 6877 2016-09-15 14:09:36Z vladimir $

% Set defaults
try options.gen_idx;   catch, options.gen_idx = 1; end
try options.ratio;     catch, options.ratio = 2; end
try options.nsubjects; catch
    if size(GCM_template,1) > 1
        options.nsubjects = size(GCM_template,1); 
    else
        options.nsubjects = 20; 
    end
end

% Load PEB
if ischar(PEB)
    PEB = load(PEB);
    PEB = PEB.PEB;
end

% Load template DCMs
if ~iscell(GCM_template)
    GCM_template = {GCM_template};
end

% Build GCM from PEB
[GCM,PEB] = peb_to_gcm(PEB, GCM_template, options);

%--------------------------------------------------------------------------
function [GCM,PEB] = peb_to_gcm(PEB, GCM_template, options)
% Generate parameters for each subject and put them in a GCM structure
%
% PEB          - PEB from which to generate DCM parameters
% GCM_template - [1 x nm] vector of template DCMs
% options      - Options structure
%
% GCM          - [subjects x nm] cell array of models
% PEB          - Updated PEB

% Unpack
nm      = size(GCM_template,2);
ns      = options.nsubjects;
Pind    = PEB.Pind;

% Get 'true' group-level parameter strengths from PEB
Ep = PEB.Ep;

if isfield(PEB,'Ce')
    % Get between-subjects (co)variance from PEB
    Ce = PEB.Ce;
else
    % Generate between-subjects (co)variance matrix
    if isfield(PEB,'Eh')
        Eh = PEB.Eh; % expected log precision
    else
        Eh = 0;
    end
        
    if ischar(GCM_template{1})
        GCM_template{1} = load(GCM_template{1});
        GCM_template{1} = GCM_template{1}.DCM;
    end    
    [i,pC] = spm_find_pC(GCM_template{1,1});
    
    pC = pC(Pind,Pind);
    Ce = exp(Eh) * (pC / PEB.beta);
    
    PEB.Ce = Ce;
end

% Sample parameters for each subject (parameters x subjects)
sEp = zeros(size(Ep,1),ns);
for s = 1:ns
    wEp      = Ep * PEB.M.X(s,:)';
    sEp(:,s) = spm_normrnd(wEp,Ce,1);
end

% Build GCM
GCM = cell(ns,nm);
for s = 1:ns
    for m = 1:nm
        % Get template DCM
        if size(GCM_template,1) > 1
            DCM = GCM_template{s,m};
        else
            DCM = GCM_template{1,m};
        end
        if ischar(DCM)
            DCM = load(DCM); 
            DCM = DCM.DCM;
        end
        
        % Get / set priors
        [i,pC,pE] = spm_find_pC(DCM);
        if ~isfield(DCM,'M') || ~isfield(DCM.M,'pC')
            DCM.M.pC = pC;
            DCM.M.pE = pE;
        end

        % Set posteriors to prior means if there are no posteriors
        if ~isfield(DCM,'Ep')
            DCM.Ep = DCM.M.pE;
        end
        
        % Set PEB-derived parameters for this subject
        ep       = spm_vec(DCM.Ep);
        ep(Pind) = sEp(:,s);    
        
        % Set any disabled parameters to their prior mean        
        i          = spm_find_pC(DCM);
        is_off     = true(length(ep),1);
        is_off(i)  = false;
        pe         = spm_vec(DCM.M.pE);
        ep(is_off) = pe(is_off);
        
        DCM.Ep   = spm_unvec(full(ep), DCM.M.pE);
        
        % Set covariance to be a fraction of the prior
        if isa(DCM.M.pC, 'struct')
            DCM.Cp = diag(spm_vec(DCM.M.pC)) / options.ratio;
        else
            DCM.Cp = diag(diag(DCM.M.pC)) / options.ratio;
        end
        
        % Set any disabled parameters to zero variance
        DCM.Cp(is_off,is_off) = 0;
        
        % Store updated DCMs
        GCM{s,m} = DCM;
    end
end