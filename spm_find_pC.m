function [i,pC,pE,Np] = spm_find_pC(varargin)
% Utility routine that finds the indices of non-zero covariance
% FORMAT [i,pC,pE,Np] = spm_find_pC(pC,pE,fields)
% FORMAT [i,pC,pE,Np] = spm_find_pC(DCM,fields)
% FORMAT [i,pC,pE,Np] = spm_find_pC(DCM)
% 
% pC     - covariance matrix or variance stucture
% pE     - parameter structure
% fields - desired fields of pE
%
% or
%
% DCM    - DCM structure
%
% i      - find(diag(pC) > TOL)
% rC     - reduced covariances
% rE     - reduced expectation
% 
%__________________________________________________________________________
% Copyright (C) 2015-2019 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_find_pC.m 7714 2019-11-26 11:25:50Z spm $

%-parse input arguments
%--------------------------------------------------------------------------
if nargin > 2
    pC     = varargin{1};
    pE     = varargin{2};
    fields = varargin{3};
elseif numel(varargin) > 1
    DCM    = varargin{1};
    fields = varargin{2};
else
    DCM    = varargin{1};
end

%-get prior density from DCM
%--------------------------------------------------------------------------
if nargin < 3
    if ischar(DCM)
        DCM = load(DCM,'DCM');
        DCM = DCM.DCM;
    end
    if any(isfield(DCM,{'options','M'}))
        try, [pC,pE] = spm_find_rC(DCM); end
    end
end

%-Deal with variance structures
%--------------------------------------------------------------------------
if isstruct(pC)
    q = spm_vec(pC);
else
    q = diag(pC);
end

%-Get indices
%--------------------------------------------------------------------------
i  = find(q > mean(q(q < 1024))/1024);
Np = numel(q);

%-subsample fields if necessary
%--------------------------------------------------------------------------
if nargin > 1
    if ischar(fields), fields = {fields}; end
    if isstruct(pE)
        j = spm_fieldindices(pE,fields{:});
        if isempty(j) && ~(~isempty(fields) && strcmp(fields{1},'none'))
            warning('%s not found. Returning all fields',...
                strjoin(cellstr(fields),','));
        else
            i = j(ismember(j,i));
        end
    end
end

return


function [pC,pE] = spm_find_rC(DCM)
% FORMAT [pC,pE] = spm_find_rC(DCM)
% model priors
%__________________________________________________________________________

% Get full priors and posteriors
%--------------------------------------------------------------------------
try
    pC  = DCM.M.pC;
    pE  = DCM.M.pE;
    return
end

% get priors from model specification
%--------------------------------------------------------------------------
if isfield(DCM.options,'spatial')
    
    % EEG or MEG
    %----------------------------------------------------------------------
    if strcmpi(DCM.options.analysis,'IND')
        [pE,dummy,pC] = spm_ind_priors(DCM.A,DCM.B,DCM.C,DCM.Nf);
    else
        [pE,  pC] = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);        
        
        try                                                     %#ok<TRYNC>
            try model   = DCM.options.model;   catch, model    = 'NMM'; end
            try spatial = DCM.options.spatial; catch, spatial  = 'LFP'; end
            DCM.M.dipfit.model = model;
            DCM.M.dipfit.type  = spatial;        
            [pE,  pC] = spm_L_priors(DCM.M.dipfit,pE,pC);
            [pE,  pC] = spm_ssr_priors(pE,pC);
        end
    end
    
else
    
    % fMRI
    %----------------------------------------------------------------------
    if ~isfield(DCM,'d')
        DCM.d = zeros(size(DCM.a,1),size(DCM.a,1),0);
    end
    [pE,pC]   = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);
end
