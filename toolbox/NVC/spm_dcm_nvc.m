function DCM  = spm_dcm_nvc(P)
% Specify and estimate a DCM for multimodal fMRI and M/EEG
% FORMAT [DCM] = spm_dcm_nvc(P)
%
% Input:
% -------------------------------------------------------------------------
% P{1} - SPM structure or location of SPM.mat
% P{2} - Cell array of VOI filenames (the same order as sources in EEG DCM)
% P{3} - Location of DCM for M/EEG .mat file or DCM structure
% P{4} - Model specification for neurovascular coupling (NVC) mechanisms
% P{5} - Which neuronal populations should drive haemodynamics 
% P{6} - Which fMRI experimental conditions to include
% P{7} - DCM options
%
% Where:
%
% P{4} - A cell array of strings with three elements:
%        
%   P{4}{1} - 'pre', 'post' or decomposed ('de') neuronal signals excite 
%              NVC. Decomposed means activity is grouped into intrinsic-
%              inhibitory, intrinisic-excitatory and extrinsic-excitatory.
%   P{4}{2} - NVC has the same ('s') or different ('d') parameters for all
%              regions. 
%   P{4}{3} - extrinsic and intrinsic ('ext') or only intrinsic ('int') 
%              neuronal activity contributes to regional BOLD 
%              (for 'post', this should be 'na').
%    
%   Supported options:
%   {'pre','d','int'},{'pre','s','int'}, {'pre','d','ext'},{'pre','s','ext'},
%   {'de','d', 'int'},{'de','d','exc'},  {'de','s','int'}, {'de','s','exc'},
%   {'post','d','na'},{'post','s','na'};
%
%   Example: P{4} = {'pre', 's', 'int'} means presynaptic neuronal drive
%   (from intrinsic connections only) inputs to a model of neurovascular
%   coupling that has the same parameters for all regions.
%
% P{5} - Which neuronal populations should drive haemodynamics, by setting
%   ones or zeros in a vector ordered:
%   [superficial pyramidal, inhibitory, excitatory, deep pyramidal]
%   (default is [1 1 1 1]).
%
%   Example: [1 0 1 1] means no NVC drive from inhibitory populations.
%
% P{6} - Binary vector indicating which experimental conditions to include.
%
% P{7} - options structure for DCM for fMRI:  
%    options.name                   % name for the DCM
%    options.maxit                  % maximum number of iterations
%    options.hE                     % expected precision of the noise
%    options.hC                     % variance of noise expectation
%    options.TE                     % echo time (default: 0.04)
%
% Evaluates:
% -------------------------------------------------------------------------
% DCM.M                              % Model structure
% DCM.Ep                             % Condition means (parameter structure)
% DCM.Cp                             % Conditional covariances
% DCM.Vp                             % Conditional variances
% DCM.Pp                             % Conditional probabilities
% DCM.H1                             % 1st order hemodynamic kernels
% DCM.H2                             % 2nd order hemodynamic kernels
% DCM.K1                             % 1st order neuronal kernels
% DCM.K2                             % 2nd order neuronal kernels
% DCM.R                              % residuals
% DCM.y                              % predicted data
% DCM.T                              % Threshold for Posterior inference
% DCM.Ce                             % Error variance for each region
% DCM.F                              % Free-energy bound on log evidence
% DCM.ID                             % Data ID
% DCM.AIC                            % Akaike Information criterion
% DCM.BIC                            % Bayesian Information criterion
%
% Notes on parameters:
% -------------------------------------------------------------------------
% This scheme estimates DCM.H (haemodynamic parameters) and DCM.J
% (neurovascular coupling parameters):
% 
% DCM.Ep.H.transit - transit time (t0)
% DCM.Ep.H.decay   - signal decay d(ds/dt)/ds)
% DCM.Ep.H.epsilon - ratio of intra- to extra-vascular components of the
%                    gradient echo signal
%
% DCM.Ep.J - neurovascular coupling parameters. The dimension depends upon 
% the requested model specification. For p populations and n regions:
%
% P{7} (DCM.model)        dim(J)   notes
% =========================================
% {'pre'  'd' 'int'}      [p n]      
% {'pre'  's' 'int'}      [p 1]
% {'pre'  'd' 'ext'}      [p n]
% {'pre'  's' 'ext'}      [p 1]
% {'de'   's' 'int}       [p 2]    dim2: intrinsic inhibitory, excitatory
% {'de'   's' 'ext'}      [p 3]    dim2: intrinsic inhibitory, excitatory, extrinsic 
% {'de'   'd' 'int}       [p 2 n]  dim2: intrinsic inhibitory, excitatory
% {'de'   'd' 'ext'}      [p 3 n]  dim2: intrinsic inhibitory, excitatory, extrinsic 
% {'post' 's' 'na'}       [p 1]
% {'post' 'd' 'na'}       [p n]
%
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%
% Friston, K.J., Preller, K.H., Mathys, C., Cagnan, H., Heinzle, J., Razi, A.
% and Zeidman, P., 2017. Dynamic causal modelling revisited. Neuroimage.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging
 
% Amirhossein Jafarian
% $Id: spm_dcm_nvc.m 7735 2019-12-02 09:15:27Z peter $

% Prepare input
%--------------------------------------------------------------------------
SPM                 =  P{1};
xY                  =  P{2};
MEEG                =  P{3};
try model           =  P{4}; catch, model        = {'pre','d','int'}; end
try n_exclude       =  P{5}; catch, n_exclude    = ones(1,4)        ; end
try sess_exclude    =  P{6}; catch, sess_exclude = 'not defined'    ; end
try options         =  P{7}; catch, options = struct()              ; end   

try options.centre;    catch, options.centre      =  1;      end
try options.hE;        catch, options.hE          =  6;      end
try options.hC;        catch, options.hC          =  1/128;  end
try options.maxit;     catch, options.maxit       =  128;    end
try options.TE;        catch, options.TE          =  0.04;   end
try name=options.name; catch, name = sprintf('DCM_%s',date); end

% Set defaults
%--------------------------------------------------------------------------
if (length(n_exclude) < 4 || isempty(n_exclude))
    n_exclude     = ones(1,4);
end
if isempty(sess_exclude)
    sess_exclude = 'not defined';
end

% Specify haemodynamic model
%--------------------------------------------------------------------------
DCM    = spm_dcm_nvc_specify(SPM,xY,MEEG,model,n_exclude,sess_exclude,options);
n      = DCM.n;         
v      = DCM.v;         
U.dt   = DCM.U.dt; 
U.u    = [];
 
% Prepare fMRI signals 
%--------------------------------------------------------------------------
Y       = DCM.Y;         
Y.y     = spm_detrend(Y.y);
scale   = max(max(Y.y)) - min(min(Y.y));
scale   = 4/max(scale,4);
Y.y     = Y.y*scale;
Y.scale = scale;
if ~isfield(Y,'X0'),Y.X0 = ones(v,1); end
if ~size(Y.X0,2),   Y.X0 = ones(v,1); end

% Set fMRI slice time sampling and echo time
%--------------------------------------------------------------------------
try M.delays = DCM.delays; catch, M.delays = ones(n,1); end
M.TE = options.TE;

% Set priors
%--------------------------------------------------------------------------
[pE,pC,x] = spm_dcm_nvc_priors(DCM);

% Prepare neuronal drive functions
%--------------------------------------------------------------------------
input     = spm_dcm_nvc_nd(DCM);

% Set hyperpriors over precision - expectation and covariance
%--------------------------------------------------------------------------
hE       = sparse(n,1) + DCM.options.hE;
hC       = speye(n,n)  * DCM.options.hC;

% Prepare model structure
%--------------------------------------------------------------------------
M.IS = @spm_nvc_gen;                       
M.x  = x;                                
M.pE = pE;                                
M.pC = pC;                                
M.hE = hE;                                
M.hC = hC;                                
M.m  = n;                            
M.n  = size(spm_vec(x),1);
M.l  = n;
M.ns = v;
M.Nmax = options.maxit ;
M.TE  = DCM.TE;
M.input = input;
M.Model = model;
M.nograph = spm('CmdLine');

% Model Identification methods
% =========================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);

% Compute the prediction, residuals and error variance
%--------------------------------------------------------------------------
yhat     = feval(M.IS,Ep,M,U);
R        = Y.y - yhat;
R        = R - Y.X0*spm_inv(Y.X0'*Y.X0)*(Y.X0'*R);
Ce       = exp(-Eh);

% Compute variance that is captured by the data
%--------------------------------------------------------------------------
PSS      = sum(yhat.^2);
RSS      = sum(R.^2);
D        = 100.*PSS./(PSS + RSS);

% Compute kernels
%--------------------------------------------------------------------------
H.f      = @spm_fx_hdm;
H.g      = @spm_gx_hdm;
H.x      = M.x;
H.m      = n; 
[H0,H1]  = spm_kernels(H,Ep.H,64,1/2);

% Compute posterior probabilities per parameter
%--------------------------------------------------------------------------
T        = full(spm_vec(pE));
sw       = warning('off','SPM:negativeVariance');
Pp       = spm_unvec(1 - spm_Ncdf(T,abs(spm_vec(Ep)),diag(Cp)),Ep);
Vp       = spm_unvec(full(diag(Cp)),Ep);
warning(sw);
try  M = rmfield(M,'nograph'); end

% Store parameter estimates
%--------------------------------------------------------------------------
DCM.M       = M;
DCM.Y       = Y;
DCM.U       = U;
DCM.Ce      = Ce;
DCM.Ep      = Ep;
DCM.Cp      = Cp;
DCM.Pp      = Pp;
DCM.Vp      = Vp;
DCM.H1      = H1;
DCM.D       = D;
DCM.R       = R;
DCM.y       = yhat;
DCM.t       = (1:size(yhat,1))*DCM.Y.dt;

% Compute data ID
%--------------------------------------------------------------------------
if isfield(M,'FS')
    try
        ID = spm_data_id(feval(M.FS,Y.y,M));
    catch
        ID = spm_data_id(feval(M.FS,Y.y));
    end
else
    ID     = spm_data_id(Y.y);
end
 
% Store evidence
%--------------------------------------------------------------------------
DCM.F      = F;
DCM.ID     = ID;

%-Save DCM
%--------------------------------------------------------------------------
DCM.name = name; 
save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
