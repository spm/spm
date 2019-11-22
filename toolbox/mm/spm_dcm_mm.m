function DCM  = spm_dcm_mm(P)
% Estimate DCM for multimodal fMRI and M/EEG
% FORMAT [DCM] = spm_dcm_mm_estimate(P)
%
%--------------------------Input-------------------------------------------
% P{1,1}  - Address of SPM or SPM.mat
% P{1,2}  - Address of VOIs(order must be the same as sources in DCM for EEG)
% P{1,3}  - Address of DCM for MEG or DCM.mat for M/EEG
% P{1,4}  - Model of neurovascular coupling (NVC) mechanisms
%  P{1,4}(1) - 'pre', 'post' or decomposed ('de') neuronal signals excite NVC.
%  P{1,4}(2) - NVC has the same ('s') or different ('d') parameters for all
%              regions. 
%  P{1,4}(3) - extrinsic ('ext') or intrinsic ('int') neuronal activity
%              contribute to regional BOLD (for `post’, this should be 'na').
%  {'pre','d','int'},{'pre','s','int'}, {'pre','d','ext'},{'pre','s','ext'},
%  {'de','d', 'int'},{'de','d','exc'}, {'de','s','int'},{'de','s','exc'},
%  {'post','d','na'},  {'post','s','na'};
%        - Example: P{1,3} = {‘pre’, ‘s’, ‘int’} means  presynaptic neuronal drive
%       (excluded extrinsic neuronal drives) inputs to a model of neurovascular
%        coupling that has the same parameters for all regions.
% P{1,5}  - Excluding some neuronal drives from DCM for fMRI, by setting
%           to zero some of the entries in the following matrix:
%           [superficial pyramidal, inhibitory, excitatory, deep pyramidal]
%           (default is [1 1 1 1]).
% - Example: [ 1 0 1 1] means that inhibition would not be considered in DCM for fMRI.
% P{1,6}  - Excluding some sessions of SPM in DCM for fMRI.
% P{1,7}  - options for DCM:  
%           P{1,7}.name                   % name of DCM
%           P{1,7}.maxit                  % maximum number of iterations
%           P{1,7}.hE                     % expected precision of the noise
%           P{1,7}.hC                     % variance of noise expectation
%
% Evaluates:
%--------------------------------------------------------------------------
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
% $Id $
%
% Input
%--------------------------------------------------------------------------
SPM                 =  P{1,1};
xY                  =  P{1,2};
MEEG                =  P{1,3};
try Model           =  P{1,4};         catch, Model        = {'pre','d','int'};end
try N_exclude       =  P{1,5};         catch, N_exclude    = ones (1,4)       ;end
try Sess_exclude    =  P{1,6};         catch, Sess_exclude = 'not defined'    ;end
try options.centre  =   P{1,7}.centre; catch, options.centre      =   1       ;end
try options.hE      =   P{1,7}.hE ;    catch, options.hE          =   6       ;end
try options.hC      =   P{1,7}.hC ;    catch, options.hC          =  1/128    ;end
try options.maxit   =   P{1,7}.maxit;  catch, options.maxit       =  128      ;end
try name            =   P{1,7}.name;   catch, name = sprintf('DCM_%s',date)   ;end

%--------------------------------------------------------------------------

if (size(P,2) <5 || isempty(P{1,5}))
    N_exclude     = ones (1,4)   ;
end
if (size(P,2) < 6 || isempty(P{1,5}))
    Sess_exclude  = 'not defined';
end

if (size(P,2) < 7 || isempty(P{1,7}))
    name                = sprintf('DCM_%s',date);
    options.centre      =   1;
    options.hE          =   6;
    options.hC          =  1/128;
    options.maxit       =  128;
end

% Specify
%--------------------------------------------------------------------------
 DCM    = spm_dcm_mm_specify(SPM,xY, MEEG,Model,N_exclude,Sess_exclude,options);
 n      = DCM.n;         
 v      = DCM.v;         
 U.dt   = DCM.U.dt; 
 U.u    = [];
% fMRI signals 
%--------------------------------------------------------------------------
Y       = DCM.Y;         
Y.y     = spm_detrend(Y.y);
scale   = max(max((Y.y))) - min(min((Y.y)));
scale   = 4/max(scale,4);
Y.y     = Y.y*scale;
Y.scale = scale;
if ~isfield(Y,'X0'),Y.X0 = ones(v,1); end
if ~size(Y.X0,2),   Y.X0 = ones(v,1); end
% fMRI slice time sampling
%--------------------------------------------------------------------------
try M.delays = DCM.delays; catch, M.delays = ones(n,1); end
try M.TE     = DCM.TE;     catch, M.TE = 0.04;          end
% prior 
%--------------------------------------------------------------------------
[pE,pC,x]       = spm_dcm_mm_priors(DCM);

% neuronal drive function
%--------------------------------------------------------------------------
input           = spm_dcm_mm_nd(DCM);

% hyperpriors over precision - expectation and covariance
%--------------------------------------------------------------------------
hE       = sparse(n,1) + DCM.options.hE;
hC       = speye(n,n)  * DCM.options.hC;
%--------------------------------------------------------------------------
M.IS = @spm_mm_gen;                       
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
M.Model = Model;
M.nograph = spm('CmdLine');

%-----------------Model Identification methods-----------------------------
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);

% Recomputing the results y
%--------------------------------------------------------------------------
yhat     = feval(M.IS,Ep,M,U);
R        = Y.y - yhat;
R        = R - Y.X0*spm_inv(Y.X0'*Y.X0)*(Y.X0'*R);
Ce       = exp(-Eh);
%Variance that is chaptured by the data
%--------------------------------------------------------------------------

PSS      = sum(yhat.^2);
RSS      = sum(R.^2);
D        = 100.*PSS./(PSS + RSS);
%Kernel
%--------------------------------------------------------------------------
H.f      = @spm_fx_hdm;
H.g      = @spm_gx_hdm;
H.x      = M.x;
H.m      = n; 
[H0,H1]  = spm_kernels(H,Ep.H,64,1/2);
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


% Data ID and log-evidence
%==========================================================================
if isfield(M,'FS')
    try
        ID = spm_data_id(feval(M.FS,Y.y,M));
    catch
        ID = spm_data_id(feval(M.FS,Y.y));
    end
else
    ID     = spm_data_id(Y.y);
end
 
% Save approximations to model evidence: negative free energy, AIC, BIC
%--------------------------------------------------------------------------
evidence   = spm_dcm_evidence(DCM);
DCM.F      = F;
DCM.ID     = ID;
DCM.AIC    = evidence.aic_overall;
DCM.BIC    = evidence.bic_overall;

%-Save DCM
%--------------------------------------------------------------------------
DCM.name = name; 
save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
end




