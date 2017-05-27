function [Ep,Cp] = spm_dcm_sparse(DCM,field)
% Bayesian model reduction of all permutations of model parameters
% FORMAT [RCM,BMR] = spm_dcm_sparse(DCM,field
%
% DCM      - A single estimated DCM (or PEB) structure:
%
%  DCM.M.pE  - prior expectation
%  DCM.M.pC  - prior covariance
%  DCM.Ep    - posterior expectation
%  DCM.Cp    - posterior covariances
%  DCM.gamma - prior variance    of reduced parameters (default: 0)
%
% field      - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%             'All' will invoke all fields (i.e. random effects)
%             If Ep is not a structure, all parameters will be considered
%
% Returns:
%  Ep    - (BMA) posterior expectation
%  Cp    - (BMA) posterior covariance
%
%--------------------------------------------------------------------------
% This routine searches over reduced (nested) models of a full model (DCM)
% using Bayesian model reduction and performs Bayesian Model Averaging.
% 'Reduced' means some free parameters (parameters with a non-
% zero prior covariance) are switched off by fixing their prior variance
% to zero.This version incorporates a sparsity  prior over models (with a
% Gaussian hyperprior). In other words, the free energy is taken to be the
% likelihood of some data under a given model. The prior on that model
% corresponds to a softmax function of the prior entropy. Finally, the
% softmax (Gibbs) parameter is equipped with a Gaussian prior. Using
% Bayesian model reduction, this routine evaluates the joint probability
% over model and softmax sparsity parameter. The marginals over model space
% are then used to form Bayesian model averaging.
%
% The greedy search in this version simply evaluates the log evidence of
% models with and without each parameter and then successively removes the
% parameters with the least evidence.
%
% See also: spm_dcm_bmr and spm_dcm_bmr_all
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: spm_dcm_sparse.m 7082 2017-05-27 19:36:36Z karl $


%-Number of parameters to consider before invoking greedy search
%--------------------------------------------------------------------------
nmax  = 8;

%-specification of null prior covariance
%--------------------------------------------------------------------------
if isfield(DCM,'beta'),  beta  = DCM.beta;  else, beta  = 0; end
if isfield(DCM,'gamma'), gamma = DCM.gamma; else, gamma = 0; end

%-Check fields of parameter stucture
%--------------------------------------------------------------------------
if nargin < 2 || isempty(field)
    field = {'A','B'};
end
if ischar(field)
    field = {field};
end

%-dela with filenames stucture
%--------------------------------------------------------------------------
if ischar(DCM)
    DCM = load(DCM,'DCM');
    DCM = DCM.DCM;
end

% Get prior covariances
%--------------------------------------------------------------------------
if isstruct(DCM.M.pC), DCM.M.pC = diag(spm_vec(DCM.M.pC)); end
if spm_length(DCM.M.pE) ~= size(DCM.M.pC,1)
    DCM.M.pC = diag(spm_vec(DCM.M.pC));
end

% Get priors and posteriors
%--------------------------------------------------------------------------
qE  = DCM.Ep;
qC  = DCM.Cp;
pE  = DCM.M.pE;
pC  = DCM.M.pC;

% Remove (a priori) null space
%--------------------------------------------------------------------------
U   = spm_svd(pC);
qE  = U'*spm_vec(qE);
pE  = U'*spm_vec(pE);
qC  = U'*qC*U;
pC  = U'*pC*U;


%-Greedy search (GS) - eliminating parameters in a top down fashion
%==========================================================================

% Accumulated reduction vector (C)
%--------------------------------------------------------------------------
q   = diag(DCM.M.pC);
if sum(q < 1024)
    C   = double(q > mean(q(q < 1024))/1024);
else
    C   = double(q > 0);
end

%-Find free coupling parameters
%----------------------------------------------------------------------
if isstruct(DCM.Ep)
    k = spm_fieldindices(DCM.Ep,field{:});
else
    k = 1:spm_length(DCM.Ep);
end
k     = k(find(C(k))); %#ok<FNDSB>

% Model search over new prior without the i-th parameter
%------------------------------------------------------------------
nparam = length(k);
for i  = 1:nparam
    
    % Identify parameters to retain r and to remove s
    %--------------------------------------------------------------
    r    = C; r(k(i)) = 0; s = 1 - r;
    
    % Create reduced priors
    %--------------------------------------------------------------
    R    = U'*diag(r + s*gamma)*U;
    rC   = R*pC*R;
    F(i) = spm_log_evidence(qE,qC,pE,pC,pE,rC);
        
end

% Find parameters with the least evidence
%--------------------------------------------------------------------------
[F,i] = sort(-F);
k     = k(i);
M     = cell(0);
for i = 1:nparam
    

    % parameters to retain (r) and to remove (s)
    %----------------------------------------------------------------------
    r    = C; r(k(1:i)) = 0; s = 1 - r;
    
    % Create reduced prior covariance matrix
    %----------------------------------------------------------------------
    R    = U'*diag(r + s*gamma)*U;
    rC   = R*pC*R;
    
    % record
    %----------------------------------------------------------------------
    M(i).F  = spm_log_evidence(qE,qC,pE,pC,pE,rC)
    M(i).H  = spm_logdet(rC);
    M(i).rC = rC;
    
end

% Sparsity hyperpriors
%--------------------------------------------------------------------------
s     = (1:64)/64;
Ps    = exp(-((1:64) - 32).^2/(2*16));
Ps    = Ps/sum(Ps);

%  model likelihood, model prior and  sparsity hyperprior
%--------------------------------------------------------------------------
Lm    = spm_softmax(spm_vec(M.F));
for i = 1:numel(s)
    Pm      = spm_softmax(-s(i)*spm_vec(M.H));
    Qm(:,i) = Lm.*Pm*Ps(i);
end

% evidence and log evidence
%--------------------------------------------------------------------------
Qm    = Qm/sum(sum(Qm));
G     = log(sum(Qm,2));

%-Bayesian model average
%==========================================================================
qE    = DCM.Ep;
qC    = DCM.Cp;
pE    = DCM.M.pE;
pC    = DCM.M.pC;
pE    = spm_vec(pE);
Gmax  = max(G);
BMA   = {};
for i = 1:length(G)
    if G(i) > (Gmax - 4)    
        [F,Ep,Cp]    = spm_log_evidence_reduce(qE,qC,pE,pC,pE,M(i).rC);
        BMA{end + 1} = struct('Ep',Ep,'Cp',Cp,'F',F);
    end
end

BMA   = spm_dcm_bma(BMA);
Ep    = BMA.Ep;
Cp    = BMA.Cp;
if isstruct(Cp) || (spm_length(Cp) == spm_length(Ep))
    Cp = diag(spm_vec(Cp));
end

% Show results
% -------------------------------------------------------------------------
GRAPHICS = 1;

if ~GRAPHICS, return, end

spm_figure('Getwin','BMR - all'); clf
subplot(3,2,1)
imagesc(log(Qm)')
title('Joint density','FontSize',16)
xlabel('model','FontSize',12)
ylabel('sparsity','FontSize',12)
axis square

subplot(3,2,2)
plot(DCM.Ep,Ep,'.',DCM.Ep,DCM.Ep,':')
title('Full and reduced expectations','FontSize',16)
xlabel('Full expectations','FontSize',12)
ylabel('Reduced expectations','FontSize',12)
axis square

subplot(3,2,3)
plot(s,sum(Qm,1))
title('Marginal over sparsity','FontSize',16)
xlabel('sparsity parameter','FontSize',12)
ylabel('probability','FontSize',12)
axis square

subplot(3,2,4)
plot(sum(Qm,2))
title('Marginal over  model','FontSize',16)
xlabel(' model','FontSize',12)
ylabel('probability','FontSize',12)
axis square
drawnow



