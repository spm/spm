function [] = spm_dcm_estimate (DCM_filename)   
% Estimate parameters of a DCM
% FORMAT [] = spm_dcm_estimate (DCM_filename)   
%
% DCM_filename  - the DCM model
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_estimate.m 1143 2008-02-07 19:33:33Z spm $

 
if nargin < 1
    
    %-display model details
    %----------------------------------------------------------------------
    Finter = spm_figure('GetWin','Interactive');
    set(Finter,'name','Dynamic Causal Modeling')
    
    %-get DCM
    %----------------------------------------------------------------------
    P = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    
    spm('Pointer','Watch')
    spm('FigName','Estimation in progress');
else
    P = DCM_filename;
end

% load and unpack
%--------------------------------------------------------------------------
load(P)
a  = DCM.a;
b  = DCM.b;
c  = DCM.c;
U  = DCM.U;
Y  = DCM.Y;
xY = DCM.xY;
n  = DCM.n;
v  = DCM.v;
X0 = DCM.Y.X0;
% check whether TE of acquisition has been defined;
% if not, default to 40 ms
if ~isfield(DCM,'TE')
    DCM.TE = 0.04;
    disp('spm_dcm_estimate: TE not defined.')
    disp('DCM.TE was set to default of 0.04 s')
else
    if (DCM.TE < 0) || (DCM.TE > 0.1)
        disp('spm_dcm_estimate: Extreme TE value found.')
        disp('Please check and adjust DCM.TE (note this value must be in seconds!).')
    end
end
TE = DCM.TE;

% priors - expectations
%--------------------------------------------------------------------------
[pE,pC,qE,qC] = spm_dcm_priors(a,b,c);

% model specification
%--------------------------------------------------------------------------
M.IS  = 'spm_int';
M.f   = 'spm_fx_dcm';
M.g   = 'spm_gx_dcm';
M.x   = sparse(n*5,1);
M.pE  = pE;
M.pC  = pC;
M.m   = size(U.u,2);
M.n   = size(M.x,1);
M.l   = n;
M.N   = 32;
M.dt  = 16/M.N;
M.ns  = size(Y.y,1);
M.TE  = TE;

try
    % extension to cover fMRI slice time sampling
    M.delays = DCM.delays;
end

% nonlinear system identification (nlsi)
%--------------------------------------------------------------------------
[Ep,Cp,Ce,H0,H1,H2,M0,M1,L1,L2,F] = spm_nlsi(M,U,Y);

% predicted responses and residuals
%--------------------------------------------------------------------------
y     = spm_int(Ep,M,U);
R     = Y.y - y;
R     = R - X0*inv(X0'*X0)*(X0'*R);

% neuronal kernels
%--------------------------------------------------------------------------
L          = sparse(1:n,[1:n] + 1,1,n,length(M0));
[K0,K1,K2] = spm_kernels(M0,M1,L,M.N,M.dt);

% Bayesian inference and reshape {default threshold T = 0}
%--------------------------------------------------------------------------
T          = 0;
warning off; % switch off NaN-related warning of spm_Ncdf
pp         = 1 - spm_Ncdf(T,abs(Ep),diag(Cp));
warning on;
[ A  B  C] = spm_dcm_reshape(Ep,M.m,n,1);
[pA pB pC] = spm_dcm_reshape(pp,M.m,n,1);

% Also record variances - this helps Bayesian inference eg. across sessions
%--------------------------------------------------------------------------
vv         = diag(Cp);
[vA vB vC] = spm_dcm_reshape(vv,M.m,n,1);

% Store parameters
%-------------------------------------------------------------------
DCM.M      = M;
DCM.Y      = Y;
DCM.U      = U;
DCM.Ep     = Ep;
DCM.Cp     = Cp;
DCM.A      = A;
DCM.B      = B;
DCM.C      = C;
DCM.pA     = pA;
DCM.pB     = pB;
DCM.pC     = pC;
DCM.vA     = vA;
DCM.vB     = vB;
DCM.vC     = vC;
DCM.H1     = H1;
DCM.H2     = H2;
DCM.K1     = K1;
DCM.K2     = K2;
DCM.R      = R;
DCM.y      = y;
DCM.xY     = xY;
DCM.T      = T;
DCM.Ce     = Ce;

% Save approximations to model evidence: negative free energy, AIC, BIC
%--------------------------------------------------------------------------
evidence   = spm_dcm_evidence (DCM);
DCM.F      = F;
DCM.AIC    = evidence.aic_overall;
DCM.BIC    = evidence.bic_overall;


%-Save and reset title
%-------------------------------------------------------------------
if spm_matlab_version_chk('7') >= 0
    save(P,'-V6','DCM');
else
    save(P,'DCM');
end;

if nargin < 1
    spm('Pointer','Arrow');
    spm_input('Thank you',1,'d');
end

return
