function [] = spm_dcm_estimate (DCM_filename)   
% Estimate parameters of a DCM model
% FORMAT [] = spm_dcm_estimate (DCM_filename)   
%
% DCM_filename  - the DCM model
%
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id$

 
if nargin < 1
    %-display model details
    Finter = spm_figure('GetWin','Interactive');
    %-------------------------------------------------------------------
    set(Finter,'name','Dynamic Causal Modeling')
    %-get results
    %-------------------------------------------------------------------
    P     = spm_get(1,'DCM*.mat',{'select DCM_???.mat'});
else
    P{1} = DCM_filename;
end

load(P{:})
if nargin < 1
    spm('Pointer','Watch')
    spm('FigName','Estimation in progress');
end

a=DCM.a;
b=DCM.b;
c=DCM.c;
U=DCM.U;
Y=DCM.Y;
xY=DCM.xY;
n=DCM.n;
v=DCM.v;
X0=DCM.Y.X0;


% priors - expectations
%-------------------------------------------------------------------
[pE,pC,qE,qC] = spm_dcm_priors(a,b,c);

% model specification and nonlinear system identification
%-------------------------------------------------------------------
M.f   = 'spm_fx_dcm';
M.g   = 'spm_lx_dcm';
M.x   = sparse(n*5,1);
M.pE  = pE;
M.pC  = pC;
M.m   = size(U.u,2);
M.n   = size(M.x,1);
M.l   = n;
M.N   = 32;
M.dt  = 16/M.N;
M.maxits=32;

[Ep,Cp,Ce,H0,H1,H2,M0,M1,L] = spm_nlsi(M,U,Y);

% predicted responses and residuals
%-------------------------------------------------------------------
y     = spm_int(Ep,M,U,v);
R     = Y.y - y;
R     = R - X0*inv(X0'*X0)*(X0'*R);

% neuronal kernels
%-------------------------------------------------------------------
L          = sparse(1:n,[1:n] + 1,1,n,length(M0));
[K0,K1,K2] = spm_kernels(M0,M1,L,M.N,M.dt);

% Bayesian inference and reshape {threshold T1/2 = log(2)/T}
%-------------------------------------------------------------------
T          = log(2)/4;			
pp         = 1 - spm_Ncdf(T,abs(Ep),diag(Cp));
[ A  B  C] = spm_dcm_reshape(Ep,M.m,n,1);
[pA pB pC] = spm_dcm_reshape(pp,M.m,n,1);

% Also record variances - this helps in doing Bayesian inference eg. across sessions
vv         = diag(Cp);
[vA vB vC] = spm_dcm_reshape(vv,M.m,n,1);


% Record parameters
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

%-Save and reset title
%-------------------------------------------------------------------
if str2num(version('-release'))>=14,
    save(P{:},'-V6','DCM');
else
    save(P{:},'DCM');
end;

if nargin < 1
    spm('Pointer','Arrow');
    spm_input('Thank you',1,'d');
end
