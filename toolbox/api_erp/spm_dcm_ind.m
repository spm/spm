function DCM = spm_dcm_ind(DCM)   
% Estimate parameters of a DCM model of spectral responses)
% FORMAT DCM = spm_dcm_ind(DCM)   
%
% DCM     
%    name: name string
%       M:  Forward model
%              M.dipfit - leadfield specification
%       xY: data   [1x1 struct]
%       xU: design [1x1 struct]
%
%   Sname: cell of source name strings
%       A: {[nr x nr double]  [nr x nr double]  [nr x nr double]}
%       B: {[nr x nr double], ...}   Connection constraints
%       C: [nr x 1 double]
%
%   options.Nmodes       - number of frequency modes
%   options.Tdcm         - [start end] time window in ms
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.h            - number of DCT drift terms (usually 1 or 2)
%   options.type         - 1 - 'ECD (EEG)'
%                          2 - 'ECD (MEG)'
%                          3 - 'Imaging'
%                          4 - 'LFP' 
%                          (see spm_erp_L)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ind.m 1040 2007-12-21 20:28:30Z karl $


% check options 
%==========================================================================
clear spm_erp_L

% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                  catch, DCM.name           = 'DCM_IND'; end
try, DCM.options.Nmodes;        catch, DCM.options.Nmodes = 4;         end
try, h     = DCM.options.h;     catch, h                  = 1;         end
try, onset = DCM.options.onset; catch, onset              = 80;        end

% Data and spatial model
%==========================================================================
DCM    = spm_dcm_erp_dipfit(DCM);
DCM    = spm_dcm_ind_data(DCM);
xY     = DCM.xY;
try
    xU = DCM.U;
catch
    xU = DCM.xU;
end

% dimensions
%--------------------------------------------------------------------------
Nt     = length(xY.xy);                 % number of trials
Nr     = size(xY.xf{1},2);              % number of sources
Ns     = size(xY.xf{1},1);              % number of samples
Nf     = size(xY.xf,2);                 % number of frequency modes
nu     = size(xU.X,2);                  % number of inputs
nx     = Nr*Nf + 1;                     % number of states


% assume noise variance is the same over modes:
%--------------------------------------------------------------------------
a      = 1/5;
xY.Q   = {spm_Q(a,Ns,1)};
xY.X0  = sparse(Ns*Nt,0);

% Inputs
%==========================================================================

% trial effects
%--------------------------------------------------------------------------
try
    if size(xU.X,2) - length(DCM.B)
        warndlg({'please ensure number of trial specific effects', ...
                 'encoded by DCM.xU.X & DCM.B are the same'})
    end
catch
    DCM.B = {};
end
try
    xU.u  = kron(xU.X,ones(Ns,1));
catch
    xU.u  = sparse(Nt*Ns,0);
end

% stimulus parameters
%--------------------------------------------------------------------------
xU.dt  = xY.dt;
xU.dur = xU.dt*(Ns - 1);

% model specification and nonlinear system identification
%==========================================================================
M      = DCM.M;
try, M = rmfield(M,'g');  end
try, M = rmfield(M,'FS'); end

% adjust onset relative to pst
%--------------------------------------------------------------------------
dur    = xU.dur;
ons    = onset - xY.pst(1);

if ons < 0; warndlg('onset time is negative; please increase'); end

% prior moments
%--------------------------------------------------------------------------
A      = DCM.A;
B      = DCM.B;
C      = kron(ones(1,length(ons)),DCM.C);

[pE,gE,pC,gC] = spm_ind_priors(A,B,C,Nf);


% likelihood model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_ind';
M.G   = 'spm_lx_ind';
M.IS  = 'spm_int_U';
M.fu  = 'spm_ind_u';
M.x   = sparse(nx,1);
M.pE  = pE;
M.pC  = pC;
M.gE  = gE;
M.gC  = gC;
M.m   = nu;
M.n   = nx;
M.l   = Nr*Nf;
M.ns  = Ns*Nt;

% and fixed parameters and functional forms
%--------------------------------------------------------------------------
M.ons = ons;
M.dur = dur;

% EM: inversion
%--------------------------------------------------------------------------
[Qp,Qg,Cp,Cg,Ce,F] = spm_nlsi_N(M,xU,xY);

% Bayesian inference {threshold = prior} NB Prior on A,B  and C = exp(0) = 1
%==========================================================================
warning off
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning on

% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
L   = feval(M.G, Qg,M);           % get gain matrix
x   = feval(M.IS,Qp,M,xU);        % prediction (source space)
y   = x*L';                       % prediction (sensor space)
r   = xY.y - y;                   % prediction error
x   = x(:,2:end);                 % remove time

% trial specific responses (in mode, channel and source space)
%--------------------------------------------------------------------------
for i = 1:Nt
    for j = 1:Nf
        t = [1:Ns] + (i - 1)*Ns;
        f = [1:Nr] + (j - 1)*Nr;
        k = [1:Nr] + (j - 1)*Nr;
        
        Hc{i,j} = y(t,f);
        Ec{i,j} = r(t,f);
        K{i,j}  = x(t,k);
    end
end

% store estimates in DCM
%--------------------------------------------------------------------------
DCM.M  = M;                    % model specification
DCM.xY = xY;                   % data structure
DCM.xU = xU;                   % input structure
DCM.Ep = Qp;                   % conditional expectation f(x,u,p)
DCM.Cp = Cp;                   % conditional covariances G(g)
DCM.Eg = Qg;                   % conditional expectation
DCM.Cg = Cg;                   % conditional covariances
DCM.Pp = Pp;                   % conditional probability
DCM.Hc = Hc;                   % conditional responses (y), channel space
DCM.K  = K;                    % conditional responses (x)
DCM.Rc = Ec;                   % conditional residuals (y), channel space
DCM.Ce = Ce;                   % ReML error covariance
DCM.F  = F;                    % Laplace log evidence


% and save
%--------------------------------------------------------------------------
if spm_matlab_version_chk('7.1') >= 0
    save(DCM.name, '-V6', 'DCM');
else
    save(DCM.name, 'DCM');
end
assignin('base','DCM',DCM)
return
