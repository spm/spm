function [] = spm_dcm_generate(DCM_filename,SNR,show_data)
% Generate data from a DCM model
% FORMAT [] = spm_dcm_generate(DCM_filename,SNR,show_data)
% 
% DCM_filename  Name of DCM data structure
% SNR           Signal to noise ratio (default=1)
% show_data     plot data 0 (no) or 1 (yes) (default=0)
%
% This routine will update the DCM.Y field as follows: 
%           Y.y     synthetic BOLD data
%           Y.secs  overall number of seconds
%           Y.Ce    Error covariance
%
% @(#)spm_dcm_generate.m	2.1 Will Penny 03/10/16


% Check parameters and load specified DCM
%-------------------------------------------------------------------
if (nargin<2) | isempty(SNR)
    SNR         = 1;
end
if (nargin<3) | isempty(show_data)
    show_data   = 0;
end

FP{1} = DCM_filename;
load(FP{:})

U   = DCM.U;
v   = DCM.v;        % number of scans
n   = DCM.n;        % number of regions
m   = size(U.u,2);  % number of inputs


% Check whether the model is stable by examining the eigenvalue 
% spectrum for the combined connectivity matrix at each time
% point (i.e. scan) - this is not very elegant, but takes into account all
% possible temporal relations of the inputs
%-------------------------------------------------------------------
stable  = 1;
s       = 2;    % index for current scan
u       = 1;    % index for current input 
TR      = DCM.Y.dt;
while (stable & s<=v)
    while (stable & u<=m)
        u_t             = round(s*TR/U.dt) + 32;    % value of current input for current time point (i.e. scan)
                                                    % NB: spm_get_ons uses an offset of 32 bins!
        eigval          = eig(DCM.A + U.u(u_t,u)*DCM.B(:,:,u));
        if max(eigval) > 0 
            disp (['Modelled system is unstable: Lyapunov exponent of combined connectivity matrix is ' num2str(max(eigval))]);
            disp ('No synthetic data generated.')
            return
        end
        u   = u + 1;
    end
    u   = 1;
    s   = s + 1;
end


% Create M matrix for spm_int
%-------------------------------------------------------------------
M.f  = 'spm_fx_dcm';
M.g  = 'spm_lx_dcm';
M.x   = sparse(n*5,1);
M.m   = size(U.u,2);
M.n   = size(M.x,1);
M.l   = n;

% Create P vector for spm_int
P=[1; DCM.A(:); DCM.B(:); DCM.C(:); DCM.H(:)];

%---------------------------------------------------------------------------
randn('state',sum(100*clock));

% Compute hemodynamic response at v sample points
%---------------------------------------------------------------------------

[y,dy]    = spm_int(P,M,U,v);


% Compute required r, the standard deviation of additive noise
r=std(y(:,1))/SNR;

p     = 1;
% AR(1) coeff
a     = 0;

a     = [1 -a];
K     = inv(spdiags(ones(v,1)*a,-[0:p],v,v));
K     = K*sqrt(v/trace(K*K'));
z     = randn(v,n);
e     = K*z;
Y     = DCM.Y;
Y.Ce  = spm_Ce(v*ones(1,n));
Y.y   = y + e*r;
Y.secs = Y.dt*v;

% Now orthogonalise data with respect to effects of no interest
% If X0 is just a vector of 1s this amounts to making the data zero mean
X0=Y.X0;
Xp=X0*inv(X0'*X0)*X0';
for i=1:n,
    Y.y(:,i)=Y.y(:,i)-Xp*Y.y(:,i);
end
DCM.Y=Y;
save(FP{:},'DCM');

if show_data,
    F = spm_figure('CreateWin','Simulated BOLD time series');
    t=Y.dt*[1:1:v];
    for i=1:n,
        subplot(n,1,i);
        plot(t,Y.y(:,i));
        title(sprintf('Region %s', Y.name{i}));
        if i<n set(gca,'XTickLabel',[]); end
    end
    xlabel('secs');
end

return
