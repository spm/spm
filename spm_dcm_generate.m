function [] = spm_dcm_generate(syn_model,source_model,SNR)

% Generate synthetic data from a DCM model specification
% FORMAT [] = spm_dcm_generate(syn_model,source_model,SNR)
% 
% syn_model     Name of synthetic DCM file
% source_model  Type of souce model specification (see spm_dcm_create)
% SNR           Signal to noise ratio (default=1)
%
% This routine will update the DCM.Y field as follows: 
%           Y.y     synthetic BOLD data
%           Y.secs  overall number of seconds
%           Y.Ce    Error covariance
%
% Will Penny & Klaas Enno Stephan $Id$

randn('state',sum(100*clock));

% Check parameters and load specified DCM
%-------------------------------------------------------------------
if (nargin<3) | isempty(SNR)
    SNR         = 1;
end

load(syn_model)

U   = DCM.U;
v   = DCM.v;        % number of scans
n   = DCM.n;        % number of regions
m   = size(U.u,2);  % number of inputs


% Check whether the model is stable by examining the eigenvalue 
% spectrum for the combined connectivity matrix at each time
% point (i.e. scan). Admittedly, this is not very elegant and over-conservative, 
% but very sensitive.
%-------------------------------------------------------------------
stable  = 1;
s       = 2;    % index for current scan
u       = 1;    % index for current input 
TR      = DCM.Y.dt;
warn    = 0;
while s<=v
    while u<=m
        % compute index for current input for current time point (i.e. scan)
        switch upper(source_model)
            case 'GUI'
                u_t             = round(s*TR/U.dt) + 32;    % NB: spm_get_ons uses an offset of 32 bins!
            otherwise   % source model name was imported or specified directly
                u_t             = round(s*TR/U.dt);         % NB: spm_dcm_ui (line 138) has already corrected for the offset!                 
        end
        eigval          = eig(DCM.A + U.u(u_t,u)*DCM.B(:,:,u));
        if max(eigval) > 0 
            warn    = 1;
        end
        u   = u + 1;
    end
    u   = 1;
    s   = s + 1;
end
% display stability warning if necessary
if warn
    disp (['Modelled system is potentially unstable: Lyapunov exponent of combined connectivity matrix is ' num2str(max(eigval))]);
    disp ('Check the output to ensure that values are in a normal range.')
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
%---------------------------------------------------------------------------
P=[1; DCM.A(:); DCM.B(:); DCM.C(:); DCM.H(:)];


% Compute hemodynamic response at v sample points
%---------------------------------------------------------------------------
[y,dy]    = spm_int(P,M,U,v);


% Add noise to all areas
%---------------------------------------------------------------------------
% Compute required r, the standard deviation of additive noise, for all areas
r   = std(y)/SNR;
% Turn r into a diagonal matrix
r   = diag(r);

% Add noise 
p       = 1;
a       = 0;    % AR(1) coeff: for the moment set to zero
a       = [1 -a];
K       = inv(spdiags(ones(v,1)*a,-[0:p],v,v));
K       = K*sqrt(v/trace(K*K'));
z       = randn(v,n);
e       = K*z;
Y       = DCM.Y;
Y.Ce    = spm_Ce(v*ones(1,n));
Y.y     = y + e*r;
Y.secs  = Y.dt*v;

% Now orthogonalise data with respect to effects of no interest
% If X0 is just a vector of 1s this amounts to making the data zero mean
X0  = Y.X0;
Xp  = X0*inv(X0'*X0)*X0';
for i = 1:n,
    Y.y(:,i) = Y.y(:,i)-Xp*Y.y(:,i);
end
DCM.Y = Y;

% Save synthetic DCM
eval (['save ' syn_model ' DCM']);

% Display the generated time series
F = spm_figure('CreateWin','Simulated BOLD time series');
t=Y.dt*[1:1:v];
for i=1:n,
    subplot(n,1,i);
    plot(t,Y.y(:,i));
    title(sprintf('Region %s', Y.name{i}));
    if i<n set(gca,'XTickLabel',[]); end
end
xlabel('secs');

return
