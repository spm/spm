function P = spm_eeg_inv_vbecd(P)
% Model inversion routine for ECDs using variational Bayesian approach
%
% FORMAT P = spm_eeg_inv_vbecd(P)
%
% Input:
% structure P with fields:
%  forward      - structure containing the forward model, i.e. the "vol"
%                 and "sens" structure in a FT compatible format
%  bad          - list of bad channels, not to use.
%  y            - data vector
%
%  Niter        - maximum number of iterations
%  priors       - priors on parameters,  as filled in (and
%                 described) in spm_eeg_inv_vbecd_gui.m.
%
% Output:
% same structure with extra fields
%  init         - initial valuse used for mu_w/s
%  dF           - successive (relative) improvement of F
%  post         - posterior value of estimated parameters and ther variance
%  Fi           - successive values of F
%  F            - Free energy final value.
%
% Reference:
% Kiebel et al., Variational Bayesian inversion of the equivalent current
% dipole model in EEG/MEG., NeuroImage, 39:728-741, 2008
% (Although this algorithm uses a function for general Bayesian inversion of
% a non-linear model - see spm_nlsi_gn)
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_inv_vbecd.m 4370 2011-06-20 14:51:19Z gareth $



% unpack model, priors, data
%--------------------------------------------------------------------------
Nd    = length(P.priors.mu_w0)/3;
mu_w0 = P.priors.mu_w0;
mu_s0 = P.priors.mu_s0;
S_w0  = P.priors.S_w0;
S_s0  = P.priors.S_s0;

% subtracting mean level from eeg data
%--------------------------------------------------------------------------
if strcmp(upper(P.modality),'EEG')
    P.y = P.y-mean(P.y);
end

% rescale data to fit into minimisation routine (same magnitude for EEG and
% MEG means same threshold and exit criteria)
%--------------------------------------------------------------------------

y    = P.y;
sc_y = 1/std(y);
y    = y*sc_y;
Y.y  = y;

U.u  = 1;
if isfield(P,'chanCov'),
    Y.Q{1}=P.chanCov*sc_y*sc_y;
end;

outsideflag=1;
while outsideflag==1, %% don't use sources which end up outside the head
    
    % set random moment, scaled by prior variances
    %----------------------------------------------------------------------
    [u,s,v] = svd(S_w0);
    mu_w   = mu_w0 + u*diag(sqrt(diag(s+eps)))*v'*randn(size(mu_w0));
    
    % a random guess for the location, based on variance of the prior
    %----------------------------------------------------------------------
    [u,s,v] = svd(S_s0);
    outside = 1;
    while(outside)
        outside = 0;
        mu_s = mu_s0 + u*diag(sqrt(diag(s+eps)))*v'*randn(size(mu_s0)); %
        for i=1:3:length(mu_s), %% check each dipole is inside the head
            pos     = mu_s(i:i+2);
            outside = outside+ ~ft_inside_vol(pos',P.forward.vol);
        end;
    end;
    
    % get lead fields
    %----------------------------------------------------------------------
    M.pE  = [mu_s;mu_w];            % prior parameter estimate
    
    
    
    M.pC  = blkdiag(S_s0,S_w0);     % prior covariance estimate
    M.hE=P.priors.hE;
    
    M.hC=P.priors.hC;
    
    
    M.IS  = 'spm_eeg_wrap_dipfit_vbecd';
    startguess=M.pE;
    M.Setup =P;             % pass volume conductor and sensor locations on
    M.sc_y =sc_y;           % pass on scaling factor
    
    %% startguess=[-0.3553  -69.8440    1.0484    0.2545    0.3428    1.8526]'
    
    [starty]=spm_eeg_wrap_dipfit_vbecd(startguess,M,U);
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);
    P.Ep = Ep;
    P.Cp = Cp;
    P.Eh = Eh;
    P.F  = F;
    [P.ypost,outsideflag,leads]=spm_eeg_wrap_dipfit_vbecd(P.Ep,M,U);
    P.ypost = P.ypost./sc_y; %%  scale it back
    
    if outsideflag
        disp('running again, one or more dipoles outside head.');
    end;
    
end; % while
P.post_mu_s = Ep(1:length(mu_s));
P.post_mu_w = Ep(length(mu_s)+1:end);
P.post_S_s  = Cp(1:length(mu_s),1:length(mu_s));
P.post_S_w  = Cp(length(mu_s)+1:end,length(mu_s)+1:end);


%% return a weight matrix to map channels to dipoles
fulldipmom=[];
fulllf=[];

for d=1:Nd,
    dippos=P.post_mu_s((d-1)*3+1:d*3);
    dipmom=P.post_mu_w((d-1)*3+1:d*3);
     
    fulllf=[fulllf squeeze(leads(d,:,:))']; % *unitmom(d,:)';
    fulldipmom=[fulldipmom ;dipmom];
end;
estdata=fulllf*fulldipmom/1000;

D1=fulldipmom*fulldipmom';
P.post_wdale=D1*fulllf'*pinv(fulllf*D1*fulllf')*1000; % 1000 to put in nAm
estdipmom=P.post_wdale*estdata; %% re-estimate moment using the dale operator
%disp(sprintf('reistimation error %3.2f percent',mean((estdipmom-fulldipmom)./fulldipmom)*100));





