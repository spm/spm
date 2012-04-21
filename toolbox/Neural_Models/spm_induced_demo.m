% Demo routine for local field potential models
%==========================================================================
%
% This is a generic demonstration of neural-mass models that illustrates
% various impulse response behaviours. It is meant to show how to specify
% a neural-mass model, examine its response properties using Volterra
% kernels and transfer functions and generate electrophysiological and
% hemodynamic responses from the same model. It is anticipated that people
% will go through the code to see how the routines relate to each other.
%
% This demo contains a linear stability analysis, which can be useful for
% identifying useful domains of parameter space (here the inhibitory time-
% constant)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_induced_demo.m 4721 2012-04-21 08:53:29Z karl $


% Model specification
%==========================================================================

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = 'CMC';
options.analysis = 'TFA';
dipfit.model = options.model;
dipfit.type  = options.spatial;
dipfit.Nc    = Nc;
dipfit.Ns    = Ns;


% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_dcm_neural_priors({0 0 0},{},1,options.model);
P       = fieldnames(pE);
[pE,pC] = spm_L_priors(dipfit,pE,pC);
[pE,pC] = spm_ssr_priors(pE,pC);
[x,f]   = spm_dcm_x_neural(pE,options.model);

% eliminate channel noise and make innovations white
%--------------------------------------------------------------------------
pE.a    = [  0; -16];                  % log amplitude ang f^(-a) exponent
pE.b    = [-32; -32];                  % log amplitude ang f^(-a) exponent
pE.c    = [-32; -32];                  % log amplitude ang f^(-a) exponent


% orders and model
%==========================================================================
np      = length(spm_vec(pE));
nx      = length(spm_vec(x ));
nu      = size(pE.C,2);
u       = sparse(1,nu);

% create LFP model
%--------------------------------------------------------------------------
M.f     = f;
M.g     = 'spm_gx_erp';
M.h     = f;
M.x     = x;
M.n     = nx;
M.pE    = pE;
M.m     = nu;
M.l     = Nc;

% Volterra Kernels and transfer functions
%==========================================================================
M.u     = u;
M.Hz    = 4:68;

% compute transfer functions for different parameters
%--------------------------------------------------------------------------
iplot = 1;
ifig  = 1;
D     = 1.5;
for k = 1:0 %length(P)
    
    
    % check paremter exists
    %----------------------------------------------------------------------
    spm_figure('GetWin',sprintf('Panel %i',ifig));
    Q = getfield(pE,P{k});
    if isnumeric(Q)
        for i = 1:size(Q,1)
            for j = 1:size(Q,2);
                
                % line search
                %----------------------------------------------------------
                dQ    = linspace(Q(i,j) - D,Q(i,j) + D,32);
                for q = 1:length(dQ)
                    qE      = pE;
                    qE      = setfield(qE,P{k},{i,j},dQ(q));
                    [G w]   = spm_csd_mtf(qE,M);
                    GW(:,q) = G{1};
                end

                if any(var(GW') > 1e-6)
                    
                    % plot
                    %------------------------------------------------------
                    subplot(4,2,2*iplot - 1)
                    plot(w,GW)
                    xlabel('frequency {Hz}')
                    title(sprintf('Param: %s(%i,%i)',P{k},i,j),'FontSize',16)
                   
                    
                    subplot(4,2,2*iplot - 0)
                    imagesc(dQ,w,log(GW))
                    title('transfer functions','FontSize',16)
                    ylabel('Frequency')
                    xlabel('Inhibitory connection','FontSize',16)
                    axis xy; drawnow
                    
                    % update graphics
                    %------------------------------------------------------
                    iplot = iplot + 1;
                    if iplot > 4
                        iplot = 1;
                        ifig  = ifig + 1;
                        spm_figure('GetWin',sprintf('Panel %i',ifig));
                    end
                    
                end
            end
        end
    end
end


% exogenous input-dependent parameters
%==========================================================================
M.f     = 'spm_fx_tfm';
pE.X    = sparse(np,nu);
pC.X    = sparse(np,nu);

% state-dependent parameters
%--------------------------------------------------------------------------
ix      = 3;
ip      = spm_fieldindices(pE,'G');jp = 2;  % gamma (High - 50 Hz)
pE.Y    = sparse(ip(jp),ix,64,np,nx);

ix      = 1;
ip      = spm_fieldindices(pE,'G');jp = 3;  % Beta (24 Hz)
pE.Y    = pE.Y + sparse(ip(jp),ix,16,np,nx);

ix      = 7;
ip      = spm_fieldindices(pE,'G');jp = 4;  % Gamms (low - 38 Hz)
pE.Y    = pE.Y + sparse(ip(jp),ix,-64,np,nx);


pC.Y    = sparse(np,nx);

% Integrate system to see response (time-frequency)
%==========================================================================

% remove M.u to invoke exogenous inputs
%--------------------------------------------------------------------------
M     = rmfield(M,'u');
N     = 128;
U.dt  = 4/1000;
t     = (1:N)'*U.dt;
U.u   = sparse(N,M.m);

% exogenous input
%--------------------------------------------------------------------------
U.u(:,1)  = exp(-(t - 128/1000).^2*1024)*32;
U.u(:,1)  = spm_conv((t > 128/1000 & t < 256/1000)*32,4);

% now integrate a generative model to simulate a time frequency response
%==========================================================================
[y,w,t,x] = spm_csd_tfm(pE,M,U);

% basline correction
%--------------------------------------------------------------------------
for i = 1:length(y)
    y{i} = log(y{i});
    y{i} = y{i} - ones(size(y{i},1),1)*y{i}(1,:);
end

% plot (in ms)
%--------------------------------------------------------------------------
spm_figure('GetWin','Simulated time frequency responses');

t     = t*1000;

subplot(4,1,1)
plot(t,U.u)
xlabel('time (ms)')
title('Exogenous input','FontSize',16)
spm_axis tight

% LFP – expectation
%--------------------------------------------------------------------------
subplot(4,1,2)
plot(t,x(:,1:2:8))
xlabel('time (ms)')
title('Hidden neuronal states','FontSize',16)
spm_axis tight

% predicted time frequency response
%--------------------------------------------------------------------------
subplot(4,1,3)
imagesc(t,w,y{1}');
title('Time-frequency response (baseline corrected)','FontSize',16)
axis  xy
xlabel('time (ms)')
ylabel('Hz')


% predicted time frequency response
%--------------------------------------------------------------------------
subplot(4,1,4)
plot(t,y{1}');
title('Time-frequency response (baseline corrected)','FontSize',16)
xlabel('time (ms)')
ylabel('Hz')


return


