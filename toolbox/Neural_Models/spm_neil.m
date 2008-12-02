% Demo routine for hemodynamic model
%==========================================================================
% For Prof Neil Burgess
% Inst of Cognitive Neuroscience (Deputy Director), and Inst of Neurology
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_neil.m 2521 2008-12-02 19:49:39Z karl $ 
 

% Model specification
%==========================================================================
clear

% create BOLD model
%--------------------------------------------------------------------------
[hE,hC] = spm_hdm_priors(1,5);
hE(end) = 1;
 
% model
%--------------------------------------------------------------------------
clear H
H.f     = 'spm_fx_hdm';
H.g     = 'spm_gx_hdm';
H.x     = [0 0 0 0]';
H.pE    = hE;    
H.m     = 1;
H.n     = 4;
H.l     = 1;
 
% create neuronal responses
%--------------------------------------------------------------------------
dt      = 1/32; % seconds
X       = 1;    % seconds
u(:,1)  = kron([1 0],kron(ones(1,8),kron([1 0 1 0], exp(-[1:(X/dt)]/32) )))';
u(:,2)  = kron([1 0],kron(ones(1,8),kron([1 0 0 0], exp(-[1:(X/dt)]/32) )))';
t       = [1:length(u)]*dt;

 
% Use response to drive a hemodynamic model
%--------------------------------------------------------------------------
U.dt      = dt;
for i = 1:size(u,2)
    u(:,i)    = u(:,i)/mean(u(:,i))*dt*16;
    U.u       = u(:,i);
    BOLD(:,i) = spm_int_L(hE,H,U);
end

subplot(2,2,1)
plot(t,u)
title('Neuronal activity')
axis square
xlabel('time (s)')

subplot(2,2,3)
plot(t,BOLD)
title('BOLD responses')
axis square
xlabel('time (s)')

subplot(2,2,2)
bar(sum(u))
title('summed activity')
axis square

subplot(2,2,4)
bar(sum(BOLD))
title('summed BOLD')
axis square

drawnow

