function [y] = spm_cornsweet(P,M,U)
% generative model for psychophysical responses
% FORMAT [y] = spm_cornsweet(P,M,U)
% P  - model parameters
% M  - model
% %
%  y{1} - matched contrast level for Cornsweet effect
%  y{2} - probability of seeing Mach bands
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_cornsweet.m 4339 2011-06-02 13:04:49Z karl $
 
 
% get variables
%--------------------------------------------------------------------------
sim    = M.sim;
nb     = 6;
Cmax   = max(sim.Con);
thresh = exp(P.thresh);
 
% get coefficients of discrete cosine basis function of contrast from sim
%==========================================================================
 
% map empirical contrast level to simulated contrasts
%--------------------------------------------------------------------------
cm     = 0.0293;
[i j]  = max(sim.Ecs);
Cm     = sim.Con(j);
scale  = cm/exp(P.sig*Cm)*exp(P.off);
 
Con_cs = log(sim.Con_cs/scale)/P.sig;
Con_mb = log(sim.Con_mb/scale)/P.sig;
 
% get coefficients of contrast functions
%--------------------------------------------------------------------------
X      = cos(pi*sim.Con(:)*((1:nb) - 1)/Cmax);
Ecs    = pinv(X)*P.con*sim.Ecs(:);
Emb    = pinv(X)*sim.Emb(:);
Vmb    = pinv(X)*log(sim.Vmb(:));
 
% Cornsweet illusion
%==========================================================================
y{1}   = cos(pi*Con_cs(:)*((1:nb) - 1)/Cmax)*Ecs;
 
% Mach band illusion
%==========================================================================
m      = cos(pi*Con_mb(:)*((1:nb) - 1)/Cmax)*Emb;
v      = cos(pi*Con_mb(:)*((1:nb) - 1)/Cmax)*Vmb;
p      = 1 - spm_Ncdf(thresh,m,exp(v));
y{2}   = (tanh((p - 1/2)*P.sen) + 1)/2;
 
if isempty(U), return, end
 
% plot predictions
%==========================================================================
 
% get empirical contrast
%--------------------------------------------------------------------------
Con    = -2:1/32:16;
 
% Cornsweet at all and simulated contrast levels
%--------------------------------------------------------------------------
css    = cos(pi*sim.Con(:)*((1:nb) - 1)/Cmax)*Ecs;
cs     = cos(pi*    Con(:)*((1:nb) - 1)/Cmax)*Ecs;
 
% Mach Band at all and simulated contrast levels
%--------------------------------------------------------------------------
m      = cos(pi*sim.Con(:)*((1:nb) - 1)/Cmax)*Emb;
v      = cos(pi*sim.Con(:)*((1:nb) - 1)/Cmax)*Vmb;
p      = 1 - spm_Ncdf(thresh,m,exp(v));
prs    = (tanh((p - 1/2)*P.sen) + 1)/2;
 
m      = cos(pi*    Con(:)*((1:nb) - 1)/Cmax)*Emb;
v      = cos(pi*    Con(:)*((1:nb) - 1)/Cmax)*Vmb;
p      = 1 - spm_Ncdf(thresh,m,exp(v));
pr     = (tanh((p - 1/2)*P.sen) + 1)/2;

 
% Graphics
%--------------------------------------------------------------------------
simcon = scale*exp(P.sig*sim.Con);
con    = scale*exp(P.sig*    Con);
p      = 0:1/64:1;
pp     = (tanh((p - 1/2)*P.sen) + 1)/2;
 
subplot(2,2,1)
plot(con,cs,simcon,css,'r.',sim.Con_cs,y{1},'b.',sim.Con_cs,U.y{1},'k.','MarkerSize',16)
xlabel('empirical contrast','Fontsize',12)
ylabel('reported contrast','Fontsize',12)
title('Cornsweet','Fontsize',16)
set(gca,'XLim',[0 .2]);
axis square
 
subplot(2,2,2)
plot(con,pr,simcon,prs,'r.',sim.Con_mb,y{2},'b.',sim.Con_mb,U.y{2},'k.','MarkerSize',16)
xlabel('empirical contrast','Fontsize',12)
ylabel('report probability','Fontsize',12)
title('Mach Bands','Fontsize',16)
set(gca,'XLim',[0 .2]);
axis square
 
subplot(2,2,3)
semilogy(Con,con)
xlabel('simulated contrast','Fontsize',12)
ylabel('empirical contrast','Fontsize',12)
title('Contrast mapping','Fontsize',16)
set(gca,'XLim',[-2 8]);
axis square
 
subplot(2,2,4)
plot(p,pp)
xlabel('posterior confidence','Fontsize',12)
ylabel('response probability','Fontsize',12)
title('Response mapping','Fontsize',16)
axis square


