function [F,p] = spm_DEM_F(DEM,ip)
% Free-energy as a function of conditional parameters
% [F,P1,P2] = spm_DEM_F(DEM))
%
% DEM    - hierarchical model
%
% F(i) - free-energy at <q(P(ip))> = p(i)
%
% where p(i) is the ip-th free-parameter. This is a bound on 
% the log-likehood (log-evidence) conditioned on the expeceted parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_F.m 4146 2010-12-23 21:01:39Z karl $
 
 
% Find paramter ranges (using prior covariance)
%--------------------------------------------------------------------------
pE   = spm_vec(DEM.M(1).pE);
p    = linspace(-6,6,16);
dp   = sqrt(DEM.M(1).pC(ip,ip))*p;
p    = dp + pE(ip);
 
% get F
%==========================================================================
DEM.M(1).E.nE = 1;
DEM.M(1).E.nN = 1;

for i = 1:length(p)
    
    % adjust paramter (through prio expecatation)
    %----------------------------------------------------------------------
    P          = pE;
    P(ip)      = p(i);
    DEM.M(1).P = spm_unvec(P,DEM.M(1).pE);
    
    % comute free-energy
    %----------------------------------------------------------------------
    DEM  = spm_DEM(DEM);
    F(i) = DEM.F(end);
    
end


% predicted F under the Laplace assumption
%==========================================================================
DEM.M(1).P = DEM.M(1).pE;
DEM        = spm_DEM(DEM);

% compute free-energy
%--------------------------------------------------------------------------
dFdp   = DEM.qP.dFdp(ip);
dFdpp  = DEM.qP.dFdpp(ip,ip);
FP     = dFdp*dp + (dFdpp*dp.^2)/2;
FP     = FP - max(FP);
    

% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Free-energy');
F    = F - max(F);
subplot(2,1,1)
plot(p,F,p,FP,':'), hold on
plot(pE(ip)*[1 1],[min(F) 0],':'),  hold on
xlabel('Parameter expecatation','FontSize',12)
ylabel('Free-energy','FontSize',12)
title('Free-energy profile','FontSize',16)
axis square, box off

