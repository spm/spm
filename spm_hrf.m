function [hrf] = spm_hrf(RT);
% returns a hemodynamic response function
% FORMAT [hrf] = spm_hrf(RT);
% RT   - scan repeat time
% hrf  - hemodynamic response function
%_______________________________________________________________________
% %W% Karl Friston %E%

% modelled hemodynamic response function - {mixture of Gammas}
%-----------------------------------------------------------------------
dt    = RT/8;
u     = 0:(32/dt);
hrf   = spm_Gpdf(u,[6 dt]) - spm_Gpdf(u,[16 dt])/6;
hrf   = spm_conv(hrf,8);
hrf   = hrf([1:32/RT]*8);
hrf   = hrf/sum(hrf);
