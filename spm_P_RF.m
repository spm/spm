function [P,p,Ec,Ek] = spm_P_RF(c,k,Z,df,STAT,R,n)
% Returns the [un]corrected P value using unifed EC theory
% FORMAT [P p Ec Ek] = spm_P_RF(c,k,z,df,STAT,R,n)
%
% c     - cluster number 
% k     - extent {RESELS}
% z     - height {minimum over n values}
% df    - [df{interest} df{error}]
% STAT  - Statistical field
%       'Z' - Gaussian field
%       'T' - T - field
%       'X' - Chi squared field
%       'F' - F - field
% R     - RESEL Count {defining search volume}
% n     - number of component SPMs in conjunction
%
% P     - corrected   P value  - P(C >= c | K >= k}
% p     - uncorrected P value
% Ec    - expected number of clusters (maxima)
% Ek    - expected number of resels per cluster
%
%__________________________________________________________________________
%
% spm_P_RF returns the probability of c or more clusters with more than
% k resels in volume process of R RESELS thresholded at u.  All p values
% can be considered special cases:
%
% spm_P_RF(1,0,z,df,STAT,1,n) = uncorrected p value
% spm_P_RF(1,0,z,df,STAT,R,n) = corrected p value {based on height z)
% spm_P_RF(1,k,u,df,STAT,R,n) = corrected p value {based on extent k at u)
% spm_P_RF(c,k,u,df,STAT,R,n) = corrected p value {based on number c at k and u)
% spm_P_RF(c,0,u,df,STAT,R,n) = omnibus   p value {based on number c at u)
%
% If n > 1 a conjunction probility over the n values of the statistic
% is returned.
%__________________________________________________________________________
%
% References:
% 
% [1] Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% [2] Friston KJ et al (1994) Assessing the Significance of Focal Activations
% Using Their Spatial Extent
% Human Brain Mapping 1:210-220
% [3] Worsley KJ et al (1996) A Unified Statistical Approach for Determining
% Significant Signals in Images of Cerebral Activation
% Human Brain Mapping 4:58-73
%__________________________________________________________________________
% Copyright (C) 1999-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_P_RF.m 5770 2013-11-27 20:12:29Z karl $


%-Get expectations
%==========================================================================

% get EC densities
%--------------------------------------------------------------------------
D   = find(R,1,'last');
R   = R(1:D);
G   = sqrt(pi)./gamma((1:D)/2);
EC  = spm_ECdensity(STAT,Z,df);
EC  = max(EC(1:D),eps);

% corrected p value
%--------------------------------------------------------------------------
P   = triu(toeplitz(EC'.*G))^n;
P   = P(1,:);
EM  = (R./G).*P;        % <maxima> over D dimensions
Ec  = sum(EM);          % <maxima>
EN  = P(1)*R(D);        % <resels>
Ek  = EN/EM(D);         % Ek = EN/EM(D);


%-Get P{n > k}
%==========================================================================

% assume a Gaussian form for P{n > k} ~ exp(-beta*k^(2/D))
% Appropriate for SPM{Z} and high d.f. SPM{T}
%--------------------------------------------------------------------------
D  = D - 1;
if ~k || ~D

    p    = 1;

elseif STAT == 'Z'

    beta = (gamma(D/2 + 1)/Ek)^(2/D);
    p    = exp(-beta*(k^(2/D)));

elseif STAT == 'T'

    beta = (gamma(D/2 + 1)/Ek)^(2/D);
    p    = exp(-beta*(k^(2/D)));

elseif STAT == 'X'

    beta = (gamma(D/2 + 1)/Ek)^(2/D);
    p    = exp(-beta*(k^(2/D)));

elseif STAT == 'F'

    beta = (gamma(D/2 + 1)/Ek)^(2/D);
    p    = exp(-beta*(k^(2/D)));

end


%-Poisson clumping heuristic {for multiple clusters}
%==========================================================================
P        = 1 - spm_Pcdf(c - 1,(Ec + eps)*p);
