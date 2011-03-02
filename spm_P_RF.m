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
% is returned
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_P_RF.m 4225 2011-03-02 15:53:05Z guillaume $


% get expectations
%==========================================================================

% get EC densities
%--------------------------------------------------------------------------
D   = find(R,1,'last');
R   = R(1:D);
G   = sqrt(pi)./gamma(([1:D])/2);
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

% get P{n > k}
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

% Poisson clumping heuristic {for multiple clusters}
%==========================================================================
P        = 1 - spm_Pcdf(c - 1,(Ec + eps)*p);


% set P and p = [] for non-implemented cases
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if k > 0 && (STAT == 'X' || STAT == 'F')
    P    = []; p = [];
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%==========================================================================
% spm_ECdensity
%==========================================================================
function [EC] = spm_ECdensity(STAT,t,df)
% Returns the EC density
%__________________________________________________________________________
%
% Reference : Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
%
%--------------------------------------------------------------------------

% EC densities (EC}
%--------------------------------------------------------------------------
t = t(:)';
if      STAT == 'Z'

    % Gaussian Field
    %----------------------------------------------------------------------
    a       = 4*log(2);
    b       = exp(-t.^2/2);

    EC(1,:) = 1 - spm_Ncdf(t);
    EC(2,:) = a^(1/2)/(2*pi)*b;
    EC(3,:) = a/((2*pi)^(3/2))*b.*t;
    EC(4,:) = a^(3/2)/((2*pi)^2)*b.*(t.^2 - 1);

elseif  STAT == 'T'

    % T - Field
    %----------------------------------------------------------------------
    v       = df(2);
    a       = 4*log(2);
    b       = exp(gammaln((v+1)/2) - gammaln(v/2));
    c       = (1+t.^2/v).^((1-v)/2);

    EC(1,:) = 1 - spm_Tcdf(t,v);
    EC(2,:) = a^(1/2)/(2*pi)*c;
    EC(3,:) = a/((2*pi)^(3/2))*c.*t/((v/2)^(1/2))*b;
    EC(4,:) = a^(3/2)/((2*pi)^2)*c.*((v-1)*(t.^2)/v - 1);

elseif  STAT == 'X'

    % X - Field
    %----------------------------------------------------------------------
    v       = df(2);
    a       = (4*log(2))/(2*pi);
    b       = t.^(1/2*(v - 1)).*exp(-t/2-gammaln(v/2))/2^((v-2)/2);

    EC(1,:) = 1 - spm_Xcdf(t,v);
    EC(2,:) = a^(1/2)*b;
    EC(3,:) = a*b.*(t-(v-1));
    EC(4,:) = a^(3/2)*b.*(t.^2-(2*v-1)*t+(v-1)*(v-2));

elseif  STAT == 'F'

    % F Field
    %----------------------------------------------------------------------
    k       = df(1);
    v       = df(2);
    a       = (4*log(2))/(2*pi);
    b       = gammaln(v/2) + gammaln(k/2);

    EC(1,:) = 1 - spm_Fcdf(t,df);
    EC(2,:) = a^(1/2)*exp(gammaln((v+k-1)/2)-b)*2^(1/2)...
              *(k*t/v).^(1/2*(k-1)).*(1+k*t/v).^(-1/2*(v+k-2));
    EC(3,:) = a*exp(gammaln((v+k-2)/2)-b)*(k*t/v).^(1/2*(k-2))...
              .*(1+k*t/v).^(-1/2*(v+k-2)).*((v-1)*k*t/v-(k-1));
    EC(4,:) = a^(3/2)*exp(gammaln((v+k-3)/2)-b)...
              *2^(-1/2)*(k*t/v).^(1/2*(k-3)).*(1+k*t/v).^(-1/2*(v+k-2))...
              .*((v-1)*(v-2)*(k*t/v).^2-(2*v*k-v-k-1)*(k*t/v)+(k-1)*(k-2));
end
