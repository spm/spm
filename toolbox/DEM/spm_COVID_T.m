function [T,R] = spm_COVID_T(x,P)
% state dependent probability transition matrices
% FORMAT [T,R] = spm_COVID_T(x,P)
% x      - probability distributions (tensor)
% P      - model parameters
% 
% T      - probability transition matrix
% R      - time varying parameters
%
% This subroutine creates a transition probability tensors as a function of
% model parameters and the joint density over four factors, each with
% several levels. Crucially the transition probabilities of any one factor
% depend only upon another factor. for example, in the factor modelling
% clinical status, the transition from acute respiratory distress (ARDS) to
% death depends upon infection status (infected or not infected) and
% location (in a critical care unit or not). This version has no absorbing
% states. States such as contributing to daily deaths or tests are modelled
% by remaining in that state for one day and then returning to another
% state.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_COVID_T.m 8038 2021-01-01 16:33:49Z karl $

% setup
%==========================================================================

% identity matrices
%--------------------------------------------------------------------------
dim   = size(x);
dim   = dim(1:4);
I     = cell(numel(dim),1);
for i = 1:numel(dim)
    I{i} = speye(dim(i));
end

% exponentiate parameters
%--------------------------------------------------------------------------
P    = spm_vecfun(P,@exp);

% daily rate
%--------------------------------------------------------------------------
Kday = exp(-1);

% divergence from endemic equilibrium (1,0)
%--------------------------------------------------------------------------
Q    = (1 + cos(2*pi*(P.t - log(P.inn)*8)/365))/2;
T    = exp(-P.t/128);

% probabilistic transitions: location
%==========================================================================

% social distancing, based on prevalence of infection
%--------------------------------------------------------------------------
b    = cell(1,dim(3));
Pout = 0;                            % fluctuating mobility
if isfield(P,'mob')
    for i = 1:numel(P.mob)
        Pout = Pout + log(P.mob(i)) * cos(i*pi*P.t/365)/8;
    end
end
Pout = P.out*exp(Pout);

% lockdown (threshold) strategy
%--------------------------------------------------------------------------
q    = spm_sum(x,[1 3 4 5]);
Prev = q(3) + P.qua*q(4);            % prevalence of infection
Psde = spm_sigma(Prev,P.sde,P.s);    % 1 - P(lockdown)
Pout = Psde*Pout;                    % P(work | asymptomatic)

% viral spread
%--------------------------------------------------------------------------
Pexp = Psde*P.exp;                   % P(exposed | asymptomatic)
Pdis = P.m*Pexp;                     % P(leaving effective population)

% bed availability
%--------------------------------------------------------------------------
Phos = P.hos;                        % P(hospital | ARDS)
Pnhs = 0.03;                         % P(hospital | NHS/care worker)
Pcap = P.ccu;                        % P(transfer to CCU | ARDS)
Piso = exp(-1/10);                   % period of self-isolation

Ph2h = (1 - Pout)*(1 - Pdis)*(1 - Pnhs);
Ph2w = Pout*(1 - Pdis)*(1 - Pnhs);
Ph2r = Pdis*(1 - Pnhs);

% marginal: location {1} | asymptomatic {3}(1)
%--------------------------------------------------------------------------
%      home       work       ccu      removed     isolated   hospital
%--------------------------------------------------------------------------
b{1} = [Ph2h       1          0          Pexp      (1 - Piso) 1;
        Ph2w       0          0          0          0         0;
        0          0          0          0          0         0;
        Ph2r       0          0          (1 - Pexp) 0         0;
        0          0          1          0          Piso      0;
        Pnhs       0          0          0          0         0];

% marginal: location {1}  | symptoms {3}(2)
%--------------------------------------------------------------------------
b{2} = [0          0          0          0          0         0;
        0          0          0          0          0         0;
        0          0          0          0          0         0;
        0          0          0          0          0         0;
        1          1          1          1          1         1;
        0          0          0          0          0         0];
    
% marginal: location {1}  | ARDS {3}(3)
%--------------------------------------------------------------------------
b{3} = [0          0          0          0          0         0;
        0          0          0          0          0         0;
        0          0          1          0          0         Pcap;
        0          0          0          0          0         0;
        (1 - Phos) (1 - Phos) 0         (1 - Phos) (1 - Phos) 0;
        Phos       Phos       0          Phos       Phos      (1 - Pcap)];


% marginal: location {1}  | deceased {3}(4)
%--------------------------------------------------------------------------
b{4} = eye(6,6);

% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{2},I{4}});
B{1} = spm_permute_kron(b,dim([1,3,2,4]),[1,3,2,4]);

% equivalent to...
%--------------------------------------------------------------------------
% for i = 1:numel(b)
%     b{i} = kron(I{2},b{i});
% end
% b    = spm_cat(spm_diag(b));
% B{1} = kron(I{4},b);

% stop isolating if asymptomatic and PCR- : third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({5,1:5,1,4},{1,1:5,1,4},dim);  B{1}(ij) = 1;
ij   = Bij({5,1:5,1,4},{5,1:5,1,4},dim);  B{1}(ij) = 0;

% isolate if positive : third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({1,1:5,1,3},{5,1:5,1,3},dim);  B{1}(ij) = 1;
ij   = Bij({1,1:5,1,3},{1,1:5,1,3},dim);  B{1}(ij) = 0;
ij   = Bij({1,1:5,1,3},{2,1:5,1,3},dim);  B{1}(ij) = 0;
ij   = Bij({1,1:5,1,3},{4,1:5,1,3},dim);  B{1}(ij) = 0;
ij   = Bij({1,1:5,1,3},{6,1:5,1,3},dim);  B{1}(ij) = 0;

% isolate if infected/ious : third order dependencies : efficacy of FTTI
%--------------------------------------------------------------------------
ij   = Bij({1,2:3,1,1},{5,2:3,1,1},dim);  B{1}(ij) = P.ttt;
ij   = Bij({1,2:3,1,1},{1,2:3,1,1},dim);  B{1}(ij) = (1 - Pout)*(1 - P.ttt);
ij   = Bij({1,2:3,1,1},{2,2:3,1,1},dim);  B{1}(ij) = Pout*(1 - P.ttt);
ij   = Bij({1,2:3,1,1},{4,2:3,1,1},dim);  B{1}(ij) = 0;
ij   = Bij({1,2:3,1,1},{6,2:3,1,1},dim);  B{1}(ij) = 0;

% probabilistic transitions: infection
%==========================================================================

% transmission probabilities
%--------------------------------------------------------------------------
b    = cell(1,dim(2));
q    = spm_sum(x,[3 4 5]);
pin  = q(1,:)/sum(q(1,:) + eps);     % infection probability at home
pou  = q(2,:)/sum(q(2,:) + eps);     % infection probability at work
phs  = q(6,:)/sum(q(6,:) + eps);     % infection probability at hospital

Pmut = 0;                            % fluctuating transmission strength
if isfield(P,'vir')
    for i = 1:numel(P.vir)
        Pmut = Pmut + log(P.vir(i)) * cos(i*pi*P.t/365)/8;
    end
end
Pmut = exp(Pmut);

Ptrn = P.trn*Q + P.trm*(1 - Q);      % transmission strength
Ptrn = erf(Ptrn*Pmut);               % modulated strength
Ptin = (1 - Ptrn*pin(3))^P.Nin;      % P(no transmission) | home
Ptou = (1 - Ptrn*pou(3))^P.Nou;      % P(no transmission) | work
Pths = (1 - Ptrn*phs(3))^P.Nou;      % P(no transmission) | hospital
Kimm = exp(-1/P.Tim);                % loss of Ab+ immunity (per day)
Kinn = exp(-1/2048);                 % loss of Ab- immunity (per day)
Kinf = exp(-1/P.Tin);                % infection rate
Kcon = exp(-1/P.Tcn);                % infectious rate
Pres = P.res;                        % resistant proportion
Pvac = P.vac;                        % vaccination rate
    
% marginal: infection {2} | home {1}(1)
%--------------------------------------------------------------------------
%    susceptible  infected           infectious     Ab+        Ab-
%--------------------------------------------------------------------------
b{1} = [Ptin       0                     0          0         (1 - Kinn);
        (1 - Ptin) Kinf                  0          0          0;
        0          (1 - Pres)*(1 - Kinf) Kcon       0          0;
        0          0                     (1 - Kcon) Kimm       0;
        0          Pres*(1 - Kinf)       0          (1 - Kimm) Kinn];
    
% marginal: infection {2} | work {1}(2)
%--------------------------------------------------------------------------
b{2} = [Ptou       0                     0          0          (1 - Kinn);
        (1 - Ptou) Kinf                  0          0          0;
        0          (1 - Pres)*(1 - Kinf) Kcon       0          0;
        0          0                     (1 - Kcon) Kimm       0;
        0          Pres*(1 - Kinf)       0          (1 - Kimm) Kinn];


% marginal: infection {2} | ccu {1}(3)
%--------------------------------------------------------------------------
b{3} = [1          0                     0          0         (1 - Kinn);
        0          Kinf                  0          0          0;
        0          (1 - Pres)*(1 - Kinf) Kcon       0          0;
        0          0                     (1 - Kcon) Kimm       0;
        0          Pres*(1 - Kinf)       0          (1 - Kimm) Kinn];

% marginal: infection {2} | removed {1}(4)
%--------------------------------------------------------------------------
b{4} = b{3};

% marginal: infection {2} | isolated {1}(5)
%--------------------------------------------------------------------------
b{5} = b{3};

% marginal: infection {2} | hospital {1}(6)
%--------------------------------------------------------------------------
b{6} = [Pths*(1 - Pvac)       0                     0          0         (1 - Kinn);
        (1 - Pths)*(1 - Pvac) Kinf                  0          0          0;
        0                     (1 - Pres)*(1 - Kinf) Kcon       0          0;
        0                     0                     (1 - Kcon) Kimm       0;
        Pvac                  Pres*(1 - Kinf)       0          (1 - Kimm) Kinn];

% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{3},I{4}});
B{2} = spm_permute_kron(b,dim([2,1,3,4]),[2,1,3,4]);

% probabilistic transitions: clinical
%==========================================================================

% probabilities of developing symptoms
%--------------------------------------------------------------------------
b    = cell(1,dim(2));

Psev = erf(P.sev*Q + P.lat*(1 - Q));     % P(ARDS | infected)
Pfat = erf(P.fat*T + P.sur*(1 - T));     % P(fatality | ARDS, CCU)
Ksym = exp(-1/P.Tsy);                    % acute symptomatic rate
Ktrd = exp(-1/P.Trd);                    % acute RDS rate
Ktic = exp(-1/P.Tic);                    % symptomatic rate


% marginal: clinical {3} | susceptible {2}(1)
%--------------------------------------------------------------------------
%  asymptomatic    symptomatic           ARDS        deceased
%--------------------------------------------------------------------------
b{1} = [1          (1 - Ksym)*(1 - Psev) 0         (1 - Kday);
        0          Ksym                  0                  0;
        0          (1 - Ksym)*Psev       Ktrd               0;
        0          0                     (1 - Ktrd)      Kday];
    
% marginal: clinical {3} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [Ktic       (1 - Ksym)*(1 - Psev) 0         (1 - Kday);
        (1 - Ktic) Ksym                  0                  0;
        0          (1 - Ksym)*Psev       Ktrd               0;
        0          0                     (1 - Ktrd)      Kday];
    
% marginal: clinical {3} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = b{2};
    
% marginal: clinical {3} | Ab+ {2}(4)
%--------------------------------------------------------------------------
b{4} = b{1};
    
% marginal: clinical {3} | Ab- {2}(5)
%--------------------------------------------------------------------------
b{5} = b{1};

% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{1}});
b    = spm_permute_kron(b,dim([3,2,1]),[3,2,1]);

% kroneckor form
%--------------------------------------------------------------------------
B{3} = spm_kron({b,I{4}});

% location dependent fatalities (in hospital): third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({3,1:5,3,1:4},{3,1:5,4,1:4},dim); B{3}(ij) = (1 - Ktrd)*Pfat;
ij   = Bij({3,1:5,3,1:4},{3,1:5,1,1:4},dim); B{3}(ij) = (1 - Ktrd)*(1 - Pfat);
ij   = Bij({6,1:5,3,1:4},{6,1:5,4,1:4},dim); B{3}(ij) = (1 - Ktrd)*Pfat;
ij   = Bij({6,1:5,3,1:4},{6,1:5,1,1:4},dim); B{3}(ij) = (1 - Ktrd)*(1 - Pfat);


% probabilistic transitions: testing
%==========================================================================

% test probabilities
%--------------------------------------------------------------------------
b    = cell(1,dim(2));

% fluctuations in testing rate (discrete cosine basis functions)
%--------------------------------------------------------------------------
Ppcr = 0;                             
if isfield(P,'pcr')
    for i = 1:numel(P.pcr)
        Ppcr = Ppcr + log(P.pcr(i)) * cos(i*pi*P.t/365)/8;
    end
end
Ppcr = exp(Ppcr);

% (pillar one, two and LFT) phases of testing
%--------------------------------------------------------------------------
Pill  = zeros(1,numel(P.lim));
for i = 1:numel(P.lim)
    p(i)    = spm_phi((P.t - P.ons(i))/P.rat(i));
    Pill(i) = P.lim(i)*p(i);
end

Psen = sum(Pill);                              % testing rate
Ptes = P.tes*(1 - p(2)) + P.tts*p(2);          % testing bias
Psen = erf(Psen*Ppcr);                         % probability of testing 
Ptes = erf(Psen*Ptes);                         % probability if infected

Pfnr = P.fnr*(1 - p(3)) + 0.25*p(3);           % false negative rate
Pfpr = P.fpr*(1 - p(3)) + 0.0032*p(3);         % false positive rate
Sens = 1 - Pfnr;                               % sensitivity
Spec = 1 - Pfpr;                               % specificity                    
Kdel = exp(-1/P.del);                          % exp(-1/waiting period)

% marginal: testing {4} | susceptible {2}(1)
%--------------------------------------------------------------------------
%    not tested  waiting                 PCR+       PCR-
%--------------------------------------------------------------------------
b{1} = [(1 - Psen) 0                    (1 - Kday) (1 - Kday);
        Psen       Kdel                  0          0;
        0          (1 - Spec)*(1 - Kdel) Kday       0;
        0          Spec*(1 - Kdel)       0          Kday];

% marginal: testing {4} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Ptes) 0                     (1 - Kday) (1 - Kday);
        Ptes       Kdel                   0          0;
        0          Sens*(1 - Kdel)        Kday       0;
        0          (1 - Sens)*(1 - Kdel)  0          Kday];
    
% marginal: testing {4} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = b{2};
    
% marginal: testing {4} | Ab+ {2}(4)
%--------------------------------------------------------------------------
b{4} = b{1};
    
% marginal: testing {4} | Ab- {2}(5)
%--------------------------------------------------------------------------
b{5} = b{1};

% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{1},I{3}});
B{4} = spm_permute_kron(b,dim([4,2,1,3]),[3,2,4,1]);

% location dependent testing (none when removed): third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({4,1:5,1:4,1},{4,1:5,1:4,1},dim); B{4}(ij) = 1;
ij   = Bij({4,1:5,1:4,1},{4,1:5,1:4,2},dim); B{4}(ij) = 0;


% probability transition matrix
%==========================================================================
T     = 1;
for i = 1:4
    T =  T*B{i};
end

% time-dependent parameters
%--------------------------------------------------------------------------
R.Pout = Pout;
R.Ptrn = Ptrn;
R.Psev = Psev;
R.Pfat = Pfat;
R.Psen = Psen;
R.Ptes = Ptes;

return


% probabilistic transitions: state of the nation
%==========================================================================
q    = spm_sum(x,[1 3 4 5]);
Prev = q(2);                           % prevalence of infection

% test probabilities
%--------------------------------------------------------------------------
Pt00 = 1 - Prev;                       % P(tier 1 | tier 0)
Pt11 = exp(-1/P.mem);                  % P(tier 0 | tier 1)

% marginal: tiers
%--------------------------------------------------------------------------
%  tier: 0         1
%--------------------------------------------------------------------------
B{5} = [Pt00  (1 - Pt11);
        (1 - Pt00) Pt11];
T    = spm_kron({T,B{5}});


% Auxiliary functions
%__________________________________________________________________________

function ij = Bij(j,i,dim)
% returns linear indices of transition matrix
% j - indices of source states
% i - indices of target states
%--------------------------------------------------------------------------
z  = zeros(dim); z(j{1},j{2},j{3},j{4}) = 1; j = find(z);
z  = zeros(dim); z(i{1},i{2},i{3},i{4}) = 1; i = find(z);
ij = i + (j-1)*prod(dim);


function p = spm_sigma(x,u,s)
% reverse sigmoid function
% FORMAT p = spm_sigma(p,u,s)
% x    - probability
% u    - threshold
% s    - sensitivity (default four)
%
% p    - probability (0 < p < 1)
%
% This function is reverse sigmoid function that scales the input argument
% by the bias and flips the (biased) input. This provides a monotonically
% decreasing sigmoid function of the input that hits 50% at the threshold
% (u). The scaling ensures the probability at x = 0 is about one, for a
% suitably large sensitivity parameter s.
%--------------------------------------------------------------------------

% default sensitivity
%--------------------------------------------------------------------------
if nargin < 3, s = 4; end

% sigmoid function
%--------------------------------------------------------------------------
p = spm_phi(s*(u - x)/u);

return

