function T = spm_COVID_SB(P,I,dim,Prev,Pcco,Pinh,Pinw)
% state dependent probability transition matrices
% FORMAT T = spm_COVID_SB(P,I,dim,Prev,Pcco,Pinh,Pinw)
% P      - model parameters
% 
% T      - probability transition matrix
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
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging




% exponentiate parameters
%--------------------------------------------------------------------------
P     = spm_vecfun(P,@exp);

% upper bound probabilities
%--------------------------------------------------------------------------
P.out = min(P.out,1);
P.trn = min(P.trn,1);
P.sev = min(P.sev,1);
P.fat = min(P.fat,1);
P.sur = min(P.sur,1);

Kday  = exp(-1);


% probabilistic transitions: location
%==========================================================================
% P.out                              % P(work | home)
% P.sde                              % social distancing threshold
% P.cap                              % bed threshold (per capita)
% social distancing, based on prevalence of infection
%--------------------------------------------------------------------------

% hard (threshold) strategy
%--------------------------------------------------------------------------
Psde = spm_sigma(Prev,P.sde);
Pout = Psde*P.out;                   % P(work | home)

% bed availability
%--------------------------------------------------------------------------
Pcca = spm_sigma(Pcco,P.cap);        % P(CCU  | home, work, ARDS)
Piso = exp(-1/7);                    % period of self-isolation
b    = cell(1,dim(3));

% viral spread
%--------------------------------------------------------------------------
Kspr = exp(-Psde*Prev/P.Tex);        % period of exemption

% marginal: location {1} | asymptomatic {3}(1)
%--------------------------------------------------------------------------
%      home       work       CCU       exempt     isolation
%--------------------------------------------------------------------------
b{1} = [(1 - Pout) 1          1          (1 - Kspr) (1 - Piso);
        Pout       0          0          0           0;
        0          0          0          0           0;
        0          0          0          Kspr        0;
        0          0          0          0           Piso];

% marginal: location {1}  | symptoms {3}(2)
%--------------------------------------------------------------------------
b{2} = [0          0          0          (1 - Kspr)  0;
        0          0          0          0           0;
        0          0          0          0           0;
        0          0          0          Kspr        0;
        1          1          1          0           1];
    
% marginal: location {1}  | ARDS {3}(3)
%--------------------------------------------------------------------------
b{3} = [0          0          0          (1 - Kspr)  0;
        0          0          0          0           0;
        Pcca       Pcca       1          0           Pcca;
        0          0          0          Kspr        0;
        (1 - Pcca) (1 - Pcca) 0          0           (1 - Pcca)];

% marginal: location {1}  | deceased {3}(4)
%--------------------------------------------------------------------------
b{4} = [0          0          0          (1 - Kspr)  0;
        0          0          0          0           0;
        0          0          0          0           0;
        1          1          1          Kspr        1;
        0          0          0          0           0];

% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{2},I{4}});
B{1} = spm_permute_kron(b,dim([1,3,2,4]),[1,3,2,4]);

% stop isolating if asymptomatic and negative : third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({5,1:5,1,4},{1,1:5,1,4},dim);  B{1}(ij) = 1;
ij   = Bij({5,1:5,1,4},{5,1:5,1,4},dim);  B{1}(ij) = 0;

% isolate if positive : third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({1,1:5,1,3},{5,1:5,1,3},dim);  B{1}(ij) = 1;
ij   = Bij({1,1:5,1,3},{1,1:5,1,3},dim);  B{1}(ij) = 0;
ij   = Bij({1,1:5,1,3},{2,1:5,1,3},dim);  B{1}(ij) = 0;

% isolate if infected : third order dependencies : efficacy of FTTI
%--------------------------------------------------------------------------
ij   = Bij({1,2,1,1},{5,2,1,1},dim);  B{1}(ij) = P.ttt;
ij   = Bij({1,2,1,1},{1,2,1,1},dim);  B{1}(ij) = (1 - Pout)*(1 - P.ttt);
ij   = Bij({1,2,1,1},{2,2,1,1},dim);  B{1}(ij) = Pout*(1 - P.ttt);


% probabilistic transitions: infection
%==========================================================================
% P.Rin                             % effective number of contacts: home
% P.Rou                             % effective number of contacts: work
% P.trn                             % P(transmission | contact)
% P.Tin                             % infected (pre-contagious) period
% P.Tcn                             % infectious (contagious) period
% transmission probabilities
%--------------------------------------------------------------------------
b    = cell(1,dim(2));
Kimm = exp(-1/P.Tim/32);            % loss of immunity (per 32 days)
Pres = P.res;                       % proportion of innate immunity
Kinf = exp(-1/P.Tin);
Kcon = exp(-1/P.Tcn);
    
% marginal: infection {2} | home {1}(1)
%--------------------------------------------------------------------------
%    susceptible  infected           infectious     immune    resistant
%--------------------------------------------------------------------------
b{1} = [Pinh       0                     0          (1 - Kimm) 0;
        (1 - Pinh) Kinf                  0          0          0;
        0          (1 - Pres)*(1 - Kinf) Kcon       0          0;
        0          0                     (1 - Kcon) Kimm       0;
        0          Pres*(1 - Kinf)       0          0          1];
    
% marginal: infection {2} | work {1}(2)
%--------------------------------------------------------------------------
b{2} = [Pinw       0                     0          (1 - Kimm) 0;
        (1 - Pinw) Kinf                  0          0          0;
        0          (1 - Pres)*(1 - Kinf) Kcon       0          0;
        0          0                     (1 - Kcon) Kimm       0;
        0          Pres*(1 - Kinf)       0          0          1];


% marginal: infection {2} | CCU {1}(3)
%--------------------------------------------------------------------------
b{3} = [1          0                     0          (1 - Kimm) 0;
        0          Kinf                  0          0          0;
        0          (1 - Pres)*(1 - Kinf) Kcon       0          0;
        0          0                     (1 - Kcon) Kimm       0;
        0          Pres*(1 - Kinf)       0          0          1];

% marginal: infection {2} | morgue {1}(4)
%--------------------------------------------------------------------------
b{4} = [0          0          0          0          0;
        0          0          0          0          0;
        0          0          0          0          0;
        0          0          0          0          0;
        1          1          1          1          1];

% marginal: infection {2} | isolation {1}(5)
%--------------------------------------------------------------------------
b{5} = b{3};


% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{3},I{4}});
B{2} = spm_permute_kron(b,dim([2,1,3,4]),[2,1,3,4]);


% probabilistic transitions: clinical
%==========================================================================
% https://en.wikipedia.org/wiki/List_of_countries_by_life_expectancy
%--------------------------------------------------------------------------
% P.dev                             % P(developing symptoms | infected)
% P.sev                             % P(severe symptoms | symptomatic)
% P.Tsy                             % symptomatic period
% P.Trd                             % acute RDS   period
% P.fat                             % P(fatality | CCU)
% P.sur                             % P(survival | home)
% probabilities of developing symptoms
%--------------------------------------------------------------------------
b    = cell(1,dim(2));
Psev = P.sev;                       % P(developing symptoms | infected)
Ksym = exp(-1/P.Tsy);               % acute symptomatic rate
Ksev = exp(-1/P.Trd);               % acute RDS rate
Kdev = exp(-1/P.Tic);               % symptomatic rate
Pfat = 1 - P.sur;                   % baseline fatality rate

% marginal: clinical {3} | susceptible {2}(1)
%--------------------------------------------------------------------------
%  asymptomatic   symptomatic            acute RDS         deceased
%--------------------------------------------------------------------------
b{1} = [1         (1 - Ksym)             (1 - Ksev)*(1 - Pfat) (1 - Kday);
        0          Ksym                   0                    0;
        0          0                      Ksev                 0;
        0          0                     (1 - Ksev)*Pfat       Kday];
    
% marginal: clinical {3} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [Kdev       (1 - Ksym)*(1 - Psev) (1 - Ksev)*(1 - Pfat) (1 - Kday);
        (1 - Kdev) Ksym                   0                    0;
        0          (1 - Ksym)*Psev        Ksev                 0;
        0          0                     (1 - Ksev)*Pfat       Kday];
    
% marginal: clinical {3} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = b{2};
    
% marginal: clinical {3} | immune {2}(4)
%--------------------------------------------------------------------------
b{4} = b{1};
    
% marginal: clinical {3} | resistant {2}(5)
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

% location dependent fatalities (P.fat in CCU): third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({3,1:5,3,1:4},{3,1:5,4,1:4},dim); B{3}(ij) = (1 - Ksev)*P.fat;
ij   = Bij({3,1:5,3,1:4},{3,1:5,1,1:4},dim); B{3}(ij) = (1 - Ksev)*(1 - P.fat);


% probabilistic transitions: testing
%==========================================================================
% P.bas                       % probability of being tested
% P.del                       % delay: testing capacity
% P.tes                       % relative probability if infected
% test probabilities
%--------------------------------------------------------------------------
b    = cell(1,dim(2));
Seni = .9;                            % PCR sensitivity (infected)
Senc = .95;                           % PCR sensitivity (contagious)
Psen = P.bas/(1 - Prev + P.tes*Prev); % probability of being tested
Ptes = Psen*P.tes;                    % probability if infected
Kdel = exp(-1/P.del);                 % exp(-1/waiting period)

% marginal: testing {4} | susceptible {2}(1)
%--------------------------------------------------------------------------
%    not tested  waiting       +ve -ve
%--------------------------------------------------------------------------
b{1} = [(1 - Psen) 0            (1 - Kday) (1 - Kday);
        Psen       Kdel         0          0;
        0          0            Kday       0;
        0          (1 - Kdel)   0          Kday];

% marginal: testing {4} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Ptes) 0                     (1 - Kday) (1 - Kday);
        Ptes       Kdel                   0          0;
        0          Seni*(1 - Kdel)        Kday       0;
        0          (1 - Seni)*(1 - Kdel)  0          Kday];
    
% marginal: testing {4} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Ptes) 0                     (1 - Kday) (1 - Kday);
        Ptes       Kdel                   0          0;
        0          Senc*(1 - Kdel)        Kday       0;
        0          (1 - Senc)*(1 - Kdel)  0          Kday];
    
% marginal: testing {4} | immune {2}(4)
%--------------------------------------------------------------------------
b{4} = b{1};
    
% marginal: testing {4} | resistant {2}(5)
%--------------------------------------------------------------------------
b{5} = b{1};

% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{1},I{3}});
B{4} = spm_permute_kron(b,dim([4,2,1,3]),[3,2,4,1]);

    
% probability transition matrix
%==========================================================================
T     = 1;
for i = 1:numel(B)
    T =  T*B{i};
end


return

% Auxiliary functions
%__________________________________________________________________________

function ij = Bij(j,i,dim)
% returns linear indices of transition matrix
% j - dimensions of source states
% i - dimensions on target states
%--------------------------------------------------------------------------
nx   =  prod(dim);
z    = zeros(dim); z(j{1},j{2},j{3},j{4}) = 1; j = find(z);
z    = zeros(dim); z(i{1},i{2},i{3},i{4}) = 1; i = find(z);
ij   = find(sparse(i,j,1,nx,nx));


function p = spm_sigma(x,u,s)
% reverse sigmoid function
% FORMAT p = spm_sigma(p,u)
% x    - probability
% u    - threshold
% u    - sensitivity (default four)
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