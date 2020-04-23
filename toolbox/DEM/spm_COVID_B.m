function T = spm_COVID_B(x,P,r)
% state dependent probability transition matrices
% FORMAT T = spm_COVID_B(p,P,r)
% x      - probability distributions (tensor)
% P      - model parameters
% R      - marginals over regions
% 
% B      - probability transition matrices
%
% this subroutine creates a transition probability tensors as a function of
% model parameters and the joint density over four factors, each with four
% levels. With one exception, the transition probabilities of anyone factor
% depend only upon another factor. The exception is the factor modelling
% clinical status, where the transition from acute respiratory distress
% (ARDS) to death depends upon infection status (infected or not
% infected) and location (in a critical careunit will not).
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_COVID_B.m 7838 2020-04-23 17:40:45Z karl $

% marginal probabilities
%==========================================================================
p     = spm_marginal(x);

% identity matrices
%--------------------------------------------------------------------------
dim   = size(x);
nx    = numel(x);
I     = cell(ndims(x),1);
for i = 1:ndims(x)
    I{i} = speye(dim(i));
end

% exponentiate parameters
%--------------------------------------------------------------------------
P     = spm_vecfun(P,@exp);

% probabilistic transitions: location
%==========================================================================
% P.out                             % P(work | home)
% P.sde                             % social distancing threshold
% P.cap                             % bed availability threshold (per capita)
% social distancing, based on prevalence of infection
%--------------------------------------------------------------------------
Ppn = p{2}(2) + p{2}(3);            % prevalence of infection
Pco = p{1}(3);                      % CCU occupancy

% hard (threshold) strategy
%--------------------------------------------------------------------------
Psd = spm_sigma(Ppn,P.sde);
Psd = spm_sigma(Pco,P.sde*P.cap*8)*Psd;

% multiregional model
%--------------------------------------------------------------------------
if nargin > 2
    
    Rpn = r{2}(2) + r{2}(3);        % marginal over regions
    Rco = r{1}(3);                  % marginal over regions
    
    % hard (threshold) strategy
    %----------------------------------------------------------------------
    Rsd = spm_sigma(Rpn,P.sde);
    Rsd = spm_sigma(Rco,P.sde*P.cap*8)*Rsd;
    
    % mixture of strategies
    %----------------------------------------------------------------------
    Psd = Psd*(1 - P.fed) + Rsd*P.fed;
    
end
Pout = 1/128 + Psd*P.out;           % P(work | home)

% bed availability
%--------------------------------------------------------------------------
Pcap = P.cap*(1 + p{1}(2));         % bed capacity threshold
Pcca = spm_sigma(p{1}(3),Pcap);     % P(CCU  | home, work, ARDS)
b    = cell(1,dim(3));


% marginal: location {1} | asymptomatic {3}(1)
%--------------------------------------------------------------------------
%      home       work       CCU       morgue
%--------------------------------------------------------------------------
b{1} = [(1 - Pout) 1          1          0;
        Pout       0          0          0;
        0          0          0          0;
        0          0          0          1];

% marginal: location {1}  | symptoms {3}(2)
%--------------------------------------------------------------------------
b{2} = [1          1          1          0;
        0          0          0          0;
        0          0          0          0;
        0          0          0          1];

% marginal: location {1}  | ARDS {3}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Pcca) (1 - Pcca) 0          0;
        0          0          0          0;
        Pcca       Pcca       1          0;
        0          0          0          1];

% marginal: location {1}  | deceased {3}(4)
%--------------------------------------------------------------------------
b{4} = [0          0          0          0;
        0          0          0          0;
        0          0          0          0;
        1          1          1          1];

% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{2},I{4}});
B{1} = spm_permute_kron(b,dim([1,3,2,4]),[1,3,2,4]);

% immunity dependent return to work: third order dependencies
%--------------------------------------------------------------------------
Pout = P.out*(Psd*(1 - P.btw) + P.btw);
z    = zeros(dim); z(1,4,1,:) = 1; j = find(z);
z    = zeros(dim); z(2,4,1,:) = 1; i = find(z);
ij   = find(sparse(i,j,1,nx,nx));
B{1}(ij) = Pout;

z    = zeros(dim); z(1,4,1,:) = 1; j = find(z);
z    = zeros(dim); z(1,4,1,:) = 1; i = find(z);
ij   = find(sparse(i,j,1,nx,nx));
B{1}(ij) = 1 - Pout;


        
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
R    = (1 - P.trn*p{2}(3));         % P(no transmission per contact)
Pinh = R^P.Rin;                     % P(no transmission) | home
Pinw = R^P.Rou;                     % P(no transmission) | work
Kimm = exp(-1/P.Tim/32);            % loss of immunity (per 32 days)
Kinf = exp(-1/P.Tin);
Kcon = exp(-1/P.Tcn);
    
% marginal: infection {2} | home {1}(1)
%--------------------------------------------------------------------------
%    susceptible  infected  infectious  immune  resistant
%--------------------------------------------------------------------------
b{1} = [Pinh       0          0          (1 - Kimm) 0;
        (1 - Pinh) Kinf       0          0          0;
        0          (1 - Kinf) Kcon       0          0;
        0          0          (1 - Kcon) Kimm       0;
        0          0          0          0          1];
    
% marginal: infection {2} | work {1}(2)
%--------------------------------------------------------------------------
b{2} = [Pinw       0          0          (1 - Kimm) 0;
        (1 - Pinw) Kinf       0          0          0;
        0          (1 - Kinf) Kcon       0          0;
        0          0          (1 - Kcon) Kimm       0;
        0          0          0          0          1];


% marginal: infection {2} | CCU {1}(3)
%--------------------------------------------------------------------------
b{3} = [1          0          0          (1 - Kimm) 0;
        0          Kinf       0          0          0;
        0          (1 - Kinf) Kcon       0          0;
        0          0          (1 - Kcon) Kimm       0;
        0          0          0          0          1];


% marginal: infection {2} | morgue {1}(4)
%--------------------------------------------------------------------------
b{4} = [0          0          0          0          0;
        0          0          0          0          0;
        0          0          0          0          0;
        1          1          1          1          0;
        0          0          0          0          1];


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
Kdev = exp(-1/P.Tic);               % acute RDS rate
Pfat = 1 - P.sur;                   % baseline fatality rate

% marginal: clinical {3} | susceptible {2}(1)
%--------------------------------------------------------------------------
%  asymptomatic symptomatic acute RDS deceased
%--------------------------------------------------------------------------
b{1} = [1          1          1         0;
        0          0          0         0;
        0          0          0         0;
        0          0          0         1];

% marginal: clinical {3} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [Kdev       (1 - Ksym)*(1 - Psev) (1 - Ksev)*(1 - Pfat) 0;
        (1 - Kdev) Ksym                   0                    0;
        0          (1 - Ksym)*Psev        Ksev                 0;
        0          0                     (1 - Ksev)*Pfat       1];
    
% marginal: clinical {3} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [Kdev       (1 - Ksym)*(1 - Psev) (1 - Ksev)*(1 - Pfat) 0;
        (1 - Kdev)  Ksym                   0                    0;
        0          (1 - Ksym)*Psev        Ksev                 0;
        0          0                     (1 - Ksev)*Pfat       1];
    
% marginal: clinical {3} | immune {2}(4)
%--------------------------------------------------------------------------
b{4} = [1          1          1         0;
        0          0          0         0;
        0          0          0         0;
        0          0          0         1];
    
% marginal: clinical {3} | resistant {2}(5)
%--------------------------------------------------------------------------
b{5} = [1          1          1         0;
        0          0          0         0;
        0          0          0         0;
        0          0          0         1];


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
z    = zeros(dim); z(3,[2,3],3,:) = 1; j = find(z);
z    = zeros(dim); z(3,[2,3],4,:) = 1; i = find(z);
ij   = find(sparse(i,j,1,nx,nx));
B{3}(ij) = (1 - Ksev)*P.fat;

z    = zeros(dim); z(3,[2,3],3,:) = 1; j = find(z);
z    = zeros(dim); z(3,[2,3],1,:) = 1; i = find(z);
ij   = find(sparse(i,j,1,nx,nx));
B{3}(ij) = (1 - Ksev)*(1 - P.fat);

% probabilistic transitions: testing
%==========================================================================
% P.tft                       % threshold:   testing capacity
% P.sen;                      % sensitivity: testing capacity
% P.del                       % delay:       testing capacity
% P.tes                       % relative probability of test if uninfected
% test availability and prevalence of symptoms
%--------------------------------------------------------------------------
b    = cell(1,dim(2));
Ptes = P.sen*spm_sigma(p{4}(2),P.tft);
Pdia = P.tes*Ptes;
Kdel = exp(-1/P.del);         % exp(-1/waiting period)

% marginal: testing {4} | susceptible {2}(1)
%--------------------------------------------------------------------------
%    not tested  waiting       +ve -ve
%--------------------------------------------------------------------------
b{1} = [(1 - Pdia) 0            0   0;
        Pdia       Kdel         0   0;
        0          0            1   0;
        0          (1 - Kdel)   0   1];

% marginal: testing {4} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Ptes) 0            0   0;
        Ptes       Kdel         0   0;
        0          (1 - Kdel)   1   0;
        0          0            0   1];
    
% marginal: testing {4} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Ptes) 0            0   0;
        Ptes       Kdel         0   0;
        0          (1 - Kdel)   1   0;
        0          0            0   1];
    
% marginal: testing {4} | immune {2}(4)
%--------------------------------------------------------------------------
b{4} = [(1 - Pdia) 0            0   0;
        Pdia       Kdel         0   0;
        0          0            1   0;
        0          (1 - Kdel)   0   1];
    
% marginal: testing {4} | resistant {2}(5)
%--------------------------------------------------------------------------
b{5} = [(1 - Pdia) 0            0   0;
        Pdia       Kdel         0   0;
        0          0            1   0;
        0          (1 - Kdel)   0   1];


% kroneckor form
%--------------------------------------------------------------------------
b     = spm_cat(spm_diag(b));
b     = spm_kron({b,I{1},I{3}});
B{4}  = spm_permute_kron(b,dim([4,2,1,3]),[3,2,4,1]);
    
% probability transition matrix
%==========================================================================
T     = 1;
for i = 1:numel(B)
    T =  T*B{i};
end


return

% Auxiliary functions
%__________________________________________________________________________

function p = spm_sigma(x,u,s)
% reverse sigmoid function
% FORMAT p = spm_sigma(p,u)
% x    - probability
% u    - threshold
% u    - sensitivity (default four)
%
% p    - probability (0 < p < 1)
%
% this function is standard sigmoid function but scales the input argument
% by the bias and flips the (biased) input. This provides a monotonically
% decreasing sigmoid function of the input that hits 50 at the threshold
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