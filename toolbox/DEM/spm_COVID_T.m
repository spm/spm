function [T,R] = spm_COVID_T(P,I)
% state dependent probability transition matrices
% FORMAT [T,R] = spm_COVID_T(P,I)
% x      - probability distributions (tensor)
% P      - model parameters
% I      - cell array of identity matrices
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
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% setup
%==========================================================================

% dimensions
%--------------------------------------------------------------------------
dim   = zeros(1,numel(I));
for i = 1:numel(I)
    dim(1,i) = length(I{i});
end

% daily rate
%--------------------------------------------------------------------------
Kday = exp(-1);

% probabilistic transitions: location
%==========================================================================

% social distancing, based on prevalence of infection
%--------------------------------------------------------------------------
b    = cell(1,dim(3));

% viral spread
%--------------------------------------------------------------------------
Pout = P.out;                        % P(leaving home | asymptomatic)
Pexp = P.out*P.exp;                  % P(exposed | asymptomatic)
Pdis = P.m*Pexp;                     % P(leaving effective population)

% bed availability
%--------------------------------------------------------------------------
Phos = erf(P.hos);                   % P(hospital | ARDS)
Pccu = erf(P.ccu);                   % P(transfer to CCU | ARDS)
Piss = erf(P.iss);                   % probability of self-isolation
Piso = exp(-1/P.iso);                % period of self-isolation

Ph2h = (1 - Pdis)*(1 - Pout);
Ph2w = (1 - Pdis)*Pout;

% marginal: location {1} | asymptomatic {3}(1)
%--------------------------------------------------------------------------
%      home       work       ccu      removed     isolated   hospital
%--------------------------------------------------------------------------
b{1} = [Ph2h       1          0          Pexp      (1 - Piso) 1;
        Ph2w       0          0          0          0         0;
        0          0          0          0          0         0;
        Pdis       0          0         (1 - Pexp)  0         0;
        0          0          1          0          Piso      0;
        0          0          0          0          0         0];

% marginal: location {1}  | symptoms {3}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Piss) (1 - Piss) (1 - Piss) 0         (1 - Piss) (1 - Piss);
        0          0          0          0          0         0;
        0          0          0          0          0         0;
        0          0          0          (1 - Piss) 0         0;
        Piss       Piss       Piss       Piss       Piss      Piss;
        0          0          0          0          0         0];
    
% marginal: location {1}  | ARDS {3}(3)
%--------------------------------------------------------------------------
b{3} = [0          0          0          0          0         0;
        0          0          0          0          0         0;
        0          0          1          0          0         Pccu;
        0          0          0          0          0         0;
       (1 - Phos) (1 - Phos)  0         (1 - Phos) (1 - Phos) 0;
        Phos       Phos       0          Phos       Phos     (1 - Pccu)];


% marginal: location {1}  | deceased {3}(4)
%--------------------------------------------------------------------------
b{4} = eye(6,6);

% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
if isa(class(Pout),'sym')
    spm_cat_T = @spm_cat_sym;
else
    spm_cat_T = @spm_cat;
end
b    = spm_cat_T(spm_diag(b));
b    = spm_kron({b,I{2},I{4}});
B{1} = spm_permute_kron(b,dim([1,3,2,4]),[1,3,2,4]);

% equivalent to...
%--------------------------------------------------------------------------
% for i = 1:numel(b)
%     b{i} = kron(I{2},b{i});
% end
% b    = spm_cat_T(spm_diag(b));
% B{1} = kron(I{4},b);

% stop isolating if asymptomatic and PCR- : third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({5,1:8,1,4},{1,1:8,1,4},dim);  B{1}(ij) = 1;
ij   = Bij({5,1:8,1,4},{5,1:8,1,4},dim);  B{1}(ij) = 0;

% isolate if asymptomatic and positive PCR/LFD+ : third order dependencies
%--------------------------------------------------------------------------
ij   = Bij({1,1:8,1,[3 5]},{5,1:8,1,[3 5]},dim);  B{1}(ij) = 1;
ij   = Bij({1,1:8,1,[3 5]},{1,1:8,1,[3 5]},dim);  B{1}(ij) = 0;
ij   = Bij({1,1:8,1,[3 5]},{2,1:8,1,[3 5]},dim);  B{1}(ij) = 0;
ij   = Bij({1,1:8,1,[3 5]},{4,1:8,1,[3 5]},dim);  B{1}(ij) = 0;
ij   = Bij({1,1:8,1,[3 5]},{6,1:8,1,[3 5]},dim);  B{1}(ij) = 0;

% isolate if infectious : third order dependencies : efficacy of FTTI
%--------------------------------------------------------------------------
ij   = Bij({1,[3 8],1,1:6},{5,[3 8],1,1:6},dim);  B{1}(ij) = P.ttt;
ij   = Bij({1,[3 8],1,1:6},{1,[3 8],1,1:6},dim);  B{1}(ij) = (1 - P.ttt)*(1 - Pout);
ij   = Bij({1,[3 8],1,1:6},{2,[3 8],1,1:6},dim);  B{1}(ij) = (1 - P.ttt)*Pout;
ij   = Bij({1,[3 8],1,1:6},{4,[3 8],1,1:6},dim);  B{1}(ij) = 0;
ij   = Bij({1,[3 8],1,1:6},{6,[3 8],1,1:6},dim);  B{1}(ij) = 0;

% vaccine-dependent transfer to CCU with ARDS: third order dependencies
%--------------------------------------------------------------------------
Pccu = Pccu*P.lnf;
ij   = Bij({6,6:8,3,1:6},{3,6:8,3,1:6},dim); B{1}(ij) = Pccu;
ij   = Bij({6,6:8,3,1:6},{6,6:8,3,1:6},dim); B{1}(ij) = (1 - Pccu);


% probabilistic transitions: infection
%==========================================================================

% transmission probabilities
%--------------------------------------------------------------------------
b    = cell(1,dim(2));

Ptin = P.tin;                              % P(no transmission) | home
Ptou = P.tou;                              % P(no transmission) | work

Kimm = exp(-1/P.Tim);                      % loss of Ab+ immunity   (per day)
Kinn = exp(-1/P.Tnn);                      % loss of Ab- immunity   (per day)
Kvac = exp(-1/P.vac);                      % Vaccine seroconversion (per day)

Kinf = exp(-1/P.Tin);                      % infection rate
Kcon = exp(-1/P.Tcn);                      % infectious rate
Pres = 1 - erf(2*P.res)/2;                 % non-resistant proportion
Pvef = erf(P.vef);                         % vaccination: Infection risk 
Prev = P.ves*Pres;                         % vaccination: Infectious risk

Pvac = P.nac;                              % 1 - vaccination rate
Rvac = 1 - Pvac;                           % vaccination rate


% marginal: infection {2} | home {1}(1)
%----------------------------------------------------------------------------------------------------------------------------------------------------------------
%    susceptible     infected            infectious      Ab+                   Ab-                      Vaccine            infected (Ab+)       infectious  (Ab+)
%----------------------------------------------------------------------------------------------------------------------------------------------------------------
Pinv = (1 - Ptin)*Pvef;
b{1} = [Ptin*Pvac        0                     0          0                    (1 - Kinn)*Ptin*Pvac        0                     0                     0; 
        (1 - Ptin)*Pvac  Kinf                  0          0                     0                          0                     0                     0;
        0               (1 - Kinf)*Pres        Kcon       0                     0                          0                     0                     0;
        0               (1 - Kinf)*(1 - Pres) (1 - Kcon) (1 - Pinv)*Kimm        Rvac                      (1 - Pinv)*(1 - Kvac)  (1 - Kinf)*(1 - Prev) (1 - Kcon) ;
        0                0                     0         (1 - Pinv)*(1 - Kimm)  Kinn*Pvac                  0                     0                     0;
        Rvac             0                     0          0                     0                         (1 - Pinv)*Kvac        0                     0;
        0                0                     0          Pinv                 (1 - Kinn)*(1 - Ptin)*Pvac  Pinv                  Kinf                  0
        0                0                     0          0                     0                          0                    (1 - Kinf)*Prev        Kcon];
    
% marginal: infection {2} | work {1}(2)
%--------------------------------------------------------------------------
Ptin = Ptou;
Pinv = (1 - Ptin)*Pvef;
b{2} = [Ptin*Pvac        0                     0          0                    (1 - Kinn)*Ptin*Pvac        0                     0                     0; 
        (1 - Ptin)*Pvac  Kinf                  0          0                     0                          0                     0                     0;
        0               (1 - Kinf)*Pres        Kcon       0                     0                          0                     0                     0;
        0               (1 - Kinf)*(1 - Pres) (1 - Kcon) (1 - Pinv)*Kimm        Rvac                      (1 - Pinv)*(1 - Kvac)  (1 - Kinf)*(1 - Prev) (1 - Kcon) ;
        0                0                     0         (1 - Pinv)*(1 - Kimm)  Kinn*Pvac                  0                     0                     0;
        Rvac             0                     0          0                     0                         (1 - Pinv)*Kvac        0                     0;
        0                0                     0          Pinv                 (1 - Kinn)*(1 - Ptin)*Pvac  Pinv                  Kinf                  0
        0                0                     0          0                     0                          0                    (1 - Kinf)*Prev        Kcon];
    
% marginal: infection {2} | ccu {1}(3)
%--------------------------------------------------------------------------
Ptin = 1;
Pinv = (1 - Ptin)*Pvef;
b{3} = [Ptin*Pvac        0                     0          0                    (1 - Kinn)*Ptin*Pvac        0                     0                     0; 
        (1 - Ptin)*Pvac  Kinf                  0          0                     0                          0                     0                     0;
        0               (1 - Kinf)*Pres        Kcon       0                     0                          0                     0                     0;
        0               (1 - Kinf)*(1 - Pres) (1 - Kcon) (1 - Pinv)*Kimm        Rvac                      (1 - Pinv)*(1 - Kvac)  (1 - Kinf)*(1 - Prev) (1 - Kcon) ;
        0                0                     0         (1 - Pinv)*(1 - Kimm)  Kinn*Pvac                  0                     0                     0;
        Rvac             0                     0          0                     0                         (1 - Pinv)*Kvac        0                     0;
        0                0                     0          Pinv                 (1 - Kinn)*(1 - Ptin)*Pvac  Pinv                  Kinf                  0
        0                0                     0          0                     0                          0                    (1 - Kinf)*Prev        Kcon];
    

% marginal: infection {2} | removed {1}(4)
%--------------------------------------------------------------------------
b{4} = b{3};

% marginal: infection {2} | isolated {1}(5)
%--------------------------------------------------------------------------
b{5} = b{3};

% marginal: infection {2} | hospital {1}(6)
%--------------------------------------------------------------------------
b{6} = b{2};
    
% kroneckor form
%--------------------------------------------------------------------------
if max(spm_vec(b)) > 1; keyboard, end

b    = spm_cat_T(spm_diag(b));
b    = spm_kron({b,I{3},I{4}});
B{2} = spm_permute_kron(b,dim([2,1,3,4]),[2,1,3,4]);

% probabilistic transitions: clinical
%==========================================================================
P.fao = 1 - (1 - P.fat)*P.oth;

% probabilities of developing symptoms
%--------------------------------------------------------------------------
b    = cell(1,dim(2));

Ktic = exp(-1/P.Tic);                        % asymptomatic rate
Ksym = exp(-1/P.Tsy);                        % symptomatic rate
Ktrd = exp(-1/P.Trd);                        % ARDS rate

Psev = P.sev;                                % P(ARDS | infected)
Pfat = P.fao;                                % P(fatality | ARDS, Other)

% marginal: clinical {3} | susceptible {2}(1)
%--------------------------------------------------------------------------
%  asymptomatic    symptomatic           ARDS                  deceased
%--------------------------------------------------------------------------
b{1} = [1          (1 - Ksym)*(1 - Psev) (1 - Ktrd)*(1 - Pfat) (1 - Kday);
        0           Ksym                  0                     0;
        0          (1 - Ksym)*Psev        Ktrd                  0;
        0           0                    (1 - Ktrd)*Pfat        Kday];
    
% marginal: clinical {3} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [Ktic       (1 - Ksym)*(1 - Psev) (1 - Ktrd)*(1 - Pfat) (1 - Kday);
       (1 - Ktic)   Ksym                  0                     0;
        0          (1 - Ksym)*Psev        Ktrd                  0;
        0           0                    (1 - Ktrd)*Pfat        Kday];
    
% marginal: clinical {3} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = b{2};

% link between pathogenicity and natural immunity
%--------------------------------------------------------------------------
Psev = P.sev*P.lnk;         % immunity: pathogenicity
Pfat = P.fao*P.lnf;         % immunity: fatality

% marginal: clinical {3} | Ab+ {2}(4)
%---------------------------------------------------------------------------
b{4} = [1          (1 - Ksym)*(1 - Psev) (1 - Ktrd)*(1 - Pfat) (1 - Kday);
        0           Ksym                  0                     0;
        0          (1 - Ksym)*Psev        Ktrd                  0;
        0           0                    (1 - Ktrd)*Pfat        Kday];
    
% marginal: clinical {3} | Ab- {2}(5)
%--------------------------------------------------------------------------
b{5} = b{4};

% marginal: clinical {3} | Vaccine {2}(6)
%--------------------------------------------------------------------------
b{6} = b{4};

% marginal: clinical {3} | infected (Ab+) {2}(7)
%--------------------------------------------------------------------------
b{7} = [Ktic       (1 - Ksym)*(1 - Psev) (1 - Ktrd)*(1 - Pfat) (1 - Kday);
       (1 - Ktic)   Ksym                  0                     0;
        0          (1 - Ksym)*Psev        Ktrd                  0;
        0           0                    (1 - Ktrd)*Pfat        Kday];
    
% marginal: clinical {3} | Infectious (Ab+) {2}(8)
%--------------------------------------------------------------------------
b{8} = b{7};


% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat_T(spm_diag(b));
b    = spm_kron({b,I{1}});
b    = spm_permute_kron(b,dim([3,2,1]),[3,2,1]);
B{3} = spm_kron({b,I{4}});

% location dependent fatalities (in hospital or CCU): third order dependencies
%--------------------------------------------------------------------------
Pfat = P.fat*P.lnf;        % vaccinated or immune
ij   = Bij({[3 6],6:8,3,1:6},{[3 6],6:8,4,1:6},dim); B{3}(ij) = (1 - Ktrd)*Pfat;
ij   = Bij({[3 6],6:8,3,1:6},{[3 6],6:8,1,1:6},dim); B{3}(ij) = (1 - Ktrd)*(1 - Pfat);
Pfat = P.fat;              % unvaccinated
ij   = Bij({[3 6],1:5,3,1:6},{[3 6],1:5,4,1:6},dim); B{3}(ij) = (1 - Ktrd)*Pfat;
ij   = Bij({[3 6],1:5,3,1:6},{[3 6],1:5,1,1:6},dim); B{3}(ij) = (1 - Ktrd)*(1 - Pfat);


% probabilistic transitions: testing
%==========================================================================

% test probabilities
%--------------------------------------------------------------------------
b    = cell(1,dim(2));

pcr1 = P.pil(1);                               % pillar one testing   PCR
pcr2 = P.pcr*P.pil(2);                         % surveillance testing PCR
Plfd = P.pcr*sum(P.pil(3) + P.pil(4));         % surveillance testing LFD
Psen = erf(pcr1 + pcr2);                       % testing rate | susceptible PCR
Ptes = erf(pcr1*P.tes(1) + pcr2*P.tes(2));     % testing rate | infection   PCR
Plen = erf(Plfd);                              % testing rate | susceptible LFD
Ples = erf(Plfd*P.tts);                        % testing rate | infection   LFD

Sens = 1 - P.fnr(1);                           % sensitivity PCR | infected
Senc = 1 - P.fnr(2);                           % sensitivity PCR | infectious
Spec = 1 - P.fpr(1);                           % specificity PCR
Speb = 1 - P.fpr(2);                           % specificity PCR | Ab +ve
Lens = 1 - P.lnr;                              % sensitivity LFD
Lpec = 1 - P.lpr;                              % specificity LFD   
Kdel = exp(-1/P.del);                          % exp(-1/waiting period) PCR

Pcon = erf(P.con);                             % LFD PCR Confirmation

% marginal: testing {4} | susceptible {2}(1)
%-------------------------------------------------------------------------------------------------------
%    not tested            waiting                PCR+       PCR-       LFD+                  LFD-
%-------------------------------------------------------------------------------------------------------
b{1} = [(1 - Psen)*(1 - Plen) 0                   (1 - Kday) (1 - Kday) (1 - Kday)*(1 - Pcon) (1 - Kday);
        Psen*(1 - Plen)       Kdel                 0          0         (1 - Kday)*Pcon        0;
        0                    (1 - Spec)*(1 - Kdel) Kday       0          0                     0;
        0                     Spec*(1 - Kdel)      0          Kday       0                     0;
        (1 - Lpec)*Plen       0                    0          0          Kday                  0;
        Lpec*Plen             0                    0          0          0                 Kday];

% marginal: testing {4} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Ptes)*(1 - Ples) 0                   (1 - Kday) (1 - Kday) (1 - Kday)*(1 - Pcon) (1 - Kday);
        Ptes*(1 - Ples)       Kdel                 0          0         (1 - Kday)*Pcon        0;
        0                     Sens*(1 - Kdel)      Kday       0          0                     0;
        0                    (1 - Sens)*(1 - Kdel) 0          Kday       0                     0;
        Lens*Ples             0                    0          0          Kday                  0;
        (1 - Lens)*Ples       0                    0          0          0                 Kday];
    
% marginal: testing {4} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Ptes)*(1 - Ples) 0                   (1 - Kday) (1 - Kday) (1 - Kday)*(1 - Pcon) (1 - Kday);
        Ptes*(1 - Ples)       Kdel                 0          0         (1 - Kday)*Pcon        0;
        0                     Senc*(1 - Kdel)      Kday       0          0                     0;
        0                    (1 - Senc)*(1 - Kdel) 0          Kday       0                     0;
        Lens*Ples             0                    0          0          Kday                  0;
        (1 - Lens)*Ples       0                    0          0          0                 Kday];

% marginal: testing {4} | Ab+ {2}(4)
%--------------------------------------------------------------------------
b{4} = [(1 - Psen)*(1 - Plen) 0                   (1 - Kday) (1 - Kday) (1 - Kday)*(1 - Pcon) (1 - Kday);
        Psen*(1 - Plen)       Kdel                 0          0         (1 - Kday)*Pcon        0;
        0                    (1 - Speb)*(1 - Kdel) Kday       0          0                     0;
        0                     Speb*(1 - Kdel)      0          Kday       0                     0;
        (1 - Lpec)*Plen       0                    0          0          Kday                  0;
        Lpec*Plen             0                    0          0          0                 Kday];
    
% marginal: testing {4} | Ab- {2}(5)
%--------------------------------------------------------------------------
b{5} = b{1};

% marginal: testing {4} | Vaccine+ {2}(6)
%--------------------------------------------------------------------------
b{6} = b{1};

% marginal: testing {4} | Infected {2}(7)
%--------------------------------------------------------------------------
b{7} = b{2};

% marginal: testing {4} | Infected {2}(7)
%--------------------------------------------------------------------------
b{8} = b{3};

% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat_T(spm_diag(b));
b    = spm_kron({b,I{1},I{3}});
B{4} = spm_permute_kron(b,dim([4,2,1,3]),[3,2,4,1]);

% location dependent PCR testing (always when in hospital)
%--------------------------------------------------------------------------
ij   = Bij({6,1:8,1:4,1},{6,1:8,1:4,1},dim); b = B{4}(ij); B{4}(ij) = 0;
ij   = Bij({6,1:8,1:4,1},{6,1:8,1:4,2},dim); B{4}(ij) = B{4}(ij) + b;


% probability transition matrix
%==========================================================================
T     = 1;
for i = 1:4
    T = T*B{i};
end

try
    if min(T(:)) < 0; keyboard, end
    if max(T(:)) > 1; keyboard, end
end


% time-dependent parameters
%--------------------------------------------------------------------------
R.Pout = Pout;
R.Pinf = Ktic;
R.Psev = P.sev;
R.Pfat = P.fat;
R.Psen = Psen;
R.Ptes = Ptes;
R.Rvac = Rvac;
R.Tin  = P.Tin;
R.Tic  = P.Tic;
R.Trd  = P.Trd;

R.Tim  = P.Tim;
R.Tnn  = P.Tnn;
R.Vac  = P.vac;

return




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

return


function [x] = spm_cat_sym(x,d)
% Convert a cell array into a matrix - a compiled routine
% FORMAT [x] = spm_cat_sym(x,d)
% x - cell array
% d - dimension over which to concatenate [default - both]
%__________________________________________________________________________
% Empty array elements are replaced by sparse zero partitions and single 0
% entries are expanded to conform to the non-empty non zero elements.
%
% e.g.:
% > x       = spm_cat_T({eye(2) []; 0 [1 1; 1 1]})
% > full(x) =
%
%     1     0     0     0
%     0     1     0     0
%     0     0     1     1
%     0     0     1     1
%
% If called with a dimension argument, a cell array is returned.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%error('spm_cat.c not compiled - see Makefile')

% check x is not already a matrix
%--------------------------------------------------------------------------
if ~iscell(x), return, end
 
% if concatenation over a specific dimension
%--------------------------------------------------------------------------
[n,m] = size(x);
if nargin > 1
 
    % concatenate over first dimension
    %----------------------------------------------------------------------
    if d == 1
        y = cell(1,m);
        for i = 1:m
            y{i} = spm_cat_sym(x(:,i));
        end
 
    % concatenate over second
    %----------------------------------------------------------------------
    elseif d == 2
 
        y = cell(n,1);
        for i = 1:n
            y{i} = spm_cat_sym(x(i,:));
        end
 
    % only viable for 2-D arrays
    %----------------------------------------------------------------------
    else
        error('uknown option')
    end
    x      = y;
    return
 
end
 
% find dimensions to fill in empty partitions
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:m
        if iscell(x{i,j})
            x{i,j} = spm_cat_sym(x{i,j});
        end
        [u,v]  = size(x{i,j});
        I(i,j) = u;
        J(i,j) = v;
    end
end
I     = max(I,[],2);
J     = max(J,[],1);
 
% sparse and empty partitions
%--------------------------------------------------------------------------
[n,m] = size(x);
for i = 1:n
    for j = 1:m
        if isempty(x{i,j})
            x{i,j} = sparse(I(i),J(j));
        end
    end
end
 
% concatenate
%--------------------------------------------------------------------------
for i = 1:n
    y{i,1} = cat(2,x{i,:});
end
try
    x = sparse(cat(1,y{:}));
catch
    x = cat(1,y{:});
end


