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
% $Id: spm_COVID_T.m 8239 2022-04-09 12:45:02Z karl $

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
Piss = erf(2);                       % probability of self-isolation
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
b{2} = [(1 - Piss) 0          (1 - Piss) 0         (1 - Piss) (1 - Piss);
        0          (1 - Piss) 0          0          0         0;
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

Kimm = exp(-1/P.Tim);                      % loss of Ab+ immunity (per day)
Kinn = exp(-1/P.Tnn);                      % loss of Ab- immunity (per day)
Kvac = exp(-1/P.vac);                      % Loss of Ab+ vaccine  (per day)

Kinf = exp(-1/P.Tin);                      % infection rate
Kcon = exp(-1/P.Tcn);                      % infectious rate
Pres = 1 - erf(2*P.res)/2;                 % non-resistant proportion
Pvef = erf(P.vef);                         % vaccination: Infection risk 
Prev = P.ves*Pres;                         % vaccination: Infectious risk

Pnac = P.nac;                              % 1 - vaccination rate
Rvac = 1 - Pnac;                           % vaccination rate

% marginal: infection {2} | home {1}(1)
%--------------------------------------------------------------------------------------------------------------------------------------------------
%    susceptible     infected            infectious      Ab+            Ab-            Vaccine+            Infected (vac)       Infectious  (vac)
%--------------------------------------------------------------------------------------------------------------------------------------------------
Pinv = (1 - Ptin)*Pvef;
b{1} = [Ptin*Pnac        0                     0          0              (1 - Kinn)*Pnac  0                     0                     0; 
        (1 - Ptin)*Pnac  Kinf                  0          0               0               0                     0                     0;
        0               (1 - Kinf)*Pres        Kcon       0               0               0                     0                     0;
        0                0                    (1 - Kcon)  Kimm*Pnac       0               0                     0                     0;
        0               (1 - Kinf)*(1 - Pres)  0         (1 - Kimm)*Pnac  Kinn*Pnac      (1 - Pinv)*(1 - Kvac)  0                     0;
        Rvac             0                     0          Rvac            Rvac           (1 - Pinv)*Kvac       (1 - Kinf)*(1 - Prev) (1 - Kcon);
        0                0                     0          0               0               Pinv                  Kinf                  0
        0                0                     0          0               0               0                    (1 - Kinf)*Prev        Kcon];
    
% marginal: infection {2} | work {1}(2)
%--------------------------------------------------------------------------
Ptin = Ptou;
Pinv = (1 - Ptin)*Pvef;
b{2} = [Ptin*Pnac        0                     0          0              (1 - Kinn)*Pnac  0                     0                     0; 
        (1 - Ptin)*Pnac  Kinf                  0          0               0               0                     0                     0;
        0               (1 - Kinf)*Pres        Kcon       0               0               0                     0                     0;
        0                0                    (1 - Kcon)  Kimm*Pnac       0               0                     0                     0;
        0               (1 - Kinf)*(1 - Pres)  0         (1 - Kimm)*Pnac  Kinn*Pnac      (1 - Pinv)*(1 - Kvac)  0                     0;
        Rvac             0                     0          Rvac            Rvac           (1 - Pinv)*Kvac       (1 - Kinf)*(1 - Prev) (1 - Kcon);
        0                0                     0          0               0               Pinv                  Kinf                  0
        0                0                     0          0               0               0                    (1 - Kinf)*Prev        Kcon];

% marginal: infection {2} | ccu {1}(3)
%--------------------------------------------------------------------------
Ptin = 1;
Pinv = (1 - Ptin)*Pvef;
b{3} = [Ptin*Pnac        0                     0          0              (1 - Kinn)*Pnac  0                     0                     0; 
        (1 - Ptin)*Pnac  Kinf                  0          0               0               0                     0                     0;
        0               (1 - Kinf)*Pres        Kcon       0               0               0                     0                     0;
        0                0                    (1 - Kcon)  Kimm*Pnac       0               0                     0                     0;
        0               (1 - Kinf)*(1 - Pres)  0         (1 - Kimm)*Pnac  Kinn*Pnac      (1 - Pinv)*(1 - Kvac)  0                     0;
        Rvac             0                     0          Rvac            Rvac           (1 - Pinv)*Kvac       (1 - Kinf)*(1 - Prev) (1 - Kcon);
        0                0                     0          0               0               Pinv                  Kinf                  0
        0                0                     0          0               0               0                    (1 - Kinf)*Prev        Kcon];

% marginal: infection {2} | removed {1}(4)
%--------------------------------------------------------------------------
b{4} = b{3};

% marginal: infection {2} | isolated {1}(5)
%--------------------------------------------------------------------------
b{5} = b{3};

% marginal: infection {2} | hospital {1}(6)
%--------------------------------------------------------------------------
b{6} = b{3};
    
% kroneckor form
%--------------------------------------------------------------------------
if max(spm_vec(b)) > 1; keyboard, end

b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{3},I{4}});
B{2} = spm_permute_kron(b,dim([2,1,3,4]),[2,1,3,4]);

% probabilistic transitions: clinical
%==========================================================================

% probabilities of developing symptoms
%--------------------------------------------------------------------------
b    = cell(1,dim(2));
Psev = P.sev;                                % P(ARDS | infected)
Pfat = 1 - (1 - P.fat)*P.oth;                % P(fatality | ARDS, Other)
Ktic = exp(-1/P.Tic);                        % asymptomatic rate
Ksym = exp(-1/P.Tsy);                        % symptomatic rate
Ktrd = exp(-1/P.Trd);                        % ARDS rate

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
b{3} = [Ktic       (1 - Ksym)*(1 - Psev) (1 - Ktrd)*(1 - Pfat) (1 - Kday);
       (1 - Ktic)   Ksym                  0                     0;
        0          (1 - Ksym)*Psev        Ktrd                  0;
        0           0                    (1 - Ktrd)*Pfat        Kday];

% marginal: clinical {3} | Ab+ {2}(4)
%--------------------------------------------------------------------------
b{4} = b{1};
    
% marginal: clinical {3} | Ab- {2}(5)
%--------------------------------------------------------------------------
b{5} = b{1};

% link between pathogenicity and vaccination
%--------------------------------------------------------------------------
Psev = Psev*P.lnk;         % vaccine efficiency: pathogenicity
Pfat = Pfat*P.lnf;         % vaccine efficiency: fatality

% marginal: clinical {3} | Vaccine+ {2}(6)
%--------------------------------------------------------------------------
b{6} = [1          (1 - Ksym)*(1 - Psev) (1 - Ktrd)*(1 - Pfat) (1 - Kday);
        0           Ksym                  0                     0;
        0          (1 - Ksym)*Psev        Ktrd                  0;
        0           0                    (1 - Ktrd)*Pfat        Kday];

% marginal: clinical {3} | infected {2}(7)
%--------------------------------------------------------------------------
b{7} = [Ktic       (1 - Ksym)*(1 - Psev) (1 - Ktrd)*(1 - Pfat) (1 - Kday);
       (1 - Ktic)   Ksym                  0                     0;
        0          (1 - Ksym)*Psev        Ktrd                  0;
        0           0                    (1 - Ktrd)*Pfat        Kday];
    
% marginal: clinical {3} | Infectious {2}(8)
%--------------------------------------------------------------------------
b{8} = b{7};
    
% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{1}});
b    = spm_permute_kron(b,dim([3,2,1]),[3,2,1]);
B{3} = spm_kron({b,I{4}});

% location dependent fatalities (in hospital or CCU): third order dependencies
%--------------------------------------------------------------------------
Pfat = P.fat*P.lnf;        % vaccinated
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

% fluctuations in testing rate: Gaussian basis functions
%--------------------------------------------------------------------------
Ppcr = 0;
dt   = 48;
if isfield(P,'pcr')
    for i = 1:numel(P.pcr)
        Ppcr = Ppcr + log(P.pcr(i)) * exp(-(P.t - i*dt).^2./((dt/2)^2));
    end
end
Ppcr = exp(erf(Ppcr)/4);

% (pillar 1, 2 and LFD) phases of testing
%--------------------------------------------------------------------------
pill   = zeros(1,numel(P.lim));
for i = 1:numel(P.lim)
    pill(i) = P.lim(i)*spm_phi((P.t - P.ons(i))/P.rat(i));
end

pcr1 = pill(1);                                % pillar one testing   PCR
pcr2 = Ppcr*pill(2);                           % surveillance testing PCR
Plfd = Ppcr*sum(pill(3) + pill(4));            % surveillance testing LFD
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
if min(spm_vec(b)) < 0; keyboard, end

b    = spm_cat(spm_diag(b));
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

% if min(spm_vec(T)) < 0; keyboard, end
% if max(spm_vec(T)) < 0; keyboard, end


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

