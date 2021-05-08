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
% $Id: spm_COVID_T.m 8101 2021-05-08 15:01:43Z karl $

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

% seasonal and monotonic fluctuations
%--------------------------------------------------------------------------
Q    = (1 + cos(2*pi*P.t/365))/2;
S    = exp(-P.t/512);

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
Pnhs = 0.03;                         % P(hospital | NHS/care worker)
Pcap = erf(P.ccu);                   % P(transfer to CCU | ARDS)
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
        Ph2r       0          0         (1 - Pexp)  0         0;
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
       (1 - Phos) (1 - Phos)  0         (1 - Phos) (1 - Phos) 0;
        Phos       Phos       0          Phos       Phos     (1 - Pcap)];


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

Ptin = P.tin;                              % P(no transmission) | home
Ptou = P.tou;                              % P(no transmission) | work
Pths = P.ths;                              % P(no transmission) | hospital

Kimm = exp(-1/P.Tim);                      % loss of Ab+ immunity (per day)
Kinn = exp(-1/256);                        % loss of Ab- immunity (per day)
Vimm = exp(-1/512);                        % Loss of Ab+ vaccine  (per day)

Kinf = exp(-1/P.Tin);                      % infection rate
Kcon = exp(-1/P.Tcn);                      % infectious rate
Pres = P.res;                              % resistant proportion

Rvac = exp(-((P.t-P.rol(2))/P.rol(3))^2);  % vaccination rollout
Rvac = Rvac*P.rol(1);                      % P(vaccination)
Rvac = Rvac*P.vac;                         % P(effective vaccination)
Pnac = 1 - Rvac;                           % 1 - P(vaccination)


% marginal: infection {2} | home {1}(1)
%--------------------------------------------------------------------------
%    susceptible     infected            infectious      Ab+             Ab-            Vaccine+
%--------------------------------------------------------------------------
b{1} = [Ptin*Pnac        0                     0          0              (1 - Kinn)*Pnac    0;
        (1 - Ptin)*Pnac  Kinf                  0          0               0                 0;
        0               (1 - Pres)*(1 - Kinf)  Kcon       0               0                 0;
        0                0                    (1 - Kcon)  Kimm            0                 0;
        0                Pres*(1 - Kinf)       0         (1 - Kimm)       Kinn*Pnac (1 - Vimm);
        Rvac             0                     0          0               Rvac           Vimm];
    
% marginal: infection {2} | work {1}(2)
%--------------------------------------------------------------------------
b{2}      = b{1};
b{2}(:,1) = [Ptou*Pnac;
            (1 - Ptou)*Pnac;
             0;
             0;
             0;
             Rvac];


% marginal: infection {2} | ccu {1}(3)
%--------------------------------------------------------------------------
b{3}      = b{1};
b{3}(:,1) = [1;
             0;
             0;
             0;
             0;
             0];

% marginal: infection {2} | removed {1}(4)
%--------------------------------------------------------------------------
b{4}      = b{1};
b{4}(:,1) = [Pnac;
             0;
             0;
             0;
             0;
             Rvac];

% marginal: infection {2} | isolated {1}(5)
%--------------------------------------------------------------------------
b{5}      = b{3};

% marginal: infection {2} | hospital {1}(6)
%--------------------------------------------------------------------------
b{6}      = b{1};
b{6}(:,1) = [Pths*Pnac;
            (1 - Pths)*Pnac;
             0;
             0;
             0;
             Rvac];
         
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
Psev = erf(P.sev*Q + P.lat*(1 - Q));         % P(ARDS | infected)
Pfat = erf(P.fat*S + P.sur*(1 - S));         % P(fatality | ARDS, CCU)
Ksym = exp(-1/P.Tsy);                        % acute symptomatic rate
Ktrd = exp(-1/P.Trd);                        % acute RDS rate
Ktic = exp(-1/P.Tic(1));                     % asymptomatic rate
Ktis = exp(-1/P.Tic(2));                     % asymptomatic rate
Ktic = (1 - Ktic);                           % symptomatic rate | infected
Ktis = (1 - Ktis);                           % symptomatic rate | infectious


% marginal: clinical {3} | susceptible {2}(1)
%--------------------------------------------------------------------------
%  asymptomatic    symptomatic           ARDS        deceased
%--------------------------------------------------------------------------
b{1} = [1          (1 - Ksym)*(1 - Psev) 0         (1 - Kday);
        0           Ksym                 0                  0;
        0          (1 - Ksym)*Psev       Ktrd               0;
        0           0                   (1 - Ktrd)      Kday];
    
% marginal: clinical {3} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Ktic) (1 - Ksym)*(1 - Psev) 0         (1 - Kday);
        Ktic        Ksym                 0                 0;
        0          (1 - Ksym)*Psev       Ktrd              0;
        0           0                   (1 - Ktrd)      Kday];
    
% marginal: clinical {3} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Ktis) (1 - Ksym)*(1 - Psev) 0         (1 - Kday);
        Ktis        Ksym                 0                 0;
        0          (1 - Ksym)*Psev       Ktrd              0;
        0           0                   (1 - Ktrd)      Kday];
    
% marginal: clinical {3} | Ab+ {2}(4)
%--------------------------------------------------------------------------
b{4} = b{1};
    
% marginal: clinical {3} | Ab- {2}(5)
%--------------------------------------------------------------------------
b{5} = b{1};

% marginal: clinical {3} | Vaccine+ {2}(5)
%--------------------------------------------------------------------------
b{6} = b{1};
    

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

% fluctuations in testing rate (Fourier basis functions)
%--------------------------------------------------------------------------
Ppcr = 0;                             
if isfield(P,'pcr')
    for i = 1:numel(P.pcr)
        if i > numel(P.pcr)/2
            Ppcr = Ppcr + log(P.pcr(i)) * cos(2*i*pi*P.t/365)/8;
        else
            Ppcr = Ppcr + log(P.pcr(i)) * sin(2*i*pi*P.t/365)/8;
        end
    end
end
Ppcr = exp(Ppcr);


% (pillar 1, 2 and LFD) phases of testing
%--------------------------------------------------------------------------
% https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/943187/S0925_Innova_Lateral_Flow_SARS-CoV-2_Antigen_test_accuracy.pdf
% https://www.bmj.com/content/bmj/369/bmj.m1808.full.pdf
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
Lens = 0.489;                                  % sensitivity LFD
Lpec = 0.9993;                                 % specificity LFD   
Kdel = exp(-1/P.del);                          % exp(-1/waiting period) PCR

% marginal: testing {4} | susceptible {2}(1)
%--------------------------------------------------------------------------
%    not tested            waiting                PCR+       PCR-       LFD+      LFD-
%--------------------------------------------------------------------------
b{1} = [(1 - Psen)*(1 - Plen) 0                   (1 - Kday) (1 - Kday)  0  (1 - Kday);
        Psen*(1 - Plen)       Kdel                 0          0         (1 - Kday)  0;
        0                    (1 - Spec)*(1 - Kdel) Kday       0          0          0;
        0                     Spec*(1 - Kdel)      0          Kday       0          0;
        (1 - Lpec)*Plen       0                    0          0          Kday       0;
        Lpec*Plen             0                    0          0          0      Kday];

% marginal: testing {4} | infected {2}(2)
%--------------------------------------------------------------------------
b{2} = [(1 - Ptes)*(1 - Ples) 0                   (1 - Kday) (1 - Kday)  0  (1 - Kday);
        Ptes*(1 - Ples)       Kdel                 0          0         (1 - Kday)  0;
        0                     Sens*(1 - Kdel)      Kday       0          0          0;
        0                    (1 - Sens)*(1 - Kdel) 0          Kday       0          0;
        Lens*Ples             0                    0          0          Kday       0;
        (1 - Lens)*Ples       0                    0          0          0       Kday];
    
% marginal: testing {4} | infectious {2}(3)
%--------------------------------------------------------------------------
b{3} = [(1 - Ptes)*(1 - Ples) 0                   (1 - Kday) (1 - Kday)  0  (1 - Kday);
        Ptes*(1 - Ples)       Kdel                 0          0         (1 - Kday)  0;
        0                     Senc*(1 - Kdel)      Kday       0          0          0;
        0                    (1 - Senc)*(1 - Kdel) 0          Kday       0          0;
        Lens*Ples             0                    0          0          Kday       0;
        (1 - Lens)*Ples       0                    0          0          0      Kday];

% marginal: testing {4} | Ab+ {2}(4)
%--------------------------------------------------------------------------
b{4} = [(1 - Psen)*(1 - Plen) 0                   (1 - Kday) (1 - Kday)  0  (1 - Kday);
        Psen*(1 - Plen)       Kdel                 0          0         (1 - Kday)  0;
        0                    (1 - Speb)*(1 - Kdel) Kday       0          0          0;
        0                     Speb*(1 - Kdel)      0          Kday       0          0;
        (1 - Lpec)*Plen       0                    0          0          Kday       0;
        Lpec*Plen             0                    0          0          0      Kday];
    
% marginal: testing {4} | Ab- {2}(5)
%--------------------------------------------------------------------------
b{5} = b{1};

% marginal: testing {4} | Vaccine+ {2}(5)
%--------------------------------------------------------------------------
b{6} = b{1};

% kroneckor form
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{1},I{3}});
B{4} = spm_permute_kron(b,dim([4,2,1,3]),[3,2,4,1]);

% location dependent PCR testing (none when removed)
%--------------------------------------------------------------------------
ij   = Bij({4,1:5,1:4,1},{4,1:5,1:4,1},dim); B{4}(ij) = 1;
ij   = Bij({4,1:5,1:4,1},{4,1:5,1:4,2},dim); B{4}(ij) = 0;

% location dependent PCR testing (always when ill in hospital)
%--------------------------------------------------------------------------
ij   = Bij({6,1:5,2:3,1},{6,1:5,2:3,1},dim); B{4}(ij) = 0;
ij   = Bij({6,1:5,2:3,1},{6,1:5,2:3,2},dim); B{4}(ij) = 1;


% probability transition matrix
%==========================================================================
T     = 1;
for i = 1:4
    T =  T*B{i};
end

% time-dependent parameters
%--------------------------------------------------------------------------
R.Pout = Pout;
R.Psev = Psev;
R.Pinf = Ktic;
R.Pfat = Pfat;
R.Psen = Psen;
R.Ptes = Ptes;

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

