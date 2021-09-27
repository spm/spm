function [y,x,z,W] = spm_SARS_gen(P,M,U,NPI,age)
% Generate predictions and hidden states of a COVID model
% FORMAT [Y,X,Z,W] = spm_SARS_gen(P,M,U,NPI,age)
% P    - model parameters
% M    - model structure (M.T - length of timeseries or data structure)
% U    - number of output variables [default: 2] or indices e.g., [4 5]
% NPI  - nonpharmaceutical intervention
%     NPI(i).period = {'dd-mm-yyyy','dd-mm-yyyy'}; % dates of epidemic
%     NPI(i).param  = {'xyz',...};                 % parameter name
%     NPI(i).Q      = (value1,...);                % parameter name
%     NPI(i).dates  = {'dd-mm-yyyy','dd-mm-yyyy'}; % dates of interevention
% age  - indices of age band (0 for average)
%
% Y(:,1)  - Daily deaths (28 days)
% Y(:,2)  - Daily confirmed cases (PCR and LFD)
% Y(:,3)  - Mechanical ventilation
% Y(:,4)  - Reproduction ratio (R)
% Y(:,5)  - Seropositive immunity (%)
% Y(:,6)  - PCR testing rate
% Y(:,7)  - Contagion risk (%)
% Y(:,8)  - Prevalence: infectious {%}
% Y(:,9)  - Daily contacts
% Y(:,10) - Daily incidence (%)
% Y(:,11) - Prevalence: infection (%)
% Y(:,12) - Number symptomatic
% Y(:,13) - Mobility (%)
% Y(:,14) - Workplace (%)
% Y(:,15) - Certified deaths
% Y(:,16) - Hospital admissions
% Y(:,17) - Hospital deaths
% Y(:,18) - Non-hospital deaths
% Y(:,19) - Daily incidence (per hundred thousand)
% Y(:,20) - Weekly confirmed cases (per hundred thousand)
% Y(:,21) - Infection fatality ratio (%)
% Y(:,22) - Percent vaccinated (%)
% Y(:,23) - PCR case positivity (%)
% Y(:,24) - Lateral flow tests
% Y(:,25) - Cumulative attack rate
% Y(:,26) - Population immunity (total)
% Y(:,27) - Hospital occupancy
% Y(:,28) - Incidence of long COVID
% Y(:,29) - population immunity (vaccine)
% Y(:,30) - cumulative admissions
%
% X       - (M.T x 4) marginal densities over four factors
% location   : {'home','out','ccu','removed','isolated','hospital'};
% infection  : {'susceptible','infected','infectious','Ab +ve','Ab -ve','Vac +ve','Vac Inf'};
% clinical   : {'asymptomatic','symptoms','ARDS','death'};
% diagnostic : {'untested','waiting','PCR +ve','PCR -ve','LFD +ve','LFD -ve'}
%
% Z{t} - joint density over hidden states at the time t
% W    - structure containing time varying parameters
%
% This function returns data Y and their latent states or causes X, given
% the parameters of a generative model. This model is a mean field
% approximation based upon population or density dynamics with certain
% conditional dependencies among the marginal densities over four factors.
% See SPM_covid_priors details. In brief, this routine transforms model
% parameters to (exponentiated) scale parameters and then generates a
% sequence of jointed densities over four factors, after assembling a state
% dependent probability transition matrix. The number in the timeseries is
% specified by M.T.
%
% Equipped with a time-dependent ensemble density, outcome measures are
% then generated as expected values. These include the rate of (new) deaths
% and cases per day. This routine can be extended to generate other
% outcomes, or indeed consider other factorisations of the probability
% transition matrices. The subroutine (spm_COVID_T) creating the
% probability transition matrices given the current states and model
% parameters defines the generative model. This model structure rests upon
% a mean field approximation to the transition probabilities that,
% crucially, depends upon (usually the marginal) densities in question.
% Working through the code below will show how this model is constructed.
%
% A more detailed description of the generative model can be found in the
% body of the script.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_SARS_gen.m 8156 2021-09-27 09:05:29Z karl $


% The generative model:
%==========================================================================
% In brief, this model generates timeseries data based on a mean field
% approximation to ensemble or population dynamics. The implicit
% probability distributions are over four latent factors, each with several
% levels or states. These factors are sufficient to generate expected
% outcomes; for example, the number of new cases or the number of people
% infected. The first factor is the location of an individual, who can be
% at home, at work, in a critical care unit (CCU), self isolated or in the
% morgue. The second factor is infection status; namely, susceptible to
% infection, infected, infectious or immune. In addition, we include a
% resistant state that does not participate in the transmission of the
% virus. This model assumes that there is a progression from a state of
% susceptibility to immunity, through a period of (pre-contagious)
% infection to an infectious (contagious) status. The third factor is the
% clinical status; namely, asymptomatic, symptomatic, acute respiratory
% distress syndrome (ARDS) or deceased. Again, there is an assumed
% progression from asymptomatic to ARDS, where people with ARDS can either
% recover to an asymptomatic state or not. Finally, the fourth factor
% represents the diagnostic or testing status of. An individual can be
% untested or waiting for the results of a test that can either be positive
% or negative. With this setup, one can be in one of five places, with any
% infectious status, expressing symptoms or not and having test results or
% not. Note that - in this construction - it is possible to be infected and
% yet be asymptomatic. However, the marginal distributions are not
% independent, in virtue of the dynamics that describe the transition among
% states within each factor. Crucially, the transitions within any factor
% depend upon the marginal distribution of other factors. For example, the
% probability of becoming infected, given that one is susceptible to
% infection, depends upon whether one is at home or at work. Similarly, the
% probability of developing symptoms depends upon whether one is infected
% or not. The probability of being tested depends upon whether one is
% symptomatic. Finally, to complete the circular dependency, the
% probability of leaving home to go to work depends upon the number of
% infected people in the population - as a result of social distancing
% (please see main text). These conditional dependencies constitute the
% mean field approximation and enable the dynamics to be solved or
% integrated over time. At any one point in time, the probability of being
% in any combination of the four factors determines what would be observed
% at the population level. For example, the occupancy of the deceased level
% of the clinical factor determines the current number of people who
% contribute to daily deaths. Similarly, the occupancy of the positive
% level of the testing factor determines the daily positive cases reported.

% References
% neutralising antibodies : https://www.nature.com/articles/s41586-020-2012-7
%--------------------------------------------------------------------------
% seropositive: https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30196-1/fulltext
%--------------------------------------------------------------------------
% For 16 patients with serum samples available 14 days or longer after
% symptom onset, rates of seropositivity were 94% for anti-NP IgG (n=15),
% 88% for anti-NP IgM (n=14), 100% for anti-RBD IgG (n=16), and 94% for
% anti-RBD IgM (n=15). Anti-SARS-CoV-2-NP or anti-SARS-CoV-2-RBD IgG levels
% correlated with virus neutralisation titre (R2>0.9). No genome mutations
% were detected on serial samples.
%--------------------------------------------------------------------------
% immunity : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2271881/
%--------------------------------------------------------------------------
% In this group, antibody concentrations started to increase 1 week after
% inoculation and reached a maximum about 1 week later. Thereafter antibody
% titres slowly declined. Although concentrations were still slightly
% raised 1 year later, this did not always prevent reinfection when
% volunteers were then challenged with the homologous virus.
%--------------------------------------------------------------------------
% long-lasting immunity : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2851497/
%--------------------------------------------------------------------------
% Among 176 patients who had had severe acute respiratory syndrome (SARS),
% SARS-specific antibodies were maintained for an average of 2 years, and
% significant reduction of immunoglobulin G-positive percentage and titers
% occurred in the third year. Thus, SARS patients might be susceptible to
% reinfection >3 years after initial exposure.


% setup and defaults (assume new deaths and cases as outcome variables)
%--------------------------------------------------------------------------
if (nargin < 3) || isempty(U), U = 1:2; end         % two outcomes
if (nargin < 4), NPI = [];              end         % interventions
if (nargin < 5), age = 0*U;             end         % age bands
try, M.T; catch, M.T = 180;             end         % over six months

% deal with data structures (asynchronous timeseries)
%--------------------------------------------------------------------------
if isstruct(M.T)        % predictions from multiple data types are required
    
    % extract data structure and specify temporal domain
    %----------------------------------------------------------------------
    D   = M.T;
    d   = spm_vec(D.date);
    if isfield(M,'date')
        try
            d0 = datenum(M.date,'dd-mm-yyyy');
        catch
            d0 = M.date;
        end
    else
        d0 = min(d);
    end
    M.T = max(d) - d0 + 1;
    d   = d0:max(d);
    U   = [D.U];
    age = [D.age];
    
end

% unpack and exponentiate parameters
%--------------------------------------------------------------------------
param = fieldnames(P);
nN    = numel(P.N);
N     = zeros(nN,1);

Q     = cell(nN,1);
R     = cell(nN,1);
x     = cell(nN,1);
p     = cell(nN,1);
for n = 1:nN
    
    % create cell array of parameters for each group
    %----------------------------------------------------------------------
    for j = 1:numel(param)
        if size(P.(param{j}),1) > 1
            R{n}.(param{j}) = P.(param{j})(n,:);
        else
            R{n}.(param{j}) = P.(param{j})(1,:);
        end
    end
    
    % unpack and exponentiate parameters
    %----------------------------------------------------------------------
    Q{n}  = spm_vecfun(R{n},@exp);
    R{n}  = Q{n};                 % baseline parameters
    
    
    % initial marginals (Dirichlet parameters)
    %----------------------------------------------------------------------
    c    = Q{n}.n;                % proportion of initial cases
    N(n) = Q{n}.N*1e6;            % population size
    r    = Q{n}.r;                % proportion of resistant cases
    h    = Q{n}.o*3/4;            % proportion at home
    w    = Q{n}.o*1/4;            % proportion at work
    m    = 1 - Q{n}.o;            % proportion of unexposed cases
    s    = 1 - c - r;             % proportion of susceptible cases

    p{n}{1} = [h w 0 m 0 0]';     % location
    p{n}{2} = [s c 0 0 r 0 0 0]'; % infection
    p{n}{3} = [1 0 0 0]';         % clinical
    p{n}{4} = [1 0 0 0 0 0]';     % testing
    
    % normalise initial marginals
    %----------------------------------------------------------------------
    Nf    = numel(p{n});
    for f = 1:Nf
        p{n}{f}  = p{n}{f}/sum(p{n}{f});
    end
    
end

% identity matrices for each latent factor
%--------------------------------------------------------------------------
for i = 1:numel(p{n})
    I{i} = speye(numel(p{n}{i}));
end

% initial ensemble density (equilibrium)
%--------------------------------------------------------------------------
for n = 1:nN
    x{n}  = spm_cross(p{n});
    P     = R{n};
    P.n   = 0;
    P.t   = 1;
    P.tin = 1;    % P(no transmission | home)
    P.tou = 1;    % P(no transmission | work)
    P.ths = 1;    % P(no transmission | hospital)
    P.nac = 1;    % P(no vaccination)
    for i = 1:8
        T    = spm_COVID_T(P,I);
        x{n} = spm_unvec(T*spm_vec(x{n}),x{n});
        x{n} = x{n}/sum(x{n}(:));
    end
    
    % normalise initial marginals
    %----------------------------------------------------------------------
    Nf    = numel(p{n});
    for f = 1:Nf
        p{n}{f}  = p{n}{f}/sum(p{n}{f});
    end
end

    
% outputs that depend upon population size
%--------------------------------------------------------------------------
uN    = [1,2,3,6,9,12,15,16,17,18,24,27,28,30];

% ensemble density tensor and solve over the specified number of days
%--------------------------------------------------------------------------
Y     = cell(nN,1);                    % outputs
X     = cell(nN,1);                    % time series of marginal densities
Z     = cell(nN,1);                    % joint densities at each time point
W     = cell(nN,1);                    % time-dependent parameters
r     = cell(nN,1);                    % probability of lockdown levels
pvac  = zeros(n,1);
for n = 1:nN
    r{n} = [1;0];
end
for i = 1:M.T
    for n = 1:nN
        
        % time-dependent parameters
        %==================================================================
        
        % nonpharmacological interventions (NPI)
        %------------------------------------------------------------------
        for j = 1:numel(NPI)
            
            % start and end dates
            %--------------------------------------------------------------
            dstart = datenum(NPI(j).dates{1},'dd-mm-yyyy') - datenum(NPI(j).period{1},'dd-mm-yyyy');
            dfinal = datenum(NPI(j).dates{2},'dd-mm-yyyy') - datenum(NPI(j).period{1},'dd-mm-yyyy');
            if (i > dstart) && (i <= dfinal)
                if ischar(NPI(j).param)
                    if isnumeric(NPI(j).Q)
                        Q{n}.(NPI(j).param) = NPI(j).Q;
                    else
                        Q{n}.(NPI(j).param) = NPI(j).Q{n};
                    end
                else
                    for k = 1:numel(NPI(j).param)
                        if isnumeric(NPI(j).Q)
                            Q{n}.(NPI(j).param{k})(1) = NPI(j).Q(k);
                        else
                            if size(NPI(j).Q{k},1) > 1 && numel(NPI(j).param) == numel(NPI(j).Q)
                                Q{n}.(NPI(j).param{k}) = NPI(j).Q{k}(n,:);
                            else
                                Q{n}.(NPI(j).param{k})(1) = NPI(j).Q{n}(k);
                            end
                        end
                    end
                end
            else
                if ischar(NPI(j).param)
                    Q{n}.(NPI(j).param) = R{n}.(NPI(j).param);
                else
                    for k = 1:numel(NPI(j).param)
                        Q{n}.(NPI(j).param{k}) = R{n}.(NPI(j).param{k});
                    end
                end
            end
        end
        

        % coupling between groups (contact rates)
        %==================================================================
        
        % fluctuations in contact rates (mobility): radial basis functions
        %------------------------------------------------------------------
        Rout = 0;
        if isfield(P,'mob')
            for j = 1:numel(Q{n}.mob)
                Rout = Rout + log(Q{n}.mob(j)) * exp(-(i - j*32).^2/256)/8;
            end
        end
        
        % probability of lockdown (a function of prevalence)
        %------------------------------------------------------------------
        q    = p{n}{2}(3)*Q{n}.sde;               % prevalence dependent
        k1   = p{n}{2}(1)*q;                      % contact rates
        k2   = exp(-1/Q{n}.qua);                  % relaxation
        r{n} = [(1 - k1) (1 - k2);
                     k1,      k2]*r{n};
        
        Pout     = r{n}(1)^Q{n}.s;                % P(going out)
        Pout     = Pout*exp(Rout);                % add fluctuations
        Q{n}.out = Pout*R{n}.out;                 % scale by baseline
       
        % seasonal variation in transmission risk
        %------------------------------------------------------------------
        S    = (1 + cos(2*pi*(i - log(Q{n}.inn)*8)/365))/2;
        
        % and fluctuations in transmissibility: cumulative error functions 
        %------------------------------------------------------------------
        Ptra = 1;
        if isfield(Q{n},'tra')
            for j = 1:numel(Q{n}.tra)
                Ptra = Ptra + erf(Q{n}.tra(j)) * (1 + erf((i - j*64)/32))/16;
            end
        end
        Ptrn = Q{n}.trn + Q{n}.trm*S;            % seasonal risk
        Ptrn = erf(Ptrn*Ptra);                   % fluctuating risk
        
        % contact rates
        %------------------------------------------------------------------
        tin  = 1;
        tou  = 1;
        ths  = 1;
        for j = 1:nN
            q   = spm_sum(x{j},[3 4]);
            pin = q(1,:)/sum(q(1,:) + eps);      % P(infection | home)
            pou = q(2,:)/sum(q(2,:) + eps);      % P(infection | work)
            phs = q(6,:)/sum(q(6,:) + eps);      % P(infection | hospital)
            
            pin = pin(3) + pin(8);               % P(infectious | home)
            pou = pou(3) + pou(8);               % P(infectious | work)
            phs = phs(3) + phs(8);               % P(infectious | hospital)

            tin = tin*(1 - Ptrn*pin)^Q{n}.Nin(j);
            tou = tou*(1 - Ptrn*pou)^Q{n}.Nou(j);
            ths = ths*(1 - Ptrn*phs)^Q{n}.Nin(j);
        end
        
        Q{n}.tin = min(tin,1);          % P(no transmission | home)
        Q{n}.tou = min(tou,1);          % P(no transmission | work)
        Q{n}.ths = min(ths,1);          % P(no transmission | hospital)
        
        % vaccination rollout: nac = 1 - vaccination rate
        %------------------------------------------------------------------
        Rvac     = Q{n}.rol(1)*(1 + erf((i - Q{n}.rol(2))/Q{n}.rol(3))) + ...
                   Q{n}.fol(1)*(1 + erf((i - Q{n}.fol(2))/Q{n}.fol(3)));
        Q{n}.nac = 1 - Rvac;     
        
        % exponential decay and time-dependent clinical susceptibility
        %------------------------------------------------------------------
        Q{n}.t   = i;                   % time (days)
        S        = exp(-i/512);
        Q{n}.sev = erf(R{n}.sev*S + R{n}.lat*(1 - S));  % P(ARDS | infected)
        Q{n}.fat = erf(R{n}.fat*S + R{n}.sur*(1 - S));  % P(fatality | ARDS, CCU)
        
        
        % update ensemble density (x)
        %==================================================================
        [T,V]  = spm_COVID_T(Q{n},I);
        x{n}   = spm_unvec(T*spm_vec(x{n}),x{n});
        x{n}   = x{n}/sum(x{n}(:));
        
        % marginal densities (p)
        %------------------------------------------------------------------
        p{n}  = spm_marginal(x{n});
        for j = 1:Nf
            X{n,j}(i,:) = p{n}{j};
        end
        
        % outcomes
        %==================================================================
        
        % probability of a positive PCR test within 28 days
        %------------------------------------------------------------------
        pcr28   = (1 - p{n}{4}(1)^28);
        pcr14   = (1 - p{n}{4}(1)^14);
        
        % time-varying parameters and other vaiables
        %------------------------------------------------------------------
        V.Ptrn  = Ptrn;
        V.pcr28 = pcr28;
        V.pcr14 = pcr14;
        W{n}(i) = V;
        
        % number of daily deaths (28 days)
        %------------------------------------------------------------------
        q         = erf(Q{n}.rel);
        Y{n}(i,1) = N(n) * p{n}{3}(4)*(q + (1 - q)*pcr28);
 
        % number of daily (positive) tests (PCR and LFD)
        %------------------------------------------------------------------
        Y{n}(i,2) = N(n) * (p{n}{4}(3) + p{n}{4}(5));
        
        % CCU bed occupancy (mechanical ventilation)
        %------------------------------------------------------------------
        Y{n}(i,3) = N(n) * p{n}{1}(3);
        
        % effective reproduction ratio (R) (based on infection prevalence)
        %------------------------------------------------------------------
        Y{n}(i,4) = p{n}{2}(3) + p{n}{2}(8);
        
        % seropositive immunity (%) Ab+ and Vaccine+
        %------------------------------------------------------------------
        Y{n}(i,5) = 100 * sum(p{n}{2}([4 6 7 8]));
        
        % total number of daily virus tests (PCR and LFD)
        %------------------------------------------------------------------
        Y{n}(i,6) = N(n) * sum(p{n}{4}(3:6));
        
        % probability of contracting virus (in a class of 15)
        %------------------------------------------------------------------
        Y{n}(i,7) = (1 - (1 - Q{n}.trn*p{n}{2}(3))^15) * 100;
        
        % prevalence of (contagious) infection (%)
        %------------------------------------------------------------------
        Y{n}(i,8) = 100 * (p{n}{2}(3) + p{n}{2}(8));
        
        % number of people at home, asymptomatic, untested but infected
        %------------------------------------------------------------------
        Y{n}(i,9) = N(n) * x{n}(1,2,1,1);
        
        % positivity based on prevalence of infection (%)
        %------------------------------------------------------------------
        q          = p{n}{2}(1)*Q{n}.fpr(1) + ...
                     p{n}{2}(2)*(1 - Q{n}.fnr(1)) + ...
                     p{n}{2}(3)*(1 - Q{n}.fnr(2)) + ...
                     p{n}{2}(4)*Q{n}.fpr(2) + ...
                     p{n}{2}(5)*Q{n}.fpr(1) + ...
                     p{n}{2}(6)*Q{n}.fpr(1) + ...
                     p{n}{2}(7)*(1 - Q{n}.fnr(1)) + ...
                     p{n}{2}(8)*(1 - Q{n}.fnr(2));
                 
        Y{n}(i,11) = 100 * q;
        
        % number of symptomatic people
        %------------------------------------------------------------------
        Y{n}(i,12) = N(n) * p{n}{3}(2);

        % mobility (% normal)
        %------------------------------------------------------------------
        q = Q{n}.out/W{n}(1).Pout;
        if isfield(Q{n},'mo')
            q = q^Q{n}.mo;
        end
        Y{n}(i,13) = 100 * q;
        
        % work (% normal)
        %------------------------------------------------------------------
        q = Q{n}.out/W{n}(1).Pout;
        if isfield(Q{n},'wo')
            q = q^Q{n}.wo;
        end
        Y{n}(i,14) = 100 * q;
        
        % gross domestic product (% 2018)
        %------------------------------------------------------------------
        q = q^Q{n}.mo;
        if isfield(Q{n},'gd')
            q = q^Q{n}.gd;
        end
        Y{n}(i,32) = 100 * q;
        
        % certified deaths per day
        %------------------------------------------------------------------
        Y{n}(i,15) = N(n) * p{n}{3}(4);
        
        % hospital occupancy (ARDS people in hospital/CCU)
        %------------------------------------------------------------------
        q          = squeeze(spm_sum(x{n},[2,4]));
        q          = sum(q([3,6],3));
        Y{n}(i,27) = N(n) * q;
        
        % hospital admissions (ARDS admitted to hospital/CCU)
        %------------------------------------------------------------------
        Y{n}(i,16) = Y{n}(i,27)/Q{n}.Trd;    
                
        % excess deaths in hospital/CCU
        %------------------------------------------------------------------
        q          = squeeze(spm_sum(x{n},[2,4]));
        Y{n}(i,17) = N(n) * sum(q([3,6],4));
        
        % excess deaths not in hospital
        %------------------------------------------------------------------
        Y{n}(i,18) = N(n) * sum(q([1,2,4,5],4));
        
        % cumulative number of people (first dose) vaccinated (%)
        %------------------------------------------------------------------
        if isfield(Q{n},'pro')
            pvac(n) = pvac(n) + (1 - erf(Q{n}.pro) - pvac(n))*V.Rvac;
        else
            pvac(n) = pvac(n) + (1 - pvac(n))*V.Rvac;
        end
        if pvac(n) < 0, keyboard, end
        Y{n}(i,22) = 100 * pvac(n);
        
        % PCR case positivity (%) (seven day rolling average)
        %------------------------------------------------------------------
        Y{n}(i,23) = 100 * (1 - (1 - p{n}{4}(3))^7)/(1 - (1 - p{n}{4}(3) - p{n}{4}(4))^7);
        
        % daily lateral flow tests (positive and negative)
        %------------------------------------------------------------------
        Y{n}(i,24) = N(n) * (p{n}{4}(5) + p{n}{4}(6));
        
        % incidence of long COVID (number): 1/400 symptomatic and 1/2 ARDS
        %------------------------------------------------------------------
        Y{n}(i,28) = N(n) * p{n}{3}(2) * (1 - exp(-1/Q{n}.Trd)) * 1/400;
        Y{n}(i,28) = N(n) * p{n}{3}(3) * (1 - exp(-1/Q{n}.Trd)) * 1/2 + Y{n}(i,28);
        
        
        % joint density
        %------------------------------------------------------------------
        Z{n,i} = x{n};
        
    end

end

for n = 1:nN
    
    % effective reproduction ratio: exp(K*Q.Tcn): K = dln(N)/dt
    %----------------------------------------------------------------------
    K          = gradient(log(Y{n}(:,4)));
    Y{n}(:,4)  = 1 + K*(Q{n}.Tin + Q{n}.Tcn) + (K.^2)*Q{n}.Tin*Q{n}.Tcn;
    
    % incidence of new infections (%)
    %----------------------------------------------------------------------
    Kinf       = 1 - exp(-1/Q{n}.Tin);
    prev       = X{n,2}(:,2) + X{n,2}(:,7);
    Y{n}(:,10) = 100 * (gradient(prev) + Kinf*prev);
    
    % Y(:,19) - Daily incidence (per hundred thousand)
    %----------------------------------------------------------------------
    Y{n}(:,19) = 1000 * Y{n}(:,10);
    
    % Y(:,20) - Weekly confirmed cases (per hundred thousand)
    %----------------------------------------------------------------------
    Y{n}(:,20) = 100000 * 7 * Y{n}(:,2) / N(n);
    
    % infection fatality ratio (%)
    %----------------------------------------------------------------------
    Psev       = [W{n}.Psev]';
    Pfat       = [W{n}.Pfat]';
    Pinf       = [W{n}.Pinf]';
    Y{n}(:,21) = 100 * (1 - Pinf.^(Q{n}.Tin + Q{n}.Tcn)).*Psev.*Pfat;

    % cumulative attack rate (%)
    %----------------------------------------------------------------------
    Y{n}(:,25) = cumsum(Y{n}(:,10));
    
    % population immunity (seropositive, seronegative) (%)
    %----------------------------------------------------------------------
    Y{n}(:,26) = 100 * sum(X{n,2}(:,4:8),2);
    
    % population immunity (vaccine) (%)
    %----------------------------------------------------------------------
    Y{n}(:,29) = 100 * sum(X{n,2}(:,6:8),2);
        
    % accommodate (3 week) delay in immunity following vaccination
    %----------------------------------------------------------------------
    Y{n}(:,22) = Y{n}(:,22) + 21*gradient(Y{n}(:,22));
    
    % cumulative admissions
    %----------------------------------------------------------------------
    Y{n}(:,30) = cumsum(Y{n}(:,16));
    
    % retain specified output variables
    %----------------------------------------------------------------------
    Y{n} = Y{n}(:,U);
    
end

% deal with mixture models
%==========================================================================

% proportion of each population
%--------------------------------------------------------------------------
n     = N/sum(N);
Nt    = M.T;
Nu    = numel(U);

if numel(age) == 1 && age > 1
    
    % return requested population
    %----------------------------------------------------------------------
    x = X(age,:);
    
else
    
    % marginalise over populations
    %----------------------------------------------------------------------
    x      = cell(1,Nf); [x{:}] = deal(0);
    z      = cell(1,Nt); [z{:}] = deal(0);
    for i = 1:nN
        for j = 1:Nf
            x{1,j} = x{1,j} + X{i,j}*n(i);
        end
    end
    for i = 1:nN
        for j = 1:Nt
            z{1,j} = z{1,j} + Z{i,j}*n(i);
        end
    end
    
end

% age-specific outcomes or averages
%--------------------------------------------------------------------------
y     = zeros(Nt,Nu);
for u = 1:Nu
    
    if age(u) > 0
        
        % age-specific outcomes
        %------------------------------------------------------------------
        y(:,u) = Y{age(u)}(:,u);

        
    else % pool over groups
        
        % absolute numbers (add)
        %------------------------------------------------------------------
        if ismember(U(u),uN)
            for i = 1:nN
                y(:,u) = y(:,u) + Y{i}(:,u);
            end
            
        else
            
            % or proportions (average)
            %--------------------------------------------------
            for i = 1:nN
                y(:,u) = y(:,u) + Y{i}(:,u)*n(i);
            end
        end
    end
end

% reporting lags (first-order approximation)
%--------------------------------------------------------------------------
try
    lag   = log(P.lag(U));
    for i = 1:Nu
        y(:,i) = y(:,i) - lag(i)*gradient(y(:,i));
    end
end

% vectorise if data are asynchronous
%--------------------------------------------------------------------------
if exist('D','var')
    for i = 1:numel(D)
        j      = ismember(d,D(i).date);
        D(i).Y = y(j,i);
    end
    y  = spm_vec(D.Y);
end

return


