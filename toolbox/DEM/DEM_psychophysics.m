function DCM = DEM_psychophysics
% FORMAT DCM = DEM_psychophysics
%
% Demonstration of psychometric curve fitting and model comparison
%__________________________________________________________________________
%
% This demonstration routine illustrates the fitting of psychometric
% functions under different models or hypotheses. The models in question
% are specified by various constraints on the model parameters; namely,
% changes in bias and sensitivity. The generative model uses a binomial
% likelihood function and a logistic function for the predicted
% psychometric curves.
% 
% A binomial likelihood model means that (under the assumption of a large
% number of trials) following a (variance stabilising) square root
% transform, the error variance of the number of correct responses is
% unity. If we scale the number of correct responses to reflect the
% proportion of correct responses, then the precision of the ensuing
% (square root transform) data feature is the total number of responses.
% This provides a simple and efficient way to specify the generative model
% as a generalised linear model, based upon a standard logistic function,
% whose parameters correspond to bias and sensitivity.
% 
% In this example, data are loaded from a CSV file and converted into the
% proportion of correct responses, under two levels or conditions of an
% experimental factor (baseline and blindspot). The generative model is
% equipped with parameters corresponding to changes in bias and
% sensitivity. Crucially, these changes have parity; namely, increases and
% decreases. This means that there are, potentially, eight models:
% corresponding to the presence or absence of an increase or decrease in
% bias and sensitivity. In the example below, three of these models are
% compared, in terms of their marginal likelihood (as approximated by a
% softmax function of the ensuing variational free energy).
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


% set up and get data
%==========================================================================
spm_figure('GetWin','SI'); clf;


%% import data: hhere, in the form of a CSV file  or spreadsheet
%--------------------------------------------------------------------------
tab   = readtable('blindspot_data_format.csv');
vname = tab.Properties.VariableNames;            % variable names

% response and explanatory variables for each subject
%--------------------------------------------------------------------------
subj  = find(ismember(vname,'ID'));              % subject identifiers
fac1  = find(ismember(vname,'foil_targ_dif'));   % psychological variable
fac2  = find(ismember(vname,'targ_in_bs'));      % conditions

% levels of response and explanatory variables for each subject
%--------------------------------------------------------------------------
subs  = unique(tab(:,subj));
lev1  = unique(tab(:,fac1));
lev2  = unique(tab(:,fac2));

% response for each entry
%--------------------------------------------------------------------------
targ  = table2array(tab(:,find(ismember(vname,'targ_chosen'))));
foil  = table2array(tab(:,find(ismember(vname,'foil_chosen'))));

% creates subject specific cell arrays for:
%--------------------------------------------------------------------------
Y     = cell(numel(subs),1);                     % proportion correct
U     = cell(numel(subs),1);                     % experimental factors
Q     = cell(numel(subs),1);                     % precision
for s = 1:numel(subs)                            % for each subject
    for f2 = 1:numel(lev2)                       % and condition factor
        for f1 = 1:numel(lev1)                   % and level                

            % find entries, j, for this level, condition and subject
            %--------------------------------------------------------------
            j = ismember(tab(:,subj),subs(s,1));
            j = j & ismember(tab(:,fac2),lev2(f2,1));
            j = j & ismember(tab(:,fac1),lev1(f1,1));
            j = find(j);

            %  record number of responses (n) and proportion correct Y{i}
            %--------------------------------------------------------------
            n               = sum(targ(j)) + sum(foil(j));
            Y{s}(end + 1,1) = sum(targ(j))/n;
            U{s}(end + 1,:) = [table2array(lev1(f1,1)) f2];
            Q{s}(end + 1,1) = n;
        end
    end
end

% set hypothesis specific priors
%==========================================================================
% For each theory, specify whether or not there is an increase or decrease
% in bias and sensitivity
%__________________________________________________________________________
theory = {'IIT','NRep','AI'};                    % theory names
Nh     = numel(theory);                          % number of theories

% theory: IT NR AI
%__________________________________________________________________________
bias   = [0  0  0;          % increase
          1  0  0];         % decrease
sens   = [1  0  0;          % increase
          1  0  1];         % decrease

% specify hypothesis or theory-specific prior constraints
%==========================================================================

% prior expectations (pE)
%--------------------------------------------------------------------------
pE.b    = mean(U{1}(:,1));                       % expected bias
pE.s    = -log(std(U{1}(:,1)));                  % expected log-sensitivity
pE.dbi  = 0;                                     % log changes
pE.dbd  = 0;
pE.dsi  = 0;
pE.dsd  = 0;

% prior [co]variance (pC)
%--------------------------------------------------------------------------
scale = 8;
for h = 1:Nh
    pc.b   = var(U{1}(:,1));
    pc.s   = scale;

    pc.dbi = bias(1,h)*scale;
    pc.dbd = bias(2,h)*scale;
    pc.dsi = sens(1,h)*scale;
    pc.dsd = sens(2,h)*scale;

    % covariance constraints for this hypothesis
    %----------------------------------------------------------------------
    pC{h}  = pc;

end

% model inversion with variational Laplace
%==========================================================================

% for each subject
%--------------------------------------------------------------------------
Ns    = numel(Y);
for s = 1:Ns

    % for each hypothesis
    %----------------------------------------------------------------------
    for h = 1:Nh

        % data structure: response variables and precision matrix
        %------------------------------------------------------------------
        xY.y = Y{s};
        xY.Q = diag(Q{s});

        % explanatory structure: levels and experimental conditions
        %------------------------------------------------------------------
        xU   = U{s};

        % model specification for this hypothesis
        %==========================================================================
        M.Nmax = 32;                   % maximum number of iterations
        M.G    = @spm_pp_gen;          % generative function
        M.FS   = @(Y)real(sqrt(Y));    % feature selection  (square root transform)
        M.pE   = pE;                   % prior expectations (parameters)
        M.pC   = pC{h};                % prior covariances  (parameters)
        M.hE   = log(1);               % prior expectation  (log-precision)
        M.hC   = 1/128;                % prior covariances  (log-precision)


        % model inversion with Variational Laplace (Gauss Newton)
        %==========================================================================
        [Ep,Cp,Eh,F] = spm_nlsi_GN(M,xU,xY);

        % save in DCM structure
        %--------------------------------------------------------------------------
        DCM(s,h).M  = M;
        DCM(s,h).Ep = Ep;
        DCM(s,h).Eh = Eh;
        DCM(s,h).Cp = Cp;
        DCM(s,h).U  = xU;
        DCM(s,h).F  = F;
        DCM(s,h).Y  = Y{s};

    end
end

% plot results
%==========================================================================

% for each subject
%--------------------------------------------------------------------------
col     = get(gca,'ColorOrder');
for s = 1:Ns

    % for each hypothesis
    %----------------------------------------------------------------------
    for h = 1:Nh

        % evaluate derivatives of outcome w.r.t. parameters
        %------------------------------------------------------------------
        [dYdP,Y] = spm_diff(@(P,M,U)spm_pp_gen(P,M,U),...
                   DCM(s,h).Ep,DCM(s,h).M,DCM(s,h).U,1);


        % and evaluate conditional covariances in data space
        %------------------------------------------------------------------
        C      = dYdP*DCM(s,h).Cp*dYdP';  % posterior predictive covariance
        F(h,s) = DCM(s,h).F;              % negative free energy or ELBO

        % plot Bayesian credible intervals and data for each subject
        %------------------------------------------------------------------
        subplot(Nh + 1, Ns, s + (h - 1)*Ns)
        set(gca,'ColorOrderIndex',1)
        for i = 1:numel(lev2)

            j = find(DCM(s,h).U(:,2) == i);
            spm_plot_ci(Y(j)',C(j,j),DCM(s,h).U(j,1)), hold on
            plot(DCM(s,h).U(j,1),DCM(s,h).Y(j),'.','color',col(i,:))
            plot(DCM(s,h).U(j,1),DCM(s,h).Y(j),'o','color',col(i,:))

            xlabel(vname(fac1))
            title(sprintf('%s: subject %i',theory{h},s),'FontSize',14)
        end

    end

    %   marginal likelihood over hypotheses (softmmax of ELBO)
    %----------------------------------------------------------------------
    subplot(Nh + 1, Ns, s + Nh*Ns)
    bar(spm_softmax(F(:,s)))
    set(gca,'XTickLabel',theory)
    title('Marginal likelihood')

end

% legends
%--------------------------------------------------------------------------
cond = table2cell(lev2);
subplot(Nh + 1, Ns, 1)
legend(cond(kron(1:numel(cond),[1,1,1,1])),'Location','southeast')



function Y = spm_pp_gen(P,M,U)
% prediction of psychometric (logistic) function
% FORMAT Y = spm_pp_gen(P,M,U)
%
% P - parameter structure
% M - model structure
% U - explanatory variables
%
% Y - predicted response variable
%__________________________________________________________________________


% prediction
%--------------------------------------------------------------------------
Y     = zeros(size(U,1),1);
for i = 1:size(U,1)

    % condition 1 (baseline)
    %----------------------------------------------------------------------
    j    = find(U(:,2) == 1);                  % indices for this condition

    b    = P.b;                                % bias
    s    = exp(P.s);                           % sensitivity

    Y(j) = spm_phi((U(j,1) - b)*s);            % logistic function

    % condition 2 (blindspot): increase or decrease to bias and sensitivity
    %----------------------------------------------------------------------
    j    = find(U(:,2) == 2);
 
    if M.pC.dbi                                % if no constraint
        b = b + exp(P.dbi);                    % increase bias
    end
    if M.pC.dbd
        b = b - exp(P.dbd);
    end
    if M.pC.dsi
        s = s + s*exp(P.dsi);
    end
    if M.pC.dsd
        s = s - s*exp(P.dsd);
    end

    Y(j) = spm_phi((U(j,1) - b)*s);

end

return


