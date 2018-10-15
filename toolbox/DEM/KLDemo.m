function KLDemo
% Illustration of information gains with Bayesian fusion
% FORMAT KLDemo)
%
%--------------------------------------------------------------------------
% This routine  illustrates the benefit of multimodal or Bayesian fusion in
% terms of conditional dependencies among parameters. In other words, it
% shows that even if one data modality contains no information about a
% particular set of parameters, it can help resolve uncertainty about
% another set and thereby disclose information contained in the other
% modality. This is illustrated here using a simple linear model with
% neuronal and haemodynamic parameters to show that EEG can provide some
% information gain, in relation to haemodynamic parameters.
% 
% comment the orthogonalisation of the fMRI design matrix below to see the
% effect of conditional dependencies on the haemodynamic information gain
% afforded by EEG data
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Setup and preliminaries
%--------------------------------------------------------------------------
m    = 64;                         % number of observations
n    = 8;                          % number of parameters
j{1} = 1:n;                        % indices of both parameters
j{2} = 1:(n/2);                    % indices of neuronal parameters
j{3} = (1 + n/2):n;                % indices of haemodynamic parameters

% design matrices
%--------------------------------------------------------------------------
XEEG      = randn(m,n);
XMRI      = randn(m,n);

% make EEG design uninformative about haemodynamic parameters
%--------------------------------------------------------------------------
XEEG(:,j{3}) = 0;

% and orthogonalised fMRI design with respect to the EEG design
%--------------------------------------------------------------------------
% XMRI      = spm_orth(XMRI);

B         = randn(n,1);            % parameters
YEEG      = XEEG*B + randn(m,1)/8; % EEG data
YMRI      = XMRI*B + randn(m,1)/4; % MRI data

% model inversion using parametric empirical Bayes
%==========================================================================
PEEG{1}.X = XEEG;
PEEG{1}.C = {eye(m,m)};
PEEG{2}.X = zeros(n,1);
PEEG{2}.C = eye(n,n);

CEEG = spm_PEB(YEEG,PEEG,1);       % inversion using EEG data

PMRI{1}.X = XMRI;
PMRI{1}.C = {eye(m,m)};
PMRI{2}.X = zeros(n,1);
PMRI{2}.C = eye(n,n);

CMRI = spm_PEB(YMRI,PMRI,1);       % inversion using MRI data

PMRI{1}.X = XMRI;
PMRI{1}.C = {eye(m,m)};
PMRI{2}.X = CEEG{2}.E;             % Bayesian belief updating
PMRI{2}.C = CEEG{2}.C;             % using posteriors from EEG inversion

CMRE = spm_PEB(YMRI,PMRI,1);       % inversion using EEG and MRI data

% evaluate fMRI posteriors using Bayesian model reduction
%==========================================================================
[F,sE,sC] = spm_log_evidence(CMRE{2}.E,CMRE{2}.C,CEEG{2}.E,CEEG{2}.C,PEEG{2}.X,PEEG{2}.C);
CMRR.E    = sE;
CMRR.C    = sC;


%  evaluate information gain is entailed divergence
%==========================================================================
%  crucially, we will use the Bayesian model reduction estimate which means
%  we never have to actually invert the fMRI data
%--------------------------------------------------------------------------
for i = 1:length(j)          %  loop over different subsets of parameters
    
    % subsets
    %----------------------------------------------------------------------
    k      = j{i};
    
    % information gains
    %----------------------------------------------------------------------
    D(1,i) = spm_kl_normal(CEEG{2}.E(k),CEEG{2}.C(k,k),PEEG{2}.X(k),PEEG{2}.C(k,k));
    D(2,i) = spm_kl_normal(   CMRR.E(k),   CMRR.C(k,k),PEEG{2}.X(k),PEEG{2}.C(k,k));
    D(3,i) = spm_kl_normal(CMRE{2}.E(k),CMRE{2}.C(k,k),CEEG{2}.E(k),CEEG{2}.C(k,k));
    D(4,i) = spm_kl_normal(CMRE{2}.E(k),CMRE{2}.C(k,k),CMRI{2}.E(k),CMRI{2}.C(k,k));
    
    %  cow divergence between reduced and direct MRI posterior
    %----------------------------------------------------------------------
    D(5,i) = spm_kl_normal(CMRR.E(k),CMRR.C(k,k),CMRI{2}.E(k),CMRI{2}.C(k,k));
    
end

% illustrate results
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf;
str{1} =  '0 -> EEG';
str{2} =  '0 -> MRI';
str{3} =  'EEG -> MRI & EEG';
str{4} =  'MRI -> MRI & EEG';


%  show for information gains
%--------------------------------------------------------------------------
for i = 1:4
    subplot(3,2,i), bar(D(i,:))
    title(str{i},'fontsize',16)
    set(gca,'XTickLabel',{'both','neuronal','haemodynamic'})
    axis square
end
a = axis;

%  show  divergence between direct and reduced posterior
%--------------------------------------------------------------------------
subplot(3,1,3),bar(D(5,:))
title('KL between estimated and direct fMRI','fontsize',12)
set(gca,'XTickLabel',{'both','neuronal','haemodynamic'})
axis (a), axis square


