function spm_api_bmc
% API to slect and compare DCMs using Bayesian model comparison
% FORMAT spm_api_bmc
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_api_bmc.m 1131 2008-02-06 11:17:09Z spm $

% get DCM.mat files
%--------------------------------------------------------------------------
[f,p]   = uigetfile('*.mat','please select DCM files',...
                     'MultiSelect','on');
if ~iscell(f)
    msgbox('please select more than one DCM')
    return
end
                 
% get Free energy approximation to log-evidence
%--------------------------------------------------------------------------
F     = [];
N     = {};
for i = 1:length(f)
    
    name = fullfile(p,f{i});
    DCM  = load(name,'-mat');
    DCM  = DCM.DCM;
    try 
        name = DCM.name;
    end
    try
        [pth, name] = fileparts(name);
    end
    try
        F(end + 1) = DCM.F;
        N{end + 1} = name;
    end
end

% compute conditional probability of DCMs under flat priors.
%--------------------------------------------------------------------------
i    = F < (max(F) - 32);
P    = F;
P(i) = max(F) - 32;
P    = P - min(P)
P    = exp(P);
P    = P/sum(P)

   
% display results
%--------------------------------------------------------------------------   
Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf

subplot(2,1,1)
bar(F)
set(gca,'XTickLabel',N)
ylabel('log-evidence (relative)')
title('Bayesian model comparison')
axis square
grid on

subplot(2,1,2)
bar(P)
set(gca,'XTickLabel',N)
ylabel('probaility')
axis square
grid on
