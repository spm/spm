function spm_api_bmc
% API to slect and compare DCMs using Bayesian model comparison
% FORMAT spm_api_bmc
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_api_bmc.m 2061 2008-09-09 18:04:42Z jean $

% get DCM.mat files
%--------------------------------------------------------------------------
f = cellstr(spm_select(Inf,'mat','please select DCM files'));
if length(f) < 2
    msgbox('please select more than one DCM')
    return
else
    for i=1:length(f)
        [p{i},f{i}] = fileparts(f{i});
    end
end

                 
% get Free energy approximation to log-evidence
%--------------------------------------------------------------------------
F     = [];
N     = {};
for i = 1:length(f)
    
    name = fullfile(p{i},f{i});
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
P    = P - min(P);
P    = exp(P);
P    = P/sum(P);

   
% display results
%--------------------------------------------------------------------------   
Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf

subplot(2,1,1)
bar(1:length(N),F)
set(gca,'XTick',1:length(N))
set(gca,'XTickLabel',N)
ylabel('log-evidence (relative)')
title('Bayesian model comparison')
axis square
grid on

subplot(2,1,2)
bar(1:length(N),P)
set(gca,'XTick',1:length(N))
set(gca,'XTickLabel',N)
ylabel('probability')
axis square
grid on
