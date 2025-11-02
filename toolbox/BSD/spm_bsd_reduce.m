function [BSDs, F, Rfq, dF] = spm_bsd_reduce(BSD)
%SPM_BSD_REDUCE Prunes frequency components from BSD model and evaluates model evidence.
%   [BSDs, F, Rfq, dF] = spm_bsd_reduce(BSD) iteratively removes the least
%   significant frequency component from the BSD structure, fits the reduced
%   model, and computes the change in model evidence (log Bayes factor).
%
%   Inputs:
%       BSD  - Structure containing BSD model parameters and frequency components.
%
%   Outputs:
%       BSDs - Cell array of BSD structures after each pruning step.
%       F    - Array of model evidence values for each model.
%       Rfq  - Cell array of removed frequency labels.
%       dF   - Array of log Bayes factors for each pruning step.

fqs = BSD.fqs; 
Nm = length(fqs);
BSDs = {}; 

% Fit full model
BSD = spm_bsd(BSD); 
F = BSD.F; 
dF = []; 
Rfq = {};

for m = 1:Nm

    % Prune least significant frequency component
    amp = BSD.Ep.a.^2;
    [~, ii] = sort(amp, 1, "ascend"); 
    kk = ii(2:end); 

    fprintf('Removing %s\n\n', fqs{ii(1)});
    BSD.name = spm_file(BSD.name, 'suffix', ['_' fqs{ii(1)}]); 
    Rfq = [Rfq, fqs{ii(1)}];

    fqs = fqs(kk); 

    BSD.fqs = fqs; 

    BSD.M.P = BSD.Ep; 

    BSD.M.P.f = BSD.M.P.f(kk);
    BSD.M.P.S = BSD.M.P.S(kk);
    BSD.M.P.a = BSD.M.P.a(kk);

    % Fit reduced model
    BSD = spm_bsd(BSD); 
    BSDs = [BSDs, BSD]; 

    F = [F, BSD.F]; 
    dF = [dF, F(end) - F(end-1)];

    fprintf('log BF: %.2e\n\n', dF(end));
end

% Plot model evidence and log Bayes factors
spm_figure('GetWin', 'BSD - Reduce'); 
subplot(2,1,1);
bar(F-F(end));
set(gca,'XTickLabel',[{'full'}, Rfq])

subplot(2,1,2);
bar(dF);
set(gca,'XTickLabel',Rfq)
end
