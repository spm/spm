spm('defaults', 'eeg');

load(spm_select(1, 'mat', 'Select DCM mat file'));

D   = spm_eeg_load(DCM.xY.Dfile);

P   = DCM.Ep;
M   = DCM.M;
xU  = DCM.xU;
Nr  = size(DCM.C,1);                    % number of sources
Nt  = length(DCM.xY.y);                % number of trials
Ns  = M.ns;
X0  = DCM.xY.X0;
pst = DCM.xY.pst;


T0     = speye(Ns) - X0*((X0'*X0)\X0');

[dummy, s] = spm_cond_units(DCM.xY.y);

% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
x0  = ones(Ns,1)*spm_vec(M.x)';        % expansion point for states
L   = feval(M.G, DCM.Eg,M);            % get gain matrix
x   = feval(M.IS,P,M,xU);              % prediction (source space)
%%

% trial-specific responses (in mode, channel and source space)
%--------------------------------------------------------------------------
j     = find(kron(DCM.M.gE.J,ones(1,Nr)));
for i = 1:Nt
    x{i} = x{i} - x0;                   % centre on expansion point
    y{i} = T0*x{i}*L'*M.U';              % prediction (sensor space)
    r{i} = T0*DCM.xY.y{i}*M.U' - y{i};   % residuals  (sensor space)
    x{i} = x{i}(:,j);                   % Depolarization in sources
    Y{i} = M.U*y{i}';
end
%%

Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
for i = 1:Nt
    subplot(Nt, 1, i)
    plot(DCM.xY.pst, Y{i});
end
%%
ftdata       = [];
ftdata.trial = Y;
ftdata.time  = repmat({1e-3*DCM.xY.pst}, 1, Nt);
ftdata.label = DCM.xY.name;
%%
DD = spm_eeg_ft2spm(ftdata, 'simulated_data.mat');
%%
DD = conditions(DD, ':', DCM.xY.code);
DD = sensors(DD, 'EEG', sensors(D, 'EEG'));
DD = fiducials(DD, D.fiducials);
DD.inv = D.inv(1);

save(DD);
