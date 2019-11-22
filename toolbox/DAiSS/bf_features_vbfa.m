function res = bf_features_vbfa(BF, S)
% Variational Bayes Factor Analysis for computing noise covariance
% Code contributed by Sri Nagarajan
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_features_vbfa.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    nl = cfg_entry;
    nl.tag = 'nl';
    nl.name = 'Factor dimensionality';
    nl.strtype = 'n';
    nl.num = [1 1];
    nl.val = {5};
    
    nem = cfg_entry;
    nem.tag = 'nem';
    nem.name = 'Number of EM iterations';
    nem.strtype = 'n';
    nem.num = [1 1];
    nem.val = {50};      
    
    vbfa      = cfg_branch;
    vbfa.tag  = 'vbfa';
    vbfa.name = 'VB Factor Analysis';
    vbfa.val  = {nl, nem};
    vbfa.help = {'This method uses code contributed by Sri Nagarajan',...
        'It should be used for computing noise covariance for Champagne'};
    
    res = vbfa;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

ntrials = length(S.trials);

nwoi            = numel(S.samples);
nsamples        = length(S.samples{1});

if nwoi ~= 2
    error('Baseline and activation windows must be specified');
end

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Reading data'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

NN    = [];
YY    = [];
for i = 1:ntrials
    for j = 1:nwoi
        Y  = squeeze(D(S.channels, S.samples{j}, S.trials(i)));
        
        Y = detrend(Y', 'constant')';
        
        if j == 1
            NN = [NN, Y];
        else
            YY = [YY Y];
        end
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

N = ntrials*nsamples*nwoi;

% In the champagne header it is advised to convert data to pT.
YY = 1e3*YY; 
NN = 1e3*NN;

Fgraph  = spm_figure('GetWin', 'VBFA'); figure(Fgraph); clf

C = vbfa_aug2015(NN,S.nl,S.nem,Fgraph);

features.C = C;
features.N = N;
features.Y = YY;

res = features;