function res = bf_features_contcov(BF, S)
% Robust covariance for continuous data
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_features_contcov.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    contcov      = cfg_branch;
    contcov.tag  = 'contcov';
    contcov.name = 'Robust covariance';
    contcov.val  = {};
    contcov.help = {'Robust covariance for continuous data'};
    
    res = contcov;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

nchans  = length(S.channels);

samples = S.samples{1};

bad  = badsamples(D, S.channels, samples, 1);

% Combine fixed 1 sec segments with change points of bad channel number
% and prune to avoid too many short segments
fixedseg = 1:round(D.fsample):D.nsamples;

chngpnt  = [1 find(any(diff(bad, [], 2)))+1];

prunetable = abs(repmat(fixedseg, length(chngpnt), 1)-repmat(chngpnt', 1, length(fixedseg)));
fixedseg(any(prunetable<D.fsample)) = [];

chngpnt = sort([fixedseg, chngpnt]);

nsegments = length(chngpnt);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nsegments, 'Computing covariance'); drawnow;
if nsegments > 100, Ibar = floor(linspace(1, nsegments,100));
else Ibar = 1:nsegments; end

YY    = zeros(nchans);

N     = YY;

id = zeros(1, nsegments);

for i = 1:nsegments
    if i<nsegments
        Y = squeeze(D(S.channels, samples(chngpnt(i)):samples((chngpnt(i+1)-1))));
    else
        Y = squeeze(D(S.channels, chngpnt(i):samples(end)));
    end
    
    goodind = ~bad(:, chngpnt(i));
    
    Y = Y(goodind, :);
    
    
    if size(Y, 2)>1
        Y = detrend(Y', 'constant')';
    end
    
    YY(goodind, goodind) = YY(goodind, goodind)+(Y*Y');
    N(goodind, goodind)  = N(goodind, goodind)+size(Y, 2);
    
    id(i) = spm_data_id(S.channels(goodind));
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

C = YY./N;

features.C = C;
features.N = N;

res = features;