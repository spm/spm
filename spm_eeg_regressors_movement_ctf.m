function res = spm_eeg_regressors_movement_ctf(S)
% Generate movement regressors for CTF MEG data
% S                     - input structure
% fields of S:
%    S.D                - M/EEG object
%
%    Additional parameters can be defined specific for each plugin
% Output:
%  res -
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided the plugin returns
%______________________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak

SVNrev = '$Rev: 6186 $';

if nargin == 0
    
    Dmov        = cfg_files;
    Dmov.tag    = 'Dmov';
    Dmov.name   = 'Movement dataset name';
    Dmov.filter = 'mat';
    Dmov.num    = [1 1];
    Dmov.help   = {'Select the M/EEG mat file containing continuous head localisation data.',...
        'This might or might not be the same as the dataset being analysed.'};
    
    
    movement_ctf = cfg_branch;
    movement_ctf.tag = 'movement_ctf';
    movement_ctf.name = 'CTF head movements';
    movement_ctf.val = {Dmov};
    
    res = movement_ctf;
    
    return
end

%read HLC-channels
%HLC0011 HLC0012 HLC0013 x, y, z coordinates of nasion-coil in m.
%HLC0021 HLC0022 HLC0023 x, y, z coordinates of lpa-coil in m.
%HLC0031 HLC0032 HLC0033 x, y, z coordinates of rpa-coil in m.
hlc_chan_label = {'HLC0011' 'HLC0012' 'HLC0013'...
    'HLC0021' 'HLC0022' 'HLC0023'...
    'HLC0031' 'HLC0032' 'HLC0033'};

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','CTF movement regressors');

if iscell(S.Dmov)
    S.Dmov = char(S.Dmov);
end

Dmov = spm_eeg_load(S.Dmov);
D    = spm_eeg_load(S.D);

hlc_chan_ind = Dmov.indchannel(hlc_chan_label);

if length(hlc_chan_ind) ~= 9
    error('Nine CTF HLC channels should be present to define movement regressors.');
end

if ~isequal(Dmov.transformtype, 'time')
    error('The movement dataset should be in the time domain.');
end

if isequal(D.type, 'continuous')
    
    if ~isequal(Dmov.type, 'continuous') || (D.time(1) < Dmov.time(1)) || (D.time(end)>Dmov.time(end))
        error('All times of the input dataset should be within the movement dataset.');
    end
    
    data = Dmov(hlc_chan_ind, :);
    
    if D.fsample ~= Dmov.fsample
        [data, alpha] = spm_timeseries_resample(data, D.fsample/Dmov.fsample);
    else
        alpha = 1;
    end
    
    start = round(alpha*Dmov.indsample(D.time(1)));
    
    data = data(:, start:(start+D.nsamples-1));
    
else
    if D.ntrials ~= Dmov.ntrials
        error('Trial numbers should be equal between input and movement dataset.');
    end
    
    data = Dmov(hlc_chan_ind, :, :);
    
    if S.summarise
        data = spm_squeeze(mean(data, 2), 2);      
    else
        data = reshape(data, size(data, 1), []);
    end

end
   
res.R     = hpi2mov(data);

res.names = {'x', 'y', 'z', 'pitch', 'roll', 'yaw'};


spm('FigName','CTF movement regressors: done');


function P = hpi2mov(data)

% John Ashburner

ref   = reshape(data(:,1),3,3)*1000;
muRef = mean(ref,2);
ref0  = ref - repmat(muRef, 1, 3);
P     = zeros(size(data,2),6);

ns   = size(data, 2);

spm_progress_bar('Init', ns, 'Computing movement parameters'); drawnow;
if ns > 100, Ibar = floor(linspace(1, ns,100));
else Ibar = 1:ns; end


for i=1:ns
    cur     = reshape(data(:,i),3,3)*1000;
    muCur   = mean(cur,2);
    cur0    = cur - repmat(muCur, 1, 3);
    [U,SS,V] = svd(cur0*ref0');
    Ri      = U*V';
    if det(Ri)<0, V(:,end)=-V(:,end); Ri=U*V'; end
    p       = spm_imatrix([Ri muCur-Ri*muRef; 0 0 0 1]);
    P(i,:)  = p(1:6);
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

P = detrend(P, 'constant');
P = P./repmat(std(P), size(P, 1), 1);

