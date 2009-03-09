function S = spm_eeg_convert2scalp(S)
% Convert epoched M/EEG data from SPM to NIfTI format by projecting
% onto the scalp surface
% FORMAT S = spm_eeg_convert2scalp(S)
%
% S         - input structure (optional)
% (optional) fields of S:
%   S.Fname            - character array of M/EEG mat-filenames
%   S.n                - size of square output image (size: n x n x ?)
%   S.interpolate_bad  - flag (0/1) whether channels should be used for 
%                        interpolation if they lie at the border of the
%                        setup [0: mask out]
% output: 
% S         - can be used to construct script (as in the history-function)
%__________________________________________________________________________
%
% spm_eeg_convert2scalp converts M/EEG data from the SPM format to the
% scalp format. The data will be in 3D format, i.e., peri-stimulus time is 
% the third dimension. The channel data is interpolated to voxel-space 
% using a linear interpolation. Each channel's data will be found in a 
% single voxel given that n is large enough. The data is written to 3D NIfTI
% images, i.e. the data of each single trial or evoked response is contained
% in one image file. The 'mask out' option interpolate_bad=0 will only have
% an effect if the bad channels are located at the edge of the setup, where
% it is assumed that there is not enough data for interpolation.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_convert2scalp.m 2843 2009-03-09 16:01:44Z guillaume $

SVNrev = '$Rev: 2843 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FnUIsetup', 'M/EEG conversion setup', 0);

%-Get parameters
%--------------------------------------------------------------------------
try
    Fname        = S.Fname;
catch
    [Fname, sts] = spm_select(Inf, 'mat', 'Select MEEG mat file(s)');
    if ~sts, S   = []; return; end
    S.Fname      = Fname;
end
Fname            = cellstr(Fname);

try
    n            = S.n;
catch
    n            = spm_input('Output image dimension', '+1', 'n', '32', 1);
    S.n          = n;
end
if length(n) > 1, error('n must be scalar'); end

try
    interpolate_bad = S.interpolate_bad;
catch
    interpolate_bad = spm_input('Interpolate bad channels or mask out?',...
        '+1', 'b', 'Interpolate|Mask out', [1,0]);
    S.interpolate_bad = interpolate_bad;
end

%-Load M/EEG data
%--------------------------------------------------------------------------
spm('Pointer', 'Watch');

D = cell(1,numel(Fname));
for i = 1:numel(Fname)
    D{i} = spm_eeg_load(Fname{i});    
end

%-For multimodal datasets, set the types of non-chosen modality to 'Other'
% This is not saved in the dataset
%--------------------------------------------------------------------------
modality = spm_eeg_modality_ui(D{1}, 1, 1);
if strcmp(modality, 'MEGPLANAR')
    error('MEG planar gradiometers are not supported yet.')
else
    for i = 1:numel(Fname)
        otherind = setdiff(1:nchannels(D{i}), strmatch(modality, chantype(D{i})));
        if ~isempty(otherind)
            D{i} = chantype(D{i}, otherind, 'Other');
        end
    end
end

%-Project M/EEG data on the scalp surface
%--------------------------------------------------------------------------
for k = 1:numel(Fname)

    [Cel, Cind, x, y] = spm_eeg_locate_channels(D{k}, n, interpolate_bad);
    
    %-Make output directory for each dataset
    %----------------------------------------------------------------------
    [P, F] = fileparts(Fname{k});
    if isempty(P), P = pwd; end
    [sts, msg] = mkdir(P, F);
    if ~sts, error(msg); end
    P  = fullfile(P, F);
    
    %-Loop over conditions
    %----------------------------------------------------------------------
    d  = spm_cond_units(D{k}(Cind, :,:));
    cl = D{k}.condlist;
    for i = 1 : D{k}.nconditions
       
        %-Make output directory for each condition
        %------------------------------------------------------------------
        dname = sprintf('type_%s', cl{i});
        [sts, msg] = mkdir(P, dname);
        if ~sts, error(msg); end
        Pi = fullfile(P, dname);
        
        %-Loop over trials
        %------------------------------------------------------------------
        Itrials = pickconditions(D{k}, cl(i), 1);
        for l = Itrials(:)'
            fname = fullfile(Pi, sprintf('trial%04d.img', l));
            
            N     = nifti;
            DIM   = [n n D{k}.nsamples];
            dat   = file_array(fname,[DIM 1],'FLOAT32-LE');
            N.dat = dat;
            V     = [136 172 100] ./ DIM;   % new voxel size
            C     = [68  100   0];          % new origin
            N.mat = [...
                V(1)  0     0                  -C(1);...
                0     V(2)  0                  -C(2);...
                0     0     1000/D{k}.fsample  time(D{k}, 1, 'ms');...
                0     0     0                  1];
            N.mat_intent = 'Aligned';
            create(N);
                        
            for j = 1 : D{k}.nsamples % time bins
                di = NaN(n,n);                
                di(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),double(d(:, j, l)),x,y,'linear');
                N.dat(:,:,j,1) = di;
            end        
            
            fprintf('Dataset %s, type %s, trial %d\n', F, cl{i}, l);    %-#
            
        end % for l = Itrials(:)'
        
    end % for i = 1 : D{k}.nconditions

end % for k = 1:numel(Fname)

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG conversion: done'); spm('Pointer','Arrow');
