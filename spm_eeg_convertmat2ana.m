function spm_eeg_convertmat2ana(S)
% Convert epoched EEG/ERP data from SPM- to analyze format by projecting
% onto the scalp surface
% FORMAT spm_eeg_convertmat2ana(S)
%
% S		    - optinal input struct
% (optional) fields of S:
% Fname		- matrix of EEG mat-files
% n         - size of quadratic output image (size: n x n x 1)
%_______________________________________________________________________
%
% spm_eeg_convertmat2ana converts EEG/MEG data from the SPM format to the
% scalp format. The channel data is interpolated to voxel-space using a
% spline interpolation. The electrodes' locations are specified by the
% channel template file. Each channel's data will be found in an individual
% voxel given that n is big enough. The data is written to 4-dim analyze
% images, i.e. the data of each single trial or ERP is contained in one
% image file.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id$

[Finter, Fgraph, CmdLine] = spm('FnUIsetup', 'EEG conversion setup',0);

% select matfiles to convert
try
    Fname = S.Fname;
catch
    Fname = spm_get(inf, '.mat', 'Select EEG mat file');
end

Nsub = size(Fname, 1);

try
    n = S.n;
catch
    n = spm_input('Output image dimension', '+1', 'n', '32', 1);
end

if length(n) > 1
    error('n must be scalar');
end

try
    interpolate_bad = S.interpolate_bad;
catch
    interpolate_bad = spm_input('Interpolate bad channels or mask out?',...
        '+1', 'b', 'Interpolate|Mask out', [1,0]);
end

spm('Pointer', 'Watch'); drawnow

% Load data set into structures
clear D
for i = 1:Nsub
    D{i} = spm_eeg_ldata(deblank(Fname(i,:)));
end

for k = 1:Nsub

    [Cel, Cind, x, y] = spm_eeg_locate_channels(D{k}, n, interpolate_bad);

    % nr of (good) channels
    Nel = length(Cel);

    % structure for spm_vol
    Vi.dim = [n n 1];
    Vi.mat = eye(4);
    Vi.pinfo = zeros(3,1);
    Vi.pinfo(1,1) = 1;
    Vi.dt(1) = 16; % float in old and new version...
    Vi.dt(2) = spm_platform('bigend');
    
    % generate data directory into which converted data goes
    [P, F] = fileparts(spm_str_manip(Fname(k, :), 'r'));
    [m, sta] = mkdir(P, spm_str_manip(Fname(k, :), 'tr'));
    cd(fullfile(P, F));
    
    d = squeeze(D{k}.data(Cind, :,:));
    
    for i = 1 : D{k}.events.Ntypes % trial types
        
        Itrials = find(D{k}.events.code == D{k}.events.types(i) & ~D{k}.events.reject);
        
        dname = sprintf('trialtype%d', D{k}.events.types(i));
        [m, sta] = mkdir(dname);
        cd(dname);
        
        for l = Itrials
            % if single trial data make new directory for single trials,
            % otherwise just write images to trialtype directory
            if D{k}.Nevents ~= D{k}.events.Ntypes
                % single trial data
                dname = sprintf('trial%d.img', l);
                fname = dname;
                [m, sta] = mkdir(dname);
                cd(dname);
            else
                fname = 'average.img';
            end
            
            Vi.fname = fname;
                        
            % remove file, if there is one
            spm_unlink([fname, '.hdr'], [fname, '.img']);
            
            for j = 1 : D{k}.Nsamples % time bins
                di = zeros(n, n);                
                di(sub2ind([n n], x, y)) = griddata(Cel(:,1), Cel(:,2), d(:, j, l),x,y, 'linear');
                
                Vi.n = j;
                Vo = spm_write_vol(Vi, di);
                
            end        
            
            if D{k}.Nevents ~= D{k}.events.Ntypes
                % single trial data
                disp(sprintf('Subject %d, type %d, trial %d', k, i, l))
                cd ..
            else
                % averaged data
                disp(sprintf('Subject %d, type %d', k, i))
            end
            
        end
        cd ..
    end
    cd ..
end

spm('Pointer', 'Arrow');
