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
% $Id: spm_eeg_convertmat2ana.m 776 2007-03-27 09:40:10Z stefan $

% [Finter, Fgraph, CmdLine] = spm('FnUIsetup', 'EEG conversion setup',0);
% 
% select matfiles to convert

try
    Fname = S.Fname;
catch
    Fname = spm_select(inf, '\.mat$', 'Select EEG mat file');
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
    
    % generate data directory into which converted data goes
    [P, F] = fileparts(spm_str_manip(Fname(k, :), 'r'));
     [m, sta] = mkdir(P, spm_str_manip(Fname(k, :), 'tr'));
    cd(fullfile(P, spm_str_manip(Fname(k, :), 'tr')));
    
    d = (D{k}.data(Cind, :,:));
    
    for i = 1 : D{k}.events.Ntypes % trial types
        
        Itrials = find(D{k}.events.code == D{k}.events.types(i) & ~D{k}.events.reject);
        
        dname = sprintf('trialtype%d', D{k}.events.types(i));
        [m, sta] = mkdir(dname);
        cd(dname);
        
        for l = Itrials
            if D{k}.Nevents ~= D{k}.events.Ntypes
                % single trial data
                fname = sprintf('trial%d.img', l);
            else
                fname = 'average.img';
            end
            
            dat = file_array(fname,[n n 1 D{k}.Nsamples],'FLOAT32');
            N = nifti;
            N.dat = dat;
            N.mat = eye(4);
            N.mat_intent = 'Aligned';
            create(N);
                        
            for j = 1 : D{k}.Nsamples % time bins
                di = ones(n, n)*NaN;                
                di(sub2ind([n n], x, y)) = griddata(Cel(:,1), Cel(:,2), d(:, j, l),x,y, 'linear');
                % griddata returns NaN for voxels outside convex hull (can
                % happen due to bad electrodes at borders of setup.)
                % Replace these by nearest neighbour interpoltation.
                tmp = find(isnan(di(sub2ind([n n], x, y))));
                di(sub2ind([n n], x(tmp), y(tmp))) =...
                    griddata(Cel(:,1), Cel(:,2), d(:, j, l),x(tmp),y(tmp), 'nearest');
                
                N.dat(:,:,1,j) = di;
            end        
            
            if D{k}.Nevents ~= D{k}.events.Ntypes
                % single trial data
                disp(sprintf('Subject %d, type %d, trial %d', k, i, l))
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
