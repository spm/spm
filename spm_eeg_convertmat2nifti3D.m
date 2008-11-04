function S = spm_eeg_convertmat2nifti3D(S)

% Convert epoched EEG/ERP data from SPM- to analyze format by projecting
% onto the scalp surface
% FORMAT spm_eeg_convertmat2nifti3D(S)
%
% S         - optinal input struct
% (optional) fields of S:
% Fname            - matrix of EEG mat-files
% n                - size of quadratic output image (size: n x n x ?)
% pixsize          - Voxel size of output image
% interpolate_bad  - flag (0/1) whether channels should be used for 
%                    interpolation if they lie at the border of the setup
%                    (0: mask out). 
% output: 
% S         - can be used to construct script (as in the history-function)
%_______________________________________________________________________
%
% spm_eeg_convertmat2nifti3D converts EEG/MEG data from the SPM format to the
% scalp format. The data will be in 3D format, i.e., peri-stimulus time is the third dimension.
% The channel data is interpolated to voxel-space using a
% linear interpolation. Each channel's data will be found in an individual
% voxel given that n is big enough. The data is written to 3-dim nifti
% images, i.e. the data of each single trial or evoked response is contained in one
% image file. The mask out option interpolate_bad=0 will only have an
% effect if the bad channels are located at the edge of the setup, where it
% is assumed that there is not enough data for interpolation.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_convertmat2nifti3D.m 2438 2008-11-04 11:21:19Z stefan $

[Finter, Fgraph, CmdLine] = spm('FnUIsetup', 'EEG conversion setup',0);


try
    Fname = S.Fname;
catch
    Fname = spm_select(inf, '\.mat$', 'Select EEG mat file(s)');
end

Nsub = size(Fname, 1);

try
    n = S.n;
catch
    n = spm_input('Output image dimension', '+1', 'n', '32', 1);
end

try
    pixsize = S.pixsize;
catch
    pixsize = spm_input('Pixel dimensions (approx)', '+1', 'n', '3', 1);
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

clear D
for i = 1:Nsub
    D{i} = spm_eeg_load(deblank(Fname(i,:)));
end

for k = 1:Nsub

    [Cel, Cind, x, y] = spm_eeg_locate_channels(D{k}, n, interpolate_bad);

    % nr of (good) channels
    Nel = length(Cel);
    
    % generate data directory into which converted data goes
    [P, F] = fileparts(spm_str_manip(Fname(k, :), 'r'));
    if isempty(P), P = '.'; end
    [m, sta] = mkdir(P, spm_str_manip(Fname(k, :), 'tr'));
    if m ~= 1
        error(sta);
        return
    end

    cd(fullfile(P, spm_str_manip(Fname(k, :), 'tr')));
    
    d = (D{k}(Cind, :,:));
    cl = D{k}.condlist;

    for i = 1 : D{k}.nconditions
        
        Itrials = pickconditions(D{k}, cl(i), 1)';

        dname = sprintf('type_%s', cl{i});
        [m, sta] = mkdir(dname);
        cd(dname);
        
        for l = Itrials(:)
            fname = sprintf('trial%d.img', l);
                               
            dat = file_array(fname,[n n D{k}.nsamples 1],'FLOAT32-LE');
            N = nifti;
            N.dat = dat;
            N.mat = [   pixsize 0 0          -n*pixsize/2;...
                0 pixsize 0      -n*pixsize/2;...
                0 0 1000/D{k}.fsample  time(D{k}, 1, 'ms');...
                0 0 0        1];
            N.mat_intent = 'Aligned';
            create(N);
                        
            for j = 1 : D{k}.nsamples % time bins
                di = ones(n, n)*NaN;                
                di(sub2ind([n n], x, y)) = griddata(Cel(:,1), Cel(:,2), double(d(:, j, l)),x,y, 'linear');
                N.dat(:,:,j,1) = di;
            end        
            
            disp(sprintf('File %d, type %d, trial %d', k, i, l))
            
        end
        cd ..
    end
    cd ..
end

spm('Pointer', 'Arrow');
