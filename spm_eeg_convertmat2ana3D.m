function spm_eeg_convertmat2ana3D(S)

% Convert epoched EEG/ERP data from SPM- to analyze format by projecting
% onto the scalp surface
% FORMAT spm_eeg_convertmat2ana3D(S)
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
% $Id: spm_eeg_convertmat2ana.m 644 2006-10-10 14:08:32Z james $

% [Finter, Fgraph, CmdLine] = spm('FnUIsetup', 'EEG conversion setup',0);
% 
% select matfiles to convert


% Changed to add time as 3rd rather than 4th dimension		Rik Henson
% Image orientation matrix includes correct dimensions and origin for time
% Pixel dimensions added as option (so coordinates could be meaningful in future? (and images "nicer" to dsiplay!!!))


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

try
    trialtypes = S.trialtypes;
catch
    S.trialtypes = [];
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
    
    if isempty(S.trialtypes)
       trialtypes = [1:D{k}.events.Ntypes];
    end
    
    for i = trialtypes
        
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
            
            dat = file_array(fname,[n n D{k}.Nsamples 1],'FLOAT32');
            N = nifti;
            N.dat = dat;
%            N.mat = eye(4);
            N.mat = [	pixsize 0 0 		 -n*pixsize/2;
			0 pixsize 0		 -n*pixsize/2;
			0 0 1000/D{k}.Radc  -D{k}.events.start*1000/D{k}.Radc;
			0 0 0		 1];
            N.mat_intent = 'Aligned';
            create(N);
                        
            for j = 1 : D{k}.Nsamples % time bins
                di = ones(n, n)*NaN;                
                di(sub2ind([n n], x, y)) = griddata(Cel(:,1), Cel(:,2), d(:, j, l),x,y, 'linear');
                N.dat(:,:,j,1) = di;
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
