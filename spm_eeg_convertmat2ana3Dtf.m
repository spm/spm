function spm_eeg_convertmat2ana3Dtf(S)

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


% Modification of spm_eeg_convertmat2ana3D to convert a 4D TF mat file with one channel to a 2D image (in 4D nifti) 			Rik Henson

try
    Fname = S.Fname;
catch
    Fname = spm_select(inf, '\.mat$', 'Select EEG mat file');
end

Nsub = size(Fname, 1);


spm('Pointer', 'Watch'); drawnow

% Load data set into structures
clear D
for i = 1:Nsub
    D{i} = spm_eeg_ldata(deblank(Fname(i,:)));

    if ~isfield(D{i},'Nfrequencies') | D{i}.Nchannels>1
	error('Only works for TF files with one channel!')
    end
end

for k = 1:Nsub
    
    % generate data directory into which converted data goes
    [P, F] = fileparts(spm_str_manip(Fname(k, :), 'r'));
    [m, sta] = mkdir(P, spm_str_manip(Fname(k, :), 'tr'));
    cd(fullfile(P, spm_str_manip(Fname(k, :), 'tr')));
    
    d = shiftdim(D{k}.data(1,:,:,:));
    
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
            
            dat = file_array(fname,[D{k}.Nsamples D{k}.Nfrequencies 1 1],'FLOAT32');
            N = nifti;
            N.dat = dat;
%            N.mat = eye(4);
%            N.mat = [  1 0 0		    min(D{k}.tf.frequencies);
%			0 1000/D{k}.Radc 0  -D{k}.events.start*1000/D{k}.Radc;
%			0 0 1		    0;
%			0 0 0		    1];
            N.mat = [	1000/D{k}.Radc 0 0  -D{k}.events.start*1000/D{k}.Radc;
			0 1 0		    min(D{k}.tf.frequencies);
			0 0 1		    0;
			0 0 0		    1];
            N.mat_intent = 'Aligned';
            create(N);
              
%	    N.dat = shiftdim(shiftdim(d(:,:,l),-1));          
            for j = 1 : D{k}.Nsamples % time bins
                N.dat(j,:,1,1) = d(:,j,l)';
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
