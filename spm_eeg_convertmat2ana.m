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
% Stefan Kiebel & Christophe Phillips $Id$

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

spm('Pointer', 'Watch');

% Load data set into structures
clear D
for i = 1:Nsub
	D{i} = spm_eeg_ldata(deblank(Fname(i,:)));
end

% load channel template file (contains location of channels)
Ctf = load(fullfile(spm('dir'), 'EEGtemplates', D{1}.channels.ctf));

Cel = [];
Nchannels = D{1}.Nchannels;

% find channel positions on 2D plane
for j = 1:Nchannels
    % don't use eog channels
	if j ~= D{1}.channels.heog & j ~= D{1}.channels.veog
		Cel = [Cel; Ctf.Cpos(:, D{1}.channels.order(j))'];
	end
end
Cel = round(Cel*n);
Nel = size(Cel, 1);

% Include all electrodes that lie within this circle (where the
% radius to the image edge is -0.5n to 0.5n).
r = 0.47*n;
[x, y] = meshgrid(1:n, 1:n);
Ip = find(((x-n/2).^2 + (y-n/2).^2) <= r^2);
x = x(Ip); y = y(Ip);

% Generation of projection matrix
%------------------
% matrix E
E = ones(Nel, 6);
E(:,2:3) = Cel;
E(:,[4 6]) = Cel.^2;
E(:, 5) = prod(Cel')';

% matrix K
K1 = repmat(Cel(:, 1), 1, Nel);
K1 = (K1 - K1').^2 + eps;
K2 = repmat(Cel(:, 2)', Nel, 1);
K2 = (K2 - K2').^2 + eps;
K  = (K1 + K2).^2 .* log(K1 + K2);

% Main matrices of the interpolation
%	[K  E] . [P] = [V]
%	[E' 0]   [Q] = [0]
% to solve for P and Q
M  = [K E; E' zeros(6)];
IM = pinv(M);

Mx = (repmat(x, 1, Nel)' - repmat(Cel(:, 1), 1, length(Ip)));
My = (repmat(y, 1, Nel)' - repmat(Cel(:, 2), 1, length(Ip)));
T = Mx.^2 + My.^2 + eps;
T = T.^2 .* log(T);

% structure for spm_vol
if isfield(D{1}, 'tf')
    Vi.dim = [n n D{1}.Nsamples 16];
    Vi.mat = eye(4);
    Vi.pinfo = zeros(3,1);
    Vi.pinfo(1,1) = 1;
else
    Vi.dim = [n n 1 16];
    Vi.mat = eye(4);
    Vi.pinfo = zeros(3,1);
    Vi.pinfo(1,1) = 1;
end

% preparing some variables for faster exexcution of main loop
g0 = zeros(n, n)+NaN;
Tt = T';
XYZ = [ones(length(Ip), 1) x y x.^2 x.*y y.^2];

for k = 1:Nsub
    
    % generate data directory into which converted data goes
    [P, F] = fileparts(spm_str_manip(Fname(k, :), 'r'));
	[m, sta] = mkdir(P, spm_str_manip(Fname(k, :), 'tr'));
	cd(fullfile(P, F));
    
    Cind = [1:Nchannels];
    
    % throw out EOG channels
    tmp = [];
    if D{k}.channels.veog
        tmp = D{k}.channels.veog;
    end
    if D{k}.channels.heog
        tmp = [tmp D{k}.channels.heog];
    end
    
    Cind(tmp) = [];
    
    % works also for time-frequency data, but only on one frequency!!
    if isfield(D{k}, 'tf')
        d = squeeze(D{k}.data(Cind, 1, :, :));
    else
        d = squeeze(D{k}.data(Cind, :,:));
    end
    
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
                dname = sprintf('trial%d', l);
                [m, sta] = mkdir(dname);
                cd(dname);
            end
            
            fname = sprintf('trial%d', l);
            Vi.fname = fname;
            
            if isfield(D{k}, 'tf')
                dg = zeros(Vi.dim(1:3));
            end
            
            for j = 1 : D{k}.Nsamples % time bins
                
                
                % Define grid
                g = g0;
                
                V = d(:, j, l);
                
                b  = [V; zeros(6, 1)];
                PQ = IM * b;
                
                P = PQ(1:Nel);
                Q = PQ(Nel+1 : Nel+6);
                
                t1 = Tt * P;
                t2 = XYZ * Q;
                C = t1 + t2;
                
                g(Ip) = C;
                g = g';
                
                if isfield(D{k}, 'tf')
                    dg(:,:, j) = g;
                else
                    Vi.n = j;
                    Vo = spm_write_vol(Vi, g);
                end
                
            end        
            
            if isfield(D{k}, 'tf')
                Vo = spm_write_vol(Vi, dg);
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
end

spm('Pointer', 'Arrow');

% test to check data conversion
% --> I found that data is correctly interpolated
%P = spm_get(inf, '*.img', 'Select img-files');
%for i = 1:size(P,1)
%	Vt(i) = spm_vol(deblank(P(i,:)));
%end
%[x, y] = meshgrid(1:Vt(1).dim(1), 1:Vt(1).dim(2));
%z = ones(size(x));
%X = zeros(Vt(1).dim(1), Vt(1).dim(2), length(Vt));
%for i = 1:length(Vt)
%	X(:,:,i) = spm_sample_vol(Vt(i),x,y,z,0);
%end
