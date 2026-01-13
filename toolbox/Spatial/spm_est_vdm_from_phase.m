function vdm_pha = spm_est_vdm_from_phase(vol1, vol1_mean, vol1_pha, total_readout_time, TE, PE_dir, outdir)

%==========================================================================
% Estimate Voxel Displacement Map (VDM) from phase images using ROMEO phase
% unwrapping tool.
% FORMAT:
% vdm_pha = spm_est_vdm_from_phase(vol1, vol1_mean, vol1_pha,
%                                         total_readout_time, TE, PE_dir, outdir)
%
% Input:
%   vol1                - cell array of file names of magnitude images with
%                         the same PE direction as the fMRI data to be corrected
%   vol1_mean          - file name of mean magnitude image with the same PE
%                         direction as the fMRI data to be corrected
%   vol1_pha           - cell array of file names of phase images with the
%                         same PE direction as the fMRI data to be corrected
%   total_readout_time - total readout time in seconds
%   TE                 - echo time in seconds
%   PE_dir             - phase-encoding direction (+1 or -1)
%   outdir             - output directory for temporary files
%
% Output:
%   vdm_pha            - Voxel Displacement Map (in mm) estimated from phase
%                         images
%
% Barbara Dymerska
% Copyright (C) 2025 Department of Imaging Neuroscience, UCL
%==========================================================================

mritools_dir = fullfile(spm('Dir'),'external', 'mritools');
romeo_bin = fullfile(mritools_dir,'bin', 'romeo') ;
if ~exist(romeo_bin, 'file')
    download_mritools(mritools_dir, '4.7.1')
end

if numel(vol1(:,1))~=numel(vol1_pha(:,1))
    error('You need the same nr of same PE direction phase volumes as magnitude volumes')
end

Nii_mag = zeros([size(nifti(vol1(1,:)).dat(:,:,:)) numel(vol1(:,1))]) ;
Nii_pha = zeros(size(Nii_mag)) ;

for f = 1:numel(vol1(:,1))
    Nii_mag_1vol = nifti(vol1(f,:)).dat(:,:,:);
    Nii_pha_1vol = nifti(vol1_pha(f,:)).dat(:,:,:);

    Nii_mag(:,:,:,f)= Nii_mag_1vol;
    Nii_pha(:,:,:,f)= Nii_pha_1vol;
end


% write 4D magnitude image for ROMEO phase unwrapping
%--------------------------------------------------------------------------
basename = spm_file(vol1(1,:),'basename');
oname    = char(spm_file(basename,'prefix','4D_mag_','ext','.nii'));
oname    = fullfile(outdir,oname);
Nio      = nifti;
Nio.dat  = file_array(oname,[size(Nii_mag)],'float32');
Nio.mat  = nifti(vol1(1,:)).mat;
create(Nio);
Nio.dat(:,:,:,:) = Nii_mag;

mag_file = oname;

% write 4D phase image for ROMEO phase unwrapping
%--------------------------------------------------------------------------
basename = spm_file(vol1_pha(1,:),'basename');
oname    = char(spm_file(basename,'prefix','4D_pha_','ext','.nii'));
oname    = fullfile(outdir,oname);
Nio      = nifti;
Nio.dat  = file_array(oname,[size(Nii_pha)],'float32');
Nio.mat  = nifti(vol1_pha(1,:)).mat;
create(Nio);
Nio.dat(:,:,:,:) = Nii_pha;

pha_file = oname;

% fix for unix library conflict for ROMEO phase unwrapping
if isunix; paths = getenv('LD_LIBRARY_PATH'); setenv('LD_LIBRARY_PATH'); end
% calling ROMEO phase unwrapping
system(sprintf('%s -p %s -m %s -o %s -t epi -v -i -g -k nomask', romeo_bin,pha_file,mag_file, outdir));
% returning to default MATLAB library environment:
if isunix; setenv('LD_LIBRARY_PATH', paths); end

pha_unwr = nifti(fullfile(outdir, 'unwrapped.nii')).dat(:,:,:,:) ;
pha_unwr = mean(pha_unwr,4) ;

vdm_pha = pha_unwr*total_readout_time/(2*pi*TE)*PE_dir ;

[~, mask] = spm_mask_from_segments(vol1_mean);

% fast VDM smoothing and extrapolation
vdm_pha = spm_smooth_extrap(vdm_pha, mask) ;


end

%% ===========================================================================
function download_mritools(outdir, version_tag)
% download_mritools  Download the correct OS-specific mritools binary
% from the GitHub release of CompileMRI.jl.
%
% Usage:
%       download_mritools('/my/download/path', '1.2.3')
% Inputs:
%   outdir      - target directory where to download and extract the files
%   version_tag - version tag of CompileMRI.jl release (e.g. '1.2.3')
%
% Note: this function requires an internet connection.
%==========================================================================

if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ============================================
%  Detect OS and derive release-year identifier
% ============================================
if ispc
    % Windows Server year detection — approximate using build number
    v = System.Environment.OSVersion.Version.Build ;
    WinYr = 2019 + (v >= 19041) * 3 + (v >= 22000) * 3;
    WinYr = min(WinYr, 2025);   % cap at latest supported
    assetKey = sprintf('mritools_windows-%d', WinYr);

elseif ismac
    % macOS version: e.g. 14.4 → major=14
    [~, verStr] = system('sw_vers -productVersion');
    parts = regexp(verStr, '(\d+)\.(\d+)', 'tokens', 'once');
    major = str2double(parts{1});
    major = max(major, 13) ; %assets only for 13+
    assetKey = sprintf('mritools_macOS-%d', major);

elseif isunix
    [~, distro] = system('lsb_release -d');
    tok = regexp(distro, 'Ubuntu\s+(\d+)', 'tokens', 'once');

    if isempty(tok)
        warning('Non-Ubuntu Linux detected; only Ubuntu assets exist hence these taken.');
        uYr = 22 ;
    else
        uYr = str2double(tok{1});
        uYr = 24 - 2 * (uYr < 24);   % map <24 → 22, else 24
    end
    assetKey = sprintf('mritools_ubuntu-%d.04', uYr);

else
    error('Unknown or unsupported operating system.');
end

fprintf('mritools (ROMEO) asset key: %s\n', assetKey);

%% ============================================
% Retrieve release page HTML
% ============================================
if ispc
    extension = 'zip' ;
else
    extension = 'tar.gz' ;
end

URL = sprintf('https://github.com/korbinian90/CompileMRI.jl/releases/download/v%s/%s_%s.%s',version_tag, assetKey, version_tag, extension);
fprintf('ROMEO binary for phase unwrapping not found. Downloading from %s \n',URL)
file_compressed = websave(sprintf('%s.%s',outdir, extension), URL);
if ispc
    unzip(file_compressed,outdir);
else
    file_tar = char(gunzip(file_compressed));
    untar(file_tar,outdir);
    movefile(fullfile(outdir,sprintf('%s_%s',assetKey, version_tag),'*'), outdir)
    rmdir(fullfile(outdir,sprintf('%s_%s',assetKey, version_tag)))
    delete(file_tar)
end
delete(file_compressed)

fprintf('mritools download complete.\n');
end
