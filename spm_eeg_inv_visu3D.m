function spm_eeg_inv_visu3D(S,Cflags)

%=======================================================================
% Enables visualisation of 3D reconstruction
% in both mesh space and voxel space
%
% FORMAT H = spm_eeg_inv_visu3D(S,Cflags)
% Input:
% S		    - input data struct (optional)
% Cflags    - [Cabs Cnorm] vector of length 2
%               - first number Cabs indicates whether absolute values
%               should be considered (1) or not (0) for display (default = 1)
%               - second number Cnorm indicates whether values should be
%               normalized to their maximum (1) or not (0) (default = 1)
%
% Note that the 4D (3D + times) raw images are stored anyway in NIFTI
% format for further statistical analysis and display
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_visu3D.m 308 2005-11-23 19:21:56Z jeremie $

spm_defaults

def_Cflags = [1 1];

if nargin == 1
    Cflags = def_Cflags
end    

    
try
    D = S;
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D = spm_eeg_ldata(D);
end

val = length(D.inv);


% Load the mesh
if ~isempty(D.inv{val}.mesh.tess_ctx)
    load(D.inv{val}.mesh.tess_ctx);
else
    error(sprintf('No specified mesh\n'));
end


% Load the interpolation matrix
if ~isempty(D.inv{val}.mesh.InterpMat)
    load(D.inv{val}.mesh.InterpMat);
else
    error(sprintf('No available interpolation matrix\n'));
end


% Load results to be displayed
if ~isempty(D.inv{val}.inverse.resfile)
    activity = D.inv{val}.inverse.activity;
    load(D.inv{val}.inverse.resfile);
else
    error(sprintf('No results available for this analysis\n'));
end


switch activity
    case 'evoked'
        MRIcurrent      = H*J';
        [pth,nam,ext]   = fileparts(D.inv{val}.inverse.resfile);
        outputimg       = fullfile(pth,[nam '.img']);
        V               = saverecons(MRIcurrent,spm_vol(D.inv{val}.mesh.sMRI),outputimg);
        
        % Display
        if Cflags(1)
            J = abs(J);
        end
        if Cflags(2)
            J = J./max(abs(J(:)));
        end

        
        
    case 'induced'
        
        
        
    case 'evoked & induced'
        
        
        
end

Jas = zeros(N_src,1);
Jas(Ifr) = Ja(:,Ts);
MRIcurrent = H*Jas';
outputname = [SubjectName '_percep_faces.img'];
spm_eeg_saverecon(MRIcurrent,outputname);

Ind = find(MRIcurrent);
Z = MRIcurrent(Ind);

%%%%%%%%%%%%%%%%%%% Change by Rik: to handle images with dif dims
%M = zeros(256,256,256);
M = zeros(IMG.dim(1:3));

M(Ind) = Z;
M = flipdim(M,3);
M = flipdim(M,2);
V = blob2img(IMG,outputname,M);


%=======================================================================
function V = saverecons(MRIcurrent,IMG,outputimg)

Ind = find(MRIcurrent);
Z = MRIcurrent(Ind);

M = zeros(IMG.dim(1:3));

M(Ind) = Z;
M = flipdim(M,3);
M = flipdim(M,2);
V = blob2img(IMG,outputname,M);

return
%=======================================================================

%=======================================================================
function V = blob2img(IMG,filename,image)

V = IMG;
[pathstr, name, ext] = fileparts(V.fname);
V = struct('fname',   filename,...
            'dim',    [V.dim(1:3), spm_type('float')],...
            'mat',    V.mat,...
            'pinfo',  [1 0 0]',...
            'descrip', 'EEG/MEG activation');
V = spm_create_vol(V);
for j=1:V.dim(3),
        V = spm_write_plane(V,image(:,:,j),j);
end
V = spm_close_vol(V);

return
%=======================================================================

