
function D = spm_eeg_make3dimage(S)
% function for converting 2D images to 3D volumes for ERPs
% FORMAT D = spm_eeg_downsample(S)
% 
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file

%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% James Kilner

try
    D = S.D;
catch
    D = spm_select(1, '\.img$', 'Select EEG image file');
end

v=spm_vol([D]);
data=zeros(v(1).dim(1),v(1).dim(1),length(v));
for n=1:length(v)
    [temp,XYZ]=spm_read_vols(v(n));
    data(:,:,n)=temp;
end
V=v(1);
[pathstr,name,ext,cersn] = fileparts(D) ;
V.fname = fullfile(pathstr,[name '3d.img']);
V.dim = [v(1).dim(1) v(1).dim(2) length(v)];
V.mat = eye(4);
%V.pinfo = [1 0 0]';
V = rmfield(V,'pinfo');

spm_write_vol(V, data); 
