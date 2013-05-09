function volumes = spm_sumvol(P)
% Compute volume of a binary mask
% FORMAT volumes = spm_sumvol(P)
% P       - list of images
% volumes - volume of each mask
%__________________________________________________________________________
%
% This is a handy-ish utility function.  Essentially, it sums up each image
% and multiplies the result by the volume of a voxel.  It is intended for
% computing ``globals'' for VBM studies.
%__________________________________________________________________________
% Copyright (C) 1998-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_sumvol.m 5485 2013-05-09 15:51:24Z john $

if nargin<1,
    P       = spm_select(Inf,'nifti','Select files');
end

volumes = zeros(size(P,1),1);                % Tissue volumes

for i=1:size(P,1),                           % Loop over images
    Nii  = nifti(deblank(P(i,:)));           % Read file information
    X    = Nii.dat(:,:,:,:,:);               % Read data
    vvol = abs(det(Nii.mat(1:3,1:3)))/100^3; % Voxel volume in litres
    volumes(i) = sum(X(isfinite(X)))*vvol;   % Compute the volume
end

