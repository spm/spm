function Vo = spm_write_filtered(Z,XYZ,DIM,M,descrip,F)
% Writes the filtered SPM as an image
% FORMAT spm_write_filtered(Z,XYZ,DIM,M,descrip)
%
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (voxel coordinates)
% DIM     - image dimensions {voxels}
% M       - voxels - > mm matrix (used for header file)
%           [default spm_matrix(-(DIM+1)/2) ]
% descrip - description string [default 'SPM-filtered']
%
%-----------------------------------------------------------------------
%
% spm_write_filtered takes a pointlist image (parallel matrixes of
% co-ordinates and voxel intensities), and writes it out into an image
% file.
%
% It is intended for writing out filtered SPM's from the results
% section of SPM, but can be used freestanding.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_write_filtered.m 3100 2009-05-06 19:00:39Z guillaume $


%-Parse arguments
%-----------------------------------------------------------------------
if nargin<3, error('Insufficient arguments'), end
if nargin<4, M = spm_matrix(-(DIM+1)/2); end
if nargin<5, descrip='SPM-filtered'; end
if nargin<6, F = spm_input('Output filename',1,'s'); end

%-Get filename
%-----------------------------------------------------------------------
Q = [spm_str_manip(F,'sd'), '.img'];
spm('Pointer','Watch')

%-Set up header information
%-----------------------------------------------------------------------
Vo      = struct(...
        'fname',    Q,...
        'dim',      DIM',...
        'dt',       [spm_type('uint8') spm_platform('bigend')],...
        'mat',      M,...
        'descrip',  descrip);

%-Reconstruct (filtered) image from XYZ & Z pointlist
%-----------------------------------------------------------------------
Y      = zeros(DIM(1:3)');
OFF    = XYZ(1,:) + DIM(1)*(XYZ(2,:)-1 + DIM(2)*(XYZ(3,:)-1));
Y(OFF) = Z.*(Z > 0);

%-Write the reconstructed volume
%-----------------------------------------------------------------------
Vo = spm_write_vol(Vo,Y);
spm('alert"',{'Written:',['    ',spm_select('CPath',Q)]},mfilename,1);

%-End
%-----------------------------------------------------------------------
spm('Pointer','Arrow');
