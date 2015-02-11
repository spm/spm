function Vo = spm_write_filtered(Z,XYZ,DIM,M,descrip,F)
% Write the filtered SPM as an image
% FORMAT Vo = spm_write_filtered(Z,XYZ,DIM,M,descrip,F)
% Z         - {1 x ?} vector point list of SPM values for MIP
% XYZ       - {3 x ?} matrix of coordinates of points (voxel coordinates)
% DIM       - image dimensions {voxels}
% M         - voxels -> mm matrix [default: spm_matrix(-(DIM+1)/2)]
% descrip   - description string [default: 'SPM-filtered']
% F         - output file's basename [default: user query]
%
% FORMAT V0 = spm_write_filtered(xSPM)
% xSPM      - SPM results structure from spm_getSPM
%
% Vo        - output image volume information
%__________________________________________________________________________
%
% spm_write_filtered takes a pointlist image (parallel matrices of
% co-ordinates and voxel intensities), and writes it out into an image
% file.
%
% It is intended for writing out filtered SPM's from the results section
% of SPM, but can be used freestanding.
%__________________________________________________________________________
% Copyright (C) 1996-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_write_filtered.m 6337 2015-02-11 18:46:30Z guillaume $


%-Parse arguments
%--------------------------------------------------------------------------
if nargin < 2
    if nargin == 0
        try
            xSPM = evalin('base', 'xSPM');
        catch
            error('Please run SPM results query first.')
        end
    else
        xSPM = Z;
    end
    Vo = spm_write_filtered(xSPM.Z, xSPM.XYZ, xSPM.DIM, xSPM.M);
    return
elseif nargin < 3
    error('Not enough input arguments.');
end
if nargin<4, M = spm_matrix(-(DIM+1)/2); end
if nargin<5, descrip = 'SPM-filtered'; end
if nargin<6 || isempty(F), F = spm_input('Output filename',1,'s'); end

%-Get filename
%--------------------------------------------------------------------------
if isempty(F), F = 'output'; end
if isempty(spm_file(F,'ext'))
    F = spm_file(F,'ext',spm_file_ext);
end
F = spm_file(F,'CPath');

%-Set up header information
%--------------------------------------------------------------------------
Vo  = struct(...
        'fname',   F,...
        'dim',     DIM',...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'pinfo',   [1 0 0]',...
        'mat',     M,...
        'descrip', descrip);
if all(Z==1) % binary map
    Vo.dt(1) = spm_type('uint8');
elseif all(ismember(Z,0:max(Z))) % n-ary map
    Vo.dt(1) = spm_type('uint16');
end

%-Reconstruct (filtered) image from XYZ & Z pointlist
%--------------------------------------------------------------------------
Y      = NaN(DIM(1:3)');
OFF    = XYZ(1,:) + DIM(1)*(XYZ(2,:)-1 + DIM(2)*(XYZ(3,:)-1));
Y(OFF) = Z.*(Z > 0);
    
%-Write the reconstructed volume
%--------------------------------------------------------------------------
Vo = spm_write_vol(Vo,Y);
