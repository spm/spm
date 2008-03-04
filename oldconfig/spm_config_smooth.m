function opts = spm_config_smooth
% Configuration file for smoothing jobs
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_config_smooth.m 1185 2008-03-04 16:31:21Z volkmar $


%_______________________________________________________________________

data.type = 'files';
data.name = 'Images to Smooth';
data.tag  = 'data';
data.filter = 'image';
data.num  = [0 Inf];
data.help = {[...
'Specify the images to smooth. ',...
'The smoothed images are written to the same subdirectories as the ',...
'original *.img and are prefixed with a ''s'' (i.e. s*.img). The prefix' ...
' can be changed by an option setting.']};
 
%------------------------------------------------------------------------

prefix.type = 'entry';
prefix.name = 'Filename Prefix';
prefix.tag  = 'prefix';
prefix.strtype = 's';
prefix.num  = [1 Inf];
prefix.val  = {'s'};
prefix.help = {[...
'Specify the string to be prepended to the filenames of the smoothed' ...
' image file(s). Default prefix is ''s''.']};

%------------------------------------------------------------------------

fwhm.type = 'entry';
fwhm.name = 'FWHM';
fwhm.tag  = 'fwhm';
fwhm.strtype = 'e';
fwhm.num  = [1 3];
fwhm.val  = {[8 8 8]};
fwhm.help = {[...
'Specify the full-width at half maximum (FWHM) of the Gaussian smoothing ',...
'kernel in mm. Three values should be entered, denoting the FWHM in the ',...
'x, y and z directions.']};

%------------------------------------------------------------------------

dtype.type = 'menu';
dtype.name = 'Data Type';
dtype.tag  = 'dtype';
dtype.labels = {'SAME','UINT8  - unsigned char','INT16 - signed short','INT32 - signed int',...
                'FLOAT - single prec. float','DOUBLE - double prec. float'};
dtype.values = {0,spm_type('uint8'),spm_type('int16'),spm_type('int32'),...
                spm_type('float32'),spm_type('float64')};
dtype.val = {0};
dtype.help = {'Data-type of output images.  SAME indicates the same datatype as the original images.'};

%------------------------------------------------------------------------

opts.type = 'branch';
opts.name = 'Smooth';
opts.tag  = 'smooth';
opts.val  = {data,fwhm,dtype,prefix};
opts.prog = @smooth;
opts.vfiles = @vfiles;
opts.help = {...
['This is for smoothing (or convolving) image volumes ',...
'with a Gaussian kernel of a specified width. ',...
'It is used as a preprocessing step to suppress noise and effects due to residual ',...
'differences in functional and gyral anatomy during inter-subject ',...
'averaging.']};

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function smooth(varargin)
job = varargin{1};

P     = strvcat(job.data);
s     = job.fwhm;
dtype = job.dtype;
n     = size(P,1);

spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
for i = 1:n
        Q = deblank(P(i,:));
        [pth,nam,ext,num] = spm_fileparts(deblank(Q));
        U = fullfile(pth,[job.prefix nam ext num]);
        spm_smooth(Q,U,s,dtype);
        spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vf = vfiles(varargin)
P  = varargin{1}.data;
vf = cell(size(P));
for i=1:numel(P),
    [pth,nam,ext,num] = spm_fileparts(P{i});
    vf{i} = fullfile(pth,[varargin{1}.prefix nam ext num]);
end;

