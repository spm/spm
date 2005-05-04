function opts = spm_config_smooth
% Configuration file for smoothing jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$


%_______________________________________________________________________

data.type = 'files';
data.name = 'Images to Smooth';
data.tag  = 'data';
data.filter = 'image';
data.num  = Inf;
data.help = {[...
'Specify the images to smooth. ',...
'The smoothed images are written to the same subdirectories as the ',...
'original *.img and are prefixed with a ''s'' (i.e. s*.img).']};
 
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

opts.type = 'branch';
opts.name = 'Smooth';
opts.tag  = 'smooth';
opts.val  = {data,fwhm};
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

P   = strvcat(job.data);
s   = job.fwhm;
n   = size(P,1);

spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
for i = 1:n
        Q = deblank(P(i,:));
        [pth,nam,xt,nm] = spm_fileparts(deblank(Q));
        U = fullfile(pth,['s' nam xt]);
        spm_smooth(Q,U,s);
        spm_progress_bar('Set',i);
end
spm_progress_bar('Clear',i);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vf = vfiles(varargin)
P  = varargin{1}.data;
vf = cell(size(P));
for i=1:numel(P),
    [pth,nam,ext,num] = spm_fileparts(P{i});
    vf{i} = fullfile(pth,['s', nam, ext, num]);
end;

