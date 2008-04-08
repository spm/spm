function fm = pm_make_fieldmap(P,flags)
% This function creates an unwrapped fieldmap (in Hz) from either 
% a single or double echo complex image volume. In the case of a 
% "single-echo" image, that will have been created by the vendor
% sequence out of two acquisitions with different echo times. 
% The complex image volume(s) may consist of either real and 
% imaginary OR phase and magnitude components 
%
% FORMAT fm = pm_make_fieldmap(P,flags);
%
% Input: 
% P           : A matrix of 2 or 4 filenames, or 
%               a struct array of 2 or 4 memory mapped image volumes.
% flags       : Struct containing parameters guiding the unwrapping.
% .iformat    : 'RI' or 'PM' 
%               'RI' - input images are Real and Imaginary. (default)
%               'PM' - input images are Phase and Magnitude
% .method     : 'Huttonish', 'Mark3D' or 'Mark2D'
%               'Huttonish': Flood-fill based unwrapping progressing 
%                from low to high uncertainty areas. 
%               'Mark3D': Region-merging based method merging 3D
%                regions starting with the big ones. (default)
%               'Mark2D': Region-merging based method merging
%                slicewise 2D regions until all connected regions
%                within slices have been merged before moving on
%                to merging the slices.
% .fwhm       : FWHM (mm) of Gaussian filter used to implement
%               a weighted (with the reciprocal of the angular
%               uncertainty) smoothing of the unwrapped maps.
%               (default: 10mm)
% .pad        : Size (in-plane voxels) of padding kernel. This
%               is an option to replace non-unwrapped voxels
%               (i.e. those that have been considered to noisy)
%               with an average of neighbouring unwrapped voxels.
%               The size defines the size of the neighbourhood.
%               (default = 0);
% .etd        : Echo time difference (ms).(default = 10) 
% .ws         : Weighted or unweighted smoothing (default = 1)
% .bmask      : Brain mask
%
% Output:
% fm          : Structure containing fieldmap information
% The elements of the fm structure are:
%   fm.upm    : unwrapped fieldmap in Hz 
%   fm.mask   : binary image used to mask fieldmap
%   fm.opm    : phase map in radians 
%   fm.jac    : Jacobian of the fieldmap
%_______________________________________________________________________
%
%     .iformat = 'RI'  (this the default mode if not specified)
%
% P(1)        : real part of complex fieldmap image
% P(2)        : imaginary part of complex fieldmap image
%       OR
% P(1)        : real part of short echo time image
% P(2)        : imaginary part of short echo time image
% P(3)        : real part of long echo time image
% P(4)        : imaginary part of long echo time image
% 
%     Mode = 'PM'
% 
% P(1)        : phase image
% P(2)        : magnitude image
%       OR
% P(1)        : phase of short echo time image
% P(2)        : magnitude of short echo time image
% P(3)        : real part of long echo time image
% P(4)        : imaginary part of long echo time image
%__________________________________________________________
% Chloe Hutton, Jesper Andersson 05/08/22
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson and Chloe Hutton
% $Id: pm_make_fieldmap.m 1317 2008-04-08 16:16:38Z chloe $

% The default values below will be used if any flags haven't been defined 
% or are empty.

def_flags = struct('iformat',     'RI',...
                   'method',      'Mark3D',...
                   'fwhm',        10,...
                   'bmask',        [],...
                   'pad',         0,...
                   'etd',         10,...
                   'ws',          1);

%
% Do input argument checking
%

if nargin < 1
   help pm_make_fieldmap
   return
end

if nargin < 2
   flags = def_flags;
end

%
% Use defaults to fill out missing or empty flags
%

def_fields = fieldnames(def_flags);
for i=1:length(def_fields)
   if ~isfield(flags,def_fields{i})
      flags = setfield(flags,def_fields{i},getfield(def_flags,def_fields{i}));
   elseif isempty(getfield(flags,def_fields{i}))
      flags = setfield(flags,def_fields{i},getfield(def_flags,def_fields{i}));
   end
end

%
% Check for unknown flags (possibly reflecting misspellings).
%

flags_fields = fieldnames(flags);
for i=1:length(flags_fields)
   if ~isfield(def_flags,flags_fields{i})
      error(sprintf('pm_make_fieldmap: Unknown flag %s passed',flags_fields{i}));
   end
end

if ~isfield(P,'dim')
   P = spm_vol(P);
end

nfiles = length(P);
if nfiles~=2 & nfiles~=4
   error('Requires 2 or 4 images');
   help pm_make_fieldmap

elseif nfiles == 2 & flags.iformat == 'RI'
   I1 = spm_read_vols(P(1));
   I2 = spm_read_vols(P(2));
   cmap = I1 + sqrt(-1)*I2;

elseif nfiles == 4 & flags.iformat == 'RI'
   cte1 = spm_read_vols(P(1)) + sqrt(-1)*spm_read_vols(P(2));
   cte2 = spm_read_vols(P(3)) + sqrt(-1)*spm_read_vols(P(4));
   cmap = abs(cte2).*exp(sqrt(-1)*angle(cte2))./exp(sqrt(-1)*angle(cte1));

elseif nfiles == 2 & flags.iformat == 'PM'
   phase = spm_read_vols(P(1));
   mag = spm_read_vols(P(2));
   cmap = mag.*exp(sqrt(-1)*phase);

elseif nfiles == 4 & flags.iformat == 'PM'
   cte1 = spm_read_vols(P(2)).*exp(sqrt(-1)*spm_read_vols(P(1)));
   cte2 = spm_read_vols(P(4)).*exp(sqrt(-1)*spm_read_vols(P(3)));
   cmap = abs(cte2).*exp(sqrt(-1)*angle(cte2))./exp(sqrt(-1)*angle(cte1));
end

% Check for any NaNs and convert them to phase=-pi
nnan=find(isnan(cmap));
cmap(nnan)=exp(sqrt(-1)*-pi);

%
% Do unwrapping...
%

pxs = sqrt(sum(P(1).mat(1:3,1:3).^2));
fm=struct('upm',[],'mask',[],'opm',[]);

[fm.upm,angvar,fm.mask,fm.opm] = pm_unwrap(cmap,pxs,flags.method);

% Chloe added these lines on 28/02/05
% If there is a brain mask then this should be used to
% mask out the outside of the brain.

if ~isempty(flags.bmask)
   fm.mask=double(flags.bmask);
end

% Do a region growing such that any voxel that hasn't been
% unwrapped and are less than flags.pad/2 voxels away 
% from an unwrapped one will be replaced by an average of 
% the surrounding unwrapped voxels.

if flags.pad ~= 0
   npxs = pxs/pxs(1);
   ks = floor((flags.pad+1)*[1 1 1]./npxs);
   if ~mod(ks(1),2) ks(1) = ks(1)+1; end
   if ~mod(ks(2),2) ks(2) = ks(2)+1; end
   if ~mod(ks(3),2) ks(3) = ks(3)+1; end
   kernel = ones(ks);
   [fm.upm,fm.mask] = pm_pad(fm.upm,fm.mask,kernel);
end

%
% Next do a weighted (by the inverse of the noise)
% gaussian smoothing of the phase-map, considering
% only the unwrapped voxels.
% Note that the resulting weighted kernel is NOT
% separable, why we have to do the whole (slow)
% 3D thingie.
% 

if flags.ws==1 
   fm.fpm = pm_smooth_phasemap(fm.upm.*fm.mask,angvar,pxs,flags.fwhm);
else
   disp('Using normal smoothing');
   fwhm = repmat(flags.fwhm,1,3)./pxs;
   fm.fpm = zeros(size(fm.upm));
   spm_smooth(fm.upm.*fm.mask,fm.fpm,fwhm);
end

fm.jac = pm_diff(fm.fpm,2);

%
% Convert from radians to Hz
%
if isfield(fm,'fpm')
   fm.fpm = fm.fpm/(2*pi*flags.etd*1e-3);
   fm.jac = fm.jac/(2*pi*flags.etd*1e-3);
end

fm.upm = fm.upm/(2*pi*flags.etd*1e-3);

return












































