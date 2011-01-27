function varargout = pm_unwrap(varargin)
% Unwrapping of phasemap. When measuring phase one
% cannot easily distinguish between e.g. a phase of 182
% degrees, and one of -178 degrees. One tries to 
% distinguish these cases by using neighbourhood
% information. So in the example above, if we find that
% a neighbouring voxel has a phase of 150 degres it
% seems much more likely that the "true" phase is
% 182 degrees than -178 degrees. It's trickier than
% it sounds.
% FORMAT: [upm,(angvar),(mask),(opm)] = pm_unwrap(ci,pxs,method)
% or
% FORMAT: [upm,(angvar),(mask),(opm)] = pm_unwrap(ci,pxs)
% or
% FORMAT: [upm,(angvar),(mask),(opm)] = pm_unwrap(P,method)
% or
% FORMAT: [upm,(angvar),(mask),(opm)] = pm_unwrap(P)
% 
% Input:
% ci          : Complex image volume corresponding
%               to  abs(te2).*exp(i*angle(te2))./exp(i*angle(te1));
%               where te1 and te2 corresponds to the complex
%               images obtained with the short and the long
%               echo-time respectively, and i denotes sqrt(-1).
% pxs         : 3x1 (or 2x1) array with pixel sizes.
% 
% or
%
% P           : File structure (from) spm_vol, containing complex
%               image volume as per above.
%
% method      : Determines which method should be used
%               for phase-unwrapping. The options are
%               'Huttonish', 'Mark2D', 'Mark3D' and 'hybrid'.
% 'Huttonish' : Loosely (hence -ish) based on method described
%               in Hutton et al. Gets an estimate of the 
%               uncertainty of the phase angle at each point
%               and unwraps in a "watershed" fashion from
%               a high certainty seed towards more uncertain
%               areas.
% 'Mark2D'    : Method suggested for high-res data in 
%               Jenkinssons MRM paper.
% 'Mark3D'    : Method suggested for low-res data in
%               Jenkinssons MRM paper.
%
% Output:
% upm         : Phasemap (corresponding to angle(ci))
%               after unwrapping of phase jumps.
% angvar      : Map of the variance of the phase-angle
%               estimates. This is used internally to
%               guide the unwrapping procedure, and
%               can also be used if one whishes to
%               do a weighted fitting of some smooth
%               basis set to the unwrapped phasemap.
% mask        : Binary mask indicating what voxels
%               have been unwrapped.
% opm         : angle(ci)
%               
% Light reading:
% 
% Examples of water-shed/flood-fill based unwrapping 
% algorithms:
%
% Hutton C, Bork A, Josephs O, Deichmann R, Ashburner J,
% Turner R. 2002. Image distortion correction in fMRI: A
% quantitative evaluation. NeuroImage 16:217-240.

% Cusack R & Papadakis N. 2002. New robust 3-D phase unwrapping
% algorithms: Application to magnetic field mapping and
% undistorting echoplanar images. NeuroImage 16:754-764.
%
% Region-merging based unwrapping algorithm.
%
% Jenkinson M. 2003. Fast, automated, N-dimensional phase-
% unwrapping algorithm. MRM 49:193-197.
%__________________________________________________________
% Jesper Andersson 30/9-03
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_unwrap.m 4175 2011-01-27 12:42:36Z chloe $

%
% The following are a set of parameters that
% guide the phase unwrapping. Some are general
% (i.e. apply to all methods) and some are specific.
%

%
% These general parameters appear fine for a wide
% variety of phase-map data.
%
% General
%
mthres = (pi^2)/6; % Angular uncertainty threshold for mask.
ndil = 1;          % No. of erode->dilate when creating mask.
rr = 0;            % Should we remove linear ramps prior to unwrapping?

%
% These parameters guide the progression of the
% unwrapping waterfront. We have found it very 
% difficult to find a set of parameters that work
% for "all" data sets. Sometimes it is better to
% "pour" quickly in the beginning (large fthres) 
% and very slowly towards the end (large nthres).
% Other times it is better to pour at an even
% pace (fthres=0).
% This difficulty of finding a general set of
% parameters have lead us to choose 'Mark3D'
% as the default method of choice.
%
% Huttonish
%
fthres = 0;             % First (safe) threshold
lthres = (pi^2)/6;      % Last threshold (when to stop)
nthres = 100;           % No. of steps in between

%
% This, single, parameter guides the unwrapping
% using Mark J's method, and has been found to be 
% useful for a large variety of data sets.
%
% Mark (2 and 3D)
nstep = 8;             % No. of angular steps for initial regions.

%
% Fiddle about with input a bit to make sure it's OK.
%
if ~isnumeric(varargin{1})
   if ischar(varargin{1}) & isstruct(spm_vol(varargin{1}))
      P = spm_vol(varargin{1});
   elseif isstruct(varargin{1}) & isfield(varargin{1},'dim')
      P = spm_vol(varargin{1});
   else
      error('');
   end
   ci = spm_read_vols(P(1));
   pxs = sqrt(sum(P(1).mat(1:3,1:3).^2));
   if nargin < 2
      method = 'Mark3D';
   else
      method = varargin{2};
   end
else
   ci = varargin{1};
   pxs = varargin{2};
   if length(size(ci)) ~= length(pxs)
      error('');
   end
   if nargin < 3
      method = 'Mark3D';
   else
      method = varargin{3};
   end
end

%
% Generate (potentially wrapped) phase-map
% 
opm = angle(ci);

%
% Try to estimate noise properties
%
angvar = pm_angvar(ci);

%
% Create mask
%
mask = pm_mask(angvar,mthres,ndil);

%
% Remove linear ramps in data prior to unwrapping
% (and restore them after unwrapping). This may or may
% not be a good idea.
%
if rr
   [ramps,opm] = pm_estimate_ramp(opm,mask);
end

%
% Here comes unwrapping per se.
%
switch lower(method)
   case 'huttonish'
      %
      % Get seed point
      %
      seed = pm_seed(angvar,mask,pxs);
      %
      % Get series of thresholds that guide 
      % evolution of unwrapping front.
      %
      thres = linspace(fthres,lthres,nthres);
      %
      % Do unwrapping
      %
      upm = opm;
      wmap = zeros(size(opm));
      wmap(seed(1),seed(2),seed(3)) = 1;
      spm_progress_bar('Init',length(thres),'Unwrapping phase','Watershed step');
      for i=1:length(thres)
         spm_progress_bar('Set',i);
         [upm,wmap] = pm_ff_unwrap(upm,angvar,wmap,mask,thres(i));
      end
      spm_progress_bar('Clear');
   case 'mark3d'
      %
      % Do initial division into connected regions with
      % limited angular span.
      %
      [irima,cn] = pm_initial_regions(opm,mask,nstep);
      %
      % Added this little bug fix which prevents pm_merge_regions crashing
      % because it has too many regions to merge. 
      while cn>1800 & nstep > 2
          nstep=nstep-1;
          [irima,cn] = pm_initial_regions(opm,mask,nstep);
      end
      % Get connectogram
      %
      [ii,jj,nn,pp] = pm_create_connectogram(irima,opm);
      %
      % Merge regions while updating connectogram
      %
      rs = histc(irima(:),[0:max(irima(:))]+0.5);
      rs = rs(1:end-1); 
      upm = pm_merge_regions(opm,irima,ii,jj,nn,pp,rs);
      wmap = mask;
   case 'mark2d'
      upm = zeros(size(opm));
      vrima = zeros(size(opm));
      for sl=1:size(opm,3)
         %
         % Make intial (slicewise) regions.
         %
         [srima,cn] = pm_initial_regions(opm(:,:,sl),mask(:,:,sl),nstep);
         %
         % Get and merge slicewise connectograms
         %
         [ii,jj,nn,pp] = pm_create_connectogram(srima,opm(:,:,sl));
         %
         % Merge regions while updating connectogram
         %
         rs = histc(srima(:),[0:max(srima(:))]+0.5);
         rs = rs(1:end-1);
         upm(:,:,sl) = pm_merge_regions(opm(:,:,sl),srima,ii,jj,nn,pp,rs);
         %
         % Get connected regions after merging (can be more than
         % one if e.g. the temporal lobes are disconnected in the
     % present slice).
         %
         [tmp,num] = spm_bwlabel(srima,6);
         if sl == 1
        vrima(:,:,sl) = tmp;
         else
        tmp(tmp(:)>0) = tmp(tmp(:)>0) + max(max(vrima(:,:,sl-1)));
        vrima(:,:,sl) = tmp;
         end
      end
      %
      % Get and merge volumewise connectogram
      %
      [ii,jj,nn,pp] = pm_create_connectogram(vrima,upm);
      rs = histc(vrima(:),[0:max(vrima(:))]+0.5);
      rs = rs(1:end-1);
      upm = pm_merge_regions(upm,vrima,ii,jj,nn,pp,rs);
      wmap = mask;
   otherwise
      error(sprintf('Unknown method %s passed to my_unwrap',method));
end

%
% Restore any ramps we may have removed.
%
if rr
   upm = pm_restore_ramp(upm,mask,ramps);
end

mask = wmap;

if nargout > 0 varargout{1} = upm; end
if nargout > 1 varargout{2} = angvar; end
if nargout > 2 varargout{3} = mask; end
if nargout > 3 varargout{4} = opm; end

return;






