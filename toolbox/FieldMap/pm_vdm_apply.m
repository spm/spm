function varargout = pm_vdm_apply(ds,flags)
%
% Applies vdm (voxel displacement map) to distortion correct images volume by volume
% FORMAT pm_vdm_apply(ds,[flags])
% or
% FORMAT P = pm_vdm_apply(ds,[flags])
%               
% ds           - a structure the containing the fields:
% 
% .P           - Images to apply correction to.
%
% .sfP         - Static field supplied by the user. It should be a 
%                filename or handle to a voxel-displacement map in
%                the same space as the first EPI image of the time-
%                series. If using the FieldMap toolbox, realignment
%                should (if necessary) have been performed as part of
%                the process of creating the VDM. Note also that the
%                VDM must be in undistorted space, i.e. if it is
%                calculated from an EPI based field-map sequence
%                it should have been inverted before using it to correct 
%                images. Again, the FieldMap toolbox will do this for you.
% .jm          - Jacobian Modulation. This option is currently disabled.
%
%            ds can also be an array of structures, each struct corresponding
%            to one sesssion and containing the corresponding images and
%            fieldmap. A preceding realignment step will ensure that
%            all images are in the space of the first image of the first session
%            and they are all resliced in the space of that image.There will
%            be a vdm file associated with each session but in practice,  
%            fieldmaps are not often measured per session. In this case the same
%            fieldmap is used for each session.
%
% flags    - a structure containing various options.  The fields are:
%
%         mask - mask output images (1 for yes, 0 for no)
%                To avoid artifactual movement-related variance the realigned
%                set of images can be internally masked, within the set (i.e.
%                if any image has a zero value at a voxel than all images have
%                zero values at that voxel).  Zero values occur when regions
%                'outside' the image are moved 'inside' the image during
%                realignment.
%
%         mean - write mean image
%                The average of all the realigned scans is written to
%                mean*.img.
%
%         interp - the interpolation method (see e.g. spm_bsplins.m).
%
%         which - Values of 0 or 1 are allowed.
%                 0   - don't create any resliced images.
%                       Useful if you only want a mean resliced image.
%                 1   - reslice all the images.
%
%         udc - Values 1 or 2 are allowed
%               1   - Do only unwarping (not correcting 
%                     for changing sampling density).
%               2   - Do both unwarping and Jacobian correction.
%
%
%             The spatially realigned and corrected images are written to the 
%             orginal subdirectory with the same filename but prefixed with an 'c'.
%             They are all aligned with the first.
%_______________________________________________________________________
% @(#)spm_vdm_apply.m	1.0 Chloe Hutton   05/02/25

global defaults

def_flags = struct('mask',       1,...
                   'mean',       1,...
                   'interp',     4,...
                   'wrap',       [0 1 0],...
                   'which',      1,...
                   'udc',        1);

defnames = fieldnames(def_flags);

%
% Replace hardcoded defaults with spm_defaults
% when exist and defined.
%
if exist('defaults','var') & isfield(defaults,'realign') & isfield(defaults.realign,'write')
   wd = defaults.realign.write;
   if isfield(wd,'interp'),    def_flags.interp = wd.interp; end
   if isfield(wd,'wrap'),      def_flags.wrap = wd.wrap; end
   if isfield(wd,'mask'),      def_flags.mask = wd.mask; end
end

if nargin < 1 | isempty(ds)
    ds   = getfield(load(spm_get(1,'*_c.mat','Select Unwarp result file')),'ds');
end

%
% Do not allow Jacobian modulation 
%
if ds(1).jm ~= 0
   ds(1).jm=0;
   def_flags.udc = 0;
end

%
% Replace defaults with user supplied values for all fields
% defined by user. Also, warn user of any invalid fields,
% probably reflecting misspellings.
%
if nargin < 2 | isempty(flags)
   flags = def_flags;
end
for i=1:length(defnames)
   if ~isfield(flags,defnames{i})
      flags = setfield(flags,defnames{i},getfield(def_flags,defnames{i}));
   end
end
flagnames = fieldnames(flags);
for i=1:length(flagnames)
   if ~isfield(def_flags,flagnames{i})
      warning(sprintf('Warning, unknown flag field %s',flagnames{i}));
   end
end

ntot = 0;
for i=1:length(ds)
   ntot = ntot + length(ds(i).P);
end

hold = [repmat(flags.interp,1,3) flags.wrap];

linfun = inline('fprintf(''%-60s%s'', x,repmat(sprintf(''\b''),1,60))');

%
% Create empty sfield for all structs.
%
[ds.sfield] = deal([]);

%
% Make space for output P-structs if required
%
if nargout > 0
   oP = cell(length(ds),1);
end

%
% First, create mask if so required.
%

if flags.mask | flags.mean,
   linfun('Computing mask..');
   spm_progress_bar('Init',ntot,'Computing available voxels',...
                    'volumes completed');
   [x,y,z] = ndgrid(1:ds(1).P(1).dim(1),1:ds(1).P(1).dim(2),1:ds(1).P(1).dim(3));
   xyz = [x(:) y(:) z(:) ones(prod(ds(1).P(1).dim(1:3)),1)]; clear x y z;
   if flags.mean
      Count    = zeros(prod(ds(1).P(1).dim(1:3)),1);
      Integral = zeros(prod(ds(1).P(1).dim(1:3)),1);
   end
   if flags.mask 
      msk = zeros(prod(ds(1).P(1).dim(1:3)),1);  
   end
   tv = 1;
   for s=1:length(ds)

      if isfield(ds(s),'sfP') & ~isempty(ds(s).sfP)
         T = ds(s).sfP.mat\ds(1).P(1).mat;
         txyz = xyz * T';
         c = spm_bsplinc(ds(s).sfP,ds(s).hold);
         ds(s).sfield = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds(s).hold);
         ds(s).sfield = ds(s).sfield(:);
         clear c txyz;
      end 

      for i = 1:prod(size(ds(s).P))
         T = inv(ds(s).P(i).mat) * ds(1).P(1).mat;
         txyz = xyz * T';
         msk = msk + real(txyz(:,1) < 1 | txyz(:,1) > ds(s).P(1).dim(1) |...
                          txyz(:,2) < 1 | txyz(:,2) > ds(s).P(1).dim(2) |...
                          txyz(:,3) < 1 | txyz(:,3) > ds(s).P(1).dim(3)); 
         spm_progress_bar('Set',tv);
         tv = tv+1;
      end

      if flags.mean, Count = Count + repmat(length(ds(s).P),prod(ds(s).P(1).dim(1:3)),1) - msk; end

      %
      % Include static field in estimation of mask.
      %
      if isfield(ds(s),'sfP') & ~isempty(ds(s).sfP)
         T = inv(ds(s).sfP.mat) * ds(1).P(1).mat;
         txyz = xyz * T';
         msk = msk + real(txyz(:,1) < 1 | txyz(:,1) > ds(s).sfP.dim(1) |...
                          txyz(:,2) < 1 | txyz(:,2) > ds(s).sfP.dim(2) |...
                          txyz(:,3) < 1 | txyz(:,3) > ds(s).sfP.dim(3)); 
      end

      if isfield(ds(s),'sfield') & ~isempty(ds(s).sfield)
         ds(s).sfield = [];
      end

   end
   if flags.mask, msk = find(msk ~= 0); end
end

linfun('Reslicing fieldmap corrected images..');
spm_progress_bar('Init',ntot,'Reslicing','volumes completed');

tiny = 5e-2; % From spm_vol_utils.c

PO = ds(1).P(1);
PO.descrip = 'spm - fieldmap corrected';
jP = ds(1).P(1);
jP = rmfield(jP,{'fname','descrip','n','private'});
jP.dim = [jP.dim(1:3) 64];
jP.pinfo = [1 0]';
tv = 1;

for s=1:length(ds)

   if isfield(ds(s),'sfP') & ~isempty(ds(s).sfP)
      T = ds(s).sfP.mat\ds(1).P(1).mat;
      txyz = xyz * T';
      c = spm_bsplinc(ds(s).sfP,ds(s).hold);
      ds(s).sfield = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds(s).hold);
      ds(s).sfield = ds(s).sfield(:);
      clear c txyz;
   end 

   for i = 1:length(ds(s).P)
      linfun(['Reslicing volume ' num2str(tv) '..']);

      %
      % Read undeformed image.
      %

      T = inv(ds(s).P(i).mat) * ds(1).P(1).mat;
      txyz = xyz * T';
      def = -ds(s).sfield;
      txyz(:,2) = txyz(:,2) + def;

      c = spm_bsplinc(ds(s).P(i),hold);
      ima = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),hold);
 
      %
      % Write it if so required.
      %
      if flags.which
         PO.fname   = prepend(ds(s).P(i).fname,'c');
         if flags.mask
            ima(msk) = NaN;
         end
         ivol = reshape(ima,PO.dim(1:3));
         tP = spm_write_vol(PO,ivol);
	 if nargout > 0
	    oP{s}(i) = tP;
	 end
      end
      %
      % Build up mean image if so required.
      %
      if flags.mean
         Integral = Integral + nan2zero(ima);
      end
      spm_progress_bar('Set',tv);
      tv = tv+1;
   end
   if isfield(ds(s),'sfield') & ~isempty(ds(s).sfield)
      ds(s).sfield = [];
   end
end

if flags.mean
   % Write integral image (16 bit signed)
   %-----------------------------------------------------------
   warning off; % Shame on me!
   Integral   = Integral./Count;
   warning on;
   PO         = ds(1).P(1);
   PO.fname   = prepend(ds(1).P(1).fname, 'meanc');
   PO.pinfo   = [max(max(max(Integral)))/32767 0 0]';
   PO.descrip = 'spm - mean undeformed image';
   PO.dim(4)  = 4;
   ivol = reshape(Integral,PO.dim(1:3));
   spm_write_vol(PO,ivol);
end

linfun(' ');
spm_figure('Clear','Interactive');

if nargout > 0
   varargout{1} = oP;
end

return;


%_______________________________________________________________________
function PO = prepend(PI,pre)
[pth,nm,xt,vr] = fileparts(deblank(PI));
PO             = fullfile(pth,[pre nm xt vr]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function vo = nan2zero(vi)
vo = vi;
vo(~isfinite(vo)) = 0;
return;
%_______________________________________________________________________

