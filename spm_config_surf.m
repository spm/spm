function opts = spm_config_surf
% Configuration file for surface extraction jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Volkmar Glauche
% $Id: spm_config_surf.m 775 2007-03-26 16:57:01Z john $

data.type = 'files';
data.name = 'Grey+white matter image';
data.tag  = 'data';
data.filter = 'image';
data.num  = [1 Inf];
data.help = {'Images to create rendering/surface from (grey and white matter segments).'};

mode.type = 'menu';
mode.name = 'Output';
mode.tag  = 'mode';
mode.labels = {'Save Rendering', 'Save Extracted Surface',...
               'Save Rendering and Surface', 'Save Surface as OBJ format'};
mode.values = {1, 2, 3, 4};
mode.val  = {3};

thresh.type    = 'entry';
thresh.name    = 'Surface isovalue(s)';
thresh.tag     = 'thresh';
thresh.num     = [1 Inf];
thresh.val     = {.5};
thresh.strtype = 'e';
thresh.help    = {['Enter one or more values at which isosurfaces through ' ...
                  'the input images will be computed.']};

opts.type = 'branch';
opts.name = 'Create Rendering/Surface';
opts.tag  = 'spm_surf';
opts.val  = {data,mode,thresh};
opts.vfiles = @filessurf;
opts.prog = @runsurf;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function runsurf(varargin)
spm_surf(strvcat(varargin{1}.data),varargin{1}.mode,varargin{1}.thresh);
return;

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vfiles=filessurf(varargin)
vfiles={};
[pth,nam,ext] = fileparts(varargin{1}.data{1});

if any(varargin{1}.mode==[1 3]),
	vfiles{1} = fullfile(pth,['render_' nam '.mat']);
end;

if any(varargin{1}.mode==[2 3 4]),
    for k=1:numel(varargin{1}.thresh)
        if numel(varargin{1}.thresh) == 1
            nam1 = nam;
        else
            nam1 = sprintf('%s-%d', nam, k);
        end;
	if any(varargin{1}.mode==[2 3]),
		vfiles{end+1} = fullfile(pth,['surf_' nam1 '.mat']);
        end;
        if any(varargin{1}.mode==[4]),
		vfiles{end+1} = fullfile(pth,[nam1 '.obj']);
        end;
    end;
end
return;
