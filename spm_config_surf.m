function opts = spm_config_surf
% Configuration file for surface extraction jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Volkmar Glauche
% $Id: spm_config_surf.m 256 2005-10-17 18:57:24Z guillaume $

data.type = 'files';
data.name = 'Grey+white matter image';
data.tag  = 'data';
data.filter = 'image';
data.num  = [2];
data.help = {'Images to create rendering/surface from (grey and white matter segments).'};

mode.type = 'menu';
mode.name = 'Output';
mode.tag  = 'mode';
mode.labels = {'Save Rendering', 'Save Extracted Surface',...
               'Save Rendering and Surface', 'Save Surface as OBJ format'};
mode.values = {1, 2, 3, 4};
mode.val  = {3};

opts.type = 'branch';
opts.name = 'Create Rendering/Surface';
opts.tag  = 'spm_surf';
opts.val  = {data,mode};
opts.vfiles = @filessurf;
opts.prog = @runsurf;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function runsurf(varargin)
spm_surf(strvcat(varargin{1}.data),varargin{1}.mode);
return;

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vfiles=filessurf(varargin)
vfiles=[];
[pth,nam,ext] = fileparts(varargin{1}.data{1});

if any(varargin{1}.mode==[1 3]),
	vfiles = fullfile(pth,['render_' nam '.mat']);
end;

if any(varargin{1}.mode==[2 3 4]),
	if any(varargin{1}.mode==[2 3]),
		vfiles = strvcat(vfiles,...
                                 fullfile(pth,['surf_' nam '.mat']));
        end;
        if any(varargin{1}.mode==[4]),
		vfiles = strvcat(vfiles,...
                                 fullfile(pth,[nam '.obj']));
        end;
end
vfiles = cellstr(vfiles);
return;
