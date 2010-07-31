function D = spm_eeg_inv_mesh_ui(varargin)
% Cortical Mesh user interface
% FORMAT D = spm_eeg_inv_mesh_ui(D, val, sMRI, Msize)
% 
% D        - input data struct (optional)
% val      - 
% sMRI     -  0 - use template (default), or string with image file name
% Msize    - 
% 
% D        - same data struct including the meshing files and variables
%__________________________________________________________________________
%
% Invokes spatial normalisation (if required) and the computation of
% the individual mesh.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_mesh_ui.m 4027 2010-07-31 12:49:19Z vladimir $

SVNrev = '$Rev: 4027 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Define head model');

%-Initialisation
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

if val == 0
    val = 1;
end

if ~isfield(D, 'inv') || ~isfield(D.inv{val}, 'comment')
    D.inv = {struct('mesh', [])};
    D.inv{val}.date    = strvcat(date,datestr(now,15));
    D.inv{val}.comment = {''};
else
    inv = struct('mesh', []);
    inv.comment = D.inv{val}.comment;
    inv.date    = D.inv{val}.date;
    D.inv{val} = inv;
end

if nargin > 2
    sMRI = varargin{3};
else
    sMRI = [];
end

if isempty(sMRI)
    template = spm_input('Select head  model', '+1','template|individual', [1 0]);
elseif ~ischar(sMRI)
    template = sMRI; % for backward compatibility
else
    template = 0;
end

if template
    sMRI = [];
elseif ~ischar(sMRI)
    % get sMRI file name
    sMRI = spm_select([0 1],'image','Select subject''s structural MRI (Press Done if none)');
    if isempty(sMRI)
        error('No structural MRI selected.');
    end
end

if nargin>3
    Msize = varargin{4};
else
    Msize = spm_input('Cortical mesh', '+1', 'coarse|normal|fine', [1 2 3]);
end

D.inv{val}.mesh = spm_eeg_inv_mesh(sMRI, Msize);

%-Check meshes and display
%--------------------------------------------------------------------------
spm_eeg_inv_checkmeshes(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','Define head model: done');
