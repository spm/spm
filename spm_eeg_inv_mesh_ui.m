function D = spm_eeg_inv_mesh_ui(varargin)
% Cortical Mesh user-interface routine
% Invokes spatial normalization (if required) and the computation of
% the proper size individual size
%
% FORMAT D = spm_eeg_inv_mesh_ui(D, val, template, Msize)
% Input:
% D        - input data struct (optional)
% Output:
% D        - same data struct including the meshing files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_mesh_ui.m 2773 2009-02-23 13:08:55Z vladimir $


% initialise
%--------------------------------------------------------------------------
[Finter] = spm('FnUIsetup','Define head model',0);

[D,val] = spm_eeg_inv_check(varargin{:});

if val == 0
    val = 1;
end

if ~isfield(D, 'inv') || ~isfield(D.inv{val}, 'comment')
    D.inv = {struct('mesh', [])};
    clck = fix(clock);
    if clck(5) < 10
        clck = [num2str(clck(4)) ':0' num2str(clck(5))];
    else
        clck = [num2str(clck(4)) ':' num2str(clck(5))];
    end
    D.inv{val}.date    = strvcat(date,clck);
    D.inv{val}.comment = {''};
else
    inv = struct('mesh', []);
    inv.comment = D.inv{val}.comment;
    inv.date    = D.inv{val}.date;
    D.inv{val} = inv;
end

if nargin>2
    template = varargin{3};
else
    template = [];
end

if isempty(template)
    template = spm_input('Select head  model', '+1','template|individual', [1 0]);
end

if template
    sMRI = [];
else
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

% check meshes and display
%--------------------------------------------------------------------------
spm_eeg_inv_checkmeshes(D);
