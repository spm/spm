function s = export(this,target)
% Export a GIfTI object into specific MATLAB struct
% FORMAT s = export(this,target)
% this   - GIfTI object
% target - string describing target output [default: Matlab]
% s      - a structure containing public fields of the object
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: export.m 2574 2008-12-18 12:56:04Z guillaume $

if numel(this) > 1, warning('Only handle scalar objects yet.'); end

if nargin <= 1, target = 'Matlab'; end

switch lower(target)
    case 'matlab'
        s = struct(this);
        
    case {'fieldtrip', 'ft'}
        s = struct('tri',[], 'pnt',[]);
        if isfield(this,'vertices')
            s.pnt = subsref(this, substruct('.', 'vertices'));
        end
        if   isfield(this,'faces')  
            s.tri = subsref(this, substruct('.', 'faces'));
        end
        
    otherwise
        error('Unknown target ''%s''.', target);
end
