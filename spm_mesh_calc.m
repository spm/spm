function Mo = spm_mesh_calc(Mi,Mo,f,varargin)
% Evaluate a function on a mesh's data
% FORMAT Mo = spm_mesh_calc(Mi,Mo,f,opts)
% Mi   - input filenames (char array or cellstr)
%        or cell array of gifti objects or patch structures
% Mo   - output filename
%        if empty, a gifti object is returned and not saved on disk
% f    - MATLAB expression to be evaluated (string or function handle)
%        (e.g., f = '(s1.*s2).^2' or f = @(s1,s2) (s1.*s2).^2)
% opts - optional list of pairs of property names and values
%        dmtx - read images into data matrix X [default: false]
%__________________________________________________________________________
% Copyright (C) 2015-2019 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_calc.m 7577 2019-04-24 08:59:56Z guillaume $


%-Check input arguments
%--------------------------------------------------------------------------
if ~isempty(Mo)
    Mo = spm_file(Mo,'ext','.gii');
    Mo = spm_file(Mo,'cpath');
end
if numel(varargin) == 1 && isstruct(varargin{1})
    varargin = reshape([fieldnames(varargin{1})';struct2cell(varargin{1})'],1,[]);
end
if mod(numel(varargin),2)
    error('Incorrect number of input arguments.');
end
opts.dmtx = false;
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case 'dmtx'
            opts.dmtx = logical(varargin{i+1});
        otherwise
            error('Unknown option %s.',varargin{i});
    end
end

%-Load data
%--------------------------------------------------------------------------
if ischar(Mi), Mi = cellstr(Mi); end
for i=1:numel(Mi)
    g = gifti(Mi{i});
    if ~isfield(g,'cdata')
        error('File %s does not contain data.',Mi{i});
    end
    D = full(g.cdata);
    if i==1
        nv = size(D);
    else
        if ~isequal(size(D),nv)
            error('Data dimension mismatch.');
        end
    end
    D = reshape(D,1,[]);
    if opts.dmtx
        X(i,:) = D;
    else
        eval(['s',num2str(i),'=D;']);
    end
end

%-Evaluate function
%--------------------------------------------------------------------------
if ischar(f)
    try
        eval(['S = ' f ';']);
    catch
        l = lasterror;
        error('%s\nCan''t evaluate "%s".',l.message,f);
    end
elseif isa(f,'function_handle')
    try
        if opts.dmtx
            S = feval(f,X);
        else
            list = sprintf('s%d,',1:numel(Mi)); list = list(1:end-1);
            eval(['S = feval(f,' list ');']);
        end
    catch
        l = lasterror;
        error('%s\nCan''t evaluate "%s".',l.message,func2str(f));
    end
else
    error('Unknown function input.');
end

%-Return or save output
%--------------------------------------------------------------------------
gs = gifti(reshape(S,nv));
g  = gifti(Mi{1});
if isfield(g,'vertices') && isfield(g,'faces')
    gs.vertices  = g.vertices;
    gs.faces     = g.faces;
elseif ~isempty(g.private.metadata)
    metadata     = g.private.metadata;
    name         = {metadata.name};
    if any(ismember(name,'SurfaceID'))
        metadata = metadata(ismember(name,'SurfaceID'));
        gs.private.metadata(1) = metadata;
    end
end
if isempty(Mo), Mo = gs; else save(gs,Mo,'ExternalFileBinary'); end
