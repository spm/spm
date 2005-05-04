function V = spm_create_vol(V,varargin)
% Create a volume
% FORMAT V = spm_create_vol(V)
% V - image volume information (see spm_vol.m)
%____________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_create_vol.m 112 2005-05-04 18:20:52Z john $


for i=1:prod(size(V)),
    if nargin>1,
        v = create_vol(V(i),varargin{:});
    else
        v = create_vol(V(i));
    end;
    f = fieldnames(v);
    for j=1:size(f,1),
        V(i).(f{j}) = v.(f{j});
    end;
end;

function V = create_vol(V,varargin)
if ~isstruct(V),        error('Not a structure.'); end;
if ~isfield(V,'fname'), error('No "fname" field'); end;

if ~isfield(V,'dim'), error('No "dim" field'); end;
if ~all(size(V.dim)==[1 3]),
    error(['"dim" field is the wrong size (' num2str(size(V.dim)) ').']);
end;

if ~isfield(V,'dt'),
    V.dt = [spm_type('float64') spm_platform('bigend')];
end;

dt{1} = spm_type(V.dt(1));
if strcmp(dt{1},'unknown'),
    error(['"' dt(1) '" is an unrecognised datatype (' num2str(V.dt(1)) ').']);
end;
if V.dt(2), dt{2} = 'BE'; else, dt{2} = 'LE'; end;

if ~isfield(V,'pinfo'), V.pinfo      = [1 0 0]'; end;
if size(V.pinfo,1)==2,  V.pinfo(3,1) = 0;        end;
V.fname       = deblank(V.fname);
[pth,nam,ext] = fileparts(V.fname);
switch ext,
case {'.img'}
    minoff = 0;
    V.pinfo(3,:) = max(V.pinfo(3,:),minoff);
case {'.nii'}
    minoff = 352;
    V.pinfo(3,:) = max(V.pinfo(3,:),minoff);
otherwise
    error(['"' ext '" is not a recognised extension.']);
end;

if ~isfield(V,'descrip'), V.descrip = '';    end;
if ~isfield(V,'n'),
    V.n       = [1 1];
else
    V.n       = [V.n(:)' 1 1];
    V.n       =  V.n(1:2);
end;
if ~isfield(V,'private'), V.private = [];    end;

dim    = [V.dim(1:3) V.n];
dat    = file_array(V.fname,dim,[dt{1} '-' dt{2}],V.pinfo(3),V.pinfo(1),V.pinfo(2));
N      = nifti;
N.dat  = dat;
N.mat  = V.mat;
N.mat0 = V.mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.descrip = V.descrip;

try,
    N0  = nifti(V.fname);
catch,
    N0  = [];
end;
if ~isempty(N0),
    tmp = [N0.dat.dim ones(1,5)];
    if any(tmp(1:3) ~= dim(1:3))
        warning(['Incompatible x,y,z dimensions in file "' V.fname '" [' num2str(tmp(1:3)) ']~=[' num2str(dim(1:3)) '].']);
    end;
    if dim(5) > tmp(5) && tmp(4) > 1,
        warning(['Incompatible 4th and 5th dimensions in file "' V.fname '" (' num2str([tmp(4:5) dim(4:5)]) ').']);
    end;
    N.dat.dim = max(dim(1:5),tmp(1:5));
    if ~strcmp(dat.dtype,N0.dat.dtype),
        warning(['Incompatible datatype in file "' V.fname '" ' N0.dat.dtype ' ~= ' dat.dtype '.']);
    end;
    if single(N.dat.scl_slope) ~= single(N0.dat.scl_slope) && (size(N0.dat,4)>1 || V.n(1)>1),
        warning(['Incompatible scalefactor in "' V.fname '" ' num2str(N0.dat.scl_slope) '~=' num2str(N.dat.scl_slope) '.']);
    end;
    if single(N.dat.scl_inter) ~= single(N0.dat.scl_inter),
        warning(['Incompatible intercept in "' V.fname '" ' num2str(N0.dat.scl_inter) '~=' num2str(N.dat.scl_inter) '.']);
    end;
    if single(N.dat.offset) ~= single(N0.dat.offset),
        warning(['Incompatible intercept in "' V.fname '" ' num2str(N0.dat.offset) '~=' num2str(N.dat.offset) '.']);
    end;

    if ~isempty(N0.extras) && isstruct(N0.extras) && isfield(N0.extras,'mat'),
        N0.extras.mat(:,:,V.n(1)) = N.mat;
        N.extras                  = N0.extras;
    end;
    if sum((V.mat(:)-N0.mat(:)).^2) > 1e-4,
        N.extras.mat(:,:,V.n(1)) = V.mat;
    end;
end;
create(N);
V.private = N;

