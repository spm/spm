function obj = subsasgn(obj,subs,dat)
% Overloaded subsasgn function for file_array objects.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id$


if isempty(subs)
    t = obj;
    return;
end;
if numel(subs)~=1, error('Expression too complicated');end;
if ~strcmp(subs.type,'()'),
if ~strcmp(subs.type,'()'),
    if strcmp(subs.type,'.'),
        %error('Attempt to reference field of non-structure array.');
        if numel(struct(obj))~=1,
            error('Can only change the fields of simple file_array objects.');
        end;
        switch(subs.subs)
        case 'fname',      obj = fname(obj,dat);
        case 'dtype',      obj = dtype(obj,dat);
        case 'offset',     obj = offset(obj,dat);
        case 'dim',        obj = dim(obj,dat);
        case 'scl_slope',  obj = scl_slope(obj,dat);
        case 'scl_inter',  obj = scl_inter(obj,dat);
        otherwise error(['Reference to non-existent field "' subs.type '".']);
        end;
        return;
    end;
    if strcmp(subs.type,'{}'), error('Cell contents reference from a non-cell array object.'); end;
end;
    if strcmp(subs.type,'{}'), error('Cell contents reference from a non-cell array object.'); end;
end;

dm   = size(obj);
sobj = struct(obj);

if length(subs.subs) < length(dm),
    l   = length(subs.subs);
    dm  = [dm(1:(l-1)) prod(dm(l:end))];
    if numel(sobj) ~= 1,
        error('Can only reshape simple file_array objects.');
    end;
    if numel(sobj.scl_slope)>1 || numel(sobj.scl_inter)>1,
        error('Can not reshape file_array objects with multiple slopes and intercepts.');
    end;
end;

do   = ones(1,16);
args = {};
for i=1:length(subs.subs),
    if ischar(subs.subs{i}),
        if ~strcmp(subs.subs{i},':'), error('This shouldn''t happen....'); end;
        args{i} = int32(1:dm(i));
    else
        args{i} = int32(subs.subs{i});
    end;
    do(i) = length(args{i});
end;
if length(sobj)==1
    sobj.dim = dm;
    if numel(dat)~=1,
        subfun(sobj,double(dat),args{:});
    else
        dat1 = double(dat) + zeros(do);
        subfun(sobj,dat1,args{:});
    end;
else
    t = zeros(do);
    for j=1:length(sobj),
        ps  = [sobj(j).pos ones(1,length(args))];
        dm  = [sobj(j).dim ones(1,length(args))];
        siz = ones(1,16);
        for i=1:length(args),
            msk      = args{i}>=ps(i) & args{i}<(ps(i)+dm(i));
            args2{i} = find(msk);
            args3{i} = int32(double(args{i}(msk))-ps(i)+1);
            siz(i)   = numel(args2{i});
        end;
        if numel(dat)~=1,
            dat1 = double(subsref(dat,struct('type','()','subs',{args2})));
        else
            dat1 = double(dat) + zeros(siz);
        end;
        subfun(sobj(j),dat1,args3{:});
    end
end
return

function sobj = subfun(sobj,dat,varargin)
va = varargin;

if ~isempty(sobj.scl_inter),
    inter = sobj.scl_inter;
    if numel(inter)>1,
        inter = resize_scales(inter,sobj.dim,varargin);
    end;
    dat = double(dat) - inter;
end;

if ~isempty(sobj.scl_slope),
    slope = sobj.scl_slope;
    if numel(slope)>1,
        slope = resize_scales(slope,sobj.dim,varargin);
        dat   = double(dat)./slope;
    else
        dat   = double(dat)/slope;
    end;
end;

dt  = datatypes;
ind = find(cat(1,dt.code)==sobj.dtype);
if isempty(ind) error('Unknown datatype'); end;
if dt(ind).isint, dat = round(dat); end;
dat   = feval(dt(ind).conv,dat);
nelem = dt(ind).nelem;
if nelem==1,
    mat2file(sobj,dat,va{:});
elseif nelem==2,
    sobj1       = sobj;
    sobj1.dim   = [2 sobj.dim];
    sobj1.dtype = dt(find(strcmp(dt(ind).prec,{dt.prec}) & (cat(2,dt.nelem)==1))).code;
    dat         = reshape(dat,[1 size(dat)]);
    dat         = [real(dat) ; imag(dat)];
    mat2file(sobj1,dat,int32([1 2]),va{:});
else
    error('Inappropriate number of elements per voxel.');
end;
return
