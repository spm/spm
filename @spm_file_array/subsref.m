function t=subsref(opt,subs)
% SUBSREF Subscripted reference
% An overloaded function...
% _________________________________________________________________________________
% %W% John Ashburner %E%

if prod(size(subs))~=1, error('Expression too complicated');end;
if ~strcmp(subs.type,'()'),
	if strcmp(subs.type,'.') error('Attempt to reference field of non-structure array.'); end;
	if strcmp(subs.type,'{}') error('Cell contents reference from a non-cell array object.'); end;
end;

if length(subs.subs) < length(opt.dim),
	l = length(subs.subs);
	opt.dim = [opt.dim(1:(l-1)) prod(opt.dim(l:end))];
end;

args = {};
for i=1:length(subs.subs),
	if ischar(subs.subs{i}),
		if ~strcmp(subs.subs{i},':'), error('This shouldn''t happen....'); end;
		args{i} = int32(1:opt.dim(i));
	else,
		args{i} = int32(subs.subs{i});
	end;
	sz(i) = length(args{i});
end;

t = spm_farray(opt,args{:});

if ~isempty(opt.scale)
	d1 = opt.scale.dim;
	d2 = setxor(opt.scale.dim, 1:length(opt.dim));
	
	% reorder, reshape and multiply
	t = permute(t, [d1 d2]);
	t = reshape(t, [prod(sz(d1)) prod(sz(d2))]);
	s = opt.scale.values(args{d1});
	s = repmat(reshape(s, prod(sz(d1)), 1), 1, prod(sz(d2)));
	t = s .* t;
	
	% reassemble
	t = reshape(t, sz([d1 d2]));
	t = ipermute(t, [d1 d2]);
end