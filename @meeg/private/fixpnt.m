function data = fixpnt(data, recurse)
% Rename point structure fields (backward compatibility)
% FORMAT data = fixpnt(data, recurse)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if nargin==1
  recurse = 1;
end

if ~isa(data, 'struct')
    return;
end

if numel(data)>1
  % loop over all individual elements
  clear tmp
  for i=1:numel(data)
    % this is to prevent an "Subscripted assignment between dissimilar structures" error
    tmp(i) = fixpnt(data(i));
  end
  data = tmp;
  clear tmp
  return
end

% replace pos by pnt
if isfield(data, 'pos')
  data.pnt = data.pos;
  data = rmfield(data, 'pos');
end

if recurse<3
  % recurse into substructures, not too deep
  fn = fieldnames(data);
  fn = setdiff(fn, {'cfg'}); % don't recurse into the cfg structure
  for i=1:length(fn)
    if isstruct(data.(fn{i}))
      data.(fn{i}) = fixpnt(data.(fn{i}), recurse+1);
    end
  end
end
