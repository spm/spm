function grad = bti2grad(hdr);

% BTI2GRAD converts the position and weights of all coils that
% compromise a gradiometer system into a structure that can be used
% by FieldTrip.

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: bti2grad.m,v $
% Revision 1.2  2006/10/09 15:39:51  roboos
% renamed BTi channels from 'MEGxxx' into 'Axxx'
%
% Revision 1.1  2006/09/18 14:52:14  roboos
% implemented bti2grad and added it to header
%

grad     = [];
grad.pnt = hdr.Meg_pos;
grad.ori = hdr.Meg_dir;
for i=1:size(grad.pnt,1)
  % grad.label{i} = sprintf('MEG%03d', i);
  grad.label{i} = sprintf('A%d', i); % according to BTi convention
end
grad.label = grad.label(:);
grad.tra = sparse(eye(size(grad.pnt,1)));

