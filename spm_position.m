function spm_position(h,u,v)
% resets graphics attribute 'Position'
% FORMAT spm_position(h,u,v)
% h    - object handle
% u    - u scaling factor
% v    - [x y] translation {units = normalized}
%____________________________________________________________________________
%
% spm_position is used to reposition or resize a graphics object
% with the appropriate attribute (usually an axis object)
%
%__________________________________________________________________________
% %W% %E%

% reset 'Position'
%----------------------------------------------------------------------------
set(h,'Units','normalized');
P  = get(h,'Position');
set(h,'Position',[(P(1:2) + v - P(3:4).*(u - 1)/2) P(3:4).*u]);
