function [samef, msg] = spm_vol_check(varargin)
% FORMAT [samef, msg] = spm_vol_check(V1, V2, ...)
% checks spm_vol structs are in same space
%
% V1, V2, etc      - arrays of spm_vol structs
%
% samef            - true if images have same dims, mats
% msg              - cell array containing helpful message if not
%_______________________________________________________________________
% %W% Matthew Brett %E%

[dims mats fnames samef msg] = deal([],[],{},1,{});

if nargin < 1
	return
end

for i = 1:prod(size(varargin)),
	vols   = varargin{i};
	dims   = cat(3,dims,vols(:).dim);
	mats   = cat(3, mats,vols(:).mat);
	fnames = {fnames{:}, vols(:).fname};
end
  
nimgs = size(dims, 3);
if nimgs < 2 
	return
end
labs = {'dimensions', 'orientation & voxel size'};

dimf = any(diff(dims(:,1:3,:),1,3));
matf = any(any(diff(mats,1,3)));
chgf = [dimf(:) matf(:)];
chgi = find(any(chgf, 2));
if ~isempty(chgi)
	samef = 0;
	e1    = chgi(1);
	msg   = {['Images don''t all have the same ' labs{chgf(e1)}],...
		'First difference between image pair:',...
		fnames{e1},...
		fnames{e1+1}};
end
return
