function V = spm_close_vol(V)
% Close image volume
% See: spm_create_vol and spm_write_plane.
%_______________________________________________________________________
% John Ashburner $Id$

for i=1:prod(size(V)),
	if isfield(V,'private') & isfield(V(i).private,'fid') & ~isempty(V(i).private.fid),
		if ~isempty(fopen(V(i).private.fid)),
			fclose(V(i).private.fid);
		end;
		V(i).private = rmfield(V(i).private,'fid');
	end;
end;
