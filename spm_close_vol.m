function V = spm_close_vol(V)
% Close image volume
% See: spm_create_vol and spm_write_plane.
%_______________________________________________________________________
% %W% John Ashburner %E%
for i=1:prod(size(V)),
	if isfield(V(i),'fid') & ~isempty(V(i).fid),
		fclose(V(i).fid);
	end;
end;
if isfield(V,'fid'), V = rmfield(V,'fid'); end;
