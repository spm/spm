function V = spm_data_hdr_write(V)
% Write data information to disk
% FORMAT V = spm_data_hdr_write(V)
% V        - a structure array (see spm_data_hdr_read)
%__________________________________________________________________________
% Copyright (C) 2012-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_data_hdr_write.m 6752 2016-03-24 16:17:25Z guillaume $


switch lower(spm_file(V(1).fname,'ext'))
    case {'nii','hdr','img'}
        if ~isfield(V(1),'dat'), V = spm_create_vol(V); end
        
    case 'gii'
        for i=1:numel(V)
            V(i).private = gifti(struct('cdata',zeros(V(i).dim)));
            V(i).private.private.data{1}.attributes.DataType = ...
                ['NIFTI_TYPE_' upper(spm_type(V(i).dt(1)))];
            if isfield(V(i),'SurfaceID')
                V(i).private.private.metadata(1).name = 'SurfaceID';
                V(i).private.private.metadata(1).value = V(i).SurfaceID;
            end
            % see also endianness, scale/offset/offset_byte, metadata
            save(V(i).private, V(i).fname, 'ExternalFileBinary');
            V(i).private = gifti(V(i).fname);
            %V(i).private.dat = file_array(spm_file(V(i).fname,'ext','dat'), ...
            %    V(i).dim, spm_type(V(i).dt(1)));
        end
        if isfield(V,'SurfaceID'), V = rmfield(V,'SurfaceID'); end
        
    otherwise
        error('File "%s" is not of a recognised type.', V(1).fname);
end
