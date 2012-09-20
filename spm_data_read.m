function Y = spm_data_read(V,varargin)
% Read data from disk [Y = V(I)]
% FORMAT Y = spm_data_read(V)
% V        - a structure array (see spm_data_hdr_read)
% Y        - an array of data values; the last dimension indexes numel(V)
%
% FORMAT Y = spm_data_read(V,'slice',S)
% V        - a structure array of image volumes (see spm_data_hdr_read)
% S        - an array of slice indices
% Y        - an array of data values with dimensions (x,y,s,v)
%
% FORMAT Y = spm_data_read(V,'xyz',XYZ)
% V        - a structure array (see spm_data_hdr_read)
% XYZ      - a [n x m] array of m coordinates {voxel (n=3 or 4)/vertex (n=1)}
% Y        - an array of data values with dimensions (v,m)
%
% FORMAT Y = spm_data_read(V,I1,I2,...)
% V        - a structure array (see spm_data_hdr_read)
% I1,I2,...- subscript arrays
% Y        - an array of data values with dimensions (v,m)
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_data_read.m 4940 2012-09-20 17:27:54Z guillaume $


if ~isstruct(V)
    V = spm_data_hdr_read(V);
end

switch class(V(1).private)
    case 'nifti'
        if isempty(varargin)
            % Y = V.private.dat(); % if numel(V)==1, is faster
            Y = spm_read_vols(V);
        elseif ischar(varargin{1}) && ~isequal(varargin{1},':')
            switch lower(varargin{1})
                case 'slice'
                    for i=1:numel(V), for p=1:numel(varargin{2})
                        Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 varargin{2}(p)]),V(i).dim(1:2),0);
                    end, end
                    if numel(V)==1, Y=Y(:,:,:,1); end
                case 'xyz'
                    Y = spm_get_data(V,varargin{2});
                otherwise
                    error('Unknown input option.');
            end
        else
            n = get_ndata(V(1).dim,varargin{:});
            Y = zeros(numel(V),prod(n));
            for i=1:numel(V)
                Y(i,:) = reshape(V(i).private.dat(varargin{:}),1,[]);
            end
        end
        
    case 'gifti'
        indices = varargin;
        if isempty(indices)
            indices = repmat({':'},1,ndims(V));
        end
        n = get_ndata(V(1).dim,indices{:});
        Y = zeros(numel(V),prod(n));
        for i=1:numel(V)
            Y(i,:) = reshape(V(i).private.cdata(indices{:}),1,[]);
        end
        if isempty(varargin), Y = Y'; end % to be coherent with spm_read_vols
        
    otherwise
        error('Unknown data type.');
end

%==========================================================================
function n = get_ndata(dim,varargin)
n = zeros(1,numel(varargin));
for i=1:numel(varargin)
    if isequal(varargin{i},':')
        if i==numel(varargin)
            n(i) = dim(i); %prod(dim(i:end));
        else
            n(i) = dim(i);
        end
    else
        n(i) = numel(varargin{i});
    end
end
