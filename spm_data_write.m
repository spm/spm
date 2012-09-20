function V = spm_data_write(V,Y,varargin)
% Write data to disk [V(I) = Y]
% FORMAT V = spm_data_write(V,Y)
% V        - a structure array (see spm_data_hdr_read)
% Y        - an array of data values
%
% FORMAT V = spm_data_write(V,Y,I)
% V        - a structure array (see spm_data_hdr_read)
% Y        - an array of data values
% I        - linear index to data values
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_data_write.m 4940 2012-09-20 17:27:54Z guillaume $


switch class(V(1).private)
    case 'nifti'
        if isempty(varargin)
            V = spm_write_vol(V,Y);
            %S = substruct('()',repmat({':'},1,numel(V.private.dat.dim)));
            %V.private.dat = subsasgn(V.private.dat,S,Y);
        else
            if numel(varargin) == 1
                try
                    V.private.dat(varargin{1}) = reshape(Y,size(varargin{1}));
                catch
                    V.private.dat(varargin{1}) = reshape(Y,size(varargin{1}))';
                end
            else
                error('not implemented yet');
            end
        end
    case 'gifti'
        if isempty(varargin)
            D = V.private.cdata;
            D = subsasgn(D,substruct('()',{':'}),Y);
            %V.private.cdata = D;
        else
            try
                %V.private.cdata(varargin{1}) = reshape(Y,size(varargin{1}));
                V.private.private.data{1}.data(varargin{1}) = reshape(Y,size(varargin{1}));
            catch
                %V.private.cdata(varargin{1}) = reshape(Y,size(varargin{1}))';
                V.private.private.data{1}.data(varargin{1}) = reshape(Y,size(varargin{1}))';
            end
        end
    otherwise
        error('Unknown data type.');
end
