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
% Copyright (C) 2012-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_data_write.m 5916 2014-03-13 13:15:02Z guillaume $


if isfield(V,'private')
    cl = class(V.private);
elseif isfield(V,'dat')
    cl = 'struct';
else
    error('Unkwown data type.');
end

switch cl
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
    case 'struct'
        if isempty(varargin)
            V.dat = Y;
        else
            if numel(varargin) == 1
                V.dat(varargin{1}) = reshape(Y,size(varargin{1}));
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
