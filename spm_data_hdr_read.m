function V = spm_data_hdr_read(P)
% Get data information from file
% FORMAT V = spm_data_hdr_read(P)
% P        - a char or cell array of filenames
%
% V        - a structure array containing data information
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_data_hdr_read.m 4940 2012-09-20 17:27:54Z guillaume $


if ~nargin
    V = default_hdr_struct;
    return;
end

if isempty(P)
    V = default_hdr_struct;
    if iscell(P), V = {V}; end
    return;
end

switch class(P)
    case 'struct'
        [V(1:numel(P),1)] = deal(default_hdr_struct);
        f = fieldnames(P);
        for i=1:numel(V), for j=1:numel(f)
            V(i).(f{j}) = P(i).(f{j});
        end, end
    
    case 'cell'
        V = cellfun(@spm_data_hdr_read, P, 'UniformOutput',false);
        
    case 'char'
        P = cellstr(P);
        switch lower(spm_file(P{1},'ext'))
            case {'nii','hdr','img'}
                V = spm_vol(char(P));
                
            case 'gii'
                [V(1:numel(P),1)] = deal(default_hdr_struct);
                for i=1:numel(P)
                    V(i).fname    = P{i};
                    V(i).private  = gifti(P{i}); % read header only
                    V(i).dim      = size(V(i).private.cdata);
                    %V(i).dt      = [datatype endianness];
                    %V(i).pinfo   = [scaling offset offset];
                end
                
            otherwise
                error('File "%s" is not of a recognised type.', P{1});
        end
        
    otherwise
        error('Don''t know what to do with input of class "%s".',class(P));
end

%==========================================================================
function V = default_hdr_struct
V = struct(...
    'fname',   '',...
    'dim',     [0 0 0],...
    'dt',      [spm_type('float64') spm_platform('bigend')],...
    'pinfo',   [1 0 0]',...
    'mat',     eye(4),...
    'n',       [1 1],...
    'descrip', '',...
    'private', []);
