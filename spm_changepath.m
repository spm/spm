function varargout = spm_changepath(Sf, oldp, newp)
% Recursively replace all occurences of a text pattern in a MATLAB variable.
% FORMAT S = spm_changepath(Sf, oldp, newp)
%
% Sf       - MATLAB variable to fix, or char array of MAT filenames
% oldp     - old string to replace
% newp     - new string replacing oldp
%
% If the pattern is found in a string, any occurence of an invalid file
% separator is replaced to match that of the current system.
%
% If MAT filenames are specified, they will be overwritten with the new
% version. A backup of the initial version is made with a ".old" suffix.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_changepath.m 3277 2009-07-15 11:47:40Z guillaume $


if ~nargin
    Sf = spm_select(Inf,'mat','Select MAT files to fix');
end
if ischar(Sf)
    S = cell(1,size(Sf,1));
    for i=1:size(Sf,1)
        try
            S{i} = load(deblank(Sf(i,:)));
        catch
            error(sprintf('Cannot load %s.',deblank(Sf(i,:))));
        end
    end
else
    S = {Sf};
end

if nargin <= 1
    oldp = spm_input('Old pattern','+1','s');
end
if nargin <= 2
    newp = spm_input('New pattern','+1','s');
end

for i=1:numel(S)
    S{i} = changepath(S{i},oldp,newp);
    if ischar(Sf)
        f   = deblank(Sf(i,:));
        tmp = S{i};
        [sts, msg] = movefile(f,[f '.old']);
        if ~sts, error(msg); end
        save(f ,'-struct','tmp');
    end
end

if numel(S) == 1
    S = S{1};
end

if nargout
    varargout = {S};
else
    varargout = {};
end

%==========================================================================
function S = changepath(S,oldp,newp)

switch class(S)
    case {'double','single','logical','int8','uint8','int16','uint16',...
            'int32','uint32','int64','uint64','function_handle'}
    
    case 'cell'
        for i=1:numel(S)
            S{i} = changepath(S{i},oldp,newp);
        end
        
    case 'struct'
        for i=1:numel(S)
            fn = fieldnames(S);
            for j=1:length(fn)
                S(i).(fn{j}) = changepath(S(i).(fn{j}),oldp,newp);
            end
        end
        
    case 'char'
        f = {'\' '/'};
        if ispc, f = fliplr(f); end
        tmp = cellstr(S);
        for i=1:numel(tmp)
            t = strrep(tmp{i},oldp,newp);
            if ~isequal(tmp{i},t)
                t = strrep(t,f{1},f{2});
                fprintf('%s\n',t);
            end
            tmp{i} = t;
        end
        S = char(tmp);
    
    case 'nifti'
        for i=1:numel(S)
            S(i).dat = changepath(S(i).dat,oldp,newp);
        end
    
    case 'file_array'
        S.fname = changepath(S.fname,oldp,newp);
        
    otherwise
        warning(sprintf('Unknown class %s.',class(S)));
end
