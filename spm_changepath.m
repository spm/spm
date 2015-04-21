function varargout = spm_changepath(Sf, oldp, newp)
% Recursively replace all occurences of a text pattern in a MATLAB variable.
% FORMAT S = spm_changepath(Sf, oldp, newp)
%
% Sf       - MATLAB variable to fix, or char array of MAT filenames,
%            or directory name (all found MAT files will be analysed)
% oldp     - old string to replace
% newp     - new string replacing oldp
%
% S        - updated MATLAB variable (only if Sf is one)
%
% If the pattern is found in a string, any occurence of an invalid file
% separator is replaced to match that of the current system.
%
% If MAT filenames are specified, they will be overwritten with the new
% version. A backup of the initial version is made with a ".old" suffix.
%__________________________________________________________________________
% Copyright (C) 2009-2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_changepath.m 6416 2015-04-21 15:34:10Z guillaume $


%-Input arguments
%--------------------------------------------------------------------------
if ~nargin
    Sf = spm_select(Inf,'mat','Select MAT files to fix');
end

if nargin <= 1
    oldp = spm_input('Old pattern','+1','s');
end

if nargin <= 2
    newp = spm_input('New pattern','+1','s');
end

%-Replace pattern in given MAT-files
%--------------------------------------------------------------------------
if ischar(Sf)
    if nargout
        error('Output argument only valid for MATLAB variable input');
    end
    if exist(Sf,'dir')
        Sf = spm_select('FPList',Sf,'^.*\.mat$'); % FPListRec for recursive
    end
    if isempty(Sf), Sf = {}; else Sf = cellstr(Sf); end
    for i=1:numel(Sf)
        f = Sf{i};
        try
            S = load(f);
        catch
            error('Cannot load "%s".',f);
        end
        tmp = changepath(S,oldp,newp);
        if ~isequalwithequalnans(tmp,S)
            fprintf('=> Fixing %s\n',f);
            [sts, msg] = movefile(f,[f '.old']);
            if ~sts, error(msg); end
            save(f ,'-struct','tmp', spm_get_defaults('mat.format'));
        end
    end
else
    varargout = { changepath(Sf,oldp,newp) };
end

%==========================================================================
function S = changepath(S,oldp,newp)

check = false;

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
                if check
                    if spm_existfile(t)
                        sts = ' (OK) ';
                    else
                        sts = ' (NOT FOUND) ';
                    end
                else
                    sts = '';
                end
                fprintf('%s%s\n',t,sts);
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
        
    case 'gifti'
        if isfield(S,'cdata') && isa(S.cdata,'file_array')
            fa = struct(S.cdata);
            fa.fname = changepath(fa.fname,oldp,newp);
            S.cdata = file_array(fa);
        end
        
    otherwise
        warning('Unknown class "%s".',class(S));
end
