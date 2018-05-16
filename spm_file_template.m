classdef spm_file_template
% Text file template engine
%
% Example:
% >> tpl = spm_file_template;
% >> tpl = file(tpl,'myfile','template.txt');
% >> tpl = var(tpl,'TITLE',spm('Ver'));
% >> tpl = var(tpl,'DATE',date);
% >> tpl = parse(tpl,'OUT','myfile');
% >> get(tpl,'OUT')
%__________________________________________________________________________
%
% If using a MATLAB version older than R2008a (7.6), run the following:
% >> D = fullfile(spm('Dir'),'@spm_file_template'); mkdir(D);
% >> movefile(fullfile(spm('Dir'),'spm_file_template.m'),D)
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_file_template.m 7312 2018-05-16 11:22:53Z guillaume $


%-Properties
%==========================================================================
properties
    root;
    unknowns = 'remove';
end

properties (SetAccess='private', GetAccess='private')
    filekeys = {};
    filenames = {};
    
    varkeys = {};
    varvals = {};
end

%-Constructor
%==========================================================================
methods
    function obj = spm_file_template(root,unknowns)
        obj.root = pwd;
        if ~nargin, return; end
        obj.root = root;
        if nargin == 1, return; end
        obj.unknowns = unknowns;
    end
end

%-Properties set and get methods
%==========================================================================
methods
    function obj = set.root(obj,root)
        if ~exist(root,'dir')
            error('"%s" not found.',root);
        end
        obj.root = root;
    end
    
    function obj = set.unknowns(obj,unknowns)
        if ~ismember(unknowns,{'comment','keep','remove'})
            error('Unknowns: ''remove'', ''comment'' or ''keep''.');
        end
        obj.unknowns = unknowns;
    end
    
    function str = get(obj,key)
        str = handleUnknowns(getvar(obj,key),obj.unknowns);
    end
end

%-Public methods
%==========================================================================
methods (Access='public')
    function obj = file(obj,key,name)
        key  = cellstr(key);
        name = cellstr(name);
        if filesep == '/',  mch = '^/'; else mch = '^.:\\'; end
        for i=1:numel(key)
            j = find(ismember(obj.filekeys,key{i}));
            if isempty(j), j = numel(obj.filekeys) + 1; end
            obj.filekeys{j} = key{i};
            obj.filenames{j} = name{i};
            if isempty(regexp(obj.filenames{j},mch,'once'))
                obj.filenames{j} = fullfile(obj.root,obj.filenames{j});
            end
        end
    end
    
    function obj = block(obj,F,B,V)
        obj = loadtpl(obj,F);
        if nargin == 3, V = B; end
        str = getvar(obj,F);
        blk = '';
        strbegin = ['<!-- BEGIN ' B ' -->'];
        strend   = ['<!-- END ' B ' -->'];
        indbegin = strfind(str,strbegin);
        indend   = strfind(str,strend);
        if ~isempty(indbegin) && ~isempty(indend)
            blk = str(indbegin+length(strbegin)+1:indend-1);
            str = [str(1:indbegin-1) ...
                '{' V '}' ...
                str(indend+length(strend)+1:end)];
        end
        obj = var(obj,B,blk);
        obj = var(obj,F,str);
    end
    
    function obj = var(obj,key,val,varargin)
        key = cellstr(key);
        if nargin > 2
            val = cellstr(val);
        else
            val = repmat({''}, size(key));
        end
        key = [key varargin{1:2:end}];
        val = [val varargin{2:2:end}];
        for i=1:numel(key)
            j = find(ismember(obj.varkeys,key{i}));
            if isempty(j), j = numel(obj.varkeys) + 1; end
            obj.varkeys{j} = key{i};
            obj.varvals{j} = val{i};
        end
    end
    
    function [obj, str] = parse(obj,target,handle,append)
        if nargin == 3, append = false; end
        if iscellstr(handle)
            for i=1:numel(handle)
                [obj, str] = subst(obj,handle{i});
                obj = var(obj,target,str);
            end
        elseif ischar(handle)
            [obj, str] = subst(obj,handle);
            if append
                obj = var(obj,target,[getvar(obj,target) str]);
            else
                obj = var(obj,target,str);
            end
        end
    end
end

%-Private methods
%==========================================================================
methods (Access='private')
    function val = getvar(obj,key)
        if nargin == 1, val = obj.varvals; return; end
        key = cellstr(key);
        val = cell(1,numel(key));
        for i=1:numel(key)
            j = find(ismember(obj.varkeys,key{i}));
            if isempty(j)
                val{i} = '';
            else
                val{i} = obj.varvals{j};
            end
        end
        if numel(val) == 1, val = char(val); end
    end
    
    function obj = loadtpl(obj,handle)
        if ~isempty(getvar(obj,handle))
            return;
        end
        ind = find(ismember(obj.filekeys,handle));
        if isempty(ind)
            error('Template handle "%s" not found.',handle);
        end
        filename = obj.filenames{ind};
        fid = fopen(filename,'rt');
        if fid == -1
            error('Cannot open template file "%s".',filename);
        end
        obj = var(obj,handle,fscanf(fid,'%c'));
        fclose(fid);
    end
    
    function [obj, str] = subst(obj,handle)
        obj = loadtpl(obj,handle);
        str = getvar(obj,handle);
        for i=1:numel(obj.varkeys)
            str = strrep(str, ['{' obj.varkeys{i} '}'],obj.varvals{i});
        end
    end
    
end

end

%-Helper functions
%==========================================================================
function str = handleUnknowns(str,unknowns)
switch lower(unknowns)
    case 'keep'
        %- do nothing
    case 'remove'
        str = regexprep(str,'{[^ \t\r\n}]+}','');
    case 'comment'
        str = regexprep(str,'{[^ \t\r\n}]+}',...
            '<!-- Template variable undefined -->');
end

end
