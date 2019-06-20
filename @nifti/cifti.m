function C = cifti(obj)
% Extract CIFTI-2 extension from a NIfTI-2 file
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging


C = [];
if isfield(obj.hdr,'ext') && ~isempty(obj.hdr.ext) && obj.hdr.ext.ecode == 32
    xml = char(obj.hdr.ext.edata(:)');
else
    return;
end
C.xml = xml;
xml = xmltree(xml);
C.hdr = xml2struct(xml,struct(),root(xml),{});
s = size(obj.dat);
try, s = s(5:7); catch, s = s(5:6); end
C.dat = reshape(obj.dat,s);


%==========================================================================
function s = xml2struct(tree,s,uid,arg)
    switch get(tree,uid,'type')
        case 'element'
            child = children(tree,uid);
            l = {};
            ll = get(tree,child(isfield(tree,child,'name')),'name');
            for i=1:length(child)
                arg2 = arg;
                if isfield(tree,child(i),'name')
                    name = get(tree,child(i),'name');
                    nboccur = sum(ismember(l,name));
                    nboccur2 = sum(ismember(ll,name));
                    l = [l, name];
                    arg2 = [arg2, name];
                    if nboccur || nboccur2 > 1
                        arg2 = [arg2, {{nboccur+1}}];
                    end
                end
                s = xml2struct(tree,s,child(i),arg2);
            end
            if isempty(child)
                s = xmlsetfield(s,arg{:},'');
            end
            %-Storing attributes (when possible)
            attrb = attributes(tree,'get',uid);
            if ~isempty(attrb)
                arg2 = [arg, 'attributes'];
                if ~isstruct(attrb), attrb = [attrb{:}]; end
                try, s = xmlsetfield(s,arg2{:},cell2struct({attrb.val},{attrb.key},2)); end
            end
        case {'chardata','cdata'}
            s = xmlsetfield(s,arg{:},get(tree,uid,'value'));
        case {'comment','pi'}
            % Processing instructions and comments are ignored
    end
    
%==========================================================================
function s = xmlsetfield(s,varargin)
% Same as setfield but using '{}' rather than '()'

subs = varargin(1:end-1);
types = repmat({'{}'},1,numel(subs));
types(cellfun(@ischar,subs)) = {'.'};
s = builtin('subsasgn', s, struct('type',types,'subs',subs), varargin{end});
