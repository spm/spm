classdef spm_provenance < handle
% Provenance using PROV Data Model
%   http://www.w3.org/TR/prov-dm/
%
% p = spm_provenance;
% p.get_default_namespace
% p.set_default_namespace(uri)
% p.add_namespace(prefix,uri)
% p.get_namespace(prefix)
% p.remove_namespace(prefix)
% p.entity(id,attributes)
% p.activity(id,startTime,endTime,attributes)
% p.agent(id,attributes)
% p.wasGeneratedBy(id,entity,activity,time,attributes)
% p.used(id,activity,entity,time,attributes)
% p.wasInformedBy(id,informed,informant,attributes)
% p.wasStartedBy(id,activity,trigger,starter,time,attributes)
% p.wasEndedBy(id,activity,trigger,ender,time,attributes)
% p.wasInvalidatedBy(id,entity,activity,time,attributes)
% p.wasDerivedFrom(id,generatedEntity,usedEntity,activity,generation,usage,attributes)
% p.revision(id,generatedEntity,usedEntity,activity,generation,usage,attributes)
% p.quotation(id,generatedEntity,usedEntity,activity,generation,usage,attributes)
% p.primarySource(id,generatedEntity,usedEntity,activity,generation,usage,attributes)
% p.wasAttributedTo(id,entity,agent,attributes)
% p.wasAssociatedWith(id,activity,agent,plan,attributes)
% p.actedOnBehalfOf(id,delegate,responsible,activity,attributes)
% p.wasInfluencedBy(id,influencee,influencer,attributes)
% p.alternateOf(alternate1,alternate2)
% p.specializationOf(specificEntity,generalEntity)
% p.collection(id,attributes)
% p.emptyCollection(id,attributes)
% p.hadMember(collection,entity)
% p.bundle(id,b)
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_provenance.m 6903 2016-10-12 11:36:41Z guillaume $


%-Properties
%==========================================================================
properties (SetAccess='private', GetAccess='private')
    namespace = struct('prefix','','uri','');
    stack = {};
end

%-Constructor
%==========================================================================
methods
    function obj = spm_provenance
        add_namespace(obj,'prov','http://www.w3.org/ns/prov#');
        add_namespace(obj,'xsd','http://www.w3.org/2001/XMLSchema#');
    end
end

%-Public methods
%==========================================================================
methods (Access='public')
    
    %-Namespaces
    %----------------------------------------------------------------------
    function uri = get_default_namespace(obj)
        uri = obj.namespace(1).uri;
    end
    
    function set_default_namespace(obj,uri)
        obj.namespace(1).uri = uri;
    end
    
    function ns = add_namespace(obj,prefix,uri)
        n = ismember({obj.namespace.prefix},prefix);
        if any(n(1:min(numel(n),3)))
            error('Namespace MUST NOT declare prefixes prov and xsd or be empty.');
        end
        if any(n), n = find(n);
        else       n = numel(obj.namespace) + 1; end
        obj.namespace(n).prefix = prefix;
        obj.namespace(n).uri = uri;
        
        if nargout
            ns = @(local_name) [prefix ':' local_name];
        end
    end
    
    function [uri,ns] = get_namespace(obj,prefix)
        if nargin == 1, uri = obj.namespace; return; end
        n = ismember({obj.namespace.prefix},prefix);
        if ~any(n), uri = '';
        else        uri = obj.namespace(n).uri; end
        if nargout > 1
            ns = @(local_name) [prefix ':' local_name];
        end
    end
    
    function remove_namespace(obj,prefix)
        n = ismember({obj.namespace.prefix},prefix);
        if any(n)
            obj.namespace(n) = [];
        else
            warning('Prefix ''%s'' not found.',prefix);
        end
    end
    
    %-Components
    %----------------------------------------------------------------------
    %function entity(obj,id,attributes)
    function entity(obj,varargin)
        parseArg(obj,'entity',varargin{:});
    end
    
    %function activity(obj,id,startTime,endTime,attributes)
    function activity(obj,varargin)
        parseArg(obj,'activity',varargin{:});
    end
    
    %function agent(obj,id,attributes)
    function agent(obj,varargin)
        parseArg(obj,'agent',varargin{:});
    end
    
    %-Relations
    %----------------------------------------------------------------------
    %function wasGeneratedBy(obj,id,entity,activity,time,attributes)
    function wasGeneratedBy(obj,varargin)
        parseArg(obj,'wasGeneratedBy',varargin{:});
    end
    
    %function used(obj,id,activity,entity,time,attributes)
    function used(obj,varargin)
        parseArg(obj,'used',varargin{:});
    end
    
    %function wasInformedBy(obj,id,informed,informant,attributes)
    function wasInformedBy(obj,varargin)
        parseArg(obj,'wasInformedBy',varargin{:});
    end
    
    %function wasStartedBy(obj,id,activity,trigger,starter,time,attributes)
    function wasStartedBy(obj,varargin)
        parseArg(obj,'wasStartedBy',varargin{:});
    end
    
    %function wasEndedBy(obj,id,activity,trigger,ender,time,attributes)
    function wasEndedBy(obj,varargin)
        parseArg(obj,'wasEndedBy',varargin{:});
    end
    
    %function wasInvalidatedBy(obj,id,entity,activity,time,attributes)
    function wasInvalidatedBy(obj,varargin)
        parseArg(obj,'wasInvalidatedBy',varargin{:});
    end
    
    %function wasDerivedFrom(obj,id,generatedEntity,usedEntity,activity,generation,usage,attributes)
    function wasDerivedFrom(obj,varargin)
        parseArg(obj,'wasDerivedFrom',varargin{:});
    end
    
    %function revision(obj,id,generatedEntity,usedEntity,activity,generation,usage,attributes)
    function revision(obj,varargin)
        attr = {'prov:type','prov:Revision'};
        [arg,attributes] = addAttr(varargin,attr);
        wasDerivedFrom(obj,arg{:},attributes);
    end
    
    %function quotation(obj,id,generatedEntity,usedEntity,activity,generation,usage,attributes)
    function quotation(obj,varargin)
        attr = {'prov:type','prov:Quotation'};
        [arg,attributes] = addAttr(varargin,attr);
        wasDerivedFrom(obj,arg{:},attributes);
    end
    
    %function primarySource(obj,id,generatedEntity,usedEntity,activity,generation,usage,attributes)
    function primarySource(obj,varargin)
        attr = {'prov:type','prov:primarySource'};
        [arg,attributes] = addAttr(varargin,attr);
        wasDerivedFrom(obj,arg{:},attributes);
    end
    
    %function wasAttributedTo(obj,id,entity,agent,attributes)
    function wasAttributedTo(obj,varargin)
        parseArg(obj,'wasAttributedTo',varargin{:});
    end
    
    %function wasAssociatedWith(obj,id,activity,agent,plan,attributes)
    function wasAssociatedWith(obj,varargin)
        parseArg(obj,'wasAssociatedWith',varargin{:});
    end
    
    %function actedOnBehalfOf(obj,id,delegate,responsible,activity,attributes)
    function actedOnBehalfOf(obj,varargin)
        parseArg(obj,'actedOnBehalfOf',varargin{:});
    end
    
    %function wasInfluencedBy(obj,id,influencee,influencer,attributes)
    function wasInfluencedBy(obj,varargin)
        parseArg(obj,'wasInfluencedBy',varargin{:});
    end
    
    %function alternateOf(obj,alternate1,alternate2)
    function alternateOf(obj,alternate1,alternate2)
        addItem(obj,'alternateOf','',alternate1,alternate2,{});
    end
    
    %function specializationOf(obj,specificEntity,generalEntity)
    function specializationOf(obj,specificEntity,generalEntity)
        addItem(obj,'specializationOf','',specificEntity,generalEntity,{});
    end
    
    %function collection(obj,id,attributes)
    function collection(obj,varargin)
        attr = {'prov:type','prov:Collection'};
        [arg,attributes] = addAttr(varargin,attr);
        entity(obj,arg{:},attributes);
    end
    
    %function emptyCollection(obj,id,attributes)
    function emptyCollection(obj,varargin)
        attr = {'prov:type','prov:emptyCollection'};
        [arg,attributes] = addAttr(varargin,attr);
        entity(obj,arg{:},attributes);
    end
    
    %function hadMember(obj,collection,entity)
    function hadMember(obj,collection,entity)
        addItem(obj,'hadMember','',collection,entity,{});
    end
    
    %function bundle(obj,id,p)
    function varargout = bundle(obj,id,p)
        if nargin < 3, p = eval(class(obj)); end
        addItem(obj,'bundle',id,p);
        if nargin < 3, varargout = {p}; end
    end
    
    %-Serialization
    %----------------------------------------------------------------------
    function varargout = serialize(obj,fmt)
        if nargin < 2, fmt = 'provn'; end
        [p,n,e] = fileparts(fmt);
        if ~isempty(e), fmt = e(2:end); end
        
        switch lower(fmt)
            case 'provn'
                %-PROV-N: the Provenance Notation
                % http://www.w3.org/TR/prov-n/
                s = sprintf('document\n');
                s = [s serialize_provn(obj)];
                s = [s sprintf('endDocument\n')];
            case 'json'
                %-PROV-JSON
                % http://www.w3.org/Submission/2013/SUBM-prov-json-20130424/
                s = sprintf('{\n');
                s = [s serialize_json(obj)];
                s = [s sprintf('}\n')];
            case 'jsonld'
                %-PROV-JSONLD
                % http://dl.acm.org/citation.cfm?id=2962062
                [s,b] = serialize_jsonld(obj);
                arr = ~isempty(b);
                while ~isempty(b)
                    s = [s sprintf(',\n')];
                    [s0,b] = serialize_jsonld(b{3},b{2});
                    s = [s s0];
                end
                if arr, s = [sprintf('[\n') s sprintf('\n]\n')]; end
            case 'ttl'
                %-Turtle
                % http://www.w3.org/TR/turtle/
                %warning('Partially implemented.');
                s = serialize_ttl(obj);
            case 'dot'
                %-GraphViz
                % http://www.graphviz.org/
                %warning('Partially implemented.');
                s = sprintf('digraph "PROV" { size="16,12"; rankdir="BT";\n');
                s = [s serialize_dot(obj)];
                s = [s sprintf('}\n')];
                %matlab.internal.strfun.dot2fig(s);
            case {'pdf','svg','png'}
                tmp = tempname;
                dotfile = [tmp '.dot'];
                if isempty(e), outfile = [tmp '.' fmt];
                else outfile = fullfile(p,[n e]); end
                serialize(obj,dotfile);
                dotexe = 'dot';
                system(['"' dotexe '" -T' fmt ' -Gdpi=350 -o "' outfile '" "' dotfile '"']);
                delete(dotfile);
                open(outfile);
                return;
            otherwise
                error('Unknown format "%s".',fmt);
        end
        
        if ~isempty(e)
            filename = fullfile(p,[n e]);
            fid = fopen(filename,'wt');
            if fid == -1, error('Cannot write "%s%".',filename); end
            fprintf(fid,'%s',s);
            fclose(fid);
        else
            varargout = {s};
        end
    end
end

%-Private methods
%==========================================================================
methods (Access='private')
    function [id,identifier,arg,attributes] = parseArg(obj,comp,varargin)
        if isempty(varargin), error('Invalid syntax.'); end
        if isstruct(varargin{1})
            id = varargin{1}.id;
            varargin = varargin(2:end);
        else
            id = '';
        end
        identifier = varargin{1};
        if iscell(varargin{end})
            attributes = attrstr(varargin{end});
            varargin = varargin(1:end-1);
        else
            attributes = {};
        end
        arg = varargin(2:end);
        
        l = list_expressions;
        i = ismember(l(:,1),comp);
        argconv = l{i,4};
        if ~ismember(comp,{'entity','activity','agent'})
            argconv = argconv(2:end);
        end
        if numel(arg) > numel(argconv)
            error('Too many input arguments.');
        end
        for j=1:numel(argconv)
            if numel(arg) < j, arg{j} = '-';
            else               arg{j} = argconv{j}(arg{j}); end
        end
        
        if ismember(comp,{'entity','activity','agent'})
            if ~isempty(id), error('Invalid syntax.'); end
            addItem(obj,comp,identifier,arg{:},attributes);
        else
            addItem(obj,comp,id,identifier,arg{:},attributes);
        end
    end
    
    function addItem(obj,varargin)
        n = numel(obj.stack) + 1;
        obj.stack{n} = varargin;
    end
    
    function str = serialize_provn(obj,step)
        if nargin < 2, step = 1; end
        o = blanks(2*step);
        str = '';
        %-Namespace
        if ~isempty(obj.namespace(1).uri)
            str = [str sprintf([o 'default <%s>\n'],obj.namespace(1).uri)];
        end
        for i=4:numel(obj.namespace)
            str = [str sprintf([o 'prefix %s <%s>\n'],...
                obj.namespace(i).prefix, obj.namespace(i).uri)];
        end
        if ~isempty(obj.namespace(1).uri) || numel(obj.namespace) > 3
            str = [str sprintf('\n')];
        end
        for i=1:numel(obj.stack)
            %-Components
            if ismember(obj.stack{i}{1},{'entity','agent'})
                str = [str sprintf([o '%s(%s'],obj.stack{i}{1:2})];
            elseif ismember(obj.stack{i}{1},{'activity'})
                if isequal(obj.stack{i}{4},'-') && isequal(obj.stack{i}{3},'-')
                    str = [str sprintf([o '%s(%s'],obj.stack{i}{1:2})];
                elseif isequal(obj.stack{i}{4},'-')
                    str = [str sprintf([o '%s(%s, %s'],obj.stack{i}{1:3})];
                else
                    str = [str sprintf([o '%s(%s, %s, %s'],obj.stack{i}{1:4})];
                end
            elseif ismember(obj.stack{i}{1},{'bundle'})
                str = [str sprintf([o 'bundle %s\n'],obj.stack{i}{2})];
                str = [str serialize_provn(obj.stack{i}{3},2)];
                str = [str sprintf([o 'endBundle\n'])];
            else
                str = [str sprintf([o '%s('],obj.stack{i}{1})];
                if ~isempty(obj.stack{i}{2})
                    str = [str sprintf('%s; ',obj.stack{i}{2})];
                end
                k = find(cellfun(@(x) ~isequal(x,'-'),obj.stack{i}(3:end-1)));
                if false %isempty(k)
                    k = 0; % remove optional '-'
                else
                    k = numel(obj.stack{i})-1;
                end
                for j=3:k
                    str = [str sprintf('%s',obj.stack{i}{j})];
                    if j~=k, str = [str sprintf(', ')]; end
                end
            end
            %-Attributes
            if ~ismember(obj.stack{i}{1},{'alternateOf','specializationOf','hadMember','bundle'})
                attr = obj.stack{i}{end};
                if ~isempty(attr)
                    str = [str sprintf([',\n' o o '['])];
                    for j=1:2:numel(attr)
                        attribute = attr{j};
                        literal = attr{j+1};
                        if iscell(literal)
                            literal = sprintf('"%s" %%%% %s',literal{:});
                        else
                            if ~isempty(parseQN(literal,'prefix'))
                                s = '''';
                            else
                                s = '"';
                            end
                            literal = sprintf([s '%s' s],literal);
                        end
                        str = [str sprintf('%s = %s',attribute,literal)];
                        if j~=numel(attr)-1, str = [str sprintf([',\n' o o])]; end
                    end
                    str = [str sprintf(']')];
                end
            end
            if ~ismember(obj.stack{i}{1},{'bundle'}), str = [str sprintf(')\n')]; end
        end
    end
    
    function str = serialize_json(obj,step)
        if nargin < 2, step = 1; end
        o = blanks(2*step);
        str = '';
        %-Namespace
        str = [str sprintf([o '"prefix": {\n'])];
        if ~isempty(obj.namespace(1).uri)
            str = [str sprintf([o o '"default": "%s"'],obj.namespace(1).uri)];
        end
        for i=4:numel(obj.namespace)
            str = [str sprintf([o o '"%s": "%s"'],...
                obj.namespace(i).prefix,obj.namespace(i).uri)];
            if i~=numel(obj.namespace), str = [str sprintf(',')]; end
            str = [str sprintf('\n')];
        end
        str = [str sprintf([o '}'])];
        %-Expressions
        s = sortprov(obj);
        for i=1:numel(s)
            if ~isempty(s(i).idx) && ~isequal(s(i).expr,'bundle')
                str = [str sprintf(',\n')];
                str = [str sprintf([o '"%s": {\n'],s(i).expr)];
                for j=s(i).idx
                    id = obj.stack{j}{2};
                    if isempty(id)
                        id = ['_:' s(i).short int2str(j)]; % change counter to start from 1
                    end
                    str = [str sprintf([o o '"%s": {\n'],id)];
                    l = find(cellfun(@(x) ~isequal(x,'-'),obj.stack{j}(3:end-1)));
                    attr = obj.stack{j}{end};
                    for k=1:numel(l)
                        str = [str sprintf([o o o '"prov:%s": "%s"'],s(i).props{k},obj.stack{j}{k+2})];
                        if k~=numel(l) || ~isempty(attr), str = [str sprintf(',')]; end
                        str = [str sprintf('\n')];
                    end
                    A = reshape(attr,2,[])';
                    while size(A,1) ~= 0
                        attribute = A{1,1};
                        l = 1;
                        for k=2:size(A,1)
                            if strcmp(A{k,1},attribute)
                                l = [l k];
                            end
                        end
                        for k=1:numel(l)
                            literal{k} = A{k,2};
                            datatype{k} = 'xsd:string';
                            if iscell(literal{k})
                                datatype{k} = literal{k}{2};
                                literal{k} = literal{k}{1};
                            else
                                if isequal(attribute,'prov:type') || strncmp(literal{k},'prov:',5)
                                    datatype{k} = 'xsd:QName';
                                end
                            end
                        end
                        str = [str sprintf([o o o '"%s":'],attribute)];
                        if numel(l) == 1
                            str = [str sprintf(' {\n')];
                        else
                            str = [str sprintf(' [\n')];
                        end
                        for k=1:numel(l)
                            if numel(l) ~= 1
                                str = [str sprintf([o o o o '{\n'])];
                            end
                            str = [str sprintf([o o o o '"$": "%s",\n'],literal{k})];
                            str = [str sprintf([o o o o '"type": "%s"\n'],datatype{k})];
                            if numel(l) ~= 1
                                str = [str sprintf([o o o o '}'])];
                                if k~=numel(l), str = [str sprintf(',')]; end
                                str = [str sprintf('\n')];
                            end
                        end
                        if numel(l) == 1
                            str = [str sprintf([o o o '}'])];
                        else
                            str = [str sprintf([o o o ']'])];
                        end
                        A(l,:) = [];
                        if size(A,1) ~= 0, str = [str sprintf(',')]; end
                        str = [str sprintf('\n')];
                    end
                    str = [str sprintf([o o '}'])];
                    if j~=s(i).idx(end), str = [str sprintf(',')]; end
                    str = [str sprintf('\n')];
                end
                str = [str sprintf([o '}'])];
            end
        end
        %-Bundles
        if ~isempty(s(end).idx) %% assumes bundle is last in the list...
            str = [str sprintf(',\n')];
            str = [str sprintf([o '"bundle": {\n'])];
            for i=1:numel(s(end).idx)
                str = [str sprintf([o o '"%s": {\n'],obj.stack{s(end).idx(i)}{2})];
                str = [str serialize_json(obj.stack{s(end).idx(i)}{3},step+1)];
                str = [str sprintf([o o '}'])];
                if i~=numel(s(end).idx), str = [str sprintf(',')]; end
                str = [str sprintf('\n')];
            end
            str = [str sprintf([o '}'])];
        end
        str = [str sprintf('\n')];
    end
    
    function [str,b] = serialize_jsonld(obj,bundle_id,step)
        if nargin < 2, bundle_id = ''; end
        if nargin < 3, step = 1; end
        o = blanks(2*step);
        str = sprintf('{\n');
        b = [];
        provo = prov_o;
        %-Namespace
        str = [str sprintf([o '"@context": [\n'])];
        str = [str sprintf([o o '"https://provenance.ecs.soton.ac.uk/prov.jsonld"'])];
        ns = obj.namespace;
        if ~isempty(bundle_id)
            ns(end+1) = struct('prefix','rdfs','uri','http://www.w3.org/2000/01/rdf-schema#');
        end
        if ~isempty(ns(1).uri) || numel(ns) > 3
            str = [str sprintf(',\n') o o sprintf('{\n')];
        else
            str = [str sprintf('\n')];
        end
        if ~isempty(ns(1).uri)
            str = [str sprintf([o o o '"@base": "%s"'],ns(1).uri)];
            if numel(ns) > 3
                str = [str sprintf(',')];
            end
            str = [str sprintf('\n')];
        end
        for i=4:numel(ns)
            str = [str sprintf([o o o '"%s": "%s"'],...
                ns(i).prefix, ns(i).uri)];
            if i~=numel(ns)
                str = [str sprintf(',')];
            end
            str = [str sprintf('\n')];
        end
        if ~isempty(ns(1).uri) || numel(ns) > 3
            str = [str o o sprintf('}\n')];
        end
        str = [str sprintf([o ']'])];
        %-Identifiant (for Bundle)
        if ~isempty(bundle_id)
            str = [str sprintf(',\n') o sprintf('"@id": "%s"',bundle_id)];
        end
        %-Graph
        str = [str sprintf(',\n') o '"@graph": [' sprintf('\n')];
        for i=1:numel(obj.stack)
            node = {};
            ont = provo(ismember(provo(:,1),obj.stack{i}{1}),:);
            if ismember(obj.stack{i}{1},{'specializationOf','alternateOf','hadMember'})
                node{1,1} = '@id';
                node{1,2} = obj.stack{i}{3};
                node{2,1} = obj.stack{i}{1};
                node{2,2} = obj.stack{i}{4};
            elseif ismember(obj.stack{i}{1},{'bundle'})
                if ~isempty(b)
                    warning('Only a single bundle is allowed.');
                end
                b = obj.stack{i};
                continue;
            else
                attr = obj.stack{i}{end};
                % identifier
                if ~isempty(obj.stack{i}{2})
                    node{end+1,1} = '@id';
                    node{end,2} = obj.stack{i}{2};
                end
                % type
                node{end+1,1} = '@type';
                node{end,2} = ont(2);
                idx = [];
                for j=1:2:numel(attr)
                    if strcmp(attr{j},'prov:type') && ~isempty(parseQN(attr{j+1}))
                        node{end,2}{end+1} = attr{j+1};
                        idx = [idx j j+1];
                    end
                end
                if numel(node{end,2}) == 1, node{end,2} = char(node{end,2}); end
                attr(idx) = [];
                % PROV attributes
                for j=1:numel(ont{3})
                    if ~isequal(obj.stack{i}{j+2},'-')
                        node{end+1,1} = ont{3}{j};
                        node{end,2} = obj.stack{i}{j+2};
                    end
                end
                % additional attributes
                for j=1:2:numel(attr)
                    if isempty(attr{j}), continue; end
                    switch attr{j}
                        case 'prov:location'
                            attr{j} = 'prov:atLocation';
                        case 'prov:label'
                            attr{j} = 'rdfs:label';
                            % remove xsd:string
                            if iscell(attr{j+1})
                                attr{j+1} = attr{j+1}{1};
                            end
                        case 'prov:hadRole'
                    end
                    node{end+1,1} = attr{j};
                    node{end,2} = attr{j+1};
                    if iscell(node{end,2})
                        node{end,2} = {struct(...
                            'type',node{end,2}{2},...
                            'value',node{end,2}{1})};
                    elseif ischar(node{end,2}) &&  ~isempty(parseQN(node{end,2}))
                        node{end,2} = {struct(...
                            'id', node{end,2})};
                    end
                    for k=(j+2):2:numel(attr)
                        if strcmp(node{end,1},attr{k})
                            attr{k} = '';
                            v = attr{k+1};
                            if iscell(v)
                                v = struct('type',v{2},'value',v{1});
                            elseif ischar(v) &&  ~isempty(parseQN(v))
                                v = struct('id',v);
                            end
                            if ischar(node{end,2})
                                node{end,2} = {node{end,2},v};
                            else
                                node{end,2}{end+1} = v;
                            end
                        end
                    end
                end
            end
            % serialize node
            str = [str o o sprintf('{\n')];
            for j=1:size(node,1)
                str = [str o o o sprintf('"%s": ',node{j,1})];
                if ischar(node{j,2})
                    str = [str sprintf('"%s"',node{j,2})];
                else
                    if iscell(node{j,2}) && numel(node{j,2})>1, str = [str '[']; end
                    for k=1:numel(node{j,2})
                        if iscell(node{j,2}) && ischar(node{j,2}{k})
                            str = [str sprintf('"%s"',node{j,2}{k})];
                        else
                            fn = fieldnames(node{j,2}{k});
                            str = [str '{'];
                            for l=1:numel(fn)
                                str = [str sprintf('"@%s": "%s"',...
                                    fn{l},node{j,2}{k}.(fn{l}))];
                                if l<numel(fn)
                                    str = [str sprintf(', ')];
                                end
                            end
                            str = [str '}'];
                        end
                        if k<numel(node{j,2})
                            str = [str sprintf(',')];
                        end
                    end
                    if iscell(node{j,2}) && numel(node{j,2})>1, str = [str ']']; end
                end
                if j<size(node,1)
                    str = [str sprintf(',')];
                end
                str = [str sprintf('\n')];
            end
            str = [str o o sprintf('}')];
            if i<numel(obj.stack) && ~(i==numel(obj.stack)-1 && strcmp(obj.stack{i+1}{1},'bundle'))
                str = [str ','];
            end
            str = [str sprintf('\n')];
        end
        str = [str o sprintf(']')];
        str = [str sprintf('\n')];
        str = [str sprintf('}')];
    end
    
    function str = serialize_ttl(obj,step)
        if nargin < 2, step = 1; end
        o = blanks(2*step);
        str = '';
        %-Namespace
        %if step == 1 % if step > 1, only write user defined namespaces
            ns = obj.namespace;
            if step == 1
                ns(end+1) = struct('prefix','rdfs','uri','http://www.w3.org/2000/01/rdf-schema#');
                if ~isempty(ns(1).uri)
                    str = [str sprintf('@prefix : <%s> .\n',ns(1).uri)];
                end
            end
            for i=2:numel(ns)
                str = [str sprintf('@prefix %s: <%s> .\n',ns(i).prefix,ns(i).uri)];
            end
            if ~isempty(ns(1).uri) || numel(ns) > 3
                str = [str sprintf('\n')];
            end
        %end
        %-Expressions
        % optional entries for activity and relations are not saved
        for i=1:numel(obj.stack)
            if ismember(obj.stack{i}{1},{'entity','activity','agent'})
                str = [str sprintf('%s\n',obj.stack{i}{2})];
                attr = obj.stack{i}{end};
                k = ismember(attr(1:2:end),'prov:type');
                a_type = [{['prov:' obj.stack{i}{1}]} attr{2*find(k)}];
                a_type{1}(6) = upper(a_type{1}(6));
                str = [str sprintf([o 'a'])];
                for j=1:numel(a_type)
                    if ~isempty(parseQN(a_type{j},'prefix'))
                        str = [str sprintf(' %s',a_type{j})];
                    else
                        str = [str sprintf(' "%s"^^xsd:string',a_type{j})];
                    end
                    if j~=numel(a_type), str = [str sprintf(',')]; end
                end
                if ~isempty(attr)
                    str = [str sprintf(' ; \n')];
                    for j=1:2:numel(attr)
                        attribute = attr{j};
                        literal = attr{j+1};
                        if strcmp(attribute,'prov:type'), continue; end
                        if strcmp(attribute,'prov:label')
                            attribute = 'rdfs:label';
                            if iscell(literal), literal = literal{1}; end
                        end
                        if strcmp(attribute,'prov:location')
                            attribute = 'prov:atLocation';
                        end
                        if iscell(literal)
                            %if ~strcmp(literal{2},'xsd:string')
                                literal = sprintf('"%s"^^%s',literal{:});
                            %else
                            %    literal = sprintf('"%s"',literal{1});
                            %end
                        else
                            if ~isempty(parseQN(literal,'prefix'))
                                s = '';
                            else
                                s = '"';
                            end
                            literal = sprintf([s '%s' s],literal);
                        end
                        str = [str sprintf([o '%s %s'],attribute,literal)];
                        if j~=numel(attr)-1, str = [str sprintf(' ;\n')]; end
                    end
                end
                str = [str sprintf(' .\n\n')];
            elseif ismember(obj.stack{i}{1},{'bundle'})
                str = [str serialize_ttl(obj.stack{i}{3},step+1)];
            else
                if strcmp(obj.stack{i}{1},'wasGeneratedBy') && ~strcmp(obj.stack{i}{5},'-')
                    % temporary hack until full PROV-N to PROV-O mapping...
                    tag = ['_:blank' int2str(i)];
                    str = [str sprintf('%s a prov:Generation .\n\n',tag)];
                    str = [str sprintf('%s prov:qualifiedGeneration %s .\n\n',obj.stack{i}{3},tag)];
                    str = [str sprintf('%s prov:atTime "%s"^^xsd:dateTime .\n\n',tag,obj.stack{i}{5})];
                else
                    str = [str sprintf('%s prov:%s %s .\n\n',obj.stack{i}{[3 1 4]})];
                end
            end
        end
    end
    
    function str = serialize_dot(obj,annn)
        s = sortprov(obj);
        str = '';
        expr = {'entity','activity','agent'};
        dot_style.entity   = {'style','filled','shape','ellipse','color','#808080','fillcolor','#FFFC87','sides','4'};
        dot_style.activity = {'style','filled','shape','polygon','color','#0000FF','fillcolor','#9FB1FC','sides','4'};
        dot_style.agent    = {'style','filled','shape','house',  'color','#000000','fillcolor','#FED37F','sides','4'};
        dot_style.default  = {'labeldistance','1.5','rotation','20','labelfontsize','8','labelangle','60.0'};
        dot_style.wasGeneratedBy    = [dot_style.default ,'color','darkgreen','fontcolor','darkgreen'];
        dot_style.used              = [dot_style.default,'color','red4','fontcolor','red'];
        dot_style.wasAttributedTo   = [dot_style.default,'color','#FED37F'];
        dot_style.wasAssociatedWith = [dot_style.default,'color','#FED37F'];
        dot_style.actedOnBehalfOf   = [dot_style.default,'color','#FED37F'];
        dot_style.wasInfluencedBy   = [dot_style.default,'color','grey'];
        dot_style.atLocation        = [dot_style.default,'color','blue','fontcolor','blue'];
        dot_style.annotationLink = {'style','dashed','color','#C0C0C0','arrowhead','none'};
        dot_style.annotation     = {'shape','note','color','gray','fontcolor','black','fontsize','10'};
        strannlab = '<<TABLE cellpadding="0" border="0">%s</TABLE>>';
        strtr = '<TR><TD align="left">%s</TD><TD align="left">%s</TD></TR>';
        if nargin < 2, annn = 0; end
        for i=1:numel(s)
            if ~isempty(s(i).idx)
                if ismember(s(i).expr,expr)
                    idx = ismember(expr,s(i).expr);
                    for j=s(i).idx
                        label = getattr(obj.stack{j}{end},'prov:label');
                        if ~isempty(label)
                            if iscell(label)
                                label = label{1};
                            end
                        else
                            label = parseQN(obj.stack{j}{2},'local');
                        end
                        url = get_url(obj,obj.stack{j}{2});
                        str = [str sprintf('n%s ',get_valid_identifier(url))];
                        str = [str dotlist([dot_style.(s(i).expr),'label',label,'URL',url])];
                        attr = obj.stack{j}{end};
                        if ~isempty(attr)
                            url_ann = sprintf('http://annot/ann%d',annn);
                            attrlist = [];
                            for k=1:2:numel(attr)
                                attribute = attr{k};
                                literal = attr{k+1};
                                if iscell(literal)
                                    literal = literal{1};
                                end
                                attrlist = [attrlist sprintf(strtr,attribute,htmlesc(literal))]; %htmlesc
                            end
                            ann_label = sprintf(strannlab,attrlist);
                            A = ['n' get_valid_identifier(url_ann)];
                            str = [str A ' ' dotlist([dot_style.annotation,'label',ann_label])];
                            B = ['n' get_valid_identifier(url)];
                            str = [str sprintf('%s -> %s ',A,B)];
                            str = [str dotlist(dot_style.annotationLink)];
                            annn = annn + 1;
                        end
                    end
                    % handle prov:location / prov:atLocation
                    if i==1
                        for j=s(i).idx
                            val = getattr(obj.stack{j}{end},'prov:location');
                            if ~isempty(val) && ~iscell(val)
                                A = ['n' get_valid_identifier(get_url(obj,val))];
                                B = ['n' get_valid_identifier(get_url(obj,obj.stack{j}{2}))];
                                str = [str sprintf('%s -> %s ',A,B)];
                                str = [str dotlist([dot_style.atLocation,'label','locationOf'])];
                            end
                        end
                    end
                elseif ~ismember(s(i).expr,{'bundle','collection','emptyCollection','alternateOf','specializationOf'})
                    for j=s(i).idx
                        if isequal(obj.stack{j}{4},'-'), continue; end
                        A = ['n' get_valid_identifier(get_url(obj,obj.stack{j}{3}))];
                        B = ['n' get_valid_identifier(get_url(obj,obj.stack{j}{4}))];
                        str = [str sprintf('%s -> %s ',A,B)];
                        try
                            str = [str dotlist([dot_style.(s(i).expr),'label',s(i).expr])];
                        catch
                            str = [str dotlist([dot_style.default,'label',s(i).expr])];
                        end
                    end
                elseif ismember(s(i).expr,{'bundle'})
                    label = obj.stack{s(i).idx}{2};
                    url = get_url(obj,obj.stack{s(i).idx}{2});
                    str = [str sprintf('subgraph clustern%s {\n',get_valid_identifier(url))];
                    str = [str sprintf('  label="%s";\n',label)];
                    str = [str sprintf('  URL="%s";\n',url)];
                    str = [str serialize_dot(obj.stack{s(i).idx}{3},annn)];
                    str = [str sprintf('}\n')];
                else
                    warning('"%s" not handled yet.',s(i).expr);
                end
            end
        end
    end
    
    function s = sortprov(obj)
        expr = list_expressions;
        l = cellfun(@(x) x{1},obj.stack,'UniformOutput',false);
        for i=1:size(expr,1)
            s(i).expr  = expr{i,1};
            s(i).short = expr{i,2};
            s(i).props = expr{i,3};
            s(i).idx   = find(ismember(l,expr{i,1}));
        end
    end
    
    function url = get_url(obj,id)
        url = [obj.get_namespace(parseQN(id,'prefix')) ...
                parseQN(id,'local')];
    end
    
end

end

%-Helper functions
%==========================================================================
function [arg,attributes] = addAttr(vararg,attr)
    if iscell(vararg{end})
        arg = vararg(1:end-1);
        attributes = [vararg{end} attr{:}];
    else
        arg = vararg;
        attributes = attr;
    end
end

function varargout = parseQN(qn,ret)
    [t,r] = strtok(qn,':');
    if isempty(r), r = t; t = ''; else r(1) = []; end
    if ~isempty(strfind(r,' ')), t = ''; r = qn; end
    if nargin == 1, ret = 'all'; end
    switch lower(ret)
        case 'all'
            varargout = {t,r};
        case 'prefix'
            varargout = {t};
        case 'local'
            varargout = {r};
        otherwise
            error('Syntax error.');
    end
end

function val = getattr(attr,key)
    val = '';
    for i=1:2:numel(attr)
        if strcmp(attr{i},key)
            val = attr{i+1};
            return;
        end
    end
end

function attr = attrstr(attr)
    for i=2:2:numel(attr)
        if isnumeric(attr{i})
            if isinteger(attr{i})
                attr{i} = intstr(attr{i});
            else
                attr{i} = floatstr(attr{i});
            end
        elseif iscell(attr{i}) && iscell(attr{i}{1})
            if numel(attr{i}) == 1
                attr{i} = jsonesc(cell2str(attr{i}{1}));
            else
                attr{i}{1} = jsonesc(cell2str(attr{i}{1}));
            end
        elseif iscell(attr{i})
            if isinteger(attr{i}{1})
                attr{i}{1} = intstr(attr{i}{1});
            else
                attr{i}{1} = floatstr(attr{i}{1});
            end
        end
    end
end

function id = esc(id)
    c = '=''(),-:;[].';
    for i=1:numel(c)
        id = strrep(id,c(i),['\' c(i)]);
    end
end

function str = htmlesc(str)
    %-Escape
    % See http://www.w3.org/TR/html4/charset.html#h-5.3.2
    str = strrep(str,'&','&amp;');
    str = strrep(str,'<','&lt;');
    str = strrep(str,'>','&gt;');
    str = strrep(str,'"','&quot;');
end

function str = jsonesc(str)
    % See http://json.org/
    str = strrep(str,'"','\"');
end

function id = get_valid_identifier(id)
    c = '/:#-.';
    for i=1:numel(c)
        id = strrep(id,c(i),'_');
    end
end

function t = timestr(t)
    if isnumeric(t)
        t = datestr(t,'yyyy-mm-ddTHH:MM:SS');
    end
end

function i = intstr(i)
    if isnumeric(i)
        i = ['[' strrep(int2str(i),'  ',', ') ']'];
    end
end

function f = floatstr(f)
    if isnumeric(f)
        if size(f,1) == 1
            f = strrep(mat2str(f),' ',', ');
        else
            ff = '[';
            for i=1:size(f,1)
                if i~=size(f,1), c=','; else c=''; end
                ff = [ff floatstr(f(i,:)) c];
            end
            ff = [ff ']'];
            f = ff;
        end
    end
end

function s = dotlist(l)
    s = '[';
    for i=1:2:numel(l)
        c = '"';
        if strncmp(l{i+1},'<<',2), c = ''; end
        s = [s sprintf('%s=%c%s%c',l{i},c,l{i+1},c)];
        if i~=numel(l)-1
            s = [s ','];
        end
    end
    s = [s ']' sprintf('\n')];
end

function s = cell2str(s)
    s = jsonesc(s);
    s = ['[' sprintf('"%s", ',s{:}) ']']; s(end-2:end-1) = [];
end

function l = list_expressions
% {expression, short_name, {property_names}, {convert_fcn}}
n = @(x) x;
l = {
    'entity',            '',      {},                                      {};...
    'activity',          '',      {'startTime','endTime'},                 {@timestr,@timestr};...
    'agent',             '',      {},                                      {};...
    'wasGeneratedBy',    'wGB',   {'entity','activity','time'},            {n,n,@timestr};...
    'used',              'u',     {'activity','entity','time'},            {n,n,@timestr};...
    'wasInformedBy',     'wInfm', {'informed','informant'},                {n,n};...
    'wasStartedBy',      'wSB',   {'activity','trigger','starter','time'}, {n,n,n,@timestr};...
    'wasEndedBy',        'wEB',   {'activity','trigger','ender','time'},   {n,n,n,@timestr};...
    'wasInvalidatedBy',  'wIB',   {'entity','activity','time'},            {n,n,@timestr};...
    'wasDerivedFrom',    'wDF',   {'generatedEntity','usedEntity','activity','generation','usage'}, {n,n,n,n,n};...
    'wasAttributedTo',   'wAT',   {'entity','agent'},                      {n,n};...
    'wasAssociatedWith', 'wAW',   {'activity','agent','plan'},             {n,n,n};...
    'actedOnBehalfOf',   'aOBO',  {'delegate','responsible','activity'},   {n,n,n};...
    'wasInfluencedBy',   'wInf',  {'influencee','influencer'},             {n,n};...
    'alternateOf',       'aO',    {'alternate1','alternate2'},             {n,n};...
    'specializationOf',  'sO',    {'specificEntity','generalEntity'},      {n,n};...
    'hadMember',         'hM',    {'collection','entity'},                 {n,n};...
    'bundle',            '',      {},                                      {};...
    };
end

function m = prov_o
m = {
    'entity',            'prov:Entity',           {};                                      ...
    'activity',          'prov:Activity',         {'startedAtTime','endedAtTime'};         ...
    'agent',             'prov:Agent',            {};                                      ...
    'wasGeneratedBy',    'prov:Generation',       {'entity_generated','activity','atTime'};            ...
    'used',              'prov:Usage',            {'activity_using','entity','atTime'};            ...
    'wasInformedBy',     'prov:Communication',    {'informed','activity'};                ...
    'wasStartedBy',      'prov:Start',            {'activity_started','entity','hadActivity','atTime'}; ...
    'wasEndedBy',        'prov:End',              {'activity_ended','entity','hadActivity','atTime'};   ...
    'wasInvalidatedBy',  'prov:Invalidation',     {'entity_invalidated','activity','atTime'};            ...
    'wasDerivedFrom',    'prov:Derivation',       {'entity_derived','entity','hadActivity','hadGeneration','hadUsage'}; ...
    'wasAttributedTo',   'prov:Attribution',      {'entity_attributed','agent'};                      ...
    'wasAssociatedWith', 'prov:Association',      {'activity_associated','agent','hadPlan'};             ...
    'actedOnBehalfOf',   'prov:Delegation',       {'delegate','agent','hadActivity'};   ...
    'wasInfluencedBy',   'prov:Influence',        {'influencee','influencer'};             ...
    'alternateOf',       'prov:alternateOf',      {'alternate1','alternate2'};             ... % no @type
    'specializationOf',  'prov:specializationOf', {'specificEntity','generalEntity'};      ... % no @type
    'hadMember',         'prov:hadMember',        {'collection','entity'};                 ... % no @type
    'bundle',            '',                      {};                                      ...
    };
end
