classdef spm_provenance < handle
% Provenance using PROV Data Model
%   http://www.w3.org/TR/prov-dm/
%
% p = spm_provenance;
% p.set_default_namespace(uri)
% p.add_namespace(prefix,uri)
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
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_provenance.m 5743 2013-11-13 18:19:47Z guillaume $

% Todo:
% - handle @xxx in attributes' literal
% - handle attributes values such as "['mm', 'mm', 'mm']"
% - escape qualified names

%-Properties
%==========================================================================
properties (SetAccess='private', GetAccess='private')
    namespace = {};
    stack = {};
end

%-Constructor
%==========================================================================
methods
    function obj = spm_provenance
    end
end

%-Public methods
%==========================================================================
methods (Access='public')
    
    %-Namespaces
    %----------------------------------------------------------------------
    function set_default_namespace(obj,uri)
        if ~isempty(obj.namespace)
            obj.namespace(cellfun(@isnumeric,obj.namespace(:,1)),:) = [];
        end
        obj.namespace = [{NaN,uri}; obj.namespace];
    end
    
    function ns = add_namespace(obj,prefix,uri)
        if ismember(prefix,{'prov','xsd'})
            error('Namespace MUST NOT declare prefix prov and xsd.');
        end
        if ~isempty(obj.namespace) && ismember(prefix,obj.namespace(:,1))
            n = find(ismember(obj.namespace(:,1),prefix));
        else
            n = size(obj.namespace,1) + 1;
        end
        obj.namespace{n,1} = prefix;
        obj.namespace{n,2} = uri;
        
        if nargout
            ns = @(x) [prefix ':' x];
        end
    end
    
    %-Components
    %----------------------------------------------------------------------
    function entity(obj,id,attributes)
        if nargin < 3 || isempty(attributes), attributes = {}; end
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'entity',id,attrstr(attributes)};
    end
     
    function activity(obj,id,varargin)
        if numel(varargin) && iscell(varargin{end})
            attributes = varargin{end};
            varargin(end) = [];
        else
            attributes = {};
        end
        if numel(varargin) < 1, startTime = '-'; else startTime = timestr(varargin{1}); end
        if numel(varargin) < 2, endTime = '-'; else endTime = timestr(varargin{2}); end
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'activity',id,startTime,endTime,attrstr(attributes)};
    end
    
    function agent(obj,id,attributes)
        if nargin < 3 || isempty(attributes), attributes = {}; end
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'agent',id,attrstr(attributes)};
    end
    
    %-Relations
    %----------------------------------------------------------------------
    %function wasGeneratedBy(obj,id,entity,activity,time,attributes)
    function wasGeneratedBy(obj,varargin)
        [id,entity,arg,attributes] = parseArg(varargin{:});
        if numel(arg) < 1, activity = '-'; else activity = arg{1}; end
        if numel(arg) < 2, time = '-'; else time = timestr(arg{2}); end
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasGeneratedBy',id,entity,activity,time,attributes};
    end
    
    %function used(obj,id,activity,entity,time,attributes)
    function used(obj,varargin)
        [id,activity,arg,attributes] = parseArg(varargin{:});
        if numel(arg) < 1, entity = '-'; else entity = arg{1}; end
        if numel(arg) < 2, time = '-'; else time = timestr(arg{2}); end
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'used',id,entity,activity,time,attributes};
    end
    
    %function wasInformedBy(obj,id,informed,informant,attributes)
    function wasInformedBy(obj,varargin)
        [id,informed,arg,attributes] = parseArg(varargin{:});
        informant = arg{1};
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasInformedBy',id,informed,informant,attributes};
    end
    
    %function wasStartedBy(obj,id,activity,trigger,starter,time,attributes)
    function wasStartedBy(obj,varargin)
        [id,activity,arg,attributes] = parseArg(varargin{:});
        if numel(arg) < 1, trigger = '-'; else trigger = arg{1}; end
        if numel(arg) < 2, starter = '-'; else starter = arg{2}; end
        if numel(arg) < 3, time = '-'; else time = timestr(arg{3}); end
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasStartedBy',id,activity,trigger,starter,time,attributes};
    end
    
    %function wasEndedBy(obj,id,activity,trigger,ender,time,attributes)
    function wasEndedBy(obj,varargin)
        [id,activity,arg,attributes] = parseArg(varargin{:});
        if numel(arg) < 1, trigger = '-'; else trigger = arg{1}; end
        if numel(arg) < 2, ender = '-'; else ender = arg{2}; end
        if numel(arg) < 3, time = '-'; else time = timestr(arg{3}); end
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasEndedBy',id,activity,trigger,ender,time,attributes};
    end
    
    %function wasInvalidatedBy(obj,id,entity,activity,time,attributes)
    function wasInvalidatedBy(obj,varargin)
        [id,entity,arg,attributes] = parseArg(varargin{:});
        if numel(arg) < 1, activity = '-'; else activity = arg{1}; end
        if numel(arg) < 2, time = '-'; else time = timestr(arg{2}); end
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasInvalidatedBy',id,entity,activity,time,attributes};
    end
    
    %function wasDerivedFrom(obj,id,generatedEntity,usedEntity,activity,generation,usage,attributes)
    function wasDerivedFrom(obj,varargin)
        [id,generatedEntity,arg,attributes] = parseArg(varargin{:});
        usedEntity = arg{1};
        if numel(arg) < 2, activity = '-'; else activity = arg{2}; end
        if numel(arg) < 3, generation = '-'; else generation = arg{3}; end
        if numel(arg) < 4, usage = '-'; else usage = arg{4}; end
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasDerivedFrom',id,generatedEntity,usedEntity,activity,generation,usage,attributes};
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
        [id,entity,arg,attributes] = parseArg(varargin{:});
        agent = arg{1};
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasAttributedTo',id,entity,agent,attributes};
    end
    
    %function wasAssociatedWith(obj,id,activity,agent,plan,attributes)
    function wasAssociatedWith(obj,varargin)
        [id,activity,arg,attributes] = parseArg(varargin{:});
        if numel(arg) < 1, agent = '-'; else agent = arg{1}; end
        if numel(arg) < 2, plan = '-'; else plan = arg{2}; end
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasAttributedTo',id,activity,agent,plan,attributes};
    end
    
    %function actedOnBehalfOf(obj,id,delegate,responsible,activity,attributes)
    function actedOnBehalfOf(obj,varargin)
        [id,delegate,arg,attributes] = parseArg(varargin{:});
        responsible = arg{1};
        if numel(arg) < 2, activity = '-'; else activity = arg{2}; end
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'actedOnBehalfOf',id,delegate,responsible,activity,attributes};
    end
    
    %function wasInfluencedBy(obj,id,influencee,influencer,attributes)
    function wasInfluencedBy(obj,varargin)
        [id,influencee,arg,attributes] = parseArg(varargin{:});
        influencer = arg{1};
        
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'wasInfluencedBy',id,influencee,influencer,attributes};
    end
    
    %function alternateOf(obj,alternate1,alternate2)
    function alternateOf(obj,alternate1,alternate2)
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'alternateOf','',alternate1,alternate2,{}};
    end
    
    %function specializationOf(obj,specificEntity,generalEntity)
    function specializationOf(obj,specificEntity,generalEntity)
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'specializationOf','',specificEntity,generalEntity,{}};
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
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'hadMember','',collection,entity,{}};
    end
    
    %function bundle(obj,id,p)
    function bundle(obj,id,p)
        n = numel(obj.stack) + 1;
        obj.stack{n} = {'bundle',id,p};
    end
    
    %-Serialization
    %----------------------------------------------------------------------
    function serialize(obj,fmt)
        if nargin < 2, fmt = 'provn'; end
        switch lower(fmt)
            case 'provn'
                %-PROV-N: the Provenance Notation
                % http://www.w3.org/TR/prov-n/
                fprintf('document\n');
                serialize_provn(obj);
                fprintf('endDocument\n');
            case 'json'
                %-PROV-JSON
                % http://www.w3.org/Submission/2013/SUBM-prov-json-20130424/
                fprintf('{\n');
                serialize_json(obj);
                fprintf('}\n');
            otherwise
                error('Unknown format "%s".',fmt);
        end
    end
end

%-Private methods
%==========================================================================
methods (Access='private')
    function serialize_provn(obj,step)
        if nargin < 2, step = 1; end
        o = blanks(2*step);
        %-Namespace
        for i=1:size(obj.namespace,1)
            if i==1 && isnumeric(obj.namespace{i,1}) && isnan(obj.namespace{i,1})
                fprintf([o 'default <%s>\n'],obj.namespace{i,2});
            else
                fprintf([o 'prefix %s <%s>\n'],...
                    obj.namespace{i,1},obj.namespace{i,2});
            end
        end
        fprintf('\n');
        for i=1:numel(obj.stack)
            %-Components
            if ismember(obj.stack{i}{1},{'entity','agent'})
                fprintf([o '%s(%s'],obj.stack{i}{1:2});
            elseif ismember(obj.stack{i}{1},{'activity'})
                fprintf([o '%s(%s, %s, %s'],obj.stack{i}{1:4});
            elseif ismember(obj.stack{i}{1},{'bundle'})
                fprintf([o 'bundle %s\n'],obj.stack{i}{2});
                serialize_provn(obj.stack{i}{3},2);
                fprintf([o 'endBundle\n']);
            else
                fprintf([o '%s('],obj.stack{i}{1});
                if ~isempty(obj.stack{i}{2})
                    fprintf('%s; ',obj.stack{i}{2});
                end
                k = find(cellfun(@(x) isequal(x,'-'),obj.stack{i}(3:end-1)));
                if isempty(k)
                    k = numel(obj.stack{i});
                else
                    k = min(k) + 2; % remove optional '-'
                end
                for j=3:k-1
                    fprintf('%s',obj.stack{i}{j});
                    if j~=k-1, fprintf(', '); end
                end
            end
            %-Attributes
            if ~ismember(obj.stack{i}{1},{'alternateOf','specializationOf','hadMember','bundle'})
                attr = obj.stack{i}{end};
                if ~isempty(attr)
                    fprintf([',\n' o o '[']);
                    for j=1:2:numel(attr)
                        attribute = attr{j};
                        literal = attr{j+1};
                        if iscell(literal)
                            literal = sprintf('"%s" %%%% %s',literal{:});
                        else
                            if isequal(attribute,'prov:type') || strncmp(literal,'prov:',5)
                                % any literal that starts with xxxx: ?
                                % or impose it with {literal}
                                s = '''';
                            else
                                s = '"';
                            end
                            literal = sprintf([s '%s' s],literal);
                        end
                        fprintf('%s = %s',attribute,literal);
                        if j~=numel(attr)-1, fprintf([',\n' o o]); end
                    end
                    fprintf(']');
                end
            end
            if ~ismember(obj.stack{i}{1},{'bundle'}), fprintf(')\n'); end
        end
    end
    
    function serialize_json(obj,step)
        if nargin < 2, step = 1; end
        o = blanks(2*step);
        %-Namespace
        fprintf([o '"prefix": {\n']);
        for i=1:size(obj.namespace,1)
            if i==1 && isnumeric(obj.namespace{i,1}) && isnan(obj.namespace{i,1})
                fprintf([o o '"default": "%s"'],obj.namespace{i,2});
            else
                fprintf([o o '"%s": "%s"'],obj.namespace{i,1:2});
            end
            if i~=size(obj.namespace,1), fprintf(','); end
            fprintf('\n');
        end
        fprintf([o '}']);
        %-Expressions
        s = sortprov(obj);
        for i=1:numel(s)
            if ~isempty(s(i).idx) && ~isequal(s(i).expr,'bundle')
                fprintf(',\n');
                fprintf([o '"%s": {\n'],s(i).expr);
                for j=s(i).idx
                    id = obj.stack{j}{2};
                    if isempty(id)
                        id = ['_:' s(i).short int2str(j)]; % change counter to start from 1
                    end
                    fprintf([o o '"%s": {\n'],id);
                    l = find(cellfun(@(x) ~isequal(x,'-'),obj.stack{j}(3:end-1)));
                    attr = obj.stack{j}{end};
                    for k=1:numel(l)
                        fprintf([o o o '"prov:%s": "%s"'],s(i).props{k},obj.stack{j}{k+2});
                        if k~=numel(l) || ~isempty(attr), fprintf(','); end
                        fprintf('\n');
                    end
                    for k=1:2:numel(attr)
                        attribute = attr{k};
                        literal = attr{k+1};
                        datatype = 'xsd:string';
                        if iscell(literal)
                            datatype = literal{2};
                            literal = literal{1};
                        else
                            if isequal(attribute,'prov:type') || strncmp(literal,'prov:',5)
                                datatype = 'xsd:QName';
                            end
                        end
                        fprintf([o o o '"%s": {\n'],attribute);
                        fprintf([o o o o '"$": "%s",\n'],literal);
                        fprintf([o o o o '"type": "%s"\n'],datatype);
                        fprintf([o o o '}']);
                        if k~=numel(attr)-1, fprintf(','); end
                        fprintf('\n');
                    end
                    fprintf([o o '}']);
                    if j~=s(i).idx(end), fprintf(','); end
                    fprintf('\n');
                end
                fprintf([o '}']);
            end
        end
        %-Bundles
        if ~isempty(s(end).idx) %% assumes bundle is last in the list...
            fprintf(',\n');
            fprintf([o '"bundle": {\n']);
            for i=1:numel(s(end).idx)
                serialize_json(obj.stack{s(end).idx(i)}{3},2);
            end
            fprintf([o '}']);
        end
        fprintf('\n');
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
    
end

end

%-Helper functions
%==========================================================================
function [id,identifier,arg,attributes] = parseArg(varargin)
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
end

function [arg,attributes] = addAttr(vararg,attr)
    if iscell(vararg{end})
        arg = vararg(1:end-1);
        attributes = [vararg{end} attr{:}];
    else
        arg = vararg;
        attributes = attr;
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
            f = strrep(mat2str(f),'  ',', ');
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

function l = list_expressions
% {expression, short_name, {property_name,...}}
l = {
    'entity',            '',      {};...
    'activity',          '',      {'startTime','endTime'};...
    'agent',             '',      {};...
    'wasGeneratedBy',    'wGB',   {'entity','activity','time'};...
    'used',              'u',     {'activity','entity','time'};...
    'wasInformedBy',     'wInfm', {'informed','informant'};...
    'wasStartedBy',      'wSB',   {'activity','trigger','starter','time'};...
    'wasEndedBy',        'wEB',   {'activity','trigger','ender','time'};...
    'wasInvalidatedBy',  'wIB',   {'entity','activity','time'};...
    'wasDerivedFrom',    'wDF',   {'generatedEntity','usedEntity','activity','generation','usage'};...
    'wasAttributedTo',   'wAT',   {'entity','agent'};...
    'wasAssociatedWith', 'wAW',   {'activity','agent','plan'};...
    'actedOnBehalfOf',   'aOBO',  {'delegate','responsible','activity'};...
    'wasInfluencedBy',   'wInf',  {'influencee','influencer'};...
    'alternateOf',       'aO',    {'alternate1','alternate2'};...
    'specializationOf',  'sO',    {'specificEntity','generalEntity'};...
    'hadMember',         'hM',    {'collection','entity'};...
    'bundle',            '',      {};...
    };
end
