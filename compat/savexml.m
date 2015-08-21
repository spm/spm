function savexml(filename, varargin)
%SAVEXML Save workspace variables to disk in XML.
%  SAVEXML FILENAME saves all workspace variables to the XML-file
%  named FILENAME.xml. The data may be retrieved with LOADXML. if
%  FILENAME has no extension, .xml is assumed.
%  
%  SAVE, by itself, creates the XML-file named 'matlab.xml'. It is
%  an error if 'matlab.xml' is not writable.
%
%  SAVE FILENAME X saves only X.
%  SAVE FILENAME X Y Z saves X, Y, and Z. The wildcard '*' can be 
%  used to save only those variables that match a pattern.
%
%  SAVE ... -APPEND adds the variables to an existing file.
%
%  Use the functional form of SAVE, such as SAVE(filename','var1','var2'),
%  when the filename or variable names are stored in strings.
%
%  See also SAVE, MAT2XML, XMLTREE.

%  Copyright 2003 Guillaume Flandin. 
%  $Revision: 6530 $  $Date: 2003/07/10 13:50 $

%  $Id: savexml.m 6530 2015-08-21 14:43:52Z guillaume $

if nargin == 0
    filename = 'matlab.xml';
    fprintf('\nSaving to: %s\n\n',filename);
else
    if ~ischar(filename)
        error('[SAVEXML] Argument must contain a string.');
    end
    [pathstr,name,ext] = fileparts(filename);
    if isempty(ext)
        filename = [filename '.xml'];
    end
end
if nargin <= 1, varargin = {'*'}; end

if nargout > 0
    error('[SAVEXML] Too many output arguments.');
end

if strcmpi(varargin{end},'-append')
    if length(varargin) > 1
        varargin = varargin(1:end-1);
    else
        varargin = {'*'};
    end
    if exist(filename,'file')
        % TODO % No need to parse the whole tree ? detect duplicate variables ?
        t = xmltree(filename);
    else
        error(sprintf(...
        '[SAVEXML] Unable to write file %s: file does not exist.',filename));
    end
else
    t = xmltree('<matfile/>');
end

for i=1:length(varargin)
    v = evalin('caller',['whos(''' varargin{i} ''')']);
    if isempty(v)
        error(['[SAVEXML] Variable ''' varargin{i} ''' not found.']);
    end
    for j=1:length(v)
        [t, uid] = add(t,root(t),'element',v(j).name);
        t = attributes(t,'add',uid,'type',v(j).class);
        t = attributes(t,'add',uid,'size',xml_num2str(v(j).size));
        t = xml_var2xml(t,evalin('caller',v(j).name),uid);
    end
end

save(t,filename);

%=======================================================================
function t = xml_var2xml(t,v,uid)

    switch class(v)
        case {'double','single','logical'}
            if ~issparse(v)
                t = add(t,uid,'chardata',xml_num2str(v));
            else % logical
                [i,j,s] = find(v);
                [t, uid2] = add(t,uid,'element','row');
                t = attributes(t,'add',uid2,'size',xml_num2str(size(i)));
                t = add(t,uid2,'chardata',xml_num2str(i));
                [t, uid2] = add(t,uid,'element','col');
                t = attributes(t,'add',uid2,'size',xml_num2str(size(j)));
                t = add(t,uid2,'chardata',xml_num2str(j));
                [t, uid2] = add(t,uid,'element','val');
                t = attributes(t,'add',uid2,'size',xml_num2str(size(s)));
                t = add(t,uid2,'chardata',xml_num2str(s));
            end
        case 'struct'
            names = fieldnames(v);
            for j=1:prod(size(v))
                for i=1:length(names)
                    [t, uid2] = add(t,uid,'element',names{i});
                    t = attributes(t,'add',uid2,'index',num2str(j));
                    t = attributes(t,'add',uid2,'type',...
                        class(getfield(v(j),names{i})));
                    t = attributes(t,'add',uid2,'size', ...
                        xml_num2str(size(getfield(v(j),names{i}))));
                    t = xml_var2xml(t,getfield(v(j),names{i}),uid2);
                end
            end
        case 'cell'
            for i=1:prod(size(v))
                [t, uid2] = add(t,uid,'element','cell'); 
                % TODO % special handling of cellstr ?
                t = attributes(t,'add',uid2,'index',num2str(i));
                t = attributes(t,'add',uid2,'type',class(v{i}));
                t = attributes(t,'add',uid2,'size',xml_num2str(size(v{i})));
                t = xml_var2xml(t,v{i},uid2);
            end
        case 'char'
            % TODO % char values should be in CData
            if size(v,1) > 1
                t = add(t,uid,'chardata',v'); % row-wise order
            else
                t = add(t,uid,'chardata',v);
            end
        case {'int8','uint8','int16','uint16','int32','uint32'}
            [t, uid] = add(t,uid,'element',class(v));
            % TODO % Handle integer formats (cannot use sprintf or num2str)
        otherwise
            if ismember('serialize',methods(class(v)))
                % TODO % is CData necessary for class output ?
                t = add(t,uid,'cdata',serialize(v));
            else
                warning(sprintf(...
                '[SAVEXML] Cannot convert from %s to XML.',class(v)));
            end
    end

%=======================================================================
function s = xml_num2str(n)
    % TODO % use format ?
    if isempty(n)
        s = '[]';
    else
        s = ['[' sprintf('%g ',n(1:end-1))];
        s = [s num2str(n(end)) ']'];
    end
