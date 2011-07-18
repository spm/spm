function varargout = loadxml(filename,varargin)
%LOADXML Load workspace variables from disk (XML file).
%  LOADXML FILENAME retrieves all variables from a file given a full 
%  pathname or a MATLABPATH relative partial pathname (see PARTIALPATH).
%  If FILENAME has no extension LOAD looks for FILENAME and FILENAME.xml 
%  and treats it as an XML file. 
%  
%  LOAD, by itself, uses the XML file named 'matlab.xml'. It is an error
%  if 'matlab.xml' is not found.
%
%  LOAD FILENAME X loads only X.
%  LOAD FILENAME X Y Z ... loads just the specified variables.  The
%  wildcard '*' loads variables that match a pattern.
%  Requested variables from FILENAME are created in the workspace.
%
%  S = LOAD(...) returns the contents of FILENAME in variable S. S is
%  a struct containing fields matching the variables retrieved.
%  
%  Use the functional form of LOAD, such as LOAD('filename'), when the
%  file name is stored in a string, when an output argument is requested,
%  or if FILENAME contains spaces.
%
%  See also LOAD, XML2MAT, XMLTREE.

%  Copyright 2003 Guillaume Flandin. 
%  $Revision: 4393 $  $Date: 2003/07/10 13:50 $

%  $Id: loadxml.m 4393 2011-07-18 14:52:32Z guillaume $

if nargin == 0
    filename = 'matlab.xml';
    fprintf('\nLoading from: %s\n\n',filename);
end

if ~ischar(filename)
    error('[LOADXML] Argument must contain a string.');
end

if ~exist(filename,'file')
    filename = [filename '.xml'];
    if ~exist(filename,'file')
        error(sprintf(...
        '[LOADXML] Unable to read file %s: file does not exist',filename));
    end
end

if nargout > 1,
    error('[LOADXML] Too many output arguments.');
end

t = xmltree(filename);

uid = children(t,root(t));

if nargout == 1
    % varargout{1} = struct([]); % Matlab 6.0 and above
end

flagfirstvar = 1;
for i=1:length(uid)
    if strcmp(get(t,uid(i),'type'),'element')
        vname = get(t,uid(i),'name');
        % TODO % No need to parse the whole tree 
        if isempty(varargin) | ismember(varargin,vname)
            v = xml_create_var(t,uid(i));
            if nargout == 1
                if flagfirstvar
                    varargout{1} = struct(vname,v);
                    flagfirstvar = 0;
                else
                    varargout{1} = setfield(varargout{1},vname,v);
                end
            else
                assignin('caller',vname,v);
            end
        end
    end
end

%=======================================================================
function v = xml_create_var(t,uid)
    type = attributes(t,'get',uid,'type');
    sz   = str2num(attributes(t,'get',uid,'size'));
    
    switch type
        case 'double'
            v = str2num(get(t,children(t,uid),'value'));
            if ~isempty(sz)
                v = reshape(v,sz);
            end
        case 'sparse'
            u = children(t,uid);
            for k=1:length(u)
                if strcmp(get(t,u(k),'name'),'row')
                    i = str2num(get(t,children(t,u(k)),'value'));
                elseif strcmp(get(t,u(k),'name'),'col')
                    j = str2num(get(t,children(t,u(k)),'value'));
                elseif strcmp(get(t,u(k),'name'),'val')
                    s = str2num(get(t,children(t,u(k)),'value'));
                end 
            end
           v = sparse(i,j,s,sz(1),sz(2));
        case 'struct'
            u = children(t,uid);
            v = []; % works with Matlab < 6.0
            for i=1:length(u)
                s(1).type = '()';
                s(1).subs = {str2num(attributes(t,'get',u(i),'index'))};
                s(2).type = '.';
                s(2).subs = get(t,u(i),'name');
                v = subsasgn(v,s,xml_create_var(t,u(i)));
            end
            if isempty(u),
                v = struct([]); % Need Matlab 6.0 and above
            end
        case 'cell'
            v = cell(sz);
            u = children(t,uid);
            for i=1:length(u)
                v{str2num(attributes(t,'get',u(i),'index'))} = ...
                    xml_create_var(t,u(i));
            end
        case 'char'
            if isempty(children(t,uid)) 
                v = '';
            else
                v = get(t,children(t,uid),'value');
            end
            try % this can fail if blank spaces are lost or entity escaping
                if ~isempty(sz)
                    if sz(1) > 1
                        v = reshape(v,fliplr(sz))'; % row-wise order
                    else
                        v = reshape(v,sz);
                    end
                end
            end
        case  {'int8','uint8','int16','uint16','int32','uint32'}
            % TODO % Handle integer formats
            warning(sprintf('%s matrices not handled.',type));
            v = 0;
        otherwise
            try,
                v = feval(class(v),get(t,uid,'value'));
            catch,
                warning(sprintf(...
                '[LOADXML] Cannot convert from XML to %s.',type));
            end
    end
