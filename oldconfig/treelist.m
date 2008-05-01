% TREELIST - Lists data in cell arrays and structs as ascii "tree"
%
% Version 1.1
%
% This functions lists the contents of structs, sub struct, cell arrays
% and sub cell array with chars: |-\ viewing the connection in the data.
% The main differents from the builtin DISP, DISPLAY is that this function
% lists all levels of sub cell arrays or struct in cells or structs.
%
% The syntax is similar to WHOS:
%
%   treelist varname1 varname2 .....
%
% By default, treelist does not show field contents and does not expand
% struct arrays. This behaviour can be changed by calling treelist as a
% function:
%
%   treelist('varname1','varname2',...,flags)
%
% where flags is a struct that may contain the following fields: 
% .dval - display values of fields? 
%         0 (default) don't display values
%         1 display values in tree structure
%         2 display values in (nearly) MATLAB script syntax, suitable for
%           output to a script file (see fname below) - see BUGs for the
%           limitations to this option
% .exps - expand struct arrays?
%         0 (default) don't expand
%         1 expand each member of struct arrays
% .fname - output file name. If not empty, output will be redirected to
%         the given file
%
% (C) Copyright 2002 Peter Rydesaeter, GNU General Public License.
%     2004 Modified by Volkmar Glauche to make struct array expansion and
%     value display optional.
%
% Se also:
%
%   WHOS, WHO, DISP, DISPLAY, CELL, STRUCT
%   
% BUGs:
% Displaying values as MATLAB code is limited by the following
% constraints: 
% * Object treatment is implemented only to a limited degree.

function treelist(varargin)
  defflags = struct('dval',0, 'exps',0, 'fid',1, 'fname','');
  if isstruct(varargin{end})
    flags = fillstruct(defflags,varargin{end});
    nvars = nargin-1;
  else
    flags = defflags;
    nvars = nargin;
  end;
  if ~isempty(flags.fname)
          flags.fid = fopen(flags.fname,'w');
  end;        
  for n=1:nvars,
      try
          v = evalin('caller',varargin{n});
          iname = varargin{n};
      catch
          v = varargin{n};
          iname = inputname(n);
      end;
    % for documentation, list things twice if creating MATLAB code and
    % writing to a file - first listing the variable structure only
    if flags.dval==2 && flags.fid > 2
            flags1=flags;
            flags1.dval=0;
            treelistsub(v,iname,'','',flags1);
    end;
    treelistsub(v,iname,'','',flags);
  end
  if flags.fid > 2
          fclose(flags.fid);
  end;
  return;
  
function treelistsub(dt,name,nameprefix,level,flags)
if nargin<4, level='';  end
if nargin<3, nameprefix='';  end
if nargin<2, name=inputname(1); end
whosdt=whos('dt');
if isobject(dt)
    dtclass='object';
else
    dtclass=whosdt.class;
end;
switch dtclass
case {'double', 'logical', 'single', 'uint8', 'uint16', 'uint32', 'uint64', ...
      'int8', 'int16', 'int32', 'int64', 'char'}
    if isempty(dt),
       if flags.dval <= 1
           dtstr={{'[]'}};
       elseif flags.dval == 2
           dtstr={{''}};
       end;
    else
        if (flags.dval == 1) || strcmp(dtclass, 'char')
            dtstr=textscan(evalc('format compact;format long g;disp(full(dt));format'), '%s', ...
                           'delimiter', char(10));
        elseif (flags.dval == 2) && ~strcmp(dtclass, 'char')
            % Display full(dt(:)) and reshape it after printing - this
            % should remedy all problems with sparse and multidim arrays
            dtstr=textscan(evalc('format compact;format long g;disp(full(dt(:)));format'),...
                           '%s', 'delimiter', char(10));
        else
            dtstr={{sprintf('%s%d %s', sprintf('%d-x-',whosdt.size(1:end-1)), ...
                            whosdt.size(end), dtclass)}};
        end;
    end
    if flags.dval < 2
        if length(level)==0,
            ss=sprintf('%s',name);
        else
            ss=sprintf('%s-%s ',level,name);
        end
        lv=length(level)+20;
        if length(ss)<lv, ss(end+1:lv)='.'; end
        idx=[1 find(ss=='-')];
        ss2=ss;
        ss2(idx(end):end)=' ';
        for k=1:numel(dtstr{1})
            if length(dtstr{1}{k})<79-length(ss),
                fprintf(flags.fid,'%% %s %s\n',ss,dtstr{1}{k});
            else
                if k==1,
                    fprintf(flags.fid,'%%%s\n', ss);
                end
                fprintf(flags.fid,'%%%s\n', dtstr{1}{k});
            end
            ss=ss2;
        end
    else
        if strcmp(dtclass,'char')
            reshead='';
            resfoot='';
        else
            reshead='reshape(';
            resfoot=sprintf(',[%s%d])',sprintf('%d,',whosdt.size(1:end-1)), ...
                            whosdt.size(end));
        end;
        if whosdt.sparse
            sphead='sparse(';
            spfoot=')';
        else
            sphead='';
            spfoot='';
        end;
        headstr=sprintf('%s%s = %s%s%s([\n', ...
                        nameprefix, name, reshead, sphead, dtclass);
        fprintf(flags.fid,headstr);
        indentstr=char(repmat(' ',1, ...
                              length(headstr)));
        if strcmp(dtclass, 'char')
            delim='''';
        else
            delim='';
        end;
        for k=1:numel(dtstr{1})
            fprintf(flags.fid,indentstr);
            if strcmp(dtclass, 'char')
                % if there is a string with single
                % quotes, replace them with double ones
                dtstr{1}{k}=strrep(dtstr{1}{k},'''', ...
                                   '''''');
            end;
            fprintf(flags.fid,'%s%s%s\n', delim,dtstr{1}{k},delim);
        end;
        fprintf(flags.fid,'%s])%s%s;\n',indentstr,spfoot,resfoot);
            
    end;
  case {'struct','object'}
    fn=fieldnames(dt);
    if numel(dt)==0,
      fprintf(flags.fid,sprintf('%% %s-%s Empty STRUCT\n',level,name));
      if flags.dval > 1
              fprintf(flags.fid,'%s%s = struct([]);\n',nameprefix,name);
      end;
      return;
    elseif numel(dt)>1,
      if flags.exps
        fprintf(flags.fid,'%% %s-%s \n',level,name);
        level(find(level=='\'))=' ';
        if flags.dval > 1
                fprintf(flags.fid,['%%==============================================' ...
                         '================\n']);
                fprintf(flags.fid,'%% %s%s\n',nameprefix,name);
                fprintf(flags.fid,['%%==============================================' ...
                         '================\n']);
        end;
        % get ndims and size to produce correct indices
        funcstr=get_index_func(dt);
        for m=1:numel(dt),
          eval(funcstr);
          newname = sprintf('%s(%s)',name,msubstr);
          if flags.dval > 1
                  fprintf(flags.fid,['%%----------------------------------------------' ...
                           '----------------\n']);
                  fprintf(flags.fid,'%% %s%s\n',nameprefix,newname);
                  fprintf(flags.fid,['%%----------------------------------------------' ...
                           '----------------\n']);
          end;
          if m==numel(dt),
            treelistsub(dt(m),newname,nameprefix,[level ' \'],flags);
          else
            treelistsub(dt(m),newname,nameprefix,[level ' |'],flags);
          end
        end
      else
        dtstr=sprintf('%s-%s %d-x-',level,name,whosdt.size(1:end-1));
        level(find(level=='\'))=' ';
        dtstr=fprintf(flags.fid,'%% %s%d struct array with fields\n',dtstr,whosdt.size(end));
        for n=1:numel(fn)
          fprintf(flags.fid,'%%%s    %s\n', level, fn{n});
        end;
      end;
      return;
    else
      if flags.dval<2
              fprintf(flags.fid,'%% %s-%s\n',level,name);
      end;
      level(find(level=='\'))=' ';
      ww=warning;                      %%HACK To remove warning msg
      warning off;
      if flags.dval==2 && strcmp(dtclass,'object')
              fprintf(flags.fid,'%s%s = %s;\n', ...
                      nameprefix,name,whosdt.class);
      end;
      for n=1:numel(fn),
        dts=getfield(dt,fn{n});
        newname=sprintf('%s',fn{n});
        newnameprefix=sprintf('%s%s.',nameprefix,name);
        if n==numel(fn),
      treelistsub(dts,newname,newnameprefix,[level ' \'],flags);
        else
      treelistsub(dts,newname,newnameprefix,[level ' |'],flags);
        end
      end;
      warning(ww);
    end
    return;
  case 'cell',
    if numel(dt)==0,
      fprintf(flags.fid,'%% %s-%s Empty CELL\n',level,name);      
      if flags.dval > 1
              fprintf(flags.fid,'%s%s = {};\n',nameprefix,name);
      end;
      return;
    end
    fprintf(flags.fid,'%% %s-%s \n',level,name);
    level(find(level=='\'))=' ';
    % get ndims and size to produce correct indices
    funcstr=get_index_func(dt);
    for m=1:numel(dt),
      eval(funcstr);
      newname=sprintf('%s{%s}',name,msubstr);
      if m==numel(dt),
    treelistsub(dt{m},newname,nameprefix,[level ' \'],flags);
      else
    treelistsub(dt{m},newname,nameprefix,[level ' |'],flags);
      end
    end
    return;
  otherwise
          fprintf('%% %s-%s Unknown item of class ''%s''\n', level, name, ...
                       whosdt.class);
  end
  return;
  
function funcstr = get_index_func(dt,varargin)
% produce code to correctly convert linear indexes to subscript indexes
% for variable dt
% Assumptions:
% - index array will be called 'ind'
% - running variable is called 'm'
% - printout of subscript index is called 'msubstr'

defflags = struct('indname','ind', 'runname','m', 'subsname','msubstr');
if nargin > 1
        flags = fillstruct(defflags,varargin{1});
else
        flags = defflags;
end;

% get ndims and size to produce correct indices
nddt = ndims(dt);
szdt = size(dt);
% omit last dimension for 1-x-N arrays (vectors), but don't do that for
% N-x-1 arrays. The former ones may be created using assignments like
% dt(k)=val_of_dt, while the latter ones must have been created with
% dt(k,1)=val_of_dt and must be assumed to be deliberately assigned to
% columns instead of rows.
if szdt(1)==1
        nddt = nddt-1;
        szdt = szdt(2:end);
end;
indstr = sprintf(sprintf('%s(%%d),',flags.indname),1:nddt);
funcstr = sprintf(['[%s] = ind2sub([%s],%s);%s=sprintf(''%%d,'',%s);' ...
'%s=%s(1:end-1);'], indstr(1:end-1), num2str(szdt), flags.runname, ...
                  flags.subsname, flags.indname, flags.subsname, flags.subsname);
