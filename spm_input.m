function varargout = spm_input(varargin)
% Input wrapper function for batch mode input
% FORMAT spm_input(gui_arg,'batch', mat-file, batch_arg)
%
% gui_arg     - Standard GUI input arguments (See spm_input_ui.m)
% 'batch'     - dummy argument delimit where batch (non GUI) arguments start
% bch_mat     - the name of a mat file containing variables that will be
%               used in spm_input
%             - NB The mat-file can be a structure obtained with 
%               a mat-file = load(mat-file).
% batch_arg   - of the form : {indexing_part, addressing_part}
%             - indexing_part should look like:
%                 {'var1',idx1,'var2',idx2, ... 'varN',idxN}
%               and all var1..N should exist in mat-file.
%             - indexing_part is specifying a variable in mat-file,
%               called last_var in this code through the indexing 
%               system. (see batch documentation: spm_bch.man).
%             - addressing_part should look like a series of 
%               'field_name' or index used to address last_var
%               (see indexing_part) to return last_var.field_name(index)
%               or last_var(index).field_name. There can be an unlimited
%               number of arguments in this addressing part. 
%
%_______________________________________________________________________
%
% To summarise, spm_input has two modes:
% 1- varargin contains the key word 'bach':
%    spm_input will run in batch mode.    
% 2- varargin does not contains the key word 'bach':
%    arguments are passed on to GUI input program for handling
%    (see spm_input_ui.m, formerly named spm_input.m)
%
% Some parts of the spm code should never be executed in batch mode.
% In these cases the calls to spm_input should have no batch marker.
%
% In case where there are no arguments for the batch mode:
%    - if global variable BCH is not empty,
%      spm_input returns immediately (with no GUI call).
%    - if global variable is empty (not in bch mode),
%      spm_input_ui is called.
%
%
% See also: spm_bch.m, spm_bch.man
%  
%_______________________________________________________________________
% %W% Jean-Baptiste Poline, Stephanie Rouquette %E%


%-Get global BCH definition, if any
% (BCH contains bch_mat, index0, and flag. See spm_bch.m for details)
%-----------------------------------------------------------------------
global BCH %- contains bch_mat, index0, and flag. see spm_bch.m


%---- ib = indice where batch starts
%-----------------------------------------------------------------------
ib = min(find(strcmp(varargin,'batch')));


%---- if no batch marker
%-----------------------------------------------------------------------
if isempty(ib)
	if nargout
		varargout = cell(1,nargout);
		[varargout{:}] = spm_input_ui(varargin{:});
	else
		spm_input_ui(varargin{:});
	end
	return
end

%---- if batch marker, not BCH mode 
%-----------------------------------------------------------------------
if isempty(BCH)
	if nargout
		varargout = cell(1,nargout);
		[varargout{:}] = spm_input_ui(varargin{1:ib-1});
	else
		spm_input_ui(varargin{1:ib-1});
	end
	return
end


%---- batch marker, BCH mode: Check arguments for BCH
%-----------------------------------------------------------------------
if ib == nargin, return, end

%----  at least one batch argument after the mat name : 
%-----------------------------------------------------------------------
bchvarin = varargin(ib+1:end);   %- contains the varargin for batch 
                                 %- after the 'batch' marker 
bchmat 	= BCH.bch_mat;	         %- mat-file (or structure in memory)
nbchin 	= length(bchvarin);      %- # of arg passed 

%-----------------------------------------------------------------------
%---------------  first construct the indices...
%-----------------------------------------------------------------------

bch_names = sf_get_var(bchmat,'bch_names');
indices   = {};
cind      = {BCH.index0{:} bchvarin{1}{:}};

if ~isempty(cind)

 if iscell(cind)

   %---  a quick check of the indexing argument.
   %--------------------------------------------------------------------
   if ~sf_check_indices(cind,bch_names), 
	cind, bch_names,
	error('parsing indices sf_check_indices'); 
   end

   %- initialise lastvar 
   %--------------------------------------------------------------------
   lastvar = sf_get_var(bchmat,cind{1});
   indices = cind(2);
   
   for k = 3:2:length(cind)
      if ~any(cind{k+1}), 
        if nbchin >= 2 %- a null indice and an argument after the
	               %- indices section => return 0
	   varargout = {0}; 
           return;
	 else	       %- null indice but no arg after the 
	               %- indice section => return empty
	   varargout = {[]}; 
           return;
        end
      else     
        try
          % Petr Janata's fix...
          % indices{1} = getfield(lastvar,indices,cind{k},cind(k+1));
          if isa(cind{k+1},'double')
            evalstr = sprintf('indices{1} = getfield(lastvar,indices,cind{k},{%d});', cind{k+1});
            eval(evalstr);
          else
            indices{1} = getfield(lastvar,indices,cind{k},cind(k+1));
          end
          lastvar = sf_get_var(bchmat,cind{k});
        catch
          indices, lastvar, cind,
	  error('parsing indices catch 1');
          return
	end  
      end 
   end %- for k = 3:2:length(cind)

 else
   error('parsing indices ~iscell(cind)');
 end

else
 error('parsing indices isempty(cind)');
end

% disp('after parsing indices'), indices
try 
  % Petr Janata's fix...
  % lastvar = lastvar(indices{:});
  if isstruct(lastvar)
    evalstr = sprintf('lastvar = lastvar(%d);', indices{:});
    eval(evalstr)
  else
    lastvar = lastvar(indices{:});
  end
catch
	lastvar
	indices{:}
	cind
	error('indices <= 0 or too large ');
end

%-----------------------------------------------------------------------
%---------------  second access the fields... 
%-----------------------------------------------------------------------
% lastvar in memory

if nbchin >= 2
  
  for i = 2:nbchin

	%------------ argument is string -------------------------------
	if isstr(bchvarin{i}), s_typ = '.'; subs = bchvarin{i};
	%------------ argument is numeric or cell ----------------------
	else
	   if iscell(lastvar), s_typ = '{}'; else, s_typ = '()'; end
	   %------------ argument is numeric ---------------------------
	   if isnumeric(bchvarin{i})
	      if prod(size(bchvarin{i})) > 1, 
                 error(' parsing : use cell to pass indices'); 
              end
	      subs = {bchvarin{i}};	
	   %------------ argument is cell ------------------------------
	   elseif iscell(bchvarin{i}) %- assumes these are indices ....
                  subs = bchvarin{i};
	   else 
                  error('parsing : what kind of argument ? ')
	   end
	end

	try,    lastvar = subsref(lastvar, substruct(s_typ,subs));
	catch
	   warning('couldn''t access the sub ref ...  ');
	   lastvar, s_typ, subs,
	   varargout = cell(1,nargout);
	   [varargout{:}] = spm_input_ui(varargin{1:ib-1});
	   return
	end

   end %- for i = 2:nbchin
	
end

% If lastvar is a string, try eval it before returning 
%-----------------------------------------------------------------------
if isstr(lastvar)
   try, varargout = {eval(lastvar)};
   catch, varargout = {lastvar};
   end
else
   varargout = {lastvar}; 
end


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================



function yn = sf_is_var_name(varname,names)
%=======================================================================
% 1 if varname is in names
% 0 otherwise
yn = any(strcmp(varname,names));



function v = sf_get_var(bchmat,varname)
%=======================================================================
% returns the variable named varname 

if ~isstr(varname), error('varname must be str'); v = []; end
if isstr(bchmat)
   try 
      eval(['v = load(''' bchmat ''',''' varname ''');']);
   catch
      error(['I''m afraid I can''t load ' varname ' in ' bchmat]);
   end
elseif isstruct(bchmat)
   v = bchmat;
else 
   bchmat
   error('What the hell is it ? should be string or struct');
end

if isempty(fieldnames(v)), 
   v = []; 
else v = getfield(v,varname);
end



function yn = sf_check_indices(cind,names)
%=======================================================================
% checks that the argument can be parsed to get the indices...
if ~length(cind) 
	yn = 1;
elseif ~rem(length(cind),2) 

	try
	  a = all(ismember(cind(1:2:end),names));
	catch %- debug
	  names
	  for i=1:length(names), isstr(names{i}), names{i}, end
	  for i=1:2:length(cind), isstr(cind{i}) , cind{i}, end
	end
	if ~a, 
	    c = cind(1:2:end);
	    warning([sprintf('''%s'' ',c{:}) ' not all members of ', ...
			sprintf('''%s'' ',names{:}) ]);
	end
	b = all(sf_are_idx(cind(2:2:end)));
	if ~b, 
	    c = cind(2:2:end)
	    warning(' arguments in indexing part are not all numerics');
	end
	yn = a & b ;
else 
	yn = 0;
	warning('length of indexing part should be even');
end



function in = sf_are_idx(n)
%=======================================================================
% checks whether these are numerics or ':' ...
in = [];
if iscell(n),
	for i=1:length(n), 
	    in = [in isnumeric(n{i}) | strcmp(n{i},':')]; 
	end
else
	in = isnumeric(n);
end
