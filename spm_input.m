function varargout = spm_input(varargin)
% FORMAT spm_input(varargin);
%---------------------------------------------------------------
% spm_input( gui_arg, 'batch', mat-file, batch_arg)
%
% all arguments are optional.
%
% gui_arg	: traditional gui arguments
% 'batch'	: batch marker to delimit where the non GUI arg
%		  start.
% mat-file	: the name of a mat file containing variables that
%		  will be used in spm_input
% 		  NB The mat-file can be a structure obtained with 
%		  a mat-file = load(mat-file).
% batch_arg	: of the form : {indexing_part}, addressing_part
%		  indexing_part should look like :
%		  {'var1',idx1,'var2',idx2, ... 'varN',idxN}
%		  and all var1..N should exist in mat-file.
%		  indexing_part is specifying a variable in mat-file,
%		  called last_var in this code through the indexing 
%		  system. (see batch documentation).
%		  addressing_part should look like a series of 
%		  'field_name' or index used to address last_var
%		  (see indexing_part) to return last_var.field_name(index)
%		  or last_var(index).field_name. There can be an unlimited
%		  number of arguments in this addressing part. 
%
% To summarise, spm_input has two modes :
% 1- varargin contains the key word 'bach', and spm_input
%    will run in batch mode. 	
% 2- GUI version if 'batch' key word does't exist in varargin
%   (cf Andrew Holmes former spm_input renamed spm_input_orig)
%
% Some part of the spm code will never be executed in batch mode.
% In these parts, the spm_input should have no batch marker.
%
% In case where there are no argument for the batch mode, and 
% batch marker is followed by a mat-file arg that is not empty,
% spm_input returns with no GUI execution. If the mat-file
% exists but is empty, then  spm_input_orig is called.
%  
%---------------------------------------------------------------
% %W% Jean-Baptiste Poline, Stephanie Rouquette %E%





%-------------------------------------------------------------------
%---- ib = indice where batch starts
%------------------------------------
ib = min(find(strcmp(varargin,'batch')));

%---- if no batch marker 
%------------------------------------
if isempty(ib), 
   if nargout
      varargout = {spm_input_orig(varargin{:})}; return; 
   else
      spm_input_orig(varargin{:}); return;
   end
end

%---- if no argument after batch marker just return, don't call GUI
%------------------------------------
if ib == nargin, warning('batch marker should be followed'); return; end  

%---- no batch mat file ?
%------------------------------------
if isempty(varargin{ib+1})
	if ib == 1, % No arg before batch marker
	   return
	else
   	   if nargout
	      varargout = {spm_input_orig(varargin{1:ib-1})}; return; 
   	   else
	      spm_input_orig(varargin{1:ib-1}); return;
	   end
	end
end;

%----  at least one batch argument after the mat name : 
%------------------------------------
bchvarin = varargin(ib+2:end);		%- contains the varargin for batch 
                                        %- after the mat-file 
bchmat 	= varargin{ib+1};		%- mat-file
nbchin 	= length(bchvarin);		%- # of arg passed 
if nbchin == 0, return, end

%-------------------------------------------------------------------
%---------------  first construct the indices...
%-------------------------------------------------------------------
curarg = 1;

%- get the bch_var_lev and bch_names : must be in bchmat
% bch_var_lev = sf_get_var(bchmat,'bch_var_lev');

bch_names = sf_get_var(bchmat,'bch_names');
indices   = {};
cind      = bchvarin{curarg};


if ~isempty(cind)

 if iscell(cind)

   %---  a quick check of the indexing argument.
   %--------------------------------------------
   if ~sf_check_indices(cind,bch_names), 
	cind, bch_names,
	error('parsing indices sf_check_indices'); 
   end

   %- initialise lastvar 
   %--------------------------------------------

   % no_ind = 1;
   lastvar = sf_get_var(bchmat,cind{1});
   indices = cind(2);
   
   for k = 3:2:length(cind)
%      no_ind = no_ind+1;
      try
          indices{1} = getfield(lastvar,indices,cind{k},cind(k+1));
          lastvar = sf_get_var(bchmat,cind{k});
      catch
          indices, lastvar, cind,
	  error('parsing indices catch 1');
          return
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
   	lastvar = lastvar(indices{:});
catch
	lastvar
	indices{:}
	cind
	error('indices <= 0 or too large ');
end

%-------------------------------------------------------------------
%---------------  second access the fields... 
%-------------------------------------------------------------------
% lastvar in memory

if nbchin >= 2

	for i = (curarg+1):nbchin

		% bchvarin{i}
		% lastvar
		%------------ argument is string ---------
		if isstr(bchvarin{i}), s_typ = '.'; subs =bchvarin{i};
		%------------ argument is numeric or cell ---------
		else
			if iscell(lastvar), s_typ = '{}'; else, s_typ = '()'; end
			%------------ argument is numeric ---------
			if isnumeric(bchvarin{i})
				if prod(size(bchvarin{i})) > 1, 
					error(' parsing : use cell to pass indices'); 
				end
				subs = {bchvarin{i}};	
			%------------ argument is cell ---------
			elseif iscell(bchvarin{i}) %- assumes these are indices ....
				subs = bchvarin{i};
			else 
				error('parsing : what kind of argument ? ')
			end
		end

		try, lastvar = subsref(lastvar, substruct(s_typ,subs));
		catch, 
			lastvar, s_typ, subs,
			varargout = {spm_input_orig(varargin{1:ib-1})}; return; 
		end

	end %- for i = 2:nbchin
	
end

% If lastvar is a string, try eval it before returning 
if isstr(lastvar)
   try, varargout = {eval(lastvar)};
   catch, varargout = {lastvar};
   end
else
   varargout = {lastvar}; 
end


%-------------------------------------------------------------------
function yn = sf_is_var_name(varname, names)
% 1 if varname is in names
% 0 ortherwise

yn = any(strcmp(varname,names));
return;
%-------------------------------------------------------------------
function v = sf_get_var(bchmat,varname)
% returns the variable named varname 
%
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
%-------------------------------------------------------------------
function yn = sf_check_indices(cind,names)
% checks that the argument can be parsed to get the indices ...
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
	    warning([sprintf('''%s'' ',c{:}) ' are not all members of  ', ...
			sprintf('''%s'' ',names{:}) ]);
	end
	b = all(sf_arenumeric(cind(2:2:end)));
	if ~b, 
	    c = cind(2:2:end)
	    warning(' arguments in indexing part are not all numerics');
	end
	yn = a & b ;
else 
	yn = 0;
	warning('length of indexing part should be even');
end
%-------------------------------------------------------------------
function in = sf_arenumeric(n)
% checks whether these are numerics ...
in = [];
if iscell(n),
	for i=1:length(n), in = [in isnumeric(n{i})]; end
else
	in = isnumeric(n);
end
%-------------------------------------------------------------------
