function bch_mat = spm_bch_bchmat(bch_mfile,typeA)
%---------------------------------------------------------------
%
% This function gets a m-file that describes the 
% bach process and returns a .mat file containing 
% the variables used during the batch.
%---------------------------------------------------------------

%---------------------------------------------------------------
% Programmers  :  when a type of analysis is added, 
% this file should be edited to define new 
% variables AND add them in the bch_names variable
% unless you wish to check the consistency of the 
% variables, which should be done here and NOT in
% function called in spm_batch


%---------------------------------------------------------------
% temporarily goes to the m file dir and comes back at the end

cwd = pwd;
try 
	cd(spm_str_manip(bch_mfile,'H'));
catch
	error([spm_str_manip(bch_mfile,'H') 'doesn''t exist']);
end

%-------------------- check the m-file exists ------------------
if ~(exist(spm_str_manip(bch_mfile,'rt')) == 2) 
   error(sprintf('%s not a m file in path', ...
	 spm_str_manip(bch_mfile,'rt')));
end

%-------------------- try to evaluate m-file -------------------
try 
  	feval(spm_str_manip(bch_mfile,'rt'));
catch
	str = sprintf('cannot eval %s',spm_str_manip(bch_mfile,'rt'));
  	error([str ' : check your paths or parameter file']);
end

cd(cwd); %- come back here 
clear cwd 


%===============================================================
switch typeA


%===============================================================
case 'analyses'

bch_names = {'analyses','type','work_dir','mfile'};


%===============================================================
case 'defaults_edit'

bch_names = {'defaults_edit', 'type_area' ...
	     'Misc','Printing','Hdr','Statistics', ...
	     'Normalisation','RealignCoreg','Reset'};

%===============================================================
case {'headers','means','realignment',...
     'normalisation','smooth'}

bch_names = {typeA};


%===============================================================
case 'model'


% missing fields added here for model       	  
%---------------------------------------------------------------
% Added fields are 
%   model(k).nsess
%   model(k).nscans
%   model(k).remain
%   model(k).interp
%---------------------------------------------------------------
%
if exist('model') == 1, 

   for k = 1:length(model)
      model(k).nsess = length(model(k).conditions);
   end

   if ~isfield(model(k),'time_sampl')
      for k = 1:length(model)
         for l = 1:model(k).nsess

   	    try %- if it can be evaluated...
   	         files = eval(model(k).files{l});
   	         e_nscans = size(files,1);
   	    catch
         	 e_nscans = size(model(k).files{l},1);
   	    end

            model(k).nscans(l) = e_nscans;
            model(k).remain{l} = 1:model(k).nscans(l);
            model(k).time_sampl{l} = [];
            model(k).interp = '';
         
          end %- for l = 1:model(k).nsess

      end %- for k = 1:length(model)
   end %- if ~isfield(model(k),'time_sampl')

   for k = 1:length(model)
      for l = 1:model(k).nsess
	   
	   try %- if it can be evaluated...
	      files = eval(model(k).files{l});
	      e_nscans = size(files,1);
	   catch
	      e_nscans = size(model(k).files{l},1);
	   end
      	   model(k).nscans(l) = e_nscans+length(model(k).time_sampl{l});
           model(k).remain{l} = setdiff(1:model(k).nscans(l),...
                                     model(k).time_sampl{l});
      
      end % for l = 1:model(k).nsess
   end % for k = 1:length(model)

end % if exist('model') == 1, 


% create regressors parametrics stochastics if necessary
%------------------------------------------------------------
if exist('model') == 1, for k = 1:length(model)

  if isempty(model(k).regressors) | ~any(model(k).regressors) ...
     | ~(exist('regressors')==1) %- not a variable in WS

	if (exist('regressors')==1) & length(regressors),
	   len = length(regressors);
	   regressors(len+1) = regressors(len);
	   regressors(len+1).number = 0;
	   model(k).regressors = (len+1)*ones(1,model(k).nsess);
	else
     	   regressors(1)  	= struct('number', 0);
	   model(k).regressors 	= ones(1,model(k).nsess);
	end
  end

  if isempty(model(k).parametrics) | ~any(model(k).parametrics) ...
     | ~(exist('parametrics')==1) %- not a variable in WS

	if (exist('parametrics')==1) & length(parametrics),
	   len = length(parametrics);
	   parametrics(len+1) = parametrics(len);
	   parametrics(len+1).type = 'none';
	   model(k).parametrics = (len+1)*ones(1,model(k).nsess);
	else
     	   parametrics(1)  	= struct('type', 'none');
	   model(k).parametrics = ones(1,model(k).nsess);
	end
  end

  if isempty(model(k).stochastics) | ~any(model(k).stochastics) ...
     | ~(exist('stochastics')==1) %- not a variable in WS

	if (exist('stochastics')==1) & length(stochastics),
	   len = length(stochastics);
	   stochastics(len+1) = stochastics(len);
	   stochastics(len+1).specify = 0;
	   model(k).stochastics = (len+1)*ones(1,model(k).nsess);
	else
     	   stochastics(1)  	= struct('specify', 0);
	   model(k).stochastics = ones(1,model(k).nsess);
	end
  end

end, end % if exist('model') == 1, for k = 1:length(model)


% missing fields added here for conditions       	  
%------------------------------------------------------------
% Added fields are 
%   conditions(k,l).number
%   conditions(k,l).fix_var_SOA
%   conditions(k,l).types
%------------------------------------------------------------

if exist('model') == 1, 

   if exist('conditions') == 1, for l = 1:length(conditions)

       conditions(l).number = length(conditions(l).names);	
       conditions(l).fix_var_SOA = 'Variable';
       conditions(l).types = cell(1,conditions(l).number);

       i_ev = find(conditions(l).bf_ev);
       i_ep = find(conditions(l).bf_ep); 
	    
       if ~isempty(i_ev)
          [conditions(l).types{i_ev}] = deal('events'); 
       end
       if ~isempty(i_ep)
          [conditions(l).types{i_ep}] = deal('epochs'); 
       end       

  end, end %- if exist('conditions') == 1, for l = 1:length(conditions) 

end % if exist('model') == 1

%- clear variables ....
clear i_ev i_ep k l str


bch_names = {'model','conditions','stochastics',...
  'regressors','parametrics','bf_ev','bf_ep'};


%===============================================================
case 'contrastes'


% missing fields added here for contrastes
%---------------------------------------------------------------
% Added fields are 
%   contrastes(k).set_action
%---------------------------------------------------------------

if exist('contrastes') == 1
   for k = 1:length(contrastes)
      if ~isfield(contrastes(k),'set_action')
       	ncont = length(contrastes(k).names);
       	if ncont
              tmp = cell(1,ncont);
              [tmp{:}] = deal('c');
              contrastes(k).set_action = tmp;
              clear tmp;
       	end
      end
   end
end %- if exist('contrastes') == 1


bch_names = {'contrastes'};

%===============================================================
case 'realignment'

   warning(sprintf('not implemented %s',typeA))


%===============================================================
otherwise
   warning(sprintf('unknown type of analyse %s',typeA))


%===============================================================
end


%-------- create bch_names with existing variables only
str = '';
names = {};

for i=1:length(bch_names), 
   if exist(bch_names{i})
       str = [str ' ' bch_names{i}];
       names{length(names)+1} = bch_names{i};
   end
end

bch_names = names;
str = [str  ' bch_names'];

% .mat saved here: here bch_mat does contain the path of the 
% bch_mfile (because we are in the working dir).
%---------------------------------------------------------------
eval(['save ' spm_str_manip(bch_mfile,'rp') ' ' str]);

bch_mat = spm_str_manip(bch_mfile,'rp'); %- here contains the path
%---------------------------------------------------------------



%===============================================================
%===============================================================
% sub function to check the consistency of the variables
%===============================================================

function ok = sf_chk_model(model)
ok=1; return;

function ok = sf_chk_conditions(conditions)

for l = 1:length(conditions)
  i_ev = find(conditions(l).bf_ev);
  i_ep = find(conditions(l).bf_ep); 
  if (any(intersect(i_ev,i_ep)) | ...
     any(setdiff(1:conditions(l).number,union(i_ev,i_ep))))
     str = sprintf('conditions(%d): bad event/epoch specif',l);
     error([str ' : check your parameter file']);
  end
end
ok = 1; return;
