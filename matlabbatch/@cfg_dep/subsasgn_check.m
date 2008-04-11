function [sts, val] = subsasgn_check(item,subs,val)

% function [sts, val] = subsasgn_check(item,subs,val)
% Do a check for proper assignments of values to fields. This routine
% will be called for derived objects from @cfg_dep/subsasgn with the
% original object as first argument and the proposed subs and val fields
% before an assignment is made.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

sts = true;
switch(subs(1).subs)
    case {'tname','sname'}
        if isempty(val)
            val = '';
        else
            sts = ischar(val);
            if ~sts
                warning('matlabbatch:cfg_dep:subsasgn_check:name', ...
                        'Dependency %s must be a string.', subs(1).subs);
            end;
        end;
    case {'tgt_exbranch','tgt_input','jtsubs','src_exbranch', ...
          'src_output'}
        if isempty(val)
            val = struct('type',{},'subs',{});
        else
            sts = isstruct(val) && numel(fieldnames(val))==2 && ...
                  all(isfield(val,{'type','subs'}));
            if ~sts
                warning('matlabbatch:cfg_dep:subsasgn_check:subs', ...
                        'Dependency reference %s must be a substruct.', ...
                        subs(1).subs);
            end;
        end;
    case {'tgt_spec'}
        if isempty(val)
            val = cfg_findspec;
        else
            sts = iscell(val);
            if sts
                for k = 1:numel(val)
                    sts = isstruct(val{k}) && numel(fieldnames(val{k}))==2 ...
                                   && all(isfield(val{k},{'name', ...
                                        'value'}));
                    if ~sts
                        break;
                    end;
                end;
            end;
            if ~sts
                warning('matlabbatch:cfg_dep:subsasgn_check:tgt_spec', ...
                        'Target specification must be a cfg_findspec.');
            end;
        end;
end;