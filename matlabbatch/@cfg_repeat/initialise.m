function item = initialise(item, val, dflag)

% function item = initialise(item, val, dflag)
% Disassemble val struct and pass its values on to child nodes. If dflags
% is set, then values are passed on to the matching values items. If
% dflags is not set, then item.val is set to the initialised copy of the
% matching values item(s).
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: initialise.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

if numel(item.values)==1 && isa(item.values{1},'cfg_branch') && ~item.forcestruct,
    if dflag
        item.values{1} = initialise(item.values{1}, val, dflag);
    else
        for k = 1:numel(val)
            item.cfg_item.val{k} = initialise(item.values{1}, val(k), dflag);
        end;
    end;
else
    if dflag
        if numel(item.values) > 1 || item.forcestruct
            % val should be either a cell array containing structs with a
            % single field (a harvested job), or a struct with multiple
            % fields (a harvested defaults tree). In the latter case,
            % convert val to a cell array before proceeding.
            if isstruct(val)
                fn = fieldnames(val);
                for k = 1:numel(fn)
                    val1{k} = struct(fn{k}, val.(fn{k}));
                end;
                val = val1;
            end;
            for l = 1:numel(val)
                vtag = fieldnames(val{l});
                for k = 1:numel(item.values)
                    if strcmp(gettag(item.values{k}), vtag{1})
                        item.values{k} = initialise(item.values{k}, ...
                            val{l}.(vtag{1}), ...
                            dflag);
                    end;
                end;
            end;
        else
            item.values{1} = initialise(item.values{1}, val{1}, dflag);
        end;
    else
        if numel(item.values) > 1 || item.forcestruct
            for l = 1:numel(val)
                % val{l} should be a struct with a single field
                vtag = fieldnames(val{l});
                for k = 1:numel(item.values)
                    if strcmp(gettag(item.values{k}), vtag{1})
                        item.cfg_item.val{l} = initialise(item.values{k}, ...
                            val{l}.(vtag{1}), ...
                            dflag);
                    end;
                end;
            end;
        else
            for l = 1:numel(val)
                item.cfg_item.val{l} = initialise(item.values{1}, ...
                    val{l}, dflag);
            end;
        end;
    end;
end;

