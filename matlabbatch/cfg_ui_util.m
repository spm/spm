function varargout = cfg_ui_util(cmd, varargin)
%CFG_UI_UTIL utility functions for displaying job, module and item values
% This function is a collection of utility functions to display a job,
% module or data summary. It also handles all value display and editing for
% a particular item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_ui_util.m 7473 2018-11-06 10:26:44Z guillaume $

rev = '$Rev: 7473 $';  %#ok<NASGU>

switch lower(cmd)
    case {'preview'}
        % cfg_ui_util('preview', ciid, dflag)
        [ciid, dflag] = deal(varargin{1:2});
        contents = cfg_ui_util('showitem', ciid, dflag);
        [tag, val] = cfg_util('harvest', ciid{:});
        try
            feval(contents{10}, val);
        end
    case {'showitemstr'}
        % [namestr datastr] = cfg_ui_util('showitemstr', contents, dflag)
        % get name and one-line data description for a single item (i.e.
        % contents{:}{k})
        [contents, dflag] = deal(varargin{1:2});
        if contents{6}-2 > 0
            indent = [' ' repmat('. ', 1, contents{6}-2)];
        else
            indent = '';
        end
        if contents{8} || (dflag && ~isempty(contents{2}))
            if any(strcmp(contents{5}, {'cfg_menu','cfg_files','cfg_entry'})) && ...
                    isa(contents{2}{1}, 'cfg_dep')
                if numel(contents{2}{1}) == 1
                    datastr = sprintf('DEP %s', contents{2}{1}.sname);
                else
                    datastr = sprintf('DEP (%d outputs)', numel(contents{2}{1}));
                end
            else
                switch contents{5}
                    case 'cfg_menu'
                        datastr = 'Unknown selection';
                        for l = 1:numel(contents{4})
                            if exist('isequalwithequalnans','builtin')
                                iseqn = isequalwithequalnans(contents{2}{1}, contents{4}{l});
                            else
                                iseqn = isequaln(contents{2}{1}, contents{4}{l});
                            end
                            if iseqn
                                datastr = contents{3}{l};
                                break;
                            end
                        end
                    case 'cfg_files'
                        if numel(contents{2}{1}) == 1
                            if isempty(contents{2}{1}{1})
                                datastr = ' ';
                            else
                                datastr = contents{2}{1}{1};
                            end
                        else
                            datastr = sprintf('%d files', numel(contents{2}{1}));
                        end
                    case 'cfg_entry'
                        csz = size(contents{2}{1});
                        % TODO use gencode like string formatting
                        if ischar(contents{2}{1}) && ...
                                numel(csz) == 2 && any(csz(1:2) == 1)
                            datastr = contents{2}{1};
                        elseif (isnumeric(contents{2}{1}) || ...
                                islogical(contents{2}{1})) && ...
                                numel(csz) == 2 && any(csz(1:2) == 1) &&...
                                numel(contents{2}{1}) <= 4
                            % always display line vector as summary
                            datastr = mat2str(contents{2}{1}(:)');
                        elseif any(csz == 0)
                            switch class(contents{2}{1})
                                case 'char'
                                    datastr = '''''';
                                case 'double'
                                    datastr = '[]';
                                otherwise
                                    datastr = sprintf('%s([])', ...
                                        class(contents{2}{1}));
                            end
                        else
                            szstr = sprintf('%dx', csz);
                            datastr = sprintf('%s %s', ...
                                szstr(1:end-1), class(contents{2}{1}));
                        end
                    otherwise
                        datastr = ' ';
                end
            end
        else
            datastr = '<-X';
        end
        namestr = sprintf('%s%s  ', indent, contents{1});
        varargout{1} = namestr;
        varargout{2} = datastr;
    case {'showitem'}
        % [contents, namestr, datastr] = cfg_ui_util('showitem', ciid, dflag)
        [ciid, dflag] = deal(varargin{1:2});
        [id, stop, contents] = ...
            cfg_util('listmod', ciid{:},...
            cfg_findspec({{'hidden',false}}), ...
            cfg_tropts({{'hidden', true}},1,1,1,1,dflag), ...
            {'name','val','labels','values','class','level', ...
            'all_set','all_set_item','num','preview'});
        contents = cellfun(@(c)subsref(c, substruct('{}',{1})), contents, 'UniformOutput', false);
        [namestr, datastr] = cfg_ui_util('showitemstr', contents, dflag);
        varargout{1} = contents;
        varargout{2} = namestr;
        varargout{3} = datastr;
    case {'showmod'}
        % [id, namestr, datastr, contents] = cfg_ui_util('showmod', cmid, dflag)
        [cmid, dflag] = deal(varargin{1:2});
        [id, stop, contents] = ...
            cfg_util('listmod', cmid{:}, [],...
            cfg_findspec({{'hidden',false}}), ...
            cfg_tropts({{'hidden', true}},1,Inf,1,Inf,dflag), ...
            {'name','val','labels','values','class','level', ...
            'all_set','all_set_item','num','preview'});
        if isempty(id) || ~cfg_util('isitem_mod_id', id{1})
            % Module not found without hidden flag
            % Try to list top level entry of module anyway, but not module items.
            [id, stop, contents] = ...
                cfg_util('listmod', cmid{:}, [],...
                cfg_findspec({}), ...
                cfg_tropts({{'hidden', true}},1,1,1,1,dflag), ...
                {'name','val','labels','values','class','level', ...
                'all_set','all_set_item','preview'});
        end
        namestr = cell(1,numel(id));
        datastr = cell(1,numel(id));
        namestr{1} = sprintf('Help on: %s',contents{1}{1});
        datastr{1} = '';
        for citem = 2:numel(id)
            [namestr{citem}, datastr{citem}] = cfg_ui_util('showitemstr', cellfun(@(c)subsref(c, substruct('{}',{citem})), contents, 'UniformOutput', false), dflag);
        end
        varargout = {id, namestr, datastr, contents};
    case {'showval'}
        % str = cfg_ui_util('showval', contents)
        % show verbose listing of contents for cfg_entry, cfg_files
        % show listing of available and selected options for cfg_menu,
        % cfg_choice
        % show listing of selected items for cfg_repeat (unused in cfg_ui)
        [contents, dflag] = deal(varargin{1:2});
        switch(contents{5})
            case {'cfg_entry','cfg_files'}
                if ~isempty(contents{2}) && isa(contents{2}{1}, 'cfg_dep')
                    str = {'Reference from'};
                    for k = 1:numel(contents{2}{1}) % we may have multiple dependencies
                        str{k+1} = contents{2}{1}(k).sname; % return something to be printed
                    end
                elseif ~isempty(contents{2})
                    if ndims(contents{2}{1}) <= 2
                        if ischar(contents{2}{1})
                            str = cellstr(contents{2}{1});
                        elseif iscellstr(contents{2}{1})
                            str = contents{2}{1};
                        elseif isnumeric(contents{2}{1}) || ...
                                islogical(contents{2}{1})
                            str = cellstr(num2str(contents{2}{1}));
                        else
                            str = gencode(contents{2}{1},'val');
                        end
                    else
                        str = gencode(contents{2}{1},'val');
                    end
                else
                    str = '';
                end
            case {'cfg_menu','cfg_choice'}
                if strcmp(contents{5},'cfg_menu') || ~dflag
                    if strcmp(contents{5},'cfg_choice')
                        % compare tag, not filled entries
                        cmpsubs = substruct('.','tag');
                    else
                        cmpsubs = struct('type',{},'subs',{});
                    end
                    valsubs = substruct('{}',{1});
                    nitem = numel(contents{4});
                    mrk = cell(1,nitem);
                    for l = 1:nitem
                        valuesubs = substruct('{}',{l});
                        if ~isempty(contents{2}) && isequal(subsref(contents{2},[valsubs cmpsubs]), subsref(contents{4},[valuesubs cmpsubs]))
                            mrk{l} = '*';
                        else
                            mrk{l} = ' ';
                        end
                    end
                    if strcmp(contents{5},'cfg_choice')
                        str = cell(1,nitem);
                        for k = 1:nitem
                            str{k} = contents{4}{k}.name;
                        end
                    else
                        str = contents{3};
                    end
                    str = strcat(mrk(:), str(:));
                end
            case {'cfg_repeat', 'cfg_mchoice'}
                if ~dflag
                    % Already selected items
                    ncitems = numel(contents{2});
                    str = cell(ncitems,1);
                    for k = 1:ncitems
                        str{k} = contents{2}{k}.name;
                    end
                end
            otherwise
                str = {};
        end
        varargout{1} = str;
    case 'showvaldeps'
        % List matching dependencies
        [job_id, mod_job_id, item_mod_id, sout] = deal(varargin{1:4});
        smatch = false(size(sout));
        % loop over sout to find whether there are dependencies that match the current item
        for k = 1:numel(sout)
            smatch(k) = cfg_util('match', job_id, mod_job_id, item_mod_id, sout(k).tgt_spec);
        end
        varargout{1} = sout(smatch);
    case 'showvaledit'
        % Fill value display boxes in a GUI, set up callbacks if necessary
        % required input args:
        % fig - figure object, created with guide. has to have guidata storing
        %       handles to (at least)
        %       handles.valshow
        %       handles.valshowLabel
        %       handles.helpbox
        % ciid  - item id
        % contents - contents for this item only
        % sout  - source output dependencies up to current module
        % dflag - defaults editing?
        % setvalcb - callback to store new value. Will be called with one
        %            argument (the new value), and should be instructed
        %            before where to store it. A default callback is
        %            provided to store the value in the referenced job
        %            item. This will be used if setvalcb is not a valid
        %            function handle.
        % updatecb - callback to redraw user interface. Must be called
        %            without any arguments.
        % GUI controls (de)activated if required
        % '^BtnVal.*'
        % '^MenuEditVal.*'
        % '^CmVal.*'
        % special menus for handling cfg_repeat editing:
        % '.*ValDelItem$'
        % '.*ValAddItem$'
        % '.*ValReplItem$'
        [fig, ciid, contents, sout, dflag, setvalcb, updatecb] = deal(varargin{1:7});
        set(findobj(fig,'-regexp', 'Tag','^BtnVal.*'), 'Visible','off');
        set(findobj(fig,'-regexp', 'Tag','^MenuEditVal.*'), 'Enable','off');
        set(findobj(fig,'-regexp', 'Tag','^CmVal.*'), 'Visible','off');
        delete(findobj(fig,'-regexp', 'Tag','^Val.*Dyn$'))
        set(findobj(fig,'-regexp', 'Tag','^valshow.*'), 'Visible','off');
        set(findobj(fig,'-regexp', 'Tag','.*Preview$'), 'Visible','off', ...
            'Enable','off');
        handles = guidata(fig);
        set(handles.valshow,'String', '','Min',0,'Max',0,'Callback',[]);
        set(handles.valshowLabel, 'String',sprintf('Current Item: %s',contents{1}));
        str = cfg_ui_util('showval', contents, dflag);
        udvalshow = local_init_udvalshow;
        if ~isempty(setvalcb) && subsasgn_check_funhandle(setvalcb)
            udvalshow.setvalcb = setvalcb;
        else
            udvalshow.setvalcb = @(nval)local_setvaledit(ciid, nval, dflag);
        end
        udvalshow.updatecb = updatecb;
        switch(contents{5})
            case {'cfg_entry','cfg_files'}
                set(findobj(fig,'-regexp', 'Tag','^valshow.*'), 'Visible','on');
                set(handles.valshow, 'Value',1, 'ListboxTop',1,'String', str, 'Userdata',udvalshow);
                if ~dflag
                    sout = cfg_ui_util('showvaldeps', ciid{:}, sout);
                    if ~isempty(sout)
                        set(findobj(fig,'-regexp','Tag','.*AddDep$'), ...
                            'Visible','on', 'Enable','on');
                    end
                end
                set(findobj(fig,'-regexp','Tag','.*EditVal$'), ...
                    'Visible','on', 'Enable','on');
                set(findobj(fig,'-regexp','Tag','.*ClearVal$'), ...
                    'Visible','on', 'Enable','on');
            case {'cfg_menu','cfg_choice'}
                if strcmp(contents{5},'cfg_menu') || ~dflag
                    cval = -1;
                    if strcmp(contents{5},'cfg_choice')
                        % compare tag, not filled entries
                        cmpsubs = substruct('.','tag');
                    else
                        cmpsubs = struct('type',{},'subs',{});
                    end
                    valsubs = substruct('{}',{1});
                    nitem = numel(contents{4});
                    for l = 1:nitem
                        valuesubs = substruct('{}',{l});
                        if ~isempty(contents{2}) && isequal(subsref(contents{2},[valsubs cmpsubs]), subsref(contents{4},[valuesubs cmpsubs]))
                            cval = l;
                        end
                    end
                    udvalshow.cval = cval;
                    if cval == -1
                        cval = 1;
                    end
                    udvalshow.cmd = num2cell(1:nitem);
                    ltop = cfg_ui_getListboxTop(handles.valshow, cval, numel(str));
                    set(findobj(fig,'-regexp', 'Tag','^valshow.*'), 'Visible','on');
                    set(handles.valshow, 'Value',cval, 'ListboxTop',ltop, 'String',str, ...
                        'Callback',@local_valedit_repeat, ...
                        'Keypressfcn',@local_valedit_key, ...
                        'Userdata',udvalshow);
                    set(findobj(fig,'-regexp','Tag','.*EditVal$'), ...
                        'Visible','on', 'Enable','on');
                    set(findobj(fig,'-regexp','Tag','.*ClearVal$'), ...
                        'Visible','on', 'Enable','on');
                end
            case {'cfg_repeat', 'cfg_mchoice'}
                if ~dflag
                    udvalshow.cval = -1;
                    % Already selected items
                    ncitems = numel(contents{2});
                    str3 = cell(ncitems,1);
                    cmd3 = cell(ncitems,1);
                    for k = 1:ncitems
                        str3{k} = sprintf('Delete: %s (%d)',...
                            contents{2}{k}.name, k);
                        cmd3{k} = [Inf k];
                        mdel = findobj(fig,'-regexp','Tag','.*ValDelItem$');
                        for cm = 1:numel(mdel)
                            uimenu(mdel(cm), ...
                                'Label',sprintf('%s (%d)', ...
                                contents{2}{k}.name, k), ...
                                'Callback',@(ob,ev)local_setvaledit(ciid, cmd3{k}, false, updatecb, ob, ev), ...
                                'Tag','ValDelItemDyn');
                        end
                    end
                    % Add/Replicate callbacks will be shown only if max number of
                    % items not yet reached
                    if (strcmp(contents{5}, 'cfg_repeat') && ncitems < contents{9}(2)) ...
                            || (strcmp(contents{5}, 'cfg_mchoice') && ncitems < numel(contents{4}))
                        % Available items
                        aitems = contents{4};
                        if strcmp(contents{5}, 'cfg_mchoice')
                            unsel = true(size(aitems));
                            for k = 1:numel(contents{2})
                                unsel(k) = ~any(cellfun(@(citem)isequal(contents{2}(k), citem), aitems));
                            end
                            aitems = aitems(unsel);
                        end
                        naitems = numel(aitems);
                        str1 = cell(naitems,1);
                        cmd1 = cell(naitems,1);
                        for k = 1:naitems
                            str1{k} = sprintf('New: %s', aitems{k}.name);
                            cmd1{k} = [k Inf];
                            madd = findobj(fig,'-regexp','Tag','.*ValAddItem$');
                            for cm = 1:numel(madd)
                                uimenu(madd(cm), ...
                                    'Label',aitems{k}.name, ...
                                    'Callback',@(ob,ev)local_setvaledit(ciid, cmd1{k}, false, updatecb, ob, ev), ...
                                    'Tag','ValAddItemDyn');
                            end
                        end
                        if strcmp(contents{5}, 'cfg_repeat')
                            str2 = cell(ncitems,1);
                            cmd2 = cell(ncitems,1);
                            for k = 1:ncitems
                                str2{k} = sprintf('Replicate: %s (%d)',...
                                    contents{2}{k}.name, k);
                                cmd2{k} = [-1 k];
                                mrepl = findobj(fig,'-regexp','Tag','.*ValReplItem$');
                                for cm = 1:numel(mrepl)
                                    uimenu(mrepl(cm), ...
                                        'Label',sprintf('%s (%d)', ...
                                        contents{2}{k}.name, k), ...
                                        'Callback',@(ob,ev)local_setvaledit(ciid, cmd2{k}, false, updatecb, ob, ev), ...
                                        'Tag','ValReplItemDyn');
                                end
                            end
                        else
                            str2 = {};
                            cmd2 = {};
                        end
                        set(findobj(fig,'-regexp','Tag','.*AddItem$'), ...
                            'Visible','on', 'Enable','on');
                        if ncitems > 0
                            set(findobj(fig,'-regexp','Tag','.*ReplItem$'), ...
                                'Visible','on', 'Enable','on');
                        end
                    else
                        str1 = {};
                        str2 = {};
                        cmd1 = {};
                        cmd2 = {};
                    end
                    str = [str1(:); str2(:); str3(:)];
                    udvalshow.cmd = [cmd1(:); cmd2(:); cmd3(:)];
                    set(findobj(fig,'-regexp', 'Tag','^valshow.*'), 'Visible','on');
                    set(handles.valshow, 'Value',1, 'ListboxTop',1, 'String', str, ...
                        'Callback',@local_valedit_repeat, ...
                        'KeyPressFcn', @local_valedit_key, ...
                        'Userdata',udvalshow);
                    set(findobj(fig,'-regexp','Tag','^Btn.*EditVal$'), ...
                        'Visible','on', 'Enable','on');
                    if ncitems > 0
                        set(findobj(fig,'-regexp','Tag','.*DelItem$'), ...
                            'Visible','on', 'Enable','on');
                    end
                    set(findobj(fig,'-regexp','Tag','.*ClearVal$'), ...
                        'Visible','on', 'Enable','on');
                end
        end
        [id, stop, help] = cfg_util('listmod', ciid{:}, cfg_findspec, ...
            cfg_tropts(cfg_findspec,1,1,1,1,false), {'showdoc'});
        set(handles.helpbox, 'Value',1, 'ListboxTop',1, 'string',cfg_justify(handles.helpbox, help{1}{1}));
        if contents{7} && ~isempty(contents{10})
            set(findobj(fig,'-regexp', 'Tag','.*Preview$'), 'Visible','on', ...
                              'Enable','on');
        end

    case 'valedit_editvalue'
        [ciid, itemname, val] = deal(varargin{1:3});
        [unused, unused, itemclass] = cfg_util('listmod', ciid{:}, cfg_findspec, ...
            cfg_tropts(cfg_findspec,1,1,1,1,false), {'class'});
        switch itemclass{1}{1}
            case {'cfg_entry'}
                [val, sts] = local_valedit_edit(ciid, itemname, val);
            case { 'cfg_files'}
                [val, sts] = local_valedit_files(ciid, itemname, val);
            case {'cfg_choice', 'cfg_mchoice', 'cfg_menu', 'cfg_repeat'}
                % does not return value - use udvalshow.updatecb inside
                % local_valedit_repeat as callback to update ui.
                sts = false;
                local_valedit_repeat(gcbf);
            otherwise
                sts = false;
        end
        if sts
            h = gcbf;
            if isempty(h) && exist('OCTAVE_VERSION', 'builtin')
                h = findall(0,'tag','cfg_ui');
            end
            handles = guidata(h);
            udvalshow = get(handles.valshow, 'Userdata');
            feval(udvalshow.setvalcb, val);
            feval(udvalshow.updatecb);
        end
    case 'setvaledit'
        local_setvaledit(varargin{:});
end

%% Callback and utility functions for showvaledit
% --------------------------------------------------------------------
function local_valedit_key(hObject, data, varargin)
if strcmpi(data.Key,'escape')
    % ESC must be checked here
    local_valedit_uiresume(hObject);
else
    % collect key info for evaluation in local_valedit_repeat
    handles = guidata(hObject);
    udvalshow = get(handles.valshow, 'Userdata');
    udvalshow.key = data;
    set(handles.valshow, 'Userdata',udvalshow);
end

% --------------------------------------------------------------------
function [val, sts] = local_valedit_edit(ciid, itemname, val)
% Normal mode. Depending on strtype, put '' or [] around entered
% input. If input has ndims > 2, isn't numeric or char, proceed with
% expert dialogue.
[id, stop, strtype] = cfg_util('listmod', ciid{:}, cfg_findspec, ...
                             cfg_tropts(cfg_findspec,1,1,1,1,false), {'strtype'});
if isempty(val) || isa(val{1}, 'cfg_dep')
    % silently clear cfg_deps
    if strtype{1}{1} == 's'
        val = {''};
    elseif strcmp(strtype{1}{1}, 's+')
        val = {{''}};
    else
        val = {[]};
    end
end
% If requested or we can't handle this, use expert mode
expmode = strcmp(cfg_get_defaults([mfilename '.ExpertEdit']), 'on') ||...
    ndims(val{1}) > 2 || ~(ischar(val{1}) || iscellstr(val{1}) || isnumeric(val{1}) || islogical(val{1}));
% Generate code for current value, if not empty
% Set dialog texts
if expmode
    if ~isequal(val, {''})
        instr = gencode(val{1},'val');
        % remove comments and put in 1-cell multiline char array
        nc = cellfun(@isempty,regexp(instr,'^\s*%'));
        instr = {char(instr(nc))};
    else
        instr = {''};
    end
    hlptxt = char({'Enter a valid MATLAB expression.', ...
        ' ', ...
        ['Strings must be enclosed in single quotes ' ...
        '(''A''), multiline string arrays in curly braces {} ' ...
        'and multiline arrays in brackets ([ ]).'], ...
        ' ', ...
        'To clear a value, enter an empty cell ''{}''.', ...
        ' ', ...
        'Leave input box with CTRL-TAB to access buttons.'});
    failtxt = {'Input could not be evaluated. Possible reasons are:',...
        '1) Input should be a vector or matrix, but is not enclosed in ''['' and '']'' brackets.',...
        '2) Input should be a character or string, but is not enclosed in '' single quotes.',...
        '3) Input should be a multiline string, but is not enclosed in ''{'' and ''}'' braces or strings not in '' single quotes.',...
        '3) Input should be a MATLAB variable, but is misspelled.',...
        '4) Input should be a MATLAB expression, but has syntax errors.'};
else
    if strtype{1}{1} == 's'
        instr = val;
        encl  = {'''' ''''};
    elseif strcmp(strtype{1}{1}, 's+')
        instr = {char(val{1})};
        encl  = {'''' ''''};
    else
        try
            instr = {num2str(val{1})};
        catch
            instr = {''};
        end
        encl  = {'[' ']'};
    end
    hlptxt = char({'Enter a value.', ...
        ' ', ...
        'To clear a value, clear the input field and accept.', ...
        ' ', ...
        'Leave input box with CTRL-TAB to access buttons.'});
    failtxt = {'Input could not be evaluated.'};
end
sts = false;
while ~sts
    % estimate size of input field based on instr
    % Maximum widthxheight 140x20, minimum 60x2
    szi = size(instr{1});
    mxwidth = 140;
    rdup = ceil(szi(2)/mxwidth)+3;
    szi = max(min([szi(1)*rdup szi(2)],[20 140]),[2,60]);
    str = inputdlg(hlptxt, ...
        itemname, ...
        szi,instr);
    if iscell(str) && isempty(str)
        % User has hit cancel button
        return;
    end
    % save instr in case of evaluation error
    instr = str;
    % str{1} is a multiline char array
    % 1) cellify it
    % 2) add newline to each string
    % 3) concatenate into one string
    cstr = cellstr(str{1});
    str = strcat(cstr, {char(10)});
    str = cat(2, str{:});
    % Evaluation is encapsulated to avoid users compromising this function
    % context - graphics handles are made invisible to avoid accidental
    % damage
    hv = cfg_ui_disable(0, 'HandleVisibility');
    [val, sts] = cfg_eval_valedit(str);
    cfg_ui_restore(hv);
    % for strtype 's', val must be a string
    sts = sts && (~strcmp(strtype{1}{1},'s') || ischar(val));
    % for strtype 's+', val must be a cellstr
    sts = sts && (~strcmp(strtype{1}{1},'s+') || iscellstr(val));
    if ~sts
        if ~expmode
            if strcmp(strtype{1}{1}, 's+')
                val = cstr;
                sts = true;
            else
                % try with matching value enclosure
                if strtype{1}{1} == 's'
                    if ishandle(val) % delete accidentally created objects
                        delete(val);
                    end
                    % escape single quotes and place the whole string in single quotes
                    str = strcat(encl(1), strrep(cstr,'''',''''''), encl(2), {char(10)});
                else
                    cestr = [encl(1); cstr(:); encl(2)]';
                    str = strcat(cestr, {char(10)});
                end
                str = cat(2, str{:});
                % Evaluation is encapsulated to avoid users compromising this function
                % context - graphics handles are made invisible to avoid accidental
                % damage
                hv = cfg_ui_disable(0, 'HandleVisibility');
                [val, sts] = cfg_eval_valedit(str);
                cfg_ui_restore(hv);
            end
        end
        if ~sts % (Still) no valid input
            uiwait(msgbox(failtxt,'Evaluation error','modal'));
        end
    end
end
% End of function will be reached with sts == true and new val

% --------------------------------------------------------------------
function [val, sts] = local_valedit_files(ciid, itemname, val)
[unused1, unused2, contents] = cfg_util('listmod', ciid{:},{},cfg_tropts({},1,1,1,1,false),{'num','filter','dir','ufilter'});
if isempty(val) || isa(val{1}, 'cfg_dep')
    inifile = '';
else
    inifile = val{1};
end
[val, sts] = cfg_getfile(contents{1}{1}, contents{2}{1}, itemname, inifile, contents{3}{1}, contents{4}{1});

% --------------------------------------------------------------------
function local_valedit_repeat(hObject,varargin)
handles = guidata(hObject);
udvalshow = get(handles.valshow, 'Userdata');
if ((isempty(udvalshow.key) || ...
        strcmpi(udvalshow.key.Key,'return')) && ...
        isequal(hObject, handles.valshow))
    % Mouse selection - no key
    % Keyboard selection - return key
    ccmd = get(hObject,'Value');
    if ccmd ~= udvalshow.cval
        feval(udvalshow.setvalcb, udvalshow.cmd{ccmd});
        feval(udvalshow.updatecb);
    end
elseif (~isempty(udvalshow.key) && strcmpi(udvalshow.key.Key,'escape') && ...
        isequal(hObject, handles.valshow))
    % callback called from handles.valshow, finish editing
    % nothing to do, if ui is not blocked
elseif ~isequal(hObject, handles.valshow)
    % callback called from elsewhere (module, menu, button)
    % if there is just one action, do it, else init editing
    if numel(udvalshow.cmd) == 1
        feval(udvalshow.setvalcb, udvalshow.cmd{1});
        feval(udvalshow.updatecb);
    else
        udvalshow.key = [];
        set(handles.valshow,'Userdata',udvalshow, 'Min',0, ...
            'Max',1);
        uicontrol(handles.valshow);
    end
else
    udvalshow.key = [];
    set(handles.valshow, 'Userdata',udvalshow);
end

% --------------------------------------------------------------------
function local_setvaledit(ciid, val, dflag, varargin)
if dflag
    cfg_util('setdef', ciid{:}, val);
else
    cfg_util('setval', ciid{:}, val);
    cfg_util('harvest', ciid{1:end-1});
end
if nargin > 3 && subsasgn_check_funhandle(varargin{1})
    % update GUI
    feval(varargin{1});
end

% --------------------------------------------------------------------
function udvalshow = local_init_udvalshow
% Initialise udvalshow to empty struct
udvalshow = struct('cval',[],'key',[],'updatecb',[],'setvalcb',[]);

