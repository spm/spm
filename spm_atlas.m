function varargout = spm_atlas(action,varargin)
% Atlas multi-function
% FORMAT xA = spm_atlas('load',atlas)
% FORMAT L = spm_atlas('list')
% FORMAT [S,sts] = spm_atlas('select',xA,label)
% FORMAT Q = spm_atlas('query',xA,XYZmm)
% FORMAT [Q,P] = spm_atlas('query',xA,xY)
% FORMAT VM = spm_atlas('mask',xA,label,opt)
% FORMAT V = spm_atlas('prob',xA,label)
% FORMAT V = spm_atlas('maxprob',xA,thresh)
% FORMAT D = spm_atlas('dir')
%
% FORMAT url = spm_atlas('weblink',XYZmm,website)
% FORMAT labels = spm_atlas('import_labels',labelfile,fmt)
% FORMAT spm_atlas('save_labels',labelfile,labels)
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_atlas.m 6656 2015-12-24 16:49:52Z guillaume $


if ~nargin, action = 'load'; end

%==========================================================================
switch lower(action), case 'dir'
%==========================================================================
    % FORMAT D = spm_atlas('dir')
    %-Return list of directories potentially containing atlas files

    d = {...
         fullfile(spm('Dir'),'tpm');...
         fullfile(spm('Dir'),'atlas');...
        };
    
    varargout = { d };
    

%==========================================================================
case 'def'
%==========================================================================
    % FORMAT def = spm_atlas('def')
    %-Return link to atlas definition file [unused]
    
    def = 'http://www.fil.ion.ucl.ac.uk/spm/ext/atlas.xml';
    
    varargout = { def };
    
    
%==========================================================================
case 'load'
%==========================================================================
    % FORMAT xA = spm_atlas('load',atlas)
    %-Load atlas

    if isempty(varargin) || isempty(varargin{1})
        [atlas,sts] = spm_atlas('select');
        if ~sts, varargout = {[]}; return; end
    else
        atlas = varargin{1};
    end
    
    if isstruct(atlas), varargout = { atlas }; return; end
    
    %-Read Description
    switch spm_file(atlas,'ext')
        case 'xml'
            T         = convert(xmltree(atlas));
            xA.info   = T.header;
            xA.info   = rmfield(xA.info,'images');
            xA.info.files.labels = atlas;
            xA.info.files.images = spm_file(T.header.images.imagefile,'path',spm_file(atlas,'fpath'));
            
            xA.VA     = spm_vol(xA.info.files.images);
            
            idx       = cellfun(@(x) str2double(x.index),T.data.label,'UniformOutput',false);
            name      = cellfun(@(x) x.name,T.data.label,'UniformOutput',false);
            xA.labels = struct('name',name,'index',idx);
            
        case {'nii','img'}
            descfile  = spm_file(atlas,'ext','xml');
            if spm_existfile(descfile)
                xA    = spm_atlas('load',descfile);
            else
                xA.info.name = spm_file(atlas,'basename');
                xA.info.files.images = atlas;
                xA.VA = spm_vol(atlas);
                % assume a single image is a label image
                if numel(xA.VA) == 1
                    xA.info.type = 'label';
                    l = unique(spm_read_vols(xA.VA));
                    % discard 0/NaN/Inf?
                else
                    xA.info.type = 'probabilistic';
                    l = 1:numel(xA.VA);
                end
                for i=1:numel(l)
                    xA.labels(i) = struct('name',sprintf('Label %04d (%d)',i,l(i)),'index',l(i));
                end
            end
        otherwise
            list      = spm_atlas('list','installed');
            idx       = find(ismember(lower({list.name}),lower(atlas)));
            if numel(idx) == 1
                xA    = preloaded(list(idx).name);
                if isempty(xA)
                    xA = spm_atlas('load',list(idx).file);
                    preloaded(list(idx).name,xA);
                end
            elseif numel(idx) > 1
                error('Two or more atlases share the same name.');
            else
                error('Unknown atlas "%s".',atlas);
            end
    end

    varargout = { xA };
    
    
%==========================================================================
case 'list'
%==========================================================================
    % FORMAT L = spm_atlas('list',{'installed','available'})
    %-Return list of installed or available atlases
    
    if isempty(varargin), varargin = {'installed'}; end
    
    switch lower(varargin{1})
        case 'installed'
            L = atlas_list_installed(varargin{2:end});
        case 'available'
            atlaslist = spm_atlas('def');
            if ismember(atlaslist(1:4),{'http','file'})
                [s, sts] = urlread(atlaslist);
                if ~sts, error('Cannot access "%s".',atlaslist); end
                atlaslist = s;
            elseif ~spm_existfile(atlaslist)
                error('Cannot open "%s".',atlaslist);
            end
            L = convert(xmltree(atlaslist));
        otherwise
            error('Unknown option.');
    end
    
    varargout = { L };
    
    
%==========================================================================
case 'select'
%==========================================================================
    % FORMAT [S,sts] = spm_atlas('select',xA,label)
    %-Select atlas or labels
    
    S = '';
    if isempty(varargin)
        d = spm_atlas('Dir');
        d = d{1};
        [S,sts] = spm_select(1, {'image','xml'},...
            'Select Atlas...', {}, d);
        if ~sts, varargout = { S, sts }; return; end
    else
        xA = spm_atlas('load',varargin{1});
        if numel(varargin) == 1
            [sel,sts] = listdlg(...
                'ListString',{xA.labels.name},...
                'SelectionMode','multiple',...
                'ListSize', [400 300],...
                'Name','Select label(s)',...
                'PromptString',sprintf('Labels from %s atlas:',xA.info.name));
            if ~sts, varargout = { S, sts }; return; end
            S  = {xA.labels(sel).name};
        else
            sts = true;
            S = filter_labels(xA,varargin{2});
        end
    end
    varargout = { S, sts };
    
    
%==========================================================================
case 'query'
%==========================================================================
    % FORMAT Q = spm_atlas('query',xA,XYZmm)
    % FORMAT [Q,P] = spm_atlas('query',xA,xY)
    %-Atlas query
    
    xA = spm_atlas('load',varargin{1});
    if nargin > 2, xY = varargin{2}; else xY = struct; end
    
    unknown = 'Unknown';
    
    if numel(xA.VA) == 1 % or xA.info.type contains type definition
        if isnumeric(xY) && size(xY,2) == 1
            %-peak
            XYZmm     = xY;
            XYZ       = xA.VA.mat\[XYZmm;1];
            vpeak     = spm_sample_vol(xA.VA,XYZ(1),XYZ(2),XYZ(3),0);
            j         = [xA.labels.index] == vpeak;
            if any(j) == 1
                Q     = xA.labels(j).name;
            else
                Q     = unknown;
            end
            
            varargout = { Q };
            if nargout > 1, varargout = { {Q}, 100 }; end
        else
            %-cluster
            v         = spm_summarise(xA.VA,xY);
            vu        = unique(v);
            vun       = histc(v,vu);
            [vun,is]  = sort(vun(:),1,'descend');
            vu        = vu(is);
            for j=1:numel(vu)
                k     = [xA.labels.index] == vu(j);
                if any(k) == 1
                    Q{j} = xA.labels(k).name;
                else
                    Q{j} = unknown;
                end
                P(j)  = 100*vun(j)/numel(v);
            end
            
            varargout = { Q, P };
        end
    else
        P             = {xA.labels.name};
        if isnumeric(xY)
            %-peak
            XYZmm     = xY;
            XYZ       = xA.VA(1).mat\[XYZmm;1];
            Q         = spm_get_data(xA.VA,XYZ); % which interp to use?
            
            varargout = { Q, P };
        else
            %-cluster
            v         = spm_summarise(xA.VA,xY);
            v         = mean(v,2);
            
            varargout = { v, P };
        end
    end
    
     
%==========================================================================
case 'mask'
%==========================================================================
    % FORMAT VM = spm_atlas('mask',xA,label,opt)
    %-Return (binary) mask for given labels

    if nargin < 2, xA = ''; else xA = varargin{1}; end
    xA = spm_atlas('load',xA);
    if nargin < 3 || isempty(varargin{2}), label = spm_atlas('select',xA);
    else label = varargin{2}; end
    label = filter_labels(xA,label);
    
    if numel(xA.VA) == 1 % or xA.info.type contains type definition
        if nargin < 4 || isempty(varargin{3}), opt = 'binary';
        else opt = lower(varargin{3}); end
        VM = struct(...
            'fname',   [xA.info.name '_mask' spm_file_ext],...
            'dim',     xA.VA(1).dim,...
            'dt',      [spm_type('uint16') spm_platform('bigend')],...
            'mat',     xA.VA(1).mat,...
            'n',       1,...
            'pinfo',   [1 0 0]',...
            'descrip', sprintf('%s mask',xA.info.name));
        if strcmp(opt,'binary')
            VM.dat = false(VM.dim);
        else
            VM.dat = uint16(zeros(VM.dim));
        end
        
        D = spm_read_vols(xA.VA);
        for i=1:numel(label)
            j = find(ismember({xA.labels.name},label{i}));
            for k=1:numel(j)
                idx = xA.labels(j(k)).index;
                if strcmp(opt,'binary')
                    VM.dat = VM.dat | (D == idx);
                elseif strcmp(opt,'atlas')
                    VM.dat(D == idx) = i;
                elseif strcmp(opt,'preserve')
                    VM.dat(D == idx) = idx;
                else
                    error('Unknown option.');
                end
            end
        end
        VM.dat = uint16(VM.dat);
    else
        if nargin < 4, thresh = 0.5; else thresh = varargin{3}; end
        VM       = spm_atlas('prob',xA,label);
        VM.dt(1) = spm_type('uint8');
        VM.dat   = uint8(VM.dat > thresh);
    end
    
    varargout = { VM };
    % The output mask can be saved to disk with:
    % VM = spm_write_vol(VM,VM.dat);
    % VM = rmfield(VM,'dat');
    
    
%==========================================================================
case 'maxprob'
%==========================================================================
    % FORMAT V = spm_atlas('maxprob',xA,thresh)
    
    if nargin < 2, xA = ''; else xA = varargin{1}; end
    xA = spm_atlas('load',xA);
    if nargin < 3 || isempty(varargin{2}), thresh = 0;
    else thresh = varargin{2}; end
    
    typ = 'int16';
    
    V = struct(...
        'fname',   [xA.info.name '_maxprob_thresh' num2str(thresh) spm_file_ext],...
        'dim',     xA.VA(1).dim,...
        'dt',      [spm_type(typ) spm_platform('bigend')],...
        'mat',     xA.VA(1).mat,...
        'n',       1,...
        'pinfo',   [1 0 0]',...
        'descrip', sprintf('%s mask',xA.info.name));
    V.dat = zeros(V.dim);
    
    for i=1:V.dim(3)
        Y = zeros(V.dim(1),V.dim(2),numel(xA.VA));
        for j=1:numel(xA.VA)
            Y(:,:,j) = spm_slice_vol(xA.VA(j),spm_matrix([0 0 i]),V.dim(1:2),0);
        end
        [Y,V.dat(:,:,i)] = max(Y,[],3);
        V.dat(:,:,i) = V.dat(:,:,i) .* (Y > thresh);
    end
    V.dat = feval(cell2mat(spm_type(typ,'conv')),V.dat);
    
    varargout = { V };
    
    
%==========================================================================
case 'prob'
%==========================================================================
    % FORMAT V = spm_atlas('prob',xA,label)
    
    if nargin < 2, xA = ''; else xA = varargin{1}; end
    xA = spm_atlas('load',xA);
    if nargin < 3, label = spm_atlas('select',xA);
    else label = varargin{2}; end
    [label,idx] = filter_labels(xA,label);
    
    if numel(idx) == 1
        descrip = label{1};
    else
        descrip = sprintf('%s prob',xA.info.name);
    end
    V = struct(...
        'fname',   [xA.info.name '_prob' spm_file_ext],...
        'dim',     xA.VA(1).dim,...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'mat',     xA.VA(1).mat,...
        'n',       1,...
        'pinfo',   [1 0 0]',...
        'descrip', descrip);
    V.dat = zeros(V.dim);
    
    for i=1:numel(idx)
        V.dat = V.dat + spm_read_vols(xA.VA(idx(i)));
    end
    V.dat = single(V.dat);
    
    varargout = { V };
    
    
%==========================================================================
case 'install'
%==========================================================================
    % FORMAT sts = spm_atlas('install',A)
    %-Install Atlas
    
    if isempty(varargin)
        h = atlas_figure('SPM Atlases');
        H = getappdata(h,'Hbrowser');
        
        if ispc, s = '/'; else s = ''; end
        localurl = @(f) sprintf(['file://' s strrep(spm('Dir'),'\','/') '/help/' f]);
        
        tpl = spm_file_template(fullfile(spm('Dir'),'help'),'keep');
        tpl = tpl.file('TPL_ATLAS','spm_atlas.tpl');
        tpl = tpl.block('TPL_ATLAS','loading','load');
        tpl = tpl.block('TPL_ATLAS','listing','list');
        tpl = tpl.block('listing','atlas','atl');
        tpl = tpl.var('SPM',spm('Ver'));
        tpl = tpl.var('SPM_CSS',localurl('spm.css'));

        tpl = tpl.var('IMG_LOADING',localurl('images/loading.gif'));
        tpl = tpl.var('TXT_LOADING','Loading description of atlases...');
        tpl = tpl.parse('load','loading',0);
        tpl = tpl.var('list','');
        tpl = tpl.parse('OUT','TPL_ATLAS');
        html = get(tpl,'OUT');
        
        spm_browser(html,H);
        
        try
            L = spm_atlas('List','Available');
        catch
            tpl = tpl.var('TXT_LOADING','Cannot access atlases description.');
            tpl = tpl.var('IMG_LOADING','');
            tpl = tpl.parse('load','loading',0);
            tpl = tpl.parse('OUT','TPL_ATLAS');
            html = get(tpl,'OUT');
            spm_browser(html,H);
            varargout = { false };
            return;
        end
        
        tpl = tpl.var('load','');
        tpl = tpl.var('tbx','');
        for i=1:numel(L.atlas)
            tpl = tpl.var('ATLAS_URL',L.atlas{i}.website);
            tpl = tpl.var('ATLAS_NAME',L.atlas{i}.name);
            tpl = tpl.var('ATLAS_ID',L.atlas{i}.name);
            if isstruct(L.atlas{i}.maintainer)
                L.atlas{i}.maintainer = {L.atlas{i}.maintainer};
            end
            for j=1:numel(L.atlas{i}.maintainer)
                tpl = tpl.var('ATLAS_AUTHOR',L.atlas{i}.maintainer{j}.name);
                tpl = tpl.var('ATLAS_EMAIL',L.atlas{i}.maintainer{j}.email);
            end
            tpl = tpl.var('ATLAS_SUMMARY',L.atlas{i}.description);
            tpl = tpl.parse('atl','atlas',1);
        end
        tpl = tpl.parse('list','listing',0);
        tpl = tpl.parse('OUT','TPL_ATLAS');
        html = get(tpl,'OUT');
        
        spm_browser(html,H);
        
        varargout = { true };
        return;
        
    else
        A  = varargin{1};
        
        if ~(nargin == 3 && strcmpi(varargin{2},'-force'))
            AI = spm_atlas('list','installed');
            if ismember(A,{AI.name})
                warning('Atlas "%s" is already installed.',A);
                varargout = { true };
                return;
            end
        end
        
        AA = spm_atlas('list','available');
        for i=1:numel(AA.atlas)
            if strcmp(A,AA.atlas{i}.name)
                url = AA.atlas{i}.download;
                if isempty(url)
                    if spm_platform('desktop')
                        str = sprintf('<a href="%s">%s</a>',url,url);
                    else
                        str = url;
                    end
                    fprintf(['Atlas "%s" is only available through website:\n' ...
                        '  ' str '\n'], A);
                    if ~spm('CmdLine'), web(url,'-browser'); end
                    varargout = { false };
                    return;
                else
                    %-Check folder permissions
                    %------------------------------------------------------
                    dest = spm_atlas('Dir');
                    dest = dest{end};
                    [sts, attrb] = fileattrib(dest);
                    if ~sts, error('"%s"\n%s',dest,attrb); end
                    if ~attrb.UserWrite
                        error('No write access to "%s".\nMaybe use "%s" instead.',...
                            dest, strrep(userpath,pathsep,''));
                    end
                    
                    %-Download atlas archive
                    %------------------------------------------------------
                    tmpfile = [tempname(dest) '.zip'];
                    try
                        F = urlwrite(url,tmpfile);
                    catch
                        l = lasterror;
                        switch l.identifier
                            case 'MATLAB:urlwrite:ConnectionFailed'
                                error('Could not access URL "%s".',url);
                            case 'MATLAB:urlwrite:InvalidOutputLocation'
                                error('Could not create output file "%s".',tmpfile);
                            otherwise
                                rethrow(l);
                        end
                    end
                    
                    %-Unzip archive in destination folder
                    %------------------------------------------------------
                    try
                        FS = unzip(F,dest);
                    catch
                        spm_unlink(F);
                        error('Error when unpackig atlas archive');
                    end
                    
                    %-Delete atlas archive
                    %------------------------------------------------------
                    spm_unlink(F);
                    
                    %-Display <ATLAS>_README.txt if present
                    %------------------------------------------------------
                    fprintf('Atlas "%s" installed.',A);
                    
                    idx = find(~cellfun('isempty',...
                        regexpi(FS,sprintf('%s_README.txt$',A))));
                    if ~isempty(idx)
                        type(FS{idx(1)});
                    end
                    
                    % Refresh list of installed atlases.
                    spm_atlas('list','installed','-refresh');
                    
                    varargout = { true };
                    return;
                end
            end
        end
        
        warning('Cannot find atlas "%s".',A);
        varargout = { false };
        return;
    end
    
    
%==========================================================================
case 'weblink'
%==========================================================================
    % FORMAT url = spm_atlas('weblink',XYZmm,website)
    %-Return URL for coordinates query
    
    XYZmm = varargin{1};
    if nargin < 3, website = ''; else website = varargin{2}; end
    
    switch lower(website)
        case ''
            url = '';
        case 'brede'
            %-Brede Database - Talairach coordinate search
            url = 'http://neuro.imm.dtu.dk/cgi-bin/brede_loc_query.pl?q=%d+%d+%d';
        case 'neurosynth'
            %-Neurosynth - structure-to-function mappings
            url = 'http://neurosynth.org/locations/%d_%d_%d';
        otherwise
            error('Unknown website "%s".',website);
    end
    
    url = sprintf(url,XYZmm);
    
    varargout = { url };
    
    
%==========================================================================
case 'import_labels'
%==========================================================================
    % FORMAT labels = spm_atlas('import_labels',labelfile,fmt)
    %-Read labels stored in other formats
    
    if isempty(varargin) || isempty(varargin{1})
        [labelfile,sts] = spm_select(1,'any','Select labels file...');
        if ~sts, varargout = { struct }; return; end
    else
        labelfile = varargin{1};
    end
    
    if nargin < 3
        fmt = 'default';
    else
        fmt = varargin{2};
    end
    
    labels = read_labels(labelfile,fmt);
    
    varargout = { labels };
    
    
%==========================================================================
case 'save_labels'
%==========================================================================
    % FORMAT spm_atlas('save_labels',labelfile,labels)
    %-Save labels to file
    
    if isempty(varargin) || isempty(varargin{1})
        labelfile = 'labels.txt';
    else
        labelfile = varargin{1};
    end
    
    if nargin < 3
        labels = spm_atlas('import_labels');
    else
        labels = varargin{2};
    end
    
    save_labels(labelfile,labels);
    
    
%==========================================================================
otherwise
%==========================================================================
    error('Unknown action.');
end


%==========================================================================
% FUNCTION labels = read_labels(labelfile,fmt)
%==========================================================================
function labels = read_labels(labelfile,fmt)
fid        = fopen(labelfile,'rt');
if fid == -1, error('Cannot open atlas labels: "%s".',labelfile); end

try
    switch lower(fmt)
        case 'default'
        % Default
        % Key Long Name
        L = textscan(fid,'%d %s','Delimiter','','ReturnOnError',false);
        labels = struct('name',L{2},'index',num2cell(L{1}));
        
        case 'aal'
        % AAL
        % ShortName Long_Name Key
        L = textscan(fid,'%s %s %d','ReturnOnError',false);
        labels = struct('name',L{2},'index',num2cell(L{3}));

        case 'freesurfer'
        % FreeSurfer
        % Key Long-Name Red Green Blue Alpha
        L = textscan(fid,'%d %s %d %d %d %d','CommentStyle','#','ReturnOnError',false);
        labels = struct('name',L{2},'index',num2cell(L{1}));

        case 'brainvisa'
        % Brainvisa
        % Key, X, Y, Z, Red, Green, Blue, Acronym, Short Name
        L = textscan(fid,'%d %f %f %f %d %d %d %s %s','Delimiter',',','HeaderLines',1,'ReturnOnError',false);
        labels = struct('name',L{9},'index',num2cell(L{1}));

        case 'talairach'
        % Talairach
        % Key <TAB> Long Name
        L = textscan(fid,'%d %s','Delimiter','\t','ReturnOnError',false);
        labels = struct('name',L{2},'index',num2cell(L{1}));

        % MRIcron JHU-WhiteMatter
        % Key <TAB> Long_Name
        % L = textscan(fid,'%d %s','Delimiter','\t','ReturnOnError',false);
        % labels = {{} L{2} L{1}};

        % MRIcron Brodmann
        % labels = {{} cellfun(@(x) sprintf('Brodmann %d',x),num2cell(1:52),'UniformOutput',0) (1:52)};

        case 'wfu_pickatlas'
        % WFU PickAtlas
        % Key <TAB> Long Name <TAB> *
        L = textscan(fid,'%d %s %*[^\n]','Delimiter','\t','CommentStyle',{'[' ']'},'ReturnOnError',false);
        labels = struct('name',L{2},'index',num2cell(L{1}));

        case 'hammers_mith'
        % Hammers_mith
        % Long_Name *
        L = textscan(fid,'%s','Delimiter',',','HeaderLines',2,'MultipleDelimsAsOne',1,'ReturnOnError',false);
        labels = struct('name',L{1}(1:end-2),'index',num2cell((1:numel(L{1})-2)'));

        case 'cerebellum'
        % Joern's Cerebellum MNIsegment-MRICroN
        % Key Long_Name ???
        L = textscan(fid,'%d %s %d','Delimiter',' ','ReturnOnError',false);
        labels = struct('name',L{2},'index',num2cell(L{1}));

        case 'fsl'
        % FSL
        % XML: <label index="Key" x="X" y="Y" z="Z">Long Name</label>
        X = xmltree(labelfile);
        I = find(X,'/atlas/data/label');
        for i=1:numel(I)
            labels(i).name = get(X,children(X,I(i)),'value');
            A = attributes(X,'get',I(i)); A = [A{:}];
            labels(i).index = str2double(A(strcmp({A.key},'index')).val) + 1; % + 1?
        end

        case 'itk-snap'
        % Colin 27: ITK-SNAP Label Description File
        % Key Red Green Blue Alpha Vis Idx "Long Name"
        L = textscan(fid,'%d %d %d %d %d %d %d "%[^"]"','CommentStyle','#','ReturnOnError',false);
        labels = struct('name',L{8},'index',num2cell(L{1}));

        otherwise
            error('Unknown label file format.');
    end
catch
    fclose(fid);
    error('Cannot read atlas labels in: "%s".',labelfile);
end
fclose(fid);


%==========================================================================
% FUNCTION save_labels(labelfile,labels)
%==========================================================================
function save_labels(labelfile,labels)
switch spm_file(labelfile,'ext')
    case 'txt'
        fid = fopen(labelfile,'wt');
        if fid == -1, error('Cannot write file "%s".',labelfile); end
        for i=1:numel(labels)
            fprintf(fid,'%d %s\n',labels(i).index,labels(i).name);
        end
        fclose(fid);
    case 'xml'
%         t = xmltree;
%         t = set(t,root(t),'name','atlas');
%         t = attributes(t,'add',root(t),'version','2.0');
%         [t,hdr] = add(t,root(t),'element','header');
%         [t,nam] = add(t,hdr,'element','name');
%         t = add(t,nam,'chardata','');
%         [t,typ] = add(t,hdr,'element','type');
%         t = add(t,typ,'chardata','');
%         [t,img] = add(t,hdr,'element','images');
%         [t,imf] = add(t,img,'element','imagefile');
%         t = add(t,imf,'chardata','');
%         %[t,ims] = add(t,img,'element','summaryimagefile');
%         %t = add(t,ims,'chardata','');
%         [t,uid] = add(t,root(t),'element','data');
%         for i=1:numel(labels)
%             [t,lab] = add(t,uid,'element','label');
%             t = add(t,lab,'chardata',labels(i).name);
%             t = attributes(t,'add',lab,'index',num2str(labels(i).index));
%             %t = attributes(t,'add',lab,'x',num2str(0));
%             %t = attributes(t,'add',lab,'y',num2str(0));
%             %t = attributes(t,'add',lab,'z',num2str(0));
%         end
%         save(t,labelfile);
        fid = fopen(labelfile,'wt');
        if fid == -1, error('Cannot write file "%s".',labelfile); end
        fprintf(fid,'<?xml version="1.0" encoding="ISO-8859-1"?>\n');
        fprintf(fid,'<atlas version="2.0">\n');
        fprintf(fid,'\t<header>\n');
        fprintf(fid,'\t\t<name>%s</name>\n',spm_file(labelfile,'basename'));
        fprintf(fid,'\t\t<version>%s</version>\n','1.0');
        fprintf(fid,'\t\t<description>%s</description>\n',spm_file(labelfile,'basename'));
        fprintf(fid,'\t\t<url>%s</url>\n','');
        fprintf(fid,'\t\t<licence>%s</licence>\n','');
        fprintf(fid,'\t\t<coordinate_system>MNI</coordinate_system>\n');
        fprintf(fid,'\t\t<type>%s</type>\n','Label'); % or 'Probabilistic'
        fprintf(fid,'\t\t<images>\n');
        fprintf(fid,'\t\t\t<imagefile>%s</imagefile>\n',spm_file(labelfile,'path','','ext',spm_file_ext));
        fprintf(fid,'\t\t</images>\n');
        fprintf(fid,'\t</header>\n');
        fprintf(fid,'\t<data>\n');
        for i=1:numel(labels)
            fprintf(fid,'\t\t<label><index>%d</index><name>%s</name></label>\n',...
                labels(i).index,labels(i).name);
        end
        fprintf(fid,'\t</data>\n');
        fprintf(fid,'</atlas>\n');
        fclose(fid);
    otherwise
        error('Unknown label file format.');
end


%==========================================================================
% FUNCTION [labels,i] = filter_labels(xA,labels)
%==========================================================================
function [labels,i] = filter_labels(xA,labels)
% calls to 'intersect' should use 'stable' option and handle repetitions
if isnumeric(labels)
    [unused,idx] = intersect([xA.labels.index],labels);
    labels  = {xA.labels(idx).name};
elseif isstruct(labels)
    labels = {labels.name};
elseif ~iscellstr(labels)
    idx = ~cellfun(@isempty,regexp({xA.labels.name},labels));
    labels = {xA.labels(idx).name};
end
[labels,i] = intersect({xA.labels.name},labels);


%==========================================================================
% FUNCTION L = atlas_list_installed(refresh)
%==========================================================================
function L = atlas_list_installed(refresh)
persistent atlas_list
if nargin && strcmpi(refresh,'-refresh'), atlas_list = []; end
if isempty(atlas_list)
    atlas_list = struct('file',{},'name',{});
    d = spm_atlas('Dir');
    for i=1:numel(d)
        L = spm_select('FPList',d{i},'^.*\.xml$');
        if isempty(L), L = {}; else L = cellstr(L); end
        for j=1:numel(L)
            A.file = L{j};
            A.name = spm_file(L{j},'basename');
            if strncmp(A.name,'labels_',7), A.name = A.name(8:end); end
            atlas_list(end+1) = A;
        end
    end
end
L = atlas_list;


%==========================================================================
% FUNCTION xA = preloaded(atlas)
%==========================================================================
function xA = preloaded(atlas,xA)
persistent pl_atlas
persistent pl_xA
if isempty(pl_atlas), pl_atlas = {}; end
if isempty(pl_xA), pl_xA = {}; end
i = find(ismember(pl_atlas,atlas));
if nargin == 1
    if isempty(i)
        xA = [];
    else
        xA = pl_xA{i};
    end
else
    if isempty(i)
        pl_atlas = [pl_atlas atlas];
        pl_xA = [pl_xA {xA}];
    else
        pl_xA{i} = xA;
    end
end

%==========================================================================
% FUNCTION atlas_figure
%==========================================================================
function h = atlas_figure(name)

if ~nargin, name = 'Atlases'; end
h = spm_figure('FindWin','SPMatlas');
if ~isempty(h), set(h,'Name',name); return; end

h = figure(...
    'MenuBar',     'none',...
    'NumberTitle', 'off',...
    'Name',        name,...
    'Resize',      'off',...
    'Toolbar',     'none',...
    'Tag',         'SPMatlas',...
    'WindowStyle', 'Normal',... %'Modal'
    'Color',       [1 1 1],...
    'Visible',     'off');
pos = get(h,'Position');
pos([3 4]) = [350 400];
set(h,'Position',pos);

[H, HC] = spm_browser('<html></html>',h,[2 2 pos(3)-4 pos(4)-4],'html');
setappdata(h,'Hbrowser',H); setappdata(h,'HCbrowser',HC);
set(h,'Resize','on','ResizeFcn',@atlas_figure_resize);

set(h,'Visible','on');


%==========================================================================
% FUNCTION atlas_resize(obj,evt,varargin)
%==========================================================================
function atlas_figure_resize(obj,evt,varargin)
old_units = get(obj,'Units');
set(obj, 'Units','pixels');
figpos    = get(obj,'Position');
H         = getappdata(obj,'HCbrowser');
set(H,   'pos',[2 2 figpos(3)-4 figpos(4)-4]);
set(obj, 'Units',old_units);
