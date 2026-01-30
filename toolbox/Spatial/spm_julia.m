function varargout = spm_julia(opt, varargin)
% Run Julia code from SPM
% FORMAT s = spm_julia('settings')
%     s - a structure containing various settings
%
% FORMAT sts = spm_julia('install')
% Install Julia within SPM. If it exists already, then do nothing.
%
% FORMAT sts = spm_julia('install','force')
% Install Julia within SPM - irrespective of whether it seems to be
% there or not.
%
% FORMAT sts = spm_julia('add-package', package1, package2, ...)
% Add Julia packages, where package1, package2, etc are names of registered
% packages.  Unregistered packages can be installed by giving their url
% within a structure. Note that fields 'rev' and 'version' can also be
% included in the structure, where 'version' could be e.g. 'v0.2.1'
% (v"0.2.1" in Julia) or '0.2' ("0.2" in Julia).  For example:
%     pkg_spec = struct('url','https://github.com/spm/PushPull.jl');
%     sts      = spm_julia('add-package', pkg_spec)
%
% FORMAT [sts,result] = spm_julia('run', cmd)
% Execute a command with julia -e
%
% FORMAT [sts,result] = spm_julia('run', cmd, package1, package2, ...)
% Checks for presence of packages, and installs them if necessary.
% Prepends use of packages before cmd.
%
% For more information about the Julia language, see https://julialang.org/
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2026 Functional Imaging Lab, UCL Institute of Neurology

    s = settings;
    setenv('JULIA_DEPOT_PATH',s.julia_depot_path)
    switch lower(opt)
        case 'add-package'
            varargout{1} = add_pkgs(s, varargin);
        case 'install'
            varargout{1} = install(s, varargin{:});
        case 'run'
            [varargout{1:nargout}] = run(s,varargin{1},varargin(2:end));
        case 'settings'
            varargout{1} = s;
        otherwise
            if nargin==1
                [varargout{1:nargout}] = run(s,opt);
            else
                [varargout{1:nargout}] = run(s,opt,varargin);
            end
    end
end

function [sts, result] = run(s, fun, pkgs)
    install(s);
    if nargin>=3
        add_pkgs(s, to_install(s,pkgs));
        expr = '';
        for i=1:length(pkgs)
            pkg  = pkg_name(pkgs{i});
            expr = [expr 'using ' pkg '; '];
        end
        expr = [expr fun];
    else
        expr = fun;
    end
    if nargout<2
        sts = run_julia_cmd(s,expr);
    else
        [sts,result] = run_julia_cmd(s,expr);
    end
end

function sts = install(s,varargin)
    sts = 0;
    if ~exist(s.cmd,'file') || (nargin>=2 && any(strcmp(varargin,'force')))
        json_file = websave(tempname,'https://julialang-s3.julialang.org/bin/versions.json');
        if isempty(json_file)
            error('Can not obtain the versions.json file from the web.')
        end
        json = spm_jsonread(json_file);
        delete(json_file);

        files = json.(['x' replace(s.version,'.','_')]).files;
        switch s.comp
            case 'GLNXA64'
                opt = findfile('x86_64-linux-gnu', files);
            case 'PCWIN64'
                opt = findfile('x86_64-w64-mingw32', files);
            case 'MACI64'
                % x86_64
                opt = findfile('x86_64-apple-darwin14', files);
            case 'MACA64'
                % Unsure
                opt = findfile('aarch64-apple-darwin14', files);
            case 'ARM'
                % Unsure
                opt = findfile('aarch64-linux-gnu', files);
            otherwise
               error(['Something went wrong! Computer is ' s.comp '.'])
        end

        if ~isempty(opt)
            mkdir_rec(s.appdir)
            fprintf('Downloading "%s" ... ', opt.url)
            compr_file = fullfile(s.appdir,['install.' opt.extension]);
            websave(compr_file, opt.url);
            fprintf('...Done\n')
            fprintf('Unpacking "%s" ... ', compr_file);
            switch opt.extension
                case {'tar.gz','tar','tgz'}
                    untar(compr_file,s.appdir)
                case 'zip'
                    unzip(compr_file,s.appdir)
                otherwise
                    error(['Something went wrong! Extension is ' opt.extension '.']);
            end
            delete(compr_file);
            fprintf(' ...Done\n')
        else
            error(['Something went wrong! Nothing suitable for ' s.comp '.']);
        end
        spm_registry = 'https://github.com/spm/SPM-registry.jl';
        [sts,result] = run_julia_cmd(s,['using Pkg; pkg"registry add General"; pkg"registry add ' spm_registry '"']);
        if sts~=0
            error('Failed to add SPM-registry!')
        end
    end
end


function sts = add_pkgs(s, pkg_list)
    sts = 0;
    if ~isempty(pkg_list)
        fprintf('\n---- Installing Julia packages ----\n')
        for i=1:length(pkg_list)
            sts = sts | add_pkg(s,pkg_list{i});
        end
        fprintf('-----------------------------------\n')
    end
end


function sts = add_pkg(s, pkg)
    [sts,result] = run_julia_cmd(s, ['import Pkg; Pkg.add(' addstr(pkg) ');']);
    if sts~=0
        fprintf('\n\n###### Installation of %s failed! ######\n\n', pkg_name(pkg))
        disp(result)
        fprintf(  '\n###### Installation of %s failed! ###### (see errors above)\n\n', pkg_name(pkg))
    end
end

function [sts,result] = run_julia_cmd(s,str)
    if s.iswin
        f = [s.cmd ' -e "' strrep(str,'"','\"') '"'];
    else
        f = [s.cmd ' -e '' ' str ' '' '];
    end
    disp(['julia -e ' str])

    if isunix
        paths = getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH');
    end
    if nargout<2
        sts = system(f);
    else
        [sts,result] = system(f);
    end
    if isunix
        setenv('LD_LIBRARY_PATH', paths);
    end
end

function mkdir_rec(dr)
    [parent,~] = fileparts(dr);
    if ~exist(dr,'dir')
        mkdir_rec(parent);
        mkdir(dr);
    end
end


function opt = findfile(arch, files)
    opt = {};
    for i=1:length(files)
        if strcmp(arch,files{i}.triplet) && ~strcmp('exe',files{i}.extension)
            opt = files{i};
            return
        end
    end
end

function s = settings
    version = '1.10.10';
    comp    = computer;
    iswin   = ispc;
    if iswin
        ext = '.exe';
    else
        ext = '';
    end
    s = struct;
    s.version = version;
    s.appdir  = fullfile(spm('dir'),'toolbox','Spatial','Apps',lower(comp));
    s.cmd     = fullfile(s.appdir,['julia-' version],'bin',['julia' ext]);
    s.julia_depot_path = fullfile(s.appdir,'julia_depot');
    s.iswin   = iswin;
    s.comp    = comp;
end

function not_installed = to_install(s,list)
    not_installed = {};
    for i = 1:length(list)
        pkg = list{i};
        if ~pkg_exists(s,pkg)
            not_installed = [not_installed, pkg];
        end
    end
end


function pkg = pkg_name(pkg)
    if isa(pkg,'struct')
        if isfield(pkg,'name')
            [~,pkg,~] = fileparts(pkg.name);
        elseif isfield(pkg,'url')
            [~,pkg,~] = fileparts(pkg.url);
        elseif isfield(pkg,'path')
            [~,pkg,~] = fileparts(pkg.path);
        else
            pkg = '';
        end
    else
        [~,pkg,~] = fileparts(pkg);
    end
end


function answer = pkg_exists(s, pkg)
    dname = fullfile(s.julia_depot_path,'packages', pkg_name(pkg));
    if exist(dname,'dir')
        answer = true;
    else
        answer = false;
    end
end


function str = addstr(pkg)
    if isa(pkg,'struct')
        if isfield(pkg,'name')
            str = ['name="' pkg.name '"'];
        elseif isfield(pkg,'url')
            str = ['url="' pkg.url '"'];
        elseif isfield(pkg,'path')
            str = ['path="' pkg.path '"'];
        end
        if isfield(pkg,'version')
            str = [str verstr(pkg.version)];
        end
        if isfield(pkg,'rev')
            str = [str ', rev="' pkg.rev '"'];
        end
    else
        str = ['"' pkg '"'];
    end
end

function str = verstr(version)
    if length(version)>1 && version(1)=='v'
        str = [', version=v"' version(2:end) '"'];
    else
        str = [', version="' version '"'];
    end
end

