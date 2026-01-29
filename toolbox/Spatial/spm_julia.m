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
% FORMAT sts = spm_julia('add-package',package1, package2, ...)
% Add Julia packages, where package1, package2, etc are names of registered
% packages.  Unregistered packages can be installed by giving their url. e.g.
% sts = spm_julia('add-package', 'https://github.com/spm/PushPull.jl')
%
% FORMAT [sts,result] = spm_julia('run', cmd)
% Execute a command with julia -e
%
% FORMAT [sts,result] = spm_julia('run', cmd, package1, package2, ...)
% Checks for presence of packages, and installs them if necessary.
% Prepends use of packages before cmd.
% Package names can be urls of unregistered packages.
%
% For more information about the Julia language, see https://julialang.org/
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2026 Functional Imaging Lab, UCL Institute of Neurology

    s = settings;
    setenv('JULIA_DEPOT_PATH',s.julia_depot_path)
    switch lower(opt)
        case 'add-package'
            varargout{1} = add_packages(s, varargin);
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

function [sts, result] = run(s, fun, packages)
    install(s);
    if nargin>=3
        add_packages(s, to_install(s,packages));
        expr = '';
        for i=1:length(packages)
            package = package_name(packages{i});
            expr    = [expr 'using ' package '; '];
        end
        expr = [expr fun];
    else
        expr = fun;
    end
    [sts,result] = run_julia_cmd(s,expr);
end

function sts = install(s,varargin)
    sts = 0;
    if ~exist(s.cmd,'file') || (nargin>=2 && any(strcmp(varargin,'force')))
        json_file = websave(tempname,'https://julialang-s3.julialang.org/bin/versions.json');
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
               error('Not ready.')
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
                    error('...Sorry!\n');
            end
            delete(compr_file);
            fprintf(' ...Done\n')
        else
            %sts = 1;
            error('Something went wrong!');
        end
    end
end

function sts = add_packages(s,pkg_list)
    sts = 0;
    if ~isempty(pkg_list)
        fprintf('\n---- Installing Julia packages ----\n')
        for i=1:length(pkg_list)
            sts = sts | add_package(s,pkg_list{i});
        end
        fprintf('-----------------------------------\n')
    end
end

function sts = add_package(s,pkg)
    if ~any(pkg == '"')
        pkg = ['"' pkg '"'];
    end
    if regexp(pkg,'https')
        pkg = ['url=' pkg];
    elseif any(pkg=='/')
        fprintf('\n###### Not sure what to do with %s! ######\n\n', pkg)
    end
    [sts,result] = run_julia_cmd(s, ['import Pkg; Pkg.add(' pkg ');']);
    if sts~=0
        fprintf('\n\n###### Installation of %s failed! ######\n\n', pkg)
        disp(result)
        fprintf(  '\n###### Installation of %s failed! ###### (see errors above)\n\n', pkg)
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
    [sts,result] = system(f);
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
        package = list{i};
        if ~package_exists(s,package)
            not_installed = [not_installed, package];
        end
    end
end

function pkg = package_name(package)
    if any(package=='/')
        [~,pkg,~] = fileparts(package);
    else
        pkg = package;
    end
end

function answer = package_exists(s, package)
    dname = fullfile(s.julia_depot_path,'packages', package_name(package));
    if exist(dname,'dir')
        answer = true;
    else
        answer = false;
    end
end

