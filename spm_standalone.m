function spm_standalone(varargin)
% Gateway function for standalone SPM
%
% References:
%
%   SPM Standalone:  https://en.wikibooks.org/wiki/SPM/Standalone
%   MATLAB Compiler: http://www.mathworks.com/products/compiler/
%
% See also: config/spm_make_standalone.m
%__________________________________________________________________________
% Copyright (C) 2010-2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_standalone.m 7372 2018-07-09 16:50:44Z guillaume $ 


if ~nargin, action = ''; else action = varargin{1}; end

if strcmpi(action,'run')
    warning('"Run" is deprecated: use "Batch".');
    action = 'batch';
end

exit_code = 0;

switch lower(action)
    
    case {'','pet','fmri','eeg','quit'}
    %----------------------------------------------------------------------
        spm(varargin{:});
    
    case {'-h','--help'}
    %----------------------------------------------------------------------
        cmd = lower(spm('Ver'));
        fprintf([...
            '%s - Statistical Parametric Mapping\n',...
            'https://www.fil.ion.ucl.ac.uk/spm/\n',...
            '\n'],...
            upper(cmd));
        fprintf([...
            'Usage: %s [ fmri | eeg | pet ]\n',...
            '       %s COMMAND [arg...]\n',...
            '       %s [ -h | --help | -v | --version ]\n',...
            '\n',...
            'Commands:\n',...
            '    batch          Run a batch job\n',...
            '    script         Execute a script\n',...
            '    function       Execute a function\n',...
            '    eval           Evaluate a MATLAB expression\n',...
            '    [NODE]         Run a specified batch node\n',...
            '\n',...
            'Options:\n',...
            '    -h, --help     Print usage statement\n',...
            '    -v, --version  Print version information\n',...
            '\n',...
            'Run ''%s [NODE] help'' for more information on a command.\n'],...
            cmd, cmd, cmd, cmd);
        
    case {'-v','--version'}
    %----------------------------------------------------------------------
        spm_banner;
        
    case 'batch'
    %----------------------------------------------------------------------
        inputs = varargin(2:end);
        flg = ismember(inputs,'--silent');
        if any(flg)
            quiet = true;
            inputs = inputs(~flg);
        else
            quiet = false;
        end
        flg = ismember(inputs,'--modality');
        if any(flg)
            idx = find(flg);
            try
                modality = inputs{idx+1};
            catch
                error('Syntax is: --modality <modality>.');
            end
            inputs([idx idx+1]) = [];
        else
            modality = 'fmri';
        end
        flg = ismember(inputs,'--cmdline');
        if any(flg)
            cmdline = true;
            inputs = inputs(~flg);
        else
            cmdline = false;
        end
        spm_banner(~quiet);
        spm('defaults',modality);
        spm_get_defaults('cmdline',cmdline);
        spm_jobman('initcfg');
        if isempty(inputs)
            h = spm_jobman;
            waitfor(h,'Visible','off');
        else
            try
                spm_jobman('run', inputs{:});
            catch
                fprintf('Execution failed: %s\n', inputs{1});
                exit_code = 1;
            end
        end
        spm('Quit');
        
    case 'script'
    %----------------------------------------------------------------------
        spm_banner;
        assignin('base','inputs',varargin(3:end));
        try
            if nargin > 1
                spm('Run',varargin(2));
            else
                spm('Run');
            end
        catch
            exit_code = 1;
        end
        
    case 'function'
    %----------------------------------------------------------------------
        spm_banner;
        if nargin == 1
            fcn = spm_input('function name','!+1','s','');
        else
            fcn = varargin{2};
        end
        try
            feval(fcn,varargin{3:end});
        catch
            exit_code = 1;
        end
    
    case 'eval'
    %----------------------------------------------------------------------
        spm_banner;
        if nargin == 1
            expr = spm_input('expression to evaluate','!+1','s','');
        else
            expr = varargin{2};
        end
        try
            eval(expr);
        catch
            exit_code = 1;
        end
        
    otherwise % cli
    %----------------------------------------------------------------------
        %spm('defaults','fmri');
        %spm_get_defaults('cmdline',true);
        spm_jobman('initcfg');
        try
            spm_cli(varargin{:});
        catch
            exit_code = 1;
        end
        
end

%-Display error message and return exit code (or use finish.m script)
%--------------------------------------------------------------------------
if exit_code ~= 0
    err = lasterror;
    msg{1} = err.message;
    if isfield(err,'stack')
        for i=1:numel(err.stack)
            if err.stack(i).line
                l = sprintf(' (line %d)',err.stack(i).line);
            else
                l = '';
            end
            msg{end+1} = sprintf('Error in %s%s',err.stack(i).name,l);
        end
    end
    fprintf('%s\n',msg{:});
    
    exit(exit_code);
end


%==========================================================================
function spm_banner(verbose)
% Display text banner
if nargin && ~verbose, return; end
[vspm,rspm] = spm('Ver');
tlkt = ver(spm_check_version);
if isdeployed, mcr = ' (standalone)'; else mcr = ''; end
fprintf('%s, version %s%s\n',vspm,rspm,mcr);
fprintf('%s, version %s\n',tlkt.Name,version);
spm('asciiwelcome');
