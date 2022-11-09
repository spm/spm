function D = spm_eeg_ctf_synth_gradiometer(S)
% Apply CTF synthetic gradiometry to MEG data
% FORMAT D = spm_opm_synth_gradiometer(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object or string to path     - Default: no Default
%   S.gradient      - Integer ranging from 0-3 defining
%                     order of gradiometry                  - Default: 3
%   S.method        - string of package to perform
%                     gradiometry correction                - Default: 'fieldtrip'
%   S.prefix        - string prefix for output MEEG object  - Default: 'g_'
% Output:
%   D               - denoised MEEG object (also written to disk)
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


%-Set default options
%--------------------------------------------------------------------------
spm('FnBanner', mfilename);
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),            error(errorMsg);            end
if ~isfield(S, 'gradient'),     S.gradient = 3;             end
if ~isfield(S, 'prefix'),       S.prefix = 'g_';            end
if ~isfield(S, 'method'),       S.method = 'fieldtrip';     end

%-Provide CTF's name of gradiometry mode
%--------------------------------------------------------------------------
switch S.gradient
    case 0
        desired = 'none';
    case 1
        desired = 'G1BR';
    case 2
        desired = 'G2BR';
    case 3
        desired = 'G3BR';
    otherwise
        error('Gradiometer order can only range from 0-3')
end

%-Sanity checks
%--------------------------------------------------------------------------
if ischar(S.D)
    S.D = spm_eeg_load(S.D);
end

errorMsg = 'MEEG object must contain CTF data!';
if isempty(strfind(S.D.sensors('meg').type,'ctf'))
    error(errorMsg);
end

errorMsg = ['Target gradiometer ' desired 'missing from MEEG data'];
if ~S.gradient && ~isfield(S.D.sensors('meg').balance,desired)
    error(errorMsg);
end

%-Gradiometry denoising
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Current order',...
    upper(S.D.sensors('meg').balance.current)); 
fprintf('%-40s: %30s\n','Target order',...
    upper(desired)); 

switch S.method
    case 'spm'
        error('gradiometry method %s not yet supported',S.method);
    case {'ft','fieldtrip'}
        
        raw             = ftraw(S.D);
        opt             = [];
        opt.gradient    = desired;
        post            = ft_denoise_synthetic(opt,raw);
        clear raw
        
        S.D = S.D.sensors('meg',post.grad);
        res = zeros(size(S.D));
        
        for ii = 1:size(S.D,3)
            res(:,:,ii) = post.trial{ii};
        end
        
        inFile  = fnamedat(S.D);
        [a,b]   = fileparts(inFile);
        outfile = fullfile(a,[S.prefix b '.dat']);
        D = clone(S.D,outfile);
        D(:,:,:) = res;
        
    otherwise
        error('gradiometry method %s not understood',S.method);
end

%-Cleanup
%--------------------------------------------------------------------------
hiSt = struct;
[a, b, ~] = fileparts(fnamedat(S.D));
hiSt.D = [a filesep b '.mat'];
hiSt.gradient = S.gradient;
hiSt.prefix = S.prefix;
hiSt.method = S.method;
D = D.history('spm_eeg_ctf_synth_gradiometer',hiSt);
save(D);

fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#