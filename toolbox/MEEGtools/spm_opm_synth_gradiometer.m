function D = spm_opm_synth_gradiometer(S)
% Denoise OPM data
% FORMAT D = spm_opm_synth_gradiometer(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                       - Default: no Default
%   S.confounds     - n x 1 cell array containing           - Default: REF
%                     channel types used for denoising.
%   S.derivative    - flag to denoise using derivatives     - Default: 1
%   S.gs            - flag to denoise using global signal   - Default: 0
%   S.prefix        - string prefix for output MEEG object  - Default 'd_'
% Output:
%   D               - denoised MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_synth_gradiometer.m 7414 2018-09-07 11:00:29Z spm $


%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),          error(errorMsg); end
if ~isfield(S, 'confounds'),  S.confounds = {'REF'}; end
if ~isfield(S, 'derivative'), S.derivative = 1; end
if ~isfield(S, 'gs'),         S.gs = 0; end
if ~isfield(S, 'prefix'),     S.prefix = 'd_'; end

% check to make sure dimensions match between reference and MEEG object
ndim = size(S.D);
nConfoundTypes = length(S.confounds);
refInd=[];

% ref, sensors and residual objects 
for i =1:nConfoundTypes
    refIndTemp = S.D.indchantype(S.confounds{i});
    refInd = [refInd,refIndTemp];
end

% need to loop over Ntrials of size winSize
Ntrials = size(S.D,3);
winSize = size(S.D,2);

% only want to denoise MEG channels (not trigs or references)
megind=S.D.indchantype('MEG');
megres = zeros(size(S.D(megind,:,:)));
ref = S.D(refInd,:,:);

% loop (inefficient due to continued memory reallocation but ...)
for i =1:Ntrials
     % add a mean column to the reference regressors
     intercept = ones(winSize,1);
     reference = ref(:,:,i)';
     reference = [reference intercept];
    
     % optionally add derivatives
     if(S.derivative)
        drefer = diff(reference);
        drefer = [drefer(1,:); drefer];
        reference = [drefer reference];
     end
     
    % optionally add global signal 
    if(S.gs)
        trial     = S.D(megind,:,i)';
        gsrefer   = mean(trial,2);
        reference = [gsrefer reference];
    end
    % reference is column major so transpose sensors
    beta = pinv(reference)*S.D(megind,:,i)';
    megres(:,:,i) = (S.D(megind,:,i)'- reference*beta)';
end

res = S.D(:,:,:); % new file will only have MEG sensors changed
res(megind,:,:) = megres;

% make sure output has the singleton dimension if matrix supplied
if ((length(size(res)))==2)
    outsize = [size(res) 1];
else
    outsize = size(res);
end


inFile  = fnamedat(S.D);
[a,b]   = fileparts(inFile);
outfile = fullfile(a,[S.prefix b '.dat']);
D = clone(S.D,outfile,outsize);
D(:,:,:) = res;
