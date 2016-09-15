function [GCM,gen] = spm_dcm_simulate(GCM, mode, noise, gen_idx)
% Populate the given group DCM array (GCM) with simulated data. If each
% subject has M models, any one of these M can be chosen to be the
% generative model, and all models for the same subject will be assigned
% the same simulated data.
%
% GCM  - subjects x model cell array where the Ep structure contains
%        connection strengths
%
% mode - zero-mean Gaussian noise is added, defined by one of:
%        'SNR_var' - signal-to-noise ratio based on the variance [default]
%        'SNR_std' - signal-to-noise ratio based on the standard deviation
%        'var'     - variance of the observation noise to be added
%
% noise - real-valued added noise (interpretation depends on mode, above)
%
% gen_idx - index of the generative model
%
% Returns:
%
% GCM  - DCM array populated with simulated data
% gen  - vector of generative models for each subject
%
% Example:
% DCM = spm_dcm_simulate(GCM, 'SNR_std', 1);
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman, Vladimir Litvak
% $Id: spm_dcm_simulate.m 6874 2016-09-15 09:29:34Z peter $

% Check parameters and load specified DCM
%--------------------------------------------------------------------------
GCM = spm_dcm_load(GCM);

if nargin < 2 || isempty(mode),    mode = 'SNR_var'; end
if nargin < 3 || isempty(noise),   noise = 1;        end
if nargin < 4 || isempty(gen_idx), gen_idx = 1;      end

model = spm_dcm_identify(GCM{1,1});

switch model
    case 'fMRI'
        [GCM, gen] = simulate_fmri(GCM, mode, noise, gen_idx);
    case 'ERP'
        [GCM, gen] = simulate_erp(GCM, mode, noise, gen_idx);
    otherwise
        error('spm_dcm_simulate not yet implemented for this type of DCM');
end
end
%--------------------------------------------------------------------------
function [GCM, gen] = simulate_fmri(GCM, mode, noise, gen_idx)
% Simulate determinstic DCM for fMRI (task-based)

[ns, nm] = size(GCM);

gen = cell(ns,1);

% Create spurious ROI information for compatibility
for i=1:GCM{1}.n
    str         = sprintf('R%d',i);
    xY(i).name  = str;
    xY(i).xyz   = [i i i]'*10;
    xY(i).XYZmm = [i i i]'*10;
    xY(i).s     = 1;
    xY(i).spec  = 1;
    xY(i).Sess  = 1;
    xY(i).u     = 1;
    xY(i).X0    = [];
end

graphics = false;

for s = 1:ns
    DCM = GCM{s,gen_idx};
    
    % Generate data
    switch lower(mode)
        case 'snr_std'
            SNR = noise;
            [Y,x,DCM_gen] = spm_dcm_generate(DCM,SNR,graphics);
        case 'snr_var'
            SNR = sqrt(noise);
            [Y,x,DCM_gen] = spm_dcm_generate(DCM,SNR,graphics);
        case 'var'
            [Y,x,DCM_gen] = spm_dcm_generate(DCM,Inf,graphics);
            
            e = sqrt(noise) .* randn(size(DCM_gen.Y.y));
            DCM_gen.Y.y = DCM_gen.y + e;
        otherwise
            error('Unknown noise definition');
    end
    
    
    for i = 1:nm
        % Store simulated timeseries
        GCM{s,i}.Y = DCM_gen.Y;
        
        % Store ROI structure
        if ~isfield(GCM{s,i},'xY')
            GCM{s,i}.xY = xY;
        end
    end
    
    % Store generative model
    gen{s} = DCM_gen;
end
end

%--------------------------------------------------------------------------
function [GCM, gen] = simulate_erp(GCM, mode, noise, gen_idx)
% Simulate determinstic DCM for ERP

[ns, nm] = size(GCM);

gen = cell(ns,1);

for i = 1:ns
    DCM = GCM{i,gen_idx};     
     
    DCM.options.DATA = 0;
    
    if isfield(DCM, 'Eg')
        Eg = DCM.Eg;
    elseif isfield(GCM{i, 1}, 'Eg')
        Eg = GCM{i, 1}.Eg;
        
        M = DCM.M;
        DCM.M = GCM{i, 1}.M;
        
        fields = fieldnames(M);
        for j = 1:numel(fields)
            DCM.M.(fields{j}) = M.(fields{j});
        end       
    else
        error('Full inverted model expected in the first column of the input.');
    end       
    
    DCM   = spm_dcm_erp_dipfit(DCM,1);
    
    % report
    %----------------------------------------------------------------------
    fprintf('Creating subject %i\n',i)
    
    
    % generate data
    %----------------------------------------------------------------------
    G   = feval(DCM.M.G, Eg, DCM.M);
    x   = feval(DCM.M.IS, DCM.Ep, DCM.M, DCM.xU);
    for c = 1:length(x)
        y{c} = x{c}*G';
        e    = spm_conv(randn(size(y{c})),8,0);
                
        switch lower(mode)
            case 'snr_std'
                e    = e*noise*mean(std(y{c}))/mean(std(e));
            case 'snr_var'
                e    = e*noise*mean(var(y{c}))/mean(var(e));                
            case 'var'
                e    = e*sqrt(noise);
            otherwise
                error('Unknown noise definition');
        end
         
        y{c} = y{c} + e;
        y{c} = DCM.M.R*y{c};
    end
    
    % specify models
    %----------------------------------------------------------------------
    for j = 1:nm
        GCM{i,j}          = rmfield(DCM,'M');
        GCM{i,j}.M.dipfit = DCM.M.dipfit;
        GCM{i,j}.xY.y     = y;
    end
    
    
    % Store generative model
    gen{i} = DCM;
end
end
