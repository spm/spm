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
%        'SNR_var' - signal-to-noise ratio based on the variance
%        'SNR_std' - signal-to-noise ratio based on the standard deviation
%        'var'     - variance of the observation noise to be added
%        'Ce'      - picks up the log noise precision from GCM{x}.Ce
%                   [default]
%
% noise - real-valued added noise (interpretation depends on mode, above)
%         if mode is set to 'hE' then this can be empty
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
% $Id: spm_dcm_simulate.m 7069 2017-04-24 10:15:20Z peter $

% Check parameters and load specified DCM
%--------------------------------------------------------------------------
GCM = spm_dcm_load(GCM);

if nargin < 2 || isempty(mode),    mode = 'ce';    end
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
        case 'ce'
            [Y,x,DCM_gen] = spm_dcm_generate(DCM,Inf,graphics);
            
            e = randn(size(DCM_gen.Y.y)) * sqrtm(diag(DCM.Ce));
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
        Ce = DCM.Ce;
    elseif isfield(GCM{i, 1}, 'Eg')
        Eg = GCM{i, 1}.Eg;
        Ce = GCM{i, 1}.Ce;
        
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
    U   = DCM.M.U;
    R   = DCM.xY.R;
    x   = feval(DCM.M.IS, DCM.Ep, DCM.M, DCM.xU);    
    for c = 1:length(x)
        y{c} = x{c}*G'*U;
                
        e    = randn(size(y{c}));
                
        switch lower(mode)
            case 'snr_std'
                e    = spm_conv(e,8,0);
                e    = e*noise*mean(std(y{c}))/mean(std(e));
            case 'snr_var'
                e    = spm_conv(e,8,0);
                e    = e*noise*mean(var(y{c}))/mean(var(e));                
            case 'var'
                e    = spm_conv(e,8,0);
                e    = e*sqrt(noise);
            case 'estimated'
                e = sqrtm(full(Ce))*e;
            otherwise
                error('Unknown noise definition');
        end
         
        y{c} = y{c} + e;
        y{c} = R*y{c}*spm_pinv(full(U));
    end
    
    % specify models
    %----------------------------------------------------------------------
    for j = 1:nm
        GCM{i,j}          = rmfield(DCM,'M');
        GCM{i,j}.M.dipfit = DCM.M.dipfit;
        GCM{i,j}.M.U      = U;
        GCM{i,j}.xY.y     = y;
    end    
    
    % Store generative model
    gen{i} = DCM;
end
end
