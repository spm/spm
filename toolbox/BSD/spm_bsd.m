function BSD = spm_bsd(BSD)
% Estimate parameters of a Bayesian Spectral Decomposition model of 
% (complex) cross-spectral density
% FORMAT BSD = spm_bsd(BSD)
%
% BSD
%    name: name string
%       xY: data   [1x1 struct]
%       xU: design [1x1 struct]
%
%   Sname: cell of source name strings
%
%   options.Nmodes       - number of spatial modes
%   options.Tdcm         - [start end] time window in ms
%   options.Fdcm         - [start end] Frequency window in Hz
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.separatenull - bool, whether to fit null model separately
%   options.fitlog       - bool, fit log power spectra
%   options.powerline    - [start end] filter with low precision
%
% Esimates:
%--------------------------------------------------------------------------
% BSD.dtf                   - directed transfer functions (source space)
% BSD.ccf                   - cross covariance functions (source space)
% BSD.coh                   - cross coherence functions (source space)
% BSD.fsd                   - specific delay functions (source space)
% BSD.pst                   - peristimulus time
% BSD.Hz                    - frequency
%
% BSD.Ep                    - conditional expectation
% BSD.Cp                    - conditional covariance
% BSD.Pp                    - conditional probability
% BSD.Hc                    - conditional responses (y), channel space
% BSD.Rc                    - conditional residuals (y), channel space
% BSD.Hs                    - conditional responses (y), source space
% BSD.Ce                    - eML error covariance
% BSD.F                     - Laplace log evidence
% BSD.ID                    -  data ID
%__________________________________________________________________________
 
% Adapted from spm_dcm_csd (Karl Friston) 
%
% Johan Medrano
% Copyright (C) 2023-2025 Wellcome Centre for Human Neuroimaging
 
 
% check options
%==========================================================================
drawnow
clear spm_erp_L
name = sprintf('BSD_%s',date);
BSD.options.analysis  = 'CSD';
 
% Filename and options
%--------------------------------------------------------------------------
try BSD.name;                        catch, BSD.name = name;        end
try model   = BSD.options.model;     catch, model    = 'DEM';       end
try spatial = BSD.options.spatial;   catch, spatial  = 'LFP';       end
try Nm      = BSD.options.Nmodes;    catch, Nm       = 8;           end
try DATA    = BSD.options.DATA;      catch, DATA     = 1;           end
try SAVE    = BSD.options.SAVE;      catch, SAVE     = 1;           end
try RESULTS = BSD.options.RESULTS;   catch, RESULTS  = 1;           end

% Model options
%--------------------------------------------------------------------------
try BSD.M.Nmax = BSD.options.Nmax;   catch, BSD.M.Nmax     = 256;   end
try BSD.M.sharefreqs = BSD.options.sharefreqs; catch, BSD.M.sharefreqs = 1; end
try BSD.M.nograph = BSD.options.nograph;  catch,  end
try BSD.M.noprint = BSD.options.noprint;  catch,  end
try BSD.M.separatenull = BSD.options.separatenull; catch, BSD.M.separatenull = 1; end
try BSD.M.fitlog = BSD.options.fitlog; catch, BSD.M.fitlog = 1; end
try BSD.M.powerline = BSD.options.powerline; catch, BSD.M.powerline = []; end
try BSD.options.Fdcm;                catch, BSD.options.Fdcm = 1:64; end
% try, BSD.M.chanavg = BSD.options.chanavg; catch, BSD.M.chanavg = 0; end


try 
    BSD.Sname = cellstr(BSD.Sname);
end

% Spatial model
%==========================================================================
BSD.options.Nmodes = Nm;
BSD.M.dipfit.model = model;
if strcmpi(spatial, 'chan')
    BSD.M.dipfit.type  = 'LFP';
else
    BSD.M.dipfit.type  = spatial;
end

if DATA
    BSD  = spm_dcm_erp_data(BSD);                   % data
    
    if strcmpi(spatial, 'chan')
        Ns = size(BSD.xY.y{1},2);  
        BSD.M.dipfit.Ns = Ns; 
        BSD.M.dipfit.Nc = Ns; 
    else 
        BSD  = spm_dcm_erp_dipfit(BSD, 1);              % spatial model
        Ns   = BSD.M.dipfit.Ns; 
    end
else 
    if isfield(BSD, 'xY') 
        if isfield(BSD.xY, 'y')
           if isnumeric(BSD.xY.y)
               BSD.xY.y = {BSD.xY.y}; 
           elseif ~iscell(BSD.xY.y)
               error('Field BSD.xY.y must be a cell array of numeric arrays [nfreqs, nchannels]');
           end
        else
            error('BSD.xY.y not specified'); 
        end

        if isfield(BSD.xY, 'Hz') 
            Nf = length(BSD.xY.Hz); 
        else
            error('BSD.xY.Hz must be specified with BSD.xY.y'); 
        end

        if ~isfield(BSD.xY, 'dt')
            BSD.xY.dt = 1; 
        end
    else
        error('BSD.xY not specified'); 
    end
    
    if size(BSD.xY.y{1}, 1) ~= Nf 
        error('Number of frequencies (BSD.xY.Hz) must match the first dimension of the data (BSD.xY.y)')
    end

    Ns = size(BSD.xY.y{1},2);  
    for i = 1:numel(BSD.xY.y)
        if ~( size(BSD.xY.y{i}, 1) == Nf && size(BSD.xY.y{i}, 2) == Ns )
            error('Shape mismatch in data for condition %d: expected (%d,%d), got (%d,%d).', i, Nf, Ns, ...
                size(BSD.xY.y{i}, 1), size(BSD.xY.y{i}, 2));
        end
    end

    BSD.M.dipfit.Ns = Ns; 
    BSD.M.dipfit.Nc = Ns; 
end

%-Feature selection using principal components (U) of lead-field
%==========================================================================
 
% Spatial modes
%--------------------------------------------------------------------------
try
    BSD.M.U = spm_dcm_eeg_channelmodes(BSD.M.dipfit,Nm);
end
 
% get data-features (in reduced eigenspace)
%==========================================================================
if DATA
    BSD  = spm_bsd_data(BSD);
end
BSD.M.Hz = BSD.options.Fdcm; 

if ~(strcmpi(spatial, 'chan') && Ns == 1)
    % scale data features (to a variance of about 8)
    %--------------------------------------------------------------------------
    ccf      = spm_csd2ccf(BSD.xY.y,BSD.xY.Hz);
    scale    = max(spm_vec(ccf))/8;
    BSD.xY.y = spm_unvec(spm_vec(BSD.xY.y)/scale,BSD.xY.y);
    BSD.xY.y = cellfun(@(x) abs(x), BSD.xY.y, 'UniformOutput', 0); 
    try BSD.xY.scale; catch BSD.xY.scale = 1; end
    BSD.xY.scale = BSD.xY.scale * scale; 
end
 
if strcmpi(spatial, 'chan') && Ns > 1
    BSD0 = BSD; 
    BSD0.options.DATA = 0; 
    BSD0.options.SAVE = 0;
    BSD0.options.RESULTS = 0;
    
    xY = BSD.xY; 
    try xY = rmfield(xY, 'name'); end
    try xY = rmfield(xY, 'Ic'); end
    try xY = rmfield(xY, 'coor2D'); end
    try xY = rmfield(xY, 'U'); end
    BSD0.xY.y = cellfun(@(x) reshape(mean(x, [2, 3]), [], 1), BSD.xY.y, ...
        'UniformOutput', false);
    
    % fit average
    BSD0 = spm_bsd(BSD0); 

    BSDi = BSD; 
    BSDi.M.P = BSD0.Ep; 
    BSDi.options.DATA = 0; 
    BSDi.options.SAVE = 0;
    BSDi.options.RESULTS = 0;

    BSDs = {}; 

    fig = spm_figure('FindWin', 'Interactive'); 
    found = ~isempty(fig);
    if ~found
        fig = spm_figure('GetWin', 'Interactive'); 
    end
    spm_progress_bar('Init', Nm, 'Fitting spectra...', 'spatial modes');

    for i = 1:Nm
        fprintf('Processing mode %d/%d.', i, Nm); 
        BSDi.xY.y = cellfun(@(x)  reshape(x(:, i, i), [], 1), BSD.xY.y, ...
            'UniformOutput', false);
        BSDs{i} = spm_bsd(BSDi); 
        spm_progress_bar('Set', i);
    end
    spm_progress_bar('Clear');

    BSDs = [BSDs{:}]; 
    BSD.models = BSDs; 
    BSD.F = sum([BSDs.F]); 

    if SAVE
        save(BSD.name, 'BSD', spm_get_defaults('mat.format'));
    end

    if RESULTS
        BSD = spm_bsd_results(BSD);
    end
    return 
end
   
% Design model and exogenous inputs
%==========================================================================
if ~isfield(BSD,'xU'),   BSD.xU.X = sparse(1 ,0); end
if ~isfield(BSD.xU,'X'), BSD.xU.X = sparse(1 ,0); end
if isempty(BSD.xU.X),    BSD.xU.X = sparse(1 ,0); end

% Bayesian spectral decomposition model
%==========================================================================

% prior moments on parameters
%--------------------------------------------------------------------------
[pE, pC, pV] = spm_bsd_priors(BSD.fqs,Ns,BSD.xU.X,BSD.M); 
P = pE; 
P.a = P.a*0+2; 
P.S = P.S + log(diff(pV.f,1,2)) + 3; 


% check to see if neuronal priors have already been specified
%--------------------------------------------------------------------------
try
    if spm_length(BSD.M.pE) == spm_length(pE)
        pE = BSD.M.pE;
        pC = BSD.M.pC;
        if ~BSD.options.noprint
            fprintf('Using existing priors\n')
        end
    end
end
try
    if spm_length(BSD.M.P) == spm_length(pE)
        P = BSD.M.P;
        if ~BSD.options.noprint
            fprintf('Using existing initial parameters\n')
        end
    end
end

% augment with priors on spatial model
%--------------------------------------------------------------------------
[pE,pC] = spm_L_priors(BSD.M.dipfit,pE,pC);
pE.J = 1; 
pC.J = 0; 

try
    if spm_length(BSD.M.pE) == spm_length(pE)
        pE = BSD.M.pE;
        pC = BSD.M.pC;
        if ~BSD.options.noprint
            fprintf('Using existing priors\n')
        end
    end
end
try
    if spm_length(BSD.M.P) == spm_length(pE)
        P = BSD.M.P;
        if ~BSD.options.noprint
            fprintf('Using existing initial parameters\n')
        end
    end
end

% check for pre-specified priors
%--------------------------------------------------------------------------
try 
    hE  = BSD.M.hE;  hC  = BSD.M.hC; 
catch
    hE = 6; hC = 1/128;
end
 
% create BSD
%--------------------------------------------------------------------------
BSD.M.IS = 'spm_bsd_int';
BSD.M.pV = pV; 
BSD.M.pE = pE;
BSD.M.pC = pC;
BSD.M.hE = hE;
BSD.M.hC = hC;
BSD.M.m  = Ns;
BSD.M.P  = P;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
BSD.M.u  = sparse(Ns,1);


% complete model specification and invert
%==========================================================================
Nm       = size(BSD.M.U,2);                    % number of spatial modes
BSD.M.l  = Nm;
BSD.M.Hz = BSD.xY.Hz;
BSD.M.dt = BSD.xY.dt;
 
% normalised precision
%--------------------------------------------------------------------------
BSD.xY.Q  = spm_dcm_csd_Q(BSD.xY.y);
if ~isempty(BSD.M.powerline)
    q = reshape(BSD.xY.Hz, [], 1); 
    q = (q > BSD.M.powerline(1)) & (q < BSD.M.powerline(2)); 
    q = logical(kron(ones(numel(BSD.xY.y), 1), q));
    BSD.xY.Q(q, q) = 0; 
end

% confounds
%--------------------------------------------------------------------------
if ~isfield(BSD.xY, 'X0')
    BSD.xY.X0 = sparse(size(BSD.xY.Q,1),0);
end

% apply  log transform
%--------------------------------------------------------------------------
if BSD.M.fitlog
    BSD.xY.y = cellfun(@(x) log(x), BSD.xY.y, 'UniformOutput', 0); 
end

if BSD.M.separatenull
    % Create and invert submodel (null model, no spectral peak) 
    BSD0 = BSD; 
    BSD0.fqs = {};
    [BSD0.M.pE, BSD0.M.pC, BSD0.M.pV] = spm_bsd_priors( ...
        BSD0.fqs,Ns,BSD0.xU.X,BSD0.M, pE, pC); 
    BSD0.M.P = BSD0.M.pE;

    
    % Variational Laplace: model inversion
    %======================================================================
    [Qp,Cp,Eh,F,L] = spm_bsd_nlsi_GN(BSD0.M,BSD0.xU,BSD0.xY);

    BSD.null    = []; 
    BSD.null.M  = BSD0.M; 
    BSD.null.Ep = Qp;
    BSD.null.Cp = Cp;
    BSD.null.F = F;
    BSD.null.L = L;
    
    BSD.M.pE.b = Qp.b;
    BSD.M.pE.L = Qp.L;
    BSD.M.pE.J = Qp.J;

    BSD.M.P.b = Qp.b;
    BSD.M.P.L = Qp.L;
    BSD.M.P.J = Qp.J;

    BSD.M.pC.b = BSD.M.pC.b / 128;
    BSD.M.hE = Eh;

    yp = spm_bsd_int(Qp,BSD0.M,BSD0.xU);
    e  = spm_unvec(spm_vec(BSD.xY.y) - spm_vec(yp),yp);  
    Hz = BSD.M.Hz;
    for i = 1:size(pV.f, 1)
        pow = 0;
        for c = 1:numel(e)
            ei = e{c};
            ei = max(ei((pV.f(i,1) <= Hz') & (Hz' <= pV.f(i, 2)))); 
            pow = pow + ei; 
        end
        BSD.M.P.a(i) = pow; % in log already
    end
 end

% Variational Laplace: model inversion
%==========================================================================
[Qp,Cp,Eh,F] = spm_bsd_nlsi_GN(BSD.M,BSD.xU,BSD.xY);

% Data ID
%--------------------------------------------------------------------------
try
    try
        ID = spm_data_id(feval(BSD.M.FS,BSD.xY.y,BSD.M));
    catch
        ID = spm_data_id(feval(BSD.M.FS,BSD.xY.y));
    end
catch
    ID = spm_data_id(BSD.xY.y);
end
 
% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');

[~,pH] = spm_bsd_int(pE,BSD.M,BSD.xU); 
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
[Hc, HC]  = spm_bsd_int(Qp,BSD.M,BSD.xU);                % prediction
Ec  = spm_unvec(spm_vec(BSD.xY.y) - spm_vec(Hc),Hc);     % prediction error

qp   = Qp; 
qp.a = qp.a .* 0 - 16;                                   % remove peaks
Hp   = spm_bsd_int(qp,BSD.M,BSD.xU);                     % prediction
Yp   = spm_unvec(spm_vec(BSD.xY.y) - spm_vec(Hp),Hp);

 
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
M             = rmfield(BSD.M,'U'); 
M.dipfit.type = 'LFP';

M.U         = 1; 
M.l         = Ns;
qp          = Qp;
qp.L        = ones(1,Ns);       % set virtual electrode gain to unity
[Hs, H]     = spm_bsd_int(qp,M,BSD.xU);
BSD.Hz      = M.Hz;
BSD.H       = H; 
BSD.Hsens   = HC; 
BSD.pH      = pH; 
BSD.xY.Yp   = Yp;

 
% store estimates in BSD
%--------------------------------------------------------------------------
BSD.Ep = Qp;                   % conditional expectation
BSD.Cp = Cp;                   % conditional covariance
BSD.Pp = Pp;                   % conditional probability
BSD.Hc = Hc;                   % conditional responses (y), channel space
BSD.Rc = Ec;                   % conditional residuals (y), channel space
BSD.Hs = Hs;                   % conditional responses (y), source space
BSD.Ce = exp(-Eh);             % ReML error covariance
BSD.F  = F;                    % Laplace log evidence
BSD.ID = ID;                   % data ID

% save
%--------------------------------------------------------------------------
BSD.options.Nmodes = Nm;
 
if SAVE
    save(BSD.name, 'BSD', spm_get_defaults('mat.format'));
end

% and show results... 
%--------------------------------------------------------------------------
if RESULTS
    try
        spm_bsd_results(BSD);
    end
end

return
end