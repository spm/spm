function res = bf_inverse_ebb(BF, S)
% Computes Empirical Bayes Beamformer filters
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% George O'Neill
% $Id: bf_inverse_ebb.m 7918 2020-08-10 10:41:19Z george $

% NOTE: this is an early developmental version so it comes with George's
% "NO RESULTS GUARENTEED (TM)" warning.

%--------------------------------------------------------------------------

if nargin == 0
    
    keeplf        = cfg_menu;
    keeplf.tag    = 'keeplf';
    keeplf.name   = 'Keep oriented leadfields';
    keeplf.labels = {'yes', 'no'};
    keeplf.values = {true, false};
    keeplf.val    = {false};
    
    corr          = cfg_menu;
    corr.tag      = 'corr';
    corr.name     = 'Correlated & homologous sources';
    corr.help     = {['Prior matrix is modified so to account for power in a location '...
        'and its correlated partner, allowing for correlated sources normally suppressed '...
        'by beamformers to be reconstructed. Correlated pairs can be definied in a matrix'...
        '(see below) or automatically guessed by looking for the mirror along the saggital plane.']};
    corr.labels   = {'yes','no'};
    corr.values   = {true, false};
    corr.val      = {false};
    
    pairs = cfg_files;
    pairs.tag = 'pairs';
    pairs.name = 'Matrix of correlated pairs';
    pairs.filter = 'mat';
    pairs.num=[1 1];
    pairs.val={''};
    pairs.help = {['[OPTIONAL] BF.mat file containing a binary adjacency matrix of correlated pairs. '...
        'If a matrix is not supplied, it will automatically look for a the homologous regions.'...
        'TIP: If you want a source to not be correlated, pair it with itself (i.e. put a 1 on the diaconal element.']};
    
    iid           = cfg_menu;
    iid.tag       = 'iid';
    iid.name      = 'Identity source covariance';
    iid.help      = {['Assumes sources are indepedent and identically distributed, equivalent to a Bayesian '...
        'minimum norm estimation. This option bypasses correlated source mode']};
    iid.labels    = {'yes','no'};
    iid.values    = {true, false};
    iid.val       = {false};
    
    noise = cfg_files;
    noise.tag = 'noise';
    noise.name = 'BF mat file containing noise matrix';
    noise.filter = 'mat';
    noise.num=[1 1];
    noise.val={''};
    noise.help = {'[OPTIONAL] BF.mat file containing empty-room noise matrix (pre-filtered)'};
    
    
    reml          = cfg_menu;
    reml.tag      = 'reml';
    reml.name     = 'Hyperprior optimisation';
    reml.help     = {['Model fitting via ReML. Select loose for (relatively) unrestricted sweep of the hyperpriors '...
        'whereas strict mode forces the them to follow a rule of the noise being ~1/100th the magnitide of the sources']};
    reml.labels   = {'Loose','Strict'};
    reml.values   = {'Loose','Strict'};
    reml.val      = {'Loose'};
    
    ebb      = cfg_branch;
    ebb.tag  = 'ebb';
    ebb.name = 'EBB';
    ebb.val  = {keeplf,iid,corr,pairs,reml,noise};
    res = ebb;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

res = [];

C       = BF.features.(S.modality).C;
invCy   = BF.features.(S.modality).Cinv;
U       = BF.features.(S.modality).U;
reduce_rank = BF.sources.reduce_rank.(S.modality(1:3));

% Nn      = 1; % The covariance has already been scaled in bf_features, so Nn=1.
Nn      = BF.features.(S.modality).N;
% C       = Nn*C;
% invCy    = pinv_plus(C);

L       = S.L;
UL      = cell(size(L));
nvert   = numel(S.L);

pow     = zeros(1,nvert);
pow_dual    = zeros(1,nvert);


% Lead field optimisation (if multiple LFs per ROI are provided).
%----------------------------------------------------------------
spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert,'Preparing lead fields'); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

if size(L{1},2)~=1
    fprintf('Optimising leadfield orentations');
    for i = 1:nvert
        if ~isnan(L{i})
            
            lf    = U'*L{i};
            
            % Robert's code - to reduce the lead fields to optimal location
            % TO DO: Add support for multiple lead fields per ROI.
            [u, ~] = svd(real(pinv_plus(lf' * invCy *lf, reduce_rank, 0)),'econ');
            eta = u(:,1);
            lf  = lf * eta;
            
            % Store the smooth, reduced lead field
            UL{i} = lf;
            
            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
        end
    end
else
    % Just apply subspace projectors.
    fprintf('Preparing lead fields\n');
    for i = 1:nvert
        if ~isnan(L{i})
            UL{i}    = U'*L{i};
            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
        end
    end
end

spm_progress_bar('Clear');

% Calculate the source covariance matrix
%------------------------------------------------------------------
if S.iid
    pow(:)     = 1;
else
    spm('Pointer', 'Watch');drawnow;
    spm_progress_bar('Init', nvert,'Generating source covariance matrix'); drawnow;
    if nvert > 100, Ibar = floor(linspace(1, nvert,100));
    else Ibar = 1:nvert; end
    
    for i = 1:nvert
        if ~isnan(UL{i})
            
            lf    = UL{i};
            
            pow_norm    = 1/(lf'*lf);
            pow_pow     = 1/(lf'*invCy*lf);
            pow(i)      = pow_pow/pow_norm;
            
            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
        end
    end
    
    spm_progress_bar('Clear');
    
    % Modifying priors for correlated sources (if needed by user)
    %------------------------------------------------------------------
    if S.corr
        
        if isempty(S.pairs)
            assert(~isempty(strmatch(BF.data.space,'MNI-aligned')),['Correlated source mode must '...
                'be used with sources in MNI-aligned space, please check you options in the data module']);
            
            pos         = BF.sources.pos;
            
            
            assert(size(pos,1)==nvert,'number of sources do not correspond with number of lead fields...');
        else
            if iscell(S.pairs)
                S.pairs = cell2mat(S.pairs);
            end
            fprintf('Loading custom pairs matrix\n');
            X = load(S.pairs);
            flds = fields(X);
            % Assume the file is the only field in the structure
            assert(length(flds) == 1,['I cannot tell which variable in the pairs mat '...
                'file is the one you want me to use, please only have one in there!']);
            eval(['pairs = X.' flds{1} ';']);
            % Check it has the correct size
            assert(length(pairs)==nvert,'size of pairs mat doesnt correspond to number of sources');
            % need to check if its symmetric - and make this fail if so.
            % (except for the condition of an indetity matrix)
            tmp = pairs - pairs';
            if sum(abs(tmp(:)))~=0
                fprintf('WARNING: pairs matrix is symmetrical, taking bottom triangle')
                pairs = tril(pairs);
            end
            % now need to check if pairs are row-wise or columnwise and fix
            % if its not what we are meant to be expecting
            if numel(unique(sum(pairs,2))) > 1
                pairs = pairs';
                assert(numel(unique(sum(pairs,2)))==1,'you can only have one correlated source per source')
            end
            
        end
        count       = zeros(1,nvert);
        spm('Pointer', 'Watch');drawnow;
        spm_progress_bar('Init', nvert,'Adding correlated sources'); drawnow;
        if nvert > 100, Ibar = floor(linspace(1, nvert,100));
        else Ibar = 1:nvert; end
        
        for i = 1:nvert
            if ~isnan(L{i})
                
                if isempty(S.pairs)
                % We are looking for sources in the mirror of the saggital
                % plane, which means flipping the position in the x-axis and
                % looking for the source which is the closest (Which isnt
                % itself!)
                target  = [-pos(i,1) pos(i,2) pos(i,3)];
                del     = bsxfun(@minus,pos,target);
                ss      = dot(del',del');
                ss(i)   = 1000;
                [~,id] = min(ss);
                else
                    id = find(pairs(i,:));
                end
                
                % pool lead fields together and calculate power
                lf2 = UL{i}+UL{id};
                pow2_norm    = 1/(lf2'*lf2);
                pow2_pow     = 1/(lf2'*invCy*lf2);
                pow2_tmp     = pow2_pow/pow2_norm;
                
                % allocate to corresponding locations in array
                % count how many times a point comes up (should be twice but
                % could be less/more depeding on assymetry in source space).
                pow_dual(i) = pow_dual(i) + pow2_tmp;
                pow_dual(id) = pow_dual(id) + pow2_tmp;
                count(i) = count(i) + 1;
                count(id) = count(id) + 1;
            end
            
            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
            
        end
        
        % correct power for the fact that a) having two sets of lead fields
        % doubles the power and b) some locations may have been visted more
        % than once,stops this skewing the covariance too much in the
        % correlated sources favour.
        count(count==0) = 1;
        pow_dual = 0.5*pow_dual./count;
        
    end
    
    spm_progress_bar('Clear');
    
end

% combine
pow = pow + pow_dual;
pow = pow./max(pow); % Scale

Qp{1} = diag(pow);
LQpL{1} = cell2mat(UL)*diag(pow)*cell2mat(UL)';

% Prepare noise covariance
%------------------------------------------------------------------
% Step 1: Default IID noise, trace normalised.
UU      = U'*U;
Qe{1}   = UU./trace(UU);

hP(1) = -5;		% assumes IID noise is 1/100th of signal.
hC(1) = 1e-64;

% Step 2: optional room noise
if ~isempty(S.noise)
    
    if iscell(S.noise)
        S.noise = cell2mat(S.noise);
    end
    
    if ~exist(S.noise,'file')
        error('Noise BF file not found')
    end
    
    try
        noise = load(S.noise,'features');
    catch
        error('Features structure not found in noise BF file');
    end
    
    tmp = noise.features.(S.modality).C;
    UQeU = U'*tmp*U;
    
    Qe{end+1} = UQeU./trace(UQeU);
    
    hP(end+1) = -3; % assumes empty noise is 1/20 of signal
    hC(end+1) = 16;
end

% FInally add the final hyperparameters for source model;
hP(end+1) = 0;
hC(end+1) = 16;

% ReML to optimise the weighted combination of covariances to
% best match the original sensor covariance. Comes in two flavours,
% STRICT: where hyperpriors and precisions have been predfined earlier
% LOOSE: hyperpriors have free reign, similar to SPM's inversions
%------------------------------------------------------------------
switch lower(S.reml)
    case 'strict'
        fprintf('Using ReML: Strict hyperprior settings\n');
        [Cy,h,~,F,Fa,Fc] = bf_spm_reml_sc(C,[],[Qe LQpL],Nn,hP,diag(hC));
    case 'loose'
        fprintf('Using ReML: Loose hyperprior settings\n');
        % Need to add a final extra term here to allow ReML to not run into
        % trouble, a fixed (co)variance componenent which is ~1/100 the
        % magnitude of the sensor covariance.
        Q0          = exp(-5)*trace(C)*Qe{1};
        [Cy,h,~,F,Fa,Fc]= bf_spm_reml_sc(C,[],[Qe LQpL],Nn,-4,16,Q0);
end

% invC_reml = pinv_plus(full(C_reml));

% Extract priors
Ne    = length(Qe);
Np    = length(Qp);
Cp    = sparse(0);
LCp   = sparse(0);
LQp{1}   = cell2mat(UL)*diag(pow);
hp    = h(Ne + (1:Np));
for j = 1:Np
    Cp  =  Cp + hp(j)*Qp{j};
    LCp = LCp + hp(j)*LQp{j};
end

% MAP estimates of instantaneous sources
%======================================================================
% This is equivalent to M = Cp*UL'*inv(Qe + UL*Cp*UL'))
% with Cp the posterior source covariance (with optimal h values)
M     = LCp'/Cy;

if S.keeplf
    L = UL;
end

sz = size(M);
res.W = mat2cell(M,ones(1,sz(1)),sz(2));
res.L = L;
res.F = F;
reml.Cy = Cy;
reml.Q = [Qe LQpL];
reml.Qtype = [repmat({'noise'},1,Ne) repmat({'source'},1,Np)];
reml.source_prior = pow;
reml.h = h;
reml.F = F;
reml.Fa = Fa;
reml.Fc = Fc;
res.reml = reml;

if S.iid
    res.inv_type = 'ebb_iid';
elseif S.corr
    res.inv_type = 'ebb_corr';
else
    res.inv_type = 'ebb';
end