function res = bf_inverse_ebb(BF, S)
% Computes Empirical Bayes Beamformer filters
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2020-2023 Wellcome Centre for Human Neuroimaging


if nargin == 0
    
    keeplf          = cfg_menu;
    keeplf.tag      = 'keeplf';
    keeplf.name     = 'Keep oriented leadfields';
    keeplf.labels   = {'yes', 'no'};
    keeplf.values   = {true, false};
    keeplf.val      = {false};
    
    corr            = cfg_menu;
    corr.tag        = 'corr';
    corr.name       = 'Correlated & homologous sources';
    corr.help       = {['Prior matrix is modified so to account for power in a location '...
        'and its correlated partner, allowing for correlated sources normally suppressed '...
        'by beamformers to be reconstructed. Correlated pairs can be definied in a matrix'...
        '(see below) or automatically guessed by looking for the mirror along the saggital plane.']};
    corr.labels     = {'yes','no'};
    corr.values     = {true, false};
    corr.val        = {false};
    
    onlycorr            = cfg_menu;
    onlycorr.tag        = 'onlycorr';
    onlycorr.name       = 'Only correlated sources';
    onlycorr.help       = {['Specifies whether the correlated prior should be considered on their '...
        'own (yes) or in combination with the uncorrelated priors (no).']};
    onlycorr.labels     = {'yes','no'};
    onlycorr.values     = {true, false};
    onlycorr.val        = {false};
    
    diags            = cfg_menu;
    diags.tag        = 'diags';
    diags.name       = '(On/Off) diagonal correlated priors';
    diags.help       = {['Sets the correlated pairs to be on the diagonal (default)'...
        'of the prior matrix or off the diagonal, which is closer '...
        'to how MSP behaves when specifying correlated priors. or both simultaneously']};
    diags.labels     = {'on','off','both'};
    diags.values     = {'on','off','both'};
    diags.val        = {'on'};
    
    mixmethod        = cfg_menu;
    mixmethod.tag    = 'mixmethod';
    mixmethod.name   = 'Prior combination method';
    mixmethod.help   = {['How should we combine the correlated and uncorrelated '...
        'priors (assuming we have both) +SUM: simple addition of the two '...
        'priors to make one matrix. +REML: Two seperate priors, where the ReML'...
        'optimisation will scale them automatically. (Default: sum)']};
    mixmethod.labels = {'sum','reml'};
    mixmethod.values = {'sum','reml'};
    mixmethod.val    = {'sum'};
    
    pairs           = cfg_files;
    pairs.tag       = 'pairs';
    pairs.name      = 'Matrix of correlated pairs';
    pairs.filter    = 'mat';
    pairs.num       = [1 1];
    pairs.val       = {''};
    pairs.help      = {['[OPTIONAL] BF.mat file containing a binary adjacency matrix of correlated pairs. '...
        'If a matrix is not supplied, it will automatically look for a the homologous regions.'...
        'TIP: If you want a source to not be correlated, pair it with itself (i.e. put a 1 on the diaconal element.']};
    
    iid             = cfg_menu;
    iid.tag         = 'iid';
    iid.name        = 'Identity source covariance';
    iid.help        = {['Assumes sources are indepedent and identically distributed, equivalent to a Bayesian '...
        'minimum norm estimation. This option bypasses correlated source mode']};
    iid.labels      = {'yes','no'};
    iid.values      = {true, false};
    iid.val         = {false};
    
    noise           = cfg_files;
    noise.tag       = 'noise';
    noise.name      = 'BF mat file containing noise matrix';
    noise.filter    = 'mat';
    noise.num       = [1 1];
    noise.val       = {''};
    noise.help      = {'[OPTIONAL] BF.mat file containing empty-room noise matrix (pre-filtered)'};
    
    reml            = cfg_menu;
    reml.tag        = 'reml';
    reml.name       = 'Hyperprior optimisation';
    reml.help       = {['Model fitting via ReML. Select loose for (relatively) unrestricted sweep of the hyperpriors '...
        'whereas strict mode forces the them to follow a rule of the noise being ~1/100th the magnitide of the sources']};
    reml.labels     = {'Loose','Strict'};
    reml.values     = {'Loose','Strict'};
    reml.val        = {'Loose'};
    
    ebb      = cfg_branch;
    ebb.tag  = 'ebb';
    ebb.name = 'EBB';
    ebb.val  = {keeplf,iid,corr,onlycorr,diags,mixmethod,pairs,reml,noise};
    res = ebb;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

res = [];

% Check for missing options and set defaults based on partial inputs
%----------------------------------------------------------------
if ~isfield(S,'keeplf'),        S.keeplf = false;       end
if ~isfield(S,'iid'),           S.iid = false;          end
if ~isfield(S,'reml'),          S.reml = 'loose';       end
if ~isfield(S,'corr'),          S.corr = false;         end
if ~isfield(S,'onlycorr'),      S.onlycorr = false;     end
if ~isfield(S,'diags'),         S.diags = 'on';         end
if ~isfield(S,'mixmethod'),     S.mixmethod = 'sum';    end
if ~isfield(S,'pairs'),         S.pairs = [];           end
if ~isfield(S,'noise'),         S.noise = [];           end

switch S.diags
    case 'on'
        ondiag = true;
        offdiag = false;
    case 'off'
        ondiag = false;
        offdiag = true;
    case 'both'
        ondiag = true;
        offdiag = true;
end

% If onlycorr is true make sure corr is a corr is also true
if S.onlycorr
    S.corr = true;
end

% Inital setup
%----------------------------------------------------------------

C       = BF.features.(S.modality).C;
invCy   = BF.features.(S.modality).Cinv;
U       = BF.features.(S.modality).U;
reduce_rank = BF.sources.reduce_rank.(S.modality(1:3));

% Look for number of samples used to generate covariance, important for
% model evidence results.
Nn      = BF.features.(S.modality).N;
% Check to see that Nn is that a single value
if numel(Nn) > 1
    warning(['multiple sample values from features detected, '...
        'checking for unique solution']);
    Nn = unique(Nn);
    if numel(Nn) ~=1
        error(['cannot determine number of samples in covariance, '...
            'please change your DAiSS Covariance method!'])
    end
end

% Check to see if tdcov has been used, warn model evidence scores may not
% behave as expected
ntrials = size(BF.data.D,3);
if Nn > ntrials*16 % maxmium possible samples using tdcov
    warning(['covariance matrix not generated using tdcov - '...
        'model evidence values may not scale as expected!'])
    Nn = 1;
end

L       = S.L;
UL      = cell(size(L));
nvert   = numel(S.L);

pow     = zeros(1,nvert);
pow_dual = sparse(nvert,nvert);

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
        
        % Correlated pairs matrix processing
        %------------------------------------------------------------------
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
            %             % need to check if its symmetric - and make this fail if so.
            %             % (except for the condition of an indetity matrix)
            %             tmp = pairs - pairs';
            %             if sum(abs(tmp(:)))~=0
            %                 fprintf('WARNING: pairs matrix is symmetrical, taking bottom triangle')
            %                 pairs = tril(pairs);
            %             end
            %             % now need to check if pairs are row-wise or columnwise and fix
            %             % if its not what we are meant to be expecting
            %             if numel(unique(sum(pairs,2))) > 1
            %                 pairs = pairs';
            %                 assert(numel(unique(sum(pairs,2)))==1,'you can only have one correlated source per source')
            %             end
            
        end
        count       = sparse(nvert,nvert);
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
                    % If using the pairs matrix, there are three possible
                    % outcomes:
                    % 1) Nothing to be found to corellate to
                    % 2) It finds itself (useful if we need a contolled mix
                    % of correlated and uncorrelated sources).
                    % 3) Finds other source(s) to correlate with.
                    id = find(pairs(i,:));
                    if numel(id) > 1
                        % check if itself has appear in the case there is
                        % more than one hit and eleminiate it
                        auto = find(id==i);
                        if ~isempty(auto)
                            id(auto) = [];
                        end
                    end
                end
                
                % only do this if there is a correlated source.
                if ~isempty(id)
                    nhits = numel(id);
                    % pool lead fields together and calculate power
                    lfpool = UL{i};
                    for j = 1:nhits
                        lfpool = lfpool + UL{id(j)};
                    end
                    pow2_norm    = 1/(lfpool'*lfpool);
                    pow2_pow     = 1/(lfpool'*invCy*lfpool);
                    pow2_tmp     = pow2_pow/pow2_norm;
                    
                    % allocate to corresponding locations in array
                    % count how many times a point comes up (should be twice but
                    % could be less/more depeding on asymmetry in source space).
                    if ondiag
                        pow_dual(i,i) = pow_dual(i,i) + pow2_tmp;
                        count(i,i) = count(i,i) + 1;
                        for j = 1:nhits
                            count(id(j),id(j)) = count(id(j),id(j)) + 1;
                            pow_dual(id(j),id(j)) = pow_dual(id(j),id(j)) + pow2_tmp;
                        end
                    end
                    if offdiag
                        for j = 1:nhits
                            pow_dual(i,id(j)) = pow_dual(i,id(j)) + pow2_tmp;
                            pow_dual(id(j),i) = pow_dual(id(j),i) + pow2_tmp;
                            count(i,id(j)) = count(i,id(j)) + 1;
                            count(id(j),i) = count(id(j),i) + 1;
                        end
                    end
                end
                
                
                if ismember(i, Ibar)
                    spm_progress_bar('Set', i); drawnow;
                end
                
            end
        end
        
        % correct power for the fact that a) having two sets of lead fields
        % doubles the power and b) some locations may have been visted more
        % than once,stops this skewing the covariance too much in the
        % correlated sources favour.
        idx = find(count);
        pow_dual(idx) = pow_dual(idx)./count(idx);
        
    end
    
    spm_progress_bar('Clear');
    
end

switch S.mixmethod
    case 'sum'
        
        if S.onlycorr
            Qpu{1} = pow_dual;
            v2 = Qpu{1};
        elseif S.corr
            Qpu{1} = diag(pow);
            Qpu{2} = pow_dual;
            v2 = Qpu{1} + Qpu{2};
        else % uncorrelated
            Qpu{1} = diag(pow);
            v2 = Qpu{1};
        end
        % scale so maximum is 1
        Qp{1} = v2./max(v2(:));
        
    case 'reml'
        % scale so maximum is 1
        Qpu{1} = diag(pow);
        Qp{1} = diag(pow)./max(pow(:));
        Qpu{2} = pow_dual;
        Qp{2} = pow_dual./max(pow_dual(:));
end

% Make sparse to save memory and disk space when saving results (estimated
% 500 MB saving per prior)
for ii = 1:numel(Qp)
    Qp{ii} = sparse(Qp{ii});
    LQpL{ii} = cell2mat(UL)*Qp{ii}*cell2mat(UL)';
end

for ii = 1:numel(Qpu)
    Qpu{ii} = sparse(Qpu{ii});
end

% Prepare noise covariance
%------------------------------------------------------------------
% Step 1: Default IID noise, trace normalised.
UU      = U'*U;
Qe{1}   = UU./trace(UU);

hP(1) = -5;     % assumes IID noise is 1/100th of signal.
hC(1) = 1e-64;

% Step 2: optional room noise: WARNING this hasnt been tested for
% compatibilty with multiple source priors yet!
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

% Finally add the final hyperparameters for source model;
for ii = 1:numel(Qp)
    hP(end+1) = 0;
    hC(end+1) = 16;
end

% ReML to optimise the weighted combination of covariances to
% best match the original sensor covariance. Comes in two flavours,
% STRICT: where hyperpriors and precisions have been predfined earlier
% LOOSE: hyperpriors have free reign, similar to SPM's inversions
%------------------------------------------------------------------
switch lower(S.reml)
    case 'strict'
        fprintf('Using ReML: Strict hyperprior settings\n');
        [Cy,h,~,F,Fa,Fc] = bf_reml_sc(C,[],[Qe LQpL],Nn,hP,diag(hC));
    case 'loose'
        fprintf('Using ReML: Loose hyperprior settings\n');
        % Need to add a final extra term here to allow ReML to not run into
        % trouble, a fixed (co)variance componenent which is ~1/100 the
        % magnitude of the sensor covariance.
        Q0          = exp(-5)*trace(C)*Qe{1};
        [Cy,h,~,F,Fa,Fc]= bf_reml_sc(C,[],[Qe LQpL],Nn,-4,16,Q0);
end

% invC_reml = pinv_plus(full(C_reml));

% Extract priors
Ne    = length(Qe);
Np    = length(Qp);
LCp   = sparse(0);
hp    = h(Ne + (1:Np));
for j = 1:Np
    LQp{j}   = cell2mat(UL)*Qp{j};
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

% Calculate variance explained

try
    R2 = calc_R2(BF.features.(S.modality).UY,UL,M);
catch
    R2 = 'Unknown';
end

sz = size(M);
res.W = mat2cell(M,ones(1,sz(1)),sz(2));
res.L = L;
res.F = F;
res.R2 = R2;
reml.Cy = Cy;
reml.Q = [Qe LQpL];
reml.Qtype = [repmat({'noise'},1,Ne) repmat({'source'},1,Np)];
reml.source_prior.scaled = Qp;
reml.source_prior.unscaled = Qpu;
reml.h = h;
reml.F = F;
reml.Fa = Fa;
reml.Fc = Fc;
res.reml = reml;
res.config = S;

if S.iid
    res.inv_type = 'ebb_iid';
elseif S.corr
    res.inv_type = 'ebb_corr';
else
    res.inv_type = 'ebb';
end

end

function R2 = calc_R2(Y,L,W)

J = W*Y;
SSR  = sum(var((Y - cell2mat(L)*J))); %% changed variance calculation
SST  = sum(var(Y));

R2  = 100*(SST - SSR)/SST;

end