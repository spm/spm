function model = fil_train(data,sett,model)
% Fit the patch-wise CCA-like model.
% FORMAT model = fil_train(data,sett,model)
% data  - a data structure encoding the images used, as well as the
%         amount of jitter etc.
% sett  - a data structure encoding settings.  Fields used are (with suggested values):
%         K       - Number of components to use in the model             [it depends]
%         nit     - Number of inner iterations for updating mu, W & Z    [5]
%         nu0     - Wishart degrees of freedom: A ~ W(I v_0 \nu_0, nu_0) [2]
%         v0      - Wishart scale parameter:    A ~ W(I v_0 \nu_0, nu_0) [6.0]
%         d1      - Patch-size (currently same in all directions)        [4]
%         r       - search radius                                        [2 voxels]
%         sd      - Standard deviation of weights within search radius   [0.75 voxels]
%         nit0    - Outer iterations                                     [8]
%         matname - filename for saving model                            [a string]
%         workers - Number of workers in parfor                          [it depends]
% model - the estimated model
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


io   = fil_io;                      % Handles to I/O functions
dat  = io.init(data{:});            % Set up data for I/O
dw   = dat.dw;                      % Image dimensions
Mw   = dat.Mw;                      % Voxel-to-world
r    = sett.r;                      % search radius
Nsub = size(dat.chan(1).fname, 1);  % Number of subjects
p0   = ones([Nsub,1],'single');     % Weights for each subject
%p0   = 1 - [0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0]'*0.5; % Full MICCAI
p    = get_weights(p0, r, sett.sd); % Weighting (when doing search)
ind  = true(numel(p), 2);           % Mask indicating which subjects have which data


% Enter some default settings
sett0 = struct('K',24,  'd1',4,  'nit',5,  'nit0',4, ...
               'nu0',NaN,  'v0',1,  'b0',1, ...
               'r',0,  'sd',3, ...
               'workers',0,  'matname','fil_blah.mat',  'verb',0);
if nargin<2 || ~isstruct(sett)
    sett = sett0;
else
    fns = fieldnames(sett0);
    for i=1:numel(fns)
        if ~isfield(sett,fns{i})
            sett.(fns{i}) = sett0.(fns{i});
        end
    end
end


if nargin<3
    % Set up the offsets defining the patches
    d1      = sett.d1;
%   offsets = {1:d1:dw(1), 1:d1:dw(2),70}; % Single slice for testing
    offsets = {1:d1:dw(1), 1:d1:dw(2), 1:d1:dw(3)};

    % Set up a data structure to hold the results
    c       = cell(cellfun(@numel,offsets));
    model   = struct('pos',c,'c',c,'mod',c,'P11',[],'P12',[]);

    for p3=1:numel(offsets{3})
        k   = offsets{3}(p3);
        pos = {[],[],k:min(k+d1-1,dw(3))};
        for p2=1:numel(offsets{2})
            j  = offsets{2}(p2);
            pos{2} = j:min(j+d1-1,dw(2));
            for p1=1:numel(offsets{1})
                i      = offsets{1}(p1);
                pos{1} = i:min(i+d1-1,dw(1));
                model(p1,p2,p3).pos = pos;
            end
        end
    end
end
dp = [size(model,1), size(model,2), size(model,3)];  % Dimensions of patch array

% Set up .mat files for latent variables
latent_names  = cell(dp(3),1);                       % .mat filenames 
[pth,nam,~] = fileparts(sett.matname);
latent = struct('Z',cell(dp(1:2)),'V',cell(dp(1:2)));
for p3=1:dp(3)
    latent_names{p3} = fullfile(pth,sprintf('%s_latent%.3d.mat', nam, p3));
    if ~exist(latent_names{p3},'file')
        latent = struct('Z',cell(dp(1:2)),'V',cell(dp(1:2)));
        [latent(:).Z] = deal(zeros([sett.K, size(ind,1)],'single'));
        [latent(:).V] = deal(zeros([sett.K, sett.K],'single'));
        save(latent_names{p3},'latent');
    end
end

% Updating uses a red-black scheme, so this function is used to identify which
% patches to update
to_update = @(p1,p2,p3,it)(((rem(p1,2)==rem(p2,2))==rem(p3,2))==rem(it,2));

for it=1:(2*sett.nit0)   % Twice the official number because of red/black scheme
    fprintf('%3d-%d:', floor((it+1)/2), rem(it-1,2)+1);

    latent3    = cell(1,3);
    tmp        = load(latent_names{1},'latent');
    latent3{3} = tmp.latent;

    for p3=1:dp(3)
        fprintf(' %d', p3);

        % Read in a new slab of data, accounting for the search radius
        cr  = ceil(r);
        dat = io.block(dat,(1-cr):(dw(1)+cr),...
                           (1-cr):(dw(2)+cr),...
                           (min(model(1,1,p3).pos{3})-cr):(max(model(1,1,p3).pos{3})+cr));

        % Move latent variables for current slab to previous, next slab to current and
        % read new latent variables for the next slab.
        latent3{1} = latent3{2};
        latent3{2} = latent3{3};
        if p3<dp(3)
            tmp        = load(latent_names{p3+1},'latent');
            latent3{3} = tmp.latent;
        else
            latent3{3} = [];
        end

        for p2=1:dp(2)
            % Collect things together for running a parfor. Can't run
            % parfor over everything because sharing large data structures
            % across nodes is expensive.
            Fs      = cell(1,size(model,1),1);
            Js      = cell(1,size(model,1),1);
            Cs      = cell(1,size(model,1),1);
            patches = cell(1,size(model,1),1);
            Zs      = cell(1,size(model,1),1);
            Vs      = cell(1,size(model,1),1);
            Z2s     = cell(1,size(model,1),1);
            V2s     = cell(1,size(model,1),1);
            for p1=1:dp(1)
                if to_update(p1,p2,p3,it)
                    patches{p1}       = model(p1,p2,p3);
                    Zs{p1}            = latent3{2}(p1,p2).Z;
                    Vs{p1}            = latent3{2}(p1,p2).V;
                    [Fs{p1},Js{p1},Cs{p1}] = io.patch(dat, patches{p1}.pos{:},r);
                    [Z2s{p1},V2s{p1}] = get_neighbours_latent(latent3,p1,p2);
                    if isempty(patches{p1}.mod), patches{p1}.c = Cs{p1}; end
                end
            end

            if sett.workers==0
                % Use a for loop instead of parfor
                for p1=1:numel(patches)
                    if to_update(p1,p2,p3,it)
                        [patches{p1},Zs{p1},Vs{p1}] = update_node(patches{p1},Zs{p1},p,Vs{p1},...
                                                                  Fs{p1},Js{p1},ind,Z2s{p1},V2s{p1},sett);
                    end
                end
            else
                % Run the parfor on the collections of stuff
                parfor(p1=1:numel(patches), sett.workers)
                    if feval(to_update,p1,p2,p3,it)
                        [patches{p1},Zs{p1},Vs{p1}] = update_node(patches{p1},Zs{p1},p,Vs{p1},...
                                                                  Fs{p1},Js{p1},ind,Z2s{p1},V2s{p1},sett);
                    end
                end
            end

            % Save the parfor results back to the model
            for p1=1:dp(1)
                if to_update(p1,p2,p3,it)
                    model(p1,p2,p3)     = patches{p1};
                    latent3{2}(p1,p2).Z = Zs{p1};
                    latent3{2}(p1,p2).V = Vs{p1};
                end
            end
            if ~rem(p2-1,8), fprintf('.'); end
        end

        latent = latent3{2};
        save(latent_names{p3},'latent');
    end
    fprintf('\n');

    if ~rem(it,4), model = fil_prune(model,sett,p); end % Prune irrelevant components
    if it>=2, model = fil_prec(model,sett,p); end       % Attach precision matrices
    save(sett.matname, 'dw','Mw','model','sett','-v7.3');
end



function p = get_weights(p0,r,sd)
% Set up weighting for augmented data
% FORMAT p = get_weights(p0,r,sd)
% p0 - original weight for each subject
% r  - search radius
% sd - standard deviation of weighting (voxels)
% p  - weighting for augmented data
if nargin<2, r  = 0;    end
if nargin<3, sd = 0.75; end 
rr  = -ceil(r):ceil(r);
[gx,gy,gz] = ndgrid(rr,rr,rr);
r2  = gx.^2+gy.^2+gz.^2;
msk = r2<=r.^2;
r2  = r2(msk);
v   = abs(sd)^2+eps;
G   = exp(-r2./(2*v));
G   = G/sum(G(:));
p   = kron(G(:),p0);


function [Z2,V2] = get_neighbours_latent(latent3,p1,p2)
% Collect latent variables of neighbouring patches
% FORMAT [Z2,V2] = get_neighbours_latent(latent3,p1,p2)
% latent3 - data structure containing latent variables
%           for previous, current and next slabs
% p1,p2   - Within slab indices of patch
% Z2      - E[Z] of neighbouring patches
% V2      - \sum Var[Z] of the neighbouting patches (with weighting) 
Z2  = [GetZ(p1  ,p2  ,latent3{3})
       GetZ(p1  ,p2  ,latent3{1})
       GetZ(p1  ,p2+1,latent3{2})
       GetZ(p1  ,p2-1,latent3{2})
       GetZ(p1+1,p2  ,latent3{2})
       GetZ(p1-1,p2  ,latent3{2})];
V2  = blkdiag(GetV(p1  ,p2  ,latent3{3}),...
              GetV(p1  ,p2  ,latent3{1}),...
              GetV(p1  ,p2+1,latent3{2}),...
              GetV(p1  ,p2-1,latent3{2}),...
              GetV(p1+1,p2  ,latent3{2}),...
              GetV(p1-1,p2  ,latent3{2}));



function Z = GetZ(p1,p2,latent)
% Extract the E[Z] values 
if p1>=1 && p1<=size(latent,1) && ...
   p2>=1 && p2<=size(latent,2) && ...
   ~isempty(latent) && ~isempty(latent(p1,p2).Z)
    Z = latent(p1,p2).Z;
else
    Z = [];
end



function V = GetV(p1,p2,latent)
% Extract the \sum Var[Z] values
if p1>=1 && p1<=size(latent,1) && ...
   p2>=1 && p2<=size(latent,2) && ...
   ~isempty(latent) && ~isempty(latent(p1,p2).Z)
    V = latent(p1,p2).V;
else
    V = [];
end



function [model,Z,V] = update_node(model,Z,p,V,F,J,ind,Z2,V2,sett)
% Update the model and latent variables for a patch
% FORMAT [model,Z,V] = update_node(model,Z,p,V,F,J,ind,Z2,V2,sett)
% model - structure array containing fields mu and W
% Z     - expectations of latent variables
% p     - column vector of weightings for augmented data
% V     - sum of variances of latent variables
% F     - cell array of 3D data, dimensions Nvox x M x N
%         Nvox - number of voxels in patch
%         M    - number of categories minus 1
%         N    - number of images (after augmentation)
% J     - Cell array of Jacobian weights, dimensions Nvox x N
% ind   - logical array of N x number of cells in F and J
% Z2    - expectations of latent variables in neighbouring patches
% V2    - Sum of variances of latent variables in neighbouring patches
% sett  - various settings

if isempty(model.mod)
    % Just initialise the fitting
    [model.mod,Z,V] = fil_fit(F,J,sett,ind,p);
else
    if isempty(Z2)
        % No information about neighbouring latent variables
        [model.mod,Z,V] = fil_fit(F,J,sett,ind,p,model.mod,Z);
    else
        % Use information about neighbours to give a prior probability
        % for Z conditional on the values of the neighbours. This is based on
        % all latent variables coming from a zero-mean multivariate Gaussian:
        %    p(z) = N(z0,inv(P11)).
        % See (e.g.) eq. 2.69 of Bishop's PRML book.

        % Expectation of Z*Z' over the central patch and the 6 neighbouring
        % patches
        EZZ = [Z *bsxfun(@times,p,Z')+V, Z *bsxfun(@times,p,Z2')
               Z2*bsxfun(@times,p,Z'),   Z2*bsxfun(@times,p,Z2')+V2];

        % Various dimensions
        K   = size(Z,1);
        K2  = size(Z2,1);
        Ns  = sum(p);

        % Expectation of precision matrix, drawn from a Wishart distribution
        nu0 = sett.nu0;
        if ~isfinite(nu0), nu0 = size(Z,1)+size(Z2,1)-0.999; end
        P   = inv(EZZ + (nu0*sett.v0)*eye(K+K2))*(Ns+nu0);

        % Determine distribution of the central patch Z ~ N(Z0,inv(P11))
        % This is used as a prior when updating the CCA-like fitting
        P11 = P(1:K,  1:K    );
        P12 = P(1:K, (1:K2)+K);
        Z0  = -P11\P12*Z2;

        % Fit the latent variable model again, and update the data structure
        [model.mod,Z,V] = fil_fit(F,J,sett,ind,p,model.mod,Z,Z0,P11);
    end
end
