function res = bf_inverse_lcmv_multicov(BF, S)
% Computes LCMV filters using spm_pca_order to constrain inverse of data
% cov matrix. Based on the paper:
% MEG beamforming using Bayesian PCA for adaptive data covariance matrix regularization.
% Woolrich M, Hunt L, Groves A, Barnes G.
% Neuroimage. 2011 Aug 15;57(4)
% 
% and allowing for multiple covariance matrices, e.g. associated with
% multiple states:
% Dynamic State Allocation for MEG Source Reconstruction
% Neuroimage. Woolrich et al. 2013.
%
% Mark Woolrich
%--------------------------------------------------------------------------
if nargin == 0      
    
    pca_order = cfg_entry;
    pca_order.tag = 'pca_order';
    pca_order.name = 'PCA order';        
    pca_order.val = {};
   
    bilateral = cfg_entry;
    bilateral.tag = 'bilateral';
    bilateral.name = 'Do bilateral beamformer';        
    bilateral.val = {};

    type        = cfg_menu;
    type.tag    = 'type';
    type.name   = 'Beamformer type';
    type.help   = {'Select Scalar or Vector'};
    type.labels = {'Vector', 'Scalar'};
    type.values = {'Vector', 'Scalar'};
    type.val    = {'Scalar'};
    
    lcmv_multicov      = cfg_branch;
    lcmv_multicov.tag  = 'lcmv_multicov';
    lcmv_multicov.name = 'LCMV use multiple covariance and PCA order';
    lcmv_multicov.val  = {pca_order,bilateral,type};      
    
    res = lcmv_multicov;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

inverse=[];

if isfield(BF.features, S.modality)

    %reduce_rank = S.reduce_rank;
    reduce_rank = BF.sources.reduce_rank.(S.modality(1:3));

    if isfield(BF.features.(S.modality),'class'),
        class = BF.features.(S.modality).class;
    else
        class={};
        class{1} = BF.features.(S.modality);
    end;

    % multiple covariances -do separate weights for each one
    for kk=1:length(class),

        if ~sum(isnan(squash(class{kk}.C))),
            [invCy, pca_order_used] = pinv_plus(class{kk}.C, S.pca_order); % rank maybe not be detected properly by just using pinv - MWW
        else
            warning(['Nans in covariance matrix for class ' num2str(kk)]);
        end;

        L = S.L;

        W = cell(size(L));

        nvert = numel(W);

        spm('Pointer', 'Watch');drawnow;
        spm_progress_bar('Init', nvert, ['Computing ' S.modality ' filters']); drawnow;
        if nvert > 100, Ibar = floor(linspace(1, nvert,100));
        else Ibar = 1:nvert; end

        for i = 1:nvert
            if ~sum(isnan(squash(L{i}))) && ~sum(isnan(squash(class{kk}.C))),
                lf = L{i};
                lf_single = lf; 

                do_lat=0;
                if S.bilateral, % beamform bilateral dipoles together
                    if(BF.sources.mni.bilateral_index(i)~=i),
                        do_lat=1;
                        lf_lat=L{BF.sources.mni.bilateral_index(i)};

                        lf=[lf_single, lf_lat];

                    end;

                end;

                switch lower(S.type)

                  case 'scalar'
                    tmp=lf' * invCy *lf;
                    nn=size(lf_single,2);

                    [u, ~] = svd(real(pinv_plus(tmp(1:nn,1:nn),reduce_rank,0)),'econ'); % this is faster,  - MWW                                                
                    eta = u(:,1);
                    lf_single  = lf_single * eta;

                    if do_lat,
                        [u, ~] = svd(real(pinv_plus(tmp(nn+1:2*nn,nn+1:2*nn),reduce_rank,0)),'econ'); % this is faster,  - MWW                                                
                        eta = u(:,1);
                        lf_lat  = lf_lat * eta;

                        lf=[lf_single, lf_lat];
                        tmp2 = pinv_plus(lf' * invCy * lf) * lf' * invCy;
                        W{i}=tmp2(1,:);

                    else
                        lf=lf_single;

                        % construct the spatial filter
                        %W{i} = pinv(lf' * invCy * lf) * lf' * invCy;
                        W{i} = lf'*invCy/(lf' * invCy * lf); % this is faster - MWW
                    end;                       


                  case 'vector'
                    W{i} = pinv_plus(lf' * invCy * lf,reduce_rank) * lf' * invCy; % - AB

                end

            else
                W{i} = NaN;
            end

             if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
        end

        inverse.class{kk}.W = W;
        inverse.class{kk}.L = L;
    end;

end

res = inverse;