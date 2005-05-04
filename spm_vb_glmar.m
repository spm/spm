function [slice] = spm_vb_glmar (Y,slice)
% Variational Bayes for GLM-AR modelling in a slice of fMRI
% FORMAT [slice] = spm_vb_glmar (Y,slice)
%
% Y             [T x N] time series with T time points, N voxels
%
% slice         data structure containing the following fields:
%
%               .X              [T x k] the design matrix       
%               .p              order of AR model               
%               .D              [N x N] spatial precision matrix
%                               (see spm_vb_set_priors.m)
%
%               The above fields are mandatory. The fields below are
%               optional or are filled in by this function.
%
%               OPTIMISIATION PARAMETERS:
%
%               .tol            termination tolerance (default = 0.01% increase in F)
%               .maxits         maximum number of iterations (default=4)
%               .verbose        '1' for description of actions (default=1)
%               .update_???     set to 1 to update parameter ??? (set to 0 to fix)
%                               eg. update_alpha=1; % update prior precision on W
%
%               ESTIMATED REGRESSION COEFFICIENTS:
%
%               .wk_mean        [k x N] VB regression coefficients
%               .wk_ols         [k x N] OLS "  "
%               .w_cov          N-element cell array with entries [k x k]
%               .w_dev          [k x N] standard deviation of regression coeffs
%
%               ESTIMATED AR COEFFICIENTS:
%
%               .ap_mean        [p x N] VB AR coefficients
%               .ap_ols         [p x N] OLS AR coefficients
%               .a_cov          N-element cell array with entries [p x p]
%
%               ESTIMATED NOISE PRECISION:
%
%               .b_lambda       [N x 1] temporal noise precisions
%               .c_lambda
%               .mean_lambda  
%
%               MODEL COMPARISON AND COEFFICIENT RESELS:
%
%               .gamma_tot      [k x 1] Coefficient RESELS 
%               .F              Negative free energy (used for model selection)
%               .F_record       [its x 1] record of F at each iteration                  
%               .elapsed_seconds  estimation time 
%               PRIORS:
%
%               .b_alpha        [k x 1] spatial prior precisions for W
%               .c_alpha    
%               .mean_alpha 
%
%               .b_beta         [p x 1] spatial prior precisions for AR
%               .c_beta    
%               .mean_beta 
%
%               HYPERPRIORS:
%
%               .b_alpha_prior   priors on alpha
%               .c_alpha_prior
%
%               .b_beta_prior    priors on beta
%               .c_beta_prior
%
%               .b_lambda_prior  priors on temporal noise precisions
%               .c_lambda_prior
%
%               There are other fields that are used internally
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_glmar.m 112 2005-05-04 18:20:52Z john $


tic;

% Check arguments and set defaults
if nargin < 2, 
   disp('Error in spm_vb_glmar: function needs two arguments'); 
   return
end
[T,N]=size(Y);
slice.T=T;
slice.N=N;
if ~isfield(slice,'X')
    disp('Error in spm_vb_glmar: mandatory field X missing');
    return
else
    X=slice.X;
end
[tmp,k]=size(X);
if ~(tmp==T)
   disp('Error in spm_vb_glmar: X is not of compatible size to Y');
   return
end
slice.k=k;
p=slice.p;

if ~isfield(slice,'Dw')
    disp('Error in spm_vb_glmar: mandatory field Dw missing');
    return
end

F = 0;
last_F=0;
slice=spm_vb_init_slice(Y,slice);

if slice.verbose
    disp(' ');
    disp('Starting VB-GLM-AR-SLICE');
end

for it = 1:slice.maxits, % Loop over iterations
    
    if slice.update_w
        slice=spm_vb_w (Y,slice);
    end
    if (slice.p>0) & (slice.update_a)
        slice=spm_vb_a (Y,slice);
    end
    if slice.update_lambda
        slice=spm_vb_lambda (Y,slice);
    end
    if slice.update_alpha
        slice=spm_vb_alpha (Y,slice);
    end
    if (slice.p>0) & (slice.update_beta)
        slice=spm_vb_beta (Y,slice);
    end
    if slice.update_F
        [F,Lav,KL]=spm_vb_F (Y,slice);
    end
    if slice.verbose
        disp(sprintf('Iteration %d, F=%1.2f',it,F));
    end
    
    if slice.update_F
        slice.F_record(it)=F;
        delta_F=F-last_F;
        if it > 2
            if delta_F < 0
                disp(sprintf('********** Warning: decrease in F of %1.4f per cent *************',100*(delta_F/F)));
                keyboard
                break;
                
                
            elseif abs(delta_F/F) < slice.tol,
                break;
            end;     
        end
        last_F=F;
    end
end;

if slice.update_F
    slice.F=F;
    slice.Lav=Lav;
    slice.KL=KL;
end
    
slice=spm_vb_gamma(Y,slice);

slice.elapsed_seconds=toc;

