function [block] = spm_vb_glmar (Y,block)
% Variational Bayes for GLM-AR modelling in a block of fMRI
% FORMAT [block] = spm_vb_glmar (Y,block)
%
% Y     -  [T x N] time series with T time points, N voxels
%
% block -  data structure containing the following fields:
%
%          .X              [T x k] the design matrix       
%          .p              order of AR model               
%          .D              [N x N] spatial precision matrix
%                          (see spm_vb_set_priors.m)
%
%          The above fields are mandatory. The fields below are
%          optional or are filled in by this function.
%
%          OPTIMISIATION PARAMETERS:
%
%          .tol            termination tolerance (default = 0.01% increase in F)
%          .maxits         maximum number of iterations (default=4)
%          .verbose        '1' for description of actions (default=1)
%          .update_???     set to 1 to update parameter ??? (set to 0 to fix)
%                          eg. update_alpha=1; % update prior precision on W
%
%          ESTIMATED REGRESSION COEFFICIENTS:
%
%          .wk_mean        [k x N] VB regression coefficients
%          .wk_ols         [k x N] OLS "  "
%          .w_cov          N-element cell array with entries [k x k]
%          .w_dev          [k x N] standard deviation of regression coeffs
%
%          ESTIMATED AR COEFFICIENTS:
%
%          .ap_mean        [p x N] VB AR coefficients
%          .ap_ols         [p x N] OLS AR coefficients
%          .a_cov          N-element cell array with entries [p x p]
%
%          ESTIMATED NOISE PRECISION:
%
%          .b_lambda       [N x 1] temporal noise precisions
%          .c_lambda
%          .mean_lambda  
%
%          MODEL COMPARISON AND COEFFICIENT RESELS:
%
%          .gamma_tot      [k x 1] Coefficient RESELS 
%          .F              Negative free energy (used for model selection)
%          .F_record       [its x 1] record of F at each iteration                  
%          .elapsed_seconds  estimation time 
%          PRIORS:
%
%          .b_alpha        [k x 1] spatial prior precisions for W
%          .c_alpha    
%          .mean_alpha 
%
%          .b_beta         [p x 1] spatial prior precisions for AR
%          .c_beta    
%          .mean_beta 
%
%          .b              [k x N] prior precision matrix
%
%          HYPERPRIORS:
%
%          .b_alpha_prior   priors on alpha
%          .c_alpha_prior
%
%          .b_beta_prior    priors on beta
%          .c_beta_prior
%
%          .b_lambda_prior  priors on temporal noise precisions
%          .c_lambda_prior
%
%          There are other fields that are used internally
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_glmar.m 4310 2011-04-18 16:07:35Z guillaume $


t0 = clock;

%-Check arguments
%-----------------------------------------------------------------------
error(nargchk(2,2,nargin));

%-Set defaults
%-----------------------------------------------------------------------
[T N]   = size(Y);
block.T = T;
block.N = N;

try
    X = block.X;
catch
    error('Mandatory field X missing.');
end

[tmp k] = size(X);
if tmp ~= T
   error('X is not of compatible size to Y.');
end

block.k = k;
p = block.p;

if ~isfield(block,'Dw')
    disp('Error in spm_vb_glmar: mandatory field Dw missing');
    return
end

%-Initialise block
%-----------------------------------------------------------------------
F      = 0;
last_F = 0;
block = spm_vb_init_block(Y,block);

%-Run Variational Bayes
%-----------------------------------------------------------------------
if block.verbose
    disp(' ');
    disp('Starting VB-GLM-AR-BLOCK');
end

for it = 1:block.maxits, % Loop over iterations
    
    if block.update_w
        block = spm_vb_w (Y,block);
    end
    if (block.p>0) && (block.update_a)
        block = spm_vb_a (Y,block);
    end
    if block.update_lambda
        block = spm_vb_lambda (Y,block);
    end
    if block.update_alpha
        block = spm_vb_alpha (Y,block);
    end
    if (block.p>0) && (block.update_beta)
        block = spm_vb_beta (Y,block);
    end
    if block.update_F
        [F, Lav, KL] = spm_vb_F (Y,block);
    end
    if block.verbose
        disp(sprintf('Iteration %d, F=%1.2f',it,F));
    end
    
    if block.update_F
        block.F_record(it)=F;
        delta_F=F-last_F;
        if it > 2
            if delta_F < 0
                disp(sprintf('********** Warning: decrease in F of %1.4f per cent *************',100*(delta_F/F)));
                break;
            elseif abs(delta_F/F) < block.tol,
                break;
            end;     
        end
        last_F = F;
    end
end

if block.update_F
    block.F   = F;
    block.Lav = Lav;
    block.KL  = KL;
end
    
block = spm_vb_gamma(Y,block);

block.elapsed_seconds = etime(clock,t0);

