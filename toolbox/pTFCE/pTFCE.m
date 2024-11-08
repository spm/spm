function [pTFCE_Z, pTFCE_p] = pTFCE(imgZ,mask, Rd, V, Nh, Zest, C, verbose)
% +-----------------------------------------------------------------------+
% |PTFCE core function for pTFCE enhancement                              | 
% +-----------------------------------------------------------------------+
% | Version 0.2.1                                                        |
% +-----------------------------------------------------------------------+
% | DESCRIPTION                                                           |
% +-----------------------------------------------------------------------+
% | The threshold-free cluster enhancement (TFCE) approach integrates cluster
% | information into voxel-wise statistical inference to enhance detectability
% | of neuroimaging signal. Despite the significantly increased sensitivity, 
% | the application of TFCE is limited by several factors: (i) generalization
% | to data structures, like brain network connectivity data is not trivial, 
% | (ii) TFCE values are in an arbitrary unit, therefore, P-values can only
% | be obtained by a computationally demanding permutation-test.
% | Here, we introduce a probabilistic approach for TFCE (pTFCE), that gives
% | a simple general framework for topology-based belief boosting.
% | The core of pTFCE is a conditional probability, calculated based on Bayes'
% | rule, from the probability of voxel intensity and the threshold-wise
% | likelihood function of the measured cluster size. We provide an estimation
% | of these distributions based on Gaussian Random Field (GRF) theory.
% | The conditional probabilities are then aggregated across cluster-forming
% | thresholds by a novel incremental aggregation method. Our approach is
% | validated on simulated and real fMRI data.
% | pTFCE is shown to be more robust to various ground truth shapes and
% | provides a stricter control over cluster "leaking" than TFCE and, in the
% | most realistic cases, further improves its sensitivity. Correction for
% | multiple comparison can be trivially performed on the enhanced P-values,
% | without the need for permutation testing, thus pTFCE is well-suitable for
% | the improvement of statistical inference in any neuroimaging workflow.
% | 
% | This matlab implementation is a port of the pTFCE R-package (v0.0.4) and
% | validated to that.
% | 
% +-----------------------------------------------------------------------+
% | PLEASE CITE                                                           |
% +-----------------------------------------------------------------------+
% | T. Spisák, Z. Spisák, M. Zunhammer, U. Bingel, S. Smith, T. Nichols, T.
% | Kincses, Probabilistic TFCE: a generalized combination of cluster size
% | and voxel intensity to increase statistical power, Neuroimage 185:12-26,
% | 2019.
% +-----------------------------------------------------------------------+
% | WEBSITE                                                               |
% | https://github.com/spisakt/pTFCE                                      |
% +-----------------------------------------------------------------------+
% | 
% +-----------------------------------------------------------------------+
% | INPUTS                                                                |
% +-----------------------------------------------------------------------+
% | - imgZ: Z-score image to enhance
% | - mask: Mask
% | - Nh: Number of thresholds
% | - Rd: Resel count (E.g. SPM.xVol.R(4). Note that the R-package uses the
% |   FSL-like Resel count, whereas the SPM GUI Toolbox uses the SPM-like
% |   (corrected) vresel count)
% | - V: Number of voxels in mask (e.g. SPM.xVol.S)
% | - ZestThr: Cluster-forming Z threshold below which P(h|c) is estimated
% | as P(h), due to limitation of GRF theory. (default: 1.3)
% | - C: voxel connectivity: 4 | 8 | 6 | 18 | 26 | 3-by-3-by- ... -by-3
% | matrix of 0s and 1s (default: 6)
% | - verbose: print out progress
% | 
% +-----------------------------------------------------------------------+
% |  DETAILS                                                              |
% +-----------------------------------------------------------------------+
% | The function takes a Z-score image and a mask image as obligatory
% | inputs.
% | Mask can be either binary or continous, in the latter case it will be
% | thresholded at 0.5.
% | Smoothness information Rd (voxels per RESEL) and number of voxels V are
% | optional and are to be interpreted as in FSL "smoothest".
% | These values can be found in the SPM.mat or estimated from the data, internally via
% | the smoothest() function of the R-package, which is a direct port of the
% | corresponding FSL function.
% | If Rd and/or V is not specified, and residual is specified, image smoothness
% | will be determined basedd on teh 4D residual data (more accurate, see ?smoothest()).
% | The default value of the parameter Nh should work with images having usual
% | ranges of Z-scores. At low values, although the processing becomes faster,
% | the estimation of the enhanced values might become inaccurate and the enhanced
% | p-value distribution becomes increasingly non-uniform. It is not recommended
% | to set it lower than 30.
% | The parameters logpmin and logpmax define the range of values the incremental
% | thresholding procedure covers. By default, these are based on the input data.
% | 
% +-----------------------------------------------------------------------+
% |  RETURNS                                                              |
% +-----------------------------------------------------------------------+
% | - pTFCE_Z: Z-score map of the pTFCE-enhanced SPM
% | - pTFCE_p: (uncorrected) log P-value map of the pTFCE enhnaced SPM
% | 
% +-----------------------------------------------------------------------+
% |  CONTATCT:                                                            |
% |  Tamas Spisak                                                         |
% |  tamas.spisak@uk-essen.de                                             |
% +-----------------------------------------------------------------------+


if nargin > 7
    error('pTFCE requires at most 3 optional inputs');
end

if nargin < 4
    error('pTFCE requires at least 4 inputs');
end

% Fill in unset optional values.
switch nargin
    case 4
        Nh = 100;
        Zest = 1.3;
        C = 6;
        verbose = 0;
    case 3
        Zest = 1.3;
        C = 6;
        verbose = 0;
    case 2
        C = 6;
        verbose = 0;
    case 1
        verbose = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START

if sum(isnan(imgZ))>0
    fprintf('%s\n','NAs detected and replaced with zero!');
    imgZ(isnan(imgZ))=0;
end

if sum(isinf(imgZ))>0
    fprintf('%s\n','Infinite values detected and replaced with zero!');
    imgZ(isinf(imgZ))=0;
end

logpmin=0;
logpmax=-log(normcdf(max(imgZ(:)), 'upper'));
logp_thres=linspace(logpmin, logpmax, Nh);
dh=logp_thres(2)-logp_thres(1);
p_thres=exp(-logp_thres);
threshs=-norminv(p_thres);
threshs(1)=-9999999999;
ndh = length(threshs);
PVC=ones([size(imgZ), ndh]);

if (verbose)
    fprintf('%s\n','* Performing pTFCE...\n');
end

spm_progress_bar('Init',ndh,'','Calculating pTFCE...','');

% do connected component analysis for the whole set of threesholds
cc = arrayfun(@(x) bwconncomp(bsxfun(@ge,imgZ,x),C), threshs);
nvox=length(imgZ);

for hi=1:ndh 
    pvc=ones(size(imgZ));
    CLUST=zeros(size(imgZ));
    spm_progress_bar('Set', hi);
    if (verbose)
        fprintf('  - P threshold: %s\n', num2str(p_thres(hi)) );
    end
    
    % assess connected components (clusters)
    ccc = cc(hi);
    voxpercc = cellfun(@numel,ccc.PixelIdxList);
    for c = 1:ccc.NumObjects
        CLUST(ind2sub(size(imgZ), ccc.PixelIdxList{c})) = voxpercc(c);
    end
    
    sizes = unique(CLUST(:));
    for i = 1:length(sizes)
        siz = sizes(i);
        if siz==0
            continue;
        end
        pvc(CLUST==siz) = pvox_clust(V, Rd, siz, threshs(hi), Zest);
        PVC(:,:,:,hi)=pvc;
        
    end 
end
% calculate pTFCE
% R-code: pTFCE=array(apply(PVC, c(1,2,3), function(x){exp( -aggregate.logpvals(-log(x), dh) )}), dim=dim(img))
% make a loop instead:
pTFCE=zeros(size(imgZ));
for i=1:size(imgZ(:))
    [x,y,z]=ind2sub(size(imgZ), i);
    pTFCE(x,y,z) = exp( -aggregate_logpvals(-log(PVC(x, y, z,:)), dh) );
end

% CONVERT BACK TO Z-SCORE
pTFCE(pTFCE==0) = realmin;
pTFCE(pTFCE==1) = 1-eps(1);
pTFCE_Z = -norminv(pTFCE); % ToDo check for overflow!

%figure;
%image(logpTFCE(:,:,30)); % display the axial slice

spm_progress_bar('Clear');
pTFCE_Z=pTFCE_Z;
pTFCE_p=pTFCE;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aggregate logpvals by the "equidistant incremental logarithmic probability aggregation".
%
% @param logpvals vector of enhnaced -logP values correspoinding tho various thresholds
% @param d delta -logP (giving the equidistant distribution in the log-space)
%
% @return aggregated -logP probability
% @export
%
% @examples
function logp = aggregate_logpvals(logpvals, d)
  % ToDo: handle overflow here if needed
  s = sum(logpvals);
  logp=0.5*(sqrt(d*(8*s+d))-d);
end

% expected value of cluster size at threshold h in imgZ of size V with RESEL count Rd
% mainly ported from the source code of FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
%
% @param h image height, that is, Z score threshold
% @param V Number of voxels
% @param Rd Rd (dLh*V)
%
% @return Expected value of cluster size
%
% @examples
function Es=Es(h, V, Rd)
  h2=h.^2;
  ret = log(V) + log(normcdf(h, 'upper'));
  
  h2=h2(h>=1.1); % GRF is inaccurate, use the estimation of FSL
  ret(h>=1.1) = ret(h>=1.1) - ( log(Rd)+log(h2-1)-h2/2+-2*log(2*pi));
  Es = exp(ret);
end


function dvox=dvox(h) % PDF of Z threshold/voxel value
  dvox = normpdf(h);
end

function pvox=pvox(h) % p-value for voxel value
  pvox = normcdf(h, 'upper');
end


function dcl=dcl(h, V, Rd, c, ZestThr) % function to integrate by dclust()
    lambda=(Es(h, V, Rd)/gamma(2.5) ).^(-2/3);
    dcl=lambda.*exp(-lambda.*c.^(2/3));
    dcl(isnan(dcl))=0; % underflow happened
    dcl(h<ZestThr)=0; % truncate distribution
    dcl;
end

function dclust=dclust(h, V, Rd, c, ZestThr) % PDF of cluster extent, given h thershold
    if nargin < 5
        ZestThr=1.3; % default value for Z threshold for GRF-based estimation
    end

    if (isnan(c))     % NaN denotes subthreshold voxel
      dclust=dvox(h);
      return
    end
    fun = @(x) dcl(x,V,Rd,c,ZestThr);
    %  patch: the normalising constant can be ignored here as the expression will be normalised again later
    % So we spare a costly numarical integration and get much faster!
    % Tests still pass for the R-package, with a reasonable tolerance
    dclust = dcl(h, V, Rd, c, ZestThr); %/ integral(fun, -Inf, Inf);
end


function dvc=dvox_clust(h, V, Rd, c, ZestThr) % PDF of Z threshold value given cluster size
    if nargin < 5
        ZestThr=1.3; % default value for Z threshold for GRF-based estimation
    end
    
    if isnan(c)
        dvc = dvox(h);
        return % return PDF of Z-threshold value in the case of subthreshold voxel
    end
        
    % ToDo: underflow handling
    %if (isnan(dclust(h, V, Rd, c)[1])) % underflow
    %    return(dvox(h)*dclust(h, V, Rd, c))
    fun = @(x) dvox(x).*dclust(x, V, Rd, c, ZestThr);
    dvc=dvox(h).*dclust(h, V, Rd, c, ZestThr) / integral(fun, -Inf, Inf);
end
 
function pvc=pvox_clust(V, Rd, c, actH, ZestThr) % p-value for Z threshold value given cluster size
    if nargin < 5
        ZestThr=1.3; % default value for Z threshold for GRF-based estimation
    end
    
    if actH<=ZestThr  % ZestThr: GRF theory might not apply at low thresholds
        pvc=pvox(actH);
        return
    end

  % ToDo: underflow handling
  %if (is.nan(dvox.clust(actH, V, Rd, c, ZestThr=ZestThr)[1])) % underflow
  %  return(exp(-745)) #TODO: make machnie independent
    fun = @(x) dvox_clust(x, V, Rd, c, ZestThr);
    pvc=integral(fun, actH, Inf);
end


