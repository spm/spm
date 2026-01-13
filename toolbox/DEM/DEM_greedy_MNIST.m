function [C,F,G,E,dom] = DEM_greedy_MNIST(OPTIONS)
% Structure learning from pixels
%__________________________________________________________________________
%
% This routine demonstrates the mapping from pixels to (discrete)
% representations of episodes. It illustrates the use of a certain kind of
% deep generative model based upon reduction and grouping operators (RG)
% found in the renormalisation group. In brief, this allows an efficient
% encoding (and generation) of high dimensional content by maintaining
% conditionally independent partitions of hidden states in a recursive
% fashion.
%
% This example uses image classification to showcase the ability of
% recursive generative models to compress images in a way that leverages
% the composition of small image components into larger arrangements. Uses
% the MNIST digit classification problem. In brief, images are discretised
% using singular value decomposition of local image patches. These are then
% recursively assembled using a fast structure learning routine. Subsequent
% training images and then assimilated using parametric learning. The
% ensuing generative model is then illustrated in terms of digit
% classification and the underlying architecture.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% This subroutine expects to find MNIST.mat in its working directory. The
% requisite MATLAB file can be downloaded from the following website
%--------------------------------------------------------------------------
% https://lucidar.me/en/matlab/load-mnist-database-of-handwritten-digits-in-matlab/



%% Greedy search
%==========================================================================
OPTIONS.N  = 10000;  % Number of training data
OPTIONS.T  = 1000;   % Number of test data
OPTIONS.q  = 16;     % concentration parameter  
OPTIONS.s  = 2;      % smoothing in pixels  
OPTIONS.S  = 24;     % histogram width

OPTIONS.Ne = 8;     % exemplars per class

OPTIONS.nd = 4;      % Diameter of tiles in pixels
OPTIONS.nb = 4;      % Number of discrete singular variates 
OPTIONS.mm = 8;      % Maximum number of singular modes

params(1).name   = 'q';
params(1).domain = fix(2.^(4:2:12));
params(2).name   = 'Ne';
params(2).domain = fix(12:32);
params(3).name   = 'nb';
params(3).domain = fix(4:2:16);
params(4).name   = 'beta';
params(4).domain = fix(2.^(0:8));
params(5).name   = 'eta';
params(5).domain = fix(2.^(6:16));
params(6).name   = 'mm';
params(6).domain = fix(4:32);

for p = 1:numel(params)
    params(p).i = 1;
    OPTIONS.(params(p).name) = params(p).domain(params(p).i);
end

C     = DEM_image_compression(OPTIONS);
clc; disp(OPTIONS); disp(C)


for t = 1:512
    for p = 1:numel(params)
        try
            OPTIONS.(params(p).name) = params(p).domain(params(p).i + 1);
            [c,f,g,e] = DEM_image_compression(OPTIONS);
            if c >= C
                params(p).i = params(p).i + 1;
                C = c;
            end
        end

        % report
        %------------------------------------------------------------------
        disp(OPTIONS); disp(C); disp(t)
    end
end

return



params = params(1:2);
for n = 1:numel(params)
    Np(n) = numel(params(n).domain);
end

dom = spm_combinations(Np);
C   = NaN(Np);
F   = NaN(Np);
G   = NaN(Np);
E   = NaN(Np);
for i = 1:size(dom,1)

    OPT   = OPTIONS;
    for j = 1:size(dom,2)
        OPT.(params(j).name) = params(j).domain(dom(i,j));
    end

    [c,f,g,e] = DEM_image_compression(OPT);
    ind       = num2cell(dom(i,:));
    C(ind{:}) = c
    E(ind{:}) = e(end)
    F(ind{:}) = f
    G(ind{:}) = mean(g)

end

spm_figure('GetWin','Model space'),clf
subplot(2,2,1)
bar3(C)
xlabel(params(2).name)
ylabel(params(1).name)
title('Classifcation accuracy')

subplot(2,2,2)
bar3(F)
xlabel(params(2).name)
ylabel(params(1).name)
title('ELBO')

subplot(2,2,3)
plot(F,C)
xlabel('ELBO')
ylabel('Classifcation accuracy'), axis square

subplot(2,2,4)
plot(F,G)
xlabel('ELBO (test)')
ylabel('ELBO (train)'), axis square


return


C = [90.3000   89.8000   88.2000   90.0000
   90.8000   90.9000   83.8000   89.8000
   89.5000   90.9000   89.7000   87.6000
   90.5000   91.5000   88.0000   91.9000
   89.6000   89.7000   89.5000   88.2000
   70.7000   90.3000   89.8000   88.6000
   91.4000   89.8000   90.7000   89.9000
   83.7000   87.3000   84.1000   89.5000
   90.7000   90.9000   90.7000   90.3000]


E = [229.0411  315.2647  371.7620  414.3841
  232.2706  316.1503  377.0265  419.4966
  235.0781  317.9290  377.9982  419.9177
  236.2536  320.5570  379.2613  424.0176
  235.7651  320.9751  381.9095  426.0738
  234.6384  321.4667  382.6702  424.8058
  237.7358  323.7133  383.5227  429.0118
  237.8516  325.0924  384.4635  428.9537
  239.3182  326.2257  385.6117  427.8780]


F =  [-13.6231  -14.4211  -14.4477  -15.0572
  -13.7034  -14.3174  -14.5218  -15.0839
  -13.4861  -14.6218  -15.1126  -14.7807
  -14.0254  -14.8677  -15.0446  -15.2778
  -13.5641  -14.8648  -15.1260  -15.3810
  -13.6049  -14.6471  -14.7732  -14.7916
  -13.7172  -14.4139  -14.6270  -15.0213
  -13.0335  -14.4133  -14.4790  -14.9801
  -13.1978  -14.5832  -14.6595  -14.7468]


G =  [-6.2552   -7.5738   -7.7339   -8.2819
   -6.6176   -7.2266   -8.0620   -8.5000
   -6.2352   -7.8203   -8.2255   -8.5121
   -6.7583   -7.7030   -8.3546   -8.6698
   -6.2430   -7.7500   -8.4243   -8.9356
   -7.2858   -7.7904   -8.0552   -8.4625
   -6.5117   -7.7055   -7.8296   -8.5658
   -6.5606   -7.9137   -8.0539   -8.3469
   -6.1798   -7.6189   -7.8633   -8.1273]






