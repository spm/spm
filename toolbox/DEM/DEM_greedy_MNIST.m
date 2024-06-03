function [C,F,G] = DEM_greedy_MNIST(OPTIONS)
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
OPTIONS.N  = 1000;   % Number of training data
OPTIONS.T  = 500;    % Number of test data
OPTIONS.q  = 16;     % concentration parameter  
OPTIONS.s  = 2;      % smoothing in pixels  
OPTIONS.S  = 24;     % histogram width

OPTIONS.Ne = 9;      % exemplars per class

OPTIONS.nd = 5;      % Diameter of tiles in pixels
OPTIONS.nb = 16;     % Number of discrete singular variates 
OPTIONS.mm = 8;      % Maximum number of singular modes


params = {'nb','q','Ne'};
d      = [2,    4,  1];
b      = [2,    4,  1];
C      = [];
F      = [];
G      = [];
for i  = 1:128
    for j = 1:numel(params)
        OPT             = OPTIONS;
        OPT.(params{j}) = OPTIONS.(params{j}) + 0;
        [c,f,g]         = DEM_image_compression(OPT);
        C(1,end + 1)    = c;
        F(1,end + 1)    = f;
        G(1,end + 1)    = mean(g);
        OPT             = OPTIONS;
        OPT.(params{j}) = OPTIONS.(params{j}) - d(j);
        [c,f,g]         = DEM_image_compression(OPT);
        C(2,end)        = c;
        F(2,end)        = f;
        G(2,end)        = mean(g);
        OPT             = OPTIONS;
        OPT.(params{j}) = OPTIONS.(params{j}) + d(j);
        [c,f,g]         = DEM_image_compression(OPT);
        C(3,end)        = c;
        F(3,end)        = f;
        G(3,end)        = mean(g);

        c = max(C(:,end)); disp(C)
        if C(2,end) == c
            OPTIONS.(params{j}) = OPTIONS.(params{j}) - d(j);
            d(j) = max(b(j),d(j) - 1);
        elseif C(3,end) == c
            OPTIONS.(params{j}) = OPTIONS.(params{j}) + d(j);
            d(j) = max(b(j),d(j) - 1);
        end

        disp(C)
        disp(OPTIONS)
        disp(d)
    end
end










