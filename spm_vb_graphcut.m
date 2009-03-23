function labels = spm_vb_graphcut(labels,index,I,W,depth,grnd_type,CUTOFF,DIM)
% Recursive bi-partition of a graph using the isoperimetric algorithm
% 
% FORMAT labels = spm_vb_graphcut(labels,index,I,W,depth,grnd_type,CUTOFF,DIM)
%
% labels     each voxel is lableled depending on whihc segment is belongs
% index      index of current node set in labels vector
% I          InMask XYZ voxel (node set) indices
% W          weight matrix i.e. adjacency matrix containing edge weights 
%            of graph
% depth      depth of recursion
% grnd_type  'random' or 'max' - ground node selected at random or the 
%            node with maximal degree
% CUTOFF     minimal number of voxels in a segment of the partition
% DIM        dimensions of data
%__________________________________________________________________________
% 
% Recursive bi-partition of a graph using the isoperimetric algorithm by 
% Grady et al. This routine is adapted from "The Graph Analysis Toolbox: 
% Image Processing on Arbitrary Graphs", available through Matlab Central
% File Exchange. See also Grady, L. Schwartz, E. L. (2006) "Isoperimetric 
% graph partitioning for image segmentation",
% IEEE Trans Pattern Anal Mach Intell, 28(3),pp469-75
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Lee Harrison
% $Id: spm_vb_graphcut.m 2923 2009-03-23 18:34:51Z guillaume $

try, grnd_type; catch, grnd_type = 'random'; end

%-Initialize
s = rand('twister'); rand('seed',0);
N           = length(index);
terminate   = 0;

%-Partition current graph
if N > CUTOFF
    reject  = 1;
    d       = sum(W,2); % "degree" of each node 
    switch grnd_type
        case 'max'
            id      = find(d==max(d));
            ground  = id(ceil(rand(1)*length(id)));
        case 'random'
            ground  = ceil(rand(1)*N);
    end
    while reject == 1,
        if N < 1e5, method = 'direct'; else, method = 'cg'; end
        parts   = bipartition(W,CUTOFF,ground,method);        
        nparts  = [length(parts{1}),length(parts{2})];
        if min(nparts) < CUTOFF, terminate = 1; break, end
        for k = 1:2, % check if partitions are contiguous
            bw                      = zeros(DIM(1),DIM(2),DIM(3));
            bw(I(parts{k}))    = 1;
            [tmp,NUM]               = spm_bwlabel(bw,6);
            if NUM > 1
                reject  = 1;
                ground  = ceil(rand(1)*N); % re-select ground node
                fprintf('depth %1.0f, partition %1.0f of 2, reject ',depth,k); fprintf('\n')
                break
            else
                reject  = 0;
                fprintf('depth %1.0f, partition %1.0f of 2, accept ',depth,k);
                fprintf('\n')
            end               
        end
    end   
else
    terminate = 1;
end

if terminate
    labels  =   labels;
    fprintf('depth %1.0f, end of branch ',depth);
    fprintf('\n')
    rand('twister',s);
else
    %Accept partition and update labels
    tmpInd                  =   find(labels>labels(index(1)));
    labels(tmpInd)          =   labels(tmpInd) + 1; %Make room for new class
    labels(index(parts{2})) =   labels(index(parts{2})) + 1; %Mark new class

    %Continue recursion on each partition
    if nparts(1) > CUTOFF
        labels = spm_vb_graphcut(labels,index(parts{1}),I(parts{1}),...
            W(parts{1},parts{1}),depth + 1,grnd_type,CUTOFF,DIM);
    end

    if nparts(2) > CUTOFF
        labels = spm_vb_graphcut(labels,index(parts{2}),I(parts{2}),...
            W(parts{2},parts{2}),depth + 1,grnd_type,CUTOFF,DIM);
    end
end

%==========================================================================
% function parts = bipartition(W,CUTOFF,ground,method) 
%==========================================================================
function parts = bipartition(W,CUTOFF,ground,method) 
% Computes bi-partition of a graph using isoperimetric algorithm.

% FORMAT parts = bipartition(W,CUTOFF,ground,method) 

% parts     1x2 cell containing indices of each partition 
% W         weight matrix 
% CUTOFF    minimal number of voxels in a segment of the partition
% ground    ground node index
% method    used to solve L0*x0=d0. Options are 'direct' or 'cg', which 
%           x0 = L0\d0 or preconditioned conjugate gradients (see Matlab 
%           rountine pcg.m)   

try, method; catch, method = 'direct'; end

%-Laplacian matrix
d   =   sum(W,2);
L   =   diag(d) - W;
N   =   length(d);

%-Compute reduced Laplacian matrix, i.e. remove ground node
index   =   [1:(ground-1),(ground+1):N];
d0      =   d(index);
L0      =   L(index,index);

%-Solve system of equations L0*x0=d0
switch method
    case 'direct'
        x0   =   L0\d0;
    case 'cg'
        x0   =   pcg(L0,d0,[],N,diag(d0));
end

%-Error catch if numerical instability occurs (due to ill-conditioned or 
%singular matrix)
minVal  =   min(min(x0)); 
if minVal < 0
    x0(find(x0 < 0))  =   max(max(x0)) + 1;
end

%-Re-insert ground point
x0   =   [x0(1:(ground-1));0;x0((ground):(N-1))];

%-Remove sparseness of output
x0   =   full(x0);

%-Determine cut point (ratio cut method)    
indicator   =   sparse(N,1);

%-Sort values
sortX       =   sortrows([x0,[1:N]'],1)';

%-Find total volume
totalVolume     =   sum(d);
halfTotalVolume =   totalVolume/2;

%-Calculate denominators
sortedDegree    =   d(sortX(2,:))';
denominators    =   cumsum(sortedDegree);
tmpIndex                =   find(denominators > halfTotalVolume);
denominators(tmpIndex)  =   totalVolume - denominators(tmpIndex);

%-Calculate numerators
L           =   L(sortX(2,:),sortX(2,:)) - diag(sortedDegree);
numerators  =   cumsum(sum((L - 2*triu(L)),2))';
if min(numerators) < 0
    %Line used to avoid negative values due to precision issues
    numerators  =   numerators - min(numerators) + eps;
end

%-Calculate ratios for Isoperimetric criteria
sw = warning('off','MATLAB:divideByZero');
[constant,minCut]   =   min(numerators(CUTOFF:(N-CUTOFF))./ ...
        denominators(CUTOFF:(N-CUTOFF)));
minCut  =   minCut + CUTOFF - 1;
warning(sw);

%-Output partitions
parts{1}   =   sortX(2,1:(minCut))';
parts{2}   =   sortX(2,(minCut+1):N);
