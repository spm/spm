function [Con] = spm_make_contrasts (k1,k2)
% Make contrasts for one/two-way ANOVAs
% FORMAT [Con] = spm_make_contrasts (k1,k2)
%
% k1    Number of levels of factor 1
% k2    Number of levels of factor 2
%
% This function computes contrasts for a generic k1-by-k2 design. 
% It is assumed that the levels of the first factor 
% change slowest. For one-way ANOVAs set k2=1. 
%
% This function generates (transposed) contrast matrices to test:
% (1) average effect of condition, (2) main effect of factor 1
% (3) main effect of factor 2, (4) interaction
%
% Con(c).c      Contrast matrix
%       .name   Name
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id$
    
% See section 5 of 
% http://www.fil.ion.ucl.ac.uk/~wpenny/publications/rik_anova.pdf
% for details of kronecker operations

C1=ones(k1,1);
C2=ones(k2,1);

D1=-1*diff(eye(k1))';
D2=-1*diff(eye(k2))';

Con(1).c=kron(C1,C2)';
Con(2).c=kron(D1,C2)';
Con(3).c=kron(C1,D2)';
Con(4).c=kron(D1,D2)';

Con(1).name='Average effect of condition';
Con(2).name='Main effect of factor 1';
Con(3).name='Main effect of factor 2';
Con(4).name='Interaction';
