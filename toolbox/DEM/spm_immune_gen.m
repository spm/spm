function [Y,X] = spm_immune_gen(P,M,U)
% Generative model of an immune response
% FORMAT [Y,X] = spm_immune_gen(P,M,U)
% Y   - timeseries data
% X   - latent states
% P   - Priors
% M   - Model
% U   - inputs (timing of measurements)
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging
 
% Thomas Parr
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% setup and defaults
%--------------------------------------------------------------------------
try, M.T; catch, M.T = 1920;       end     % 40 days
try, U;   catch, U   = [1:M.T]/24; end     % Days on which measurements are taken


% exponentiate parameters
%--------------------------------------------------------------------------
Q    = spm_vecfun(P,@exp);
n    = Q.n;                 % Pre-existing cross-reactive memory B-cells
m    = Q.m;                 % Amount of (intracellular) virus at start

% initial marginals
%--------------------------------------------------------------------------

p{1} = [1 0 0]';                           % IgM antibodies (absent, present, neutralising)
p{2} = [1-n*exp(-4) 0 0 n*exp(-4) 0]';     % B-cells (naive, IgM plasma, inactive memory, active memory, IgG plasma)
p{3} = [1 0 0]';                           % T-cells (naive, CD4+, CD8+)
p{4} = [1-m*exp(-8) m*exp(-8) 0]';         % Virus (absent, extracellular, intracellular)
p{5} = [1 0 0]';                           % IgG antibodies (absent, present, neutralising)


% normalise initial marginals
%--------------------------------------------------------------------------
Nf    = numel(p);
for f = 1:Nf
    p{f}  = p{f}/sum(p{f});
end

x     = spm_cross(p);

for i = 1:M.T
    % update ensemble density, with probability dependent transitions
    %----------------------------------------------------------------------
    B     = spm_immune_B(x,P);
    x     = spm_unvec(B*spm_vec(x),x);
    x     = x/sum(x(:));
    
    % probabilistic mappings: outcomes based on marginal densities (p)
    %======================================================================
    p     = spm_marginal(x);
    for j = 1:Nf
        X{j}(i,:) = p{j};
    end
    
    Y(i,1) = 128*Q.Abm*(p{1}(2) + p{1}(3)); % IgM antibody
    Y(i,2) = 128*Q.Abm*(p{5}(2) + p{5}(3)); % IgG antibody
    Y(i,3) = 10000*Q.VLm*p{4}(2);           % Viral load (extracellular)
    Y(i,4) = Q.INF*p{3}(2);                 % IFN gamma CD4+
    Y(i,5) = Q.INF*p{3}(3);                 % IFN gamma CD8+
end
Y = Y(U*24,:);


function T = spm_immune_B(x,P)


% marginal probabilities
%==========================================================================
p     = spm_marginal(x);

% identity matrices
%--------------------------------------------------------------------------
dim   = size(x);
I     = cell(ndims(x),1);
for i = 1:ndims(x)
    I{i} = speye(dim(i));
end

% exponentiate parameters
%--------------------------------------------------------------------------
P     = spm_vecfun(P,@exp);

% upper bound probabilities
%--------------------------------------------------------------------------

% probabilistic transitions: antibodies
%==========================================================================

% marginal: IgM antibodies {1} | B-cells {2}(1) naive
%--------------------------------------------------------------------------
%       absent           present       neutralising  
%--------------------------------------------------------------------------
b{1} = [1                P.dAb         P.dAb; 
        0                1 - P.dAb     0;
        0                0             1 - P.dAb];

% marginal: IgM antibodies {1} | B-cells {2}(2) plasma (IgM)
%--------------------------------------------------------------------------
%       absent           present       neutralising  
%--------------------------------------------------------------------------
b{2} = [1-P.pAb          P.dAb         P.dAb; 
        P.pAb*(1-P.nAb)  1 - P.dAb     0;
        P.pAb*P.nAb      0             1 - P.dAb];
    
% marginal: IgM antibodies {1} | B-cells {2}(3) inactive memory
%--------------------------------------------------------------------------
%       absent           present       neutralising  
%--------------------------------------------------------------------------
b{3} = [1                P.dAb         P.dAb; 
        0                1 - P.dAb     0;
        0                0             1 - P.dAb];
    
% marginal: IgM antibodies {1} | B-cells {2}(4) active memory
%--------------------------------------------------------------------------
b{4} = b{3};

% marginal: IgM antibodies {1} | B-cells {2}(5) plasma (IgG)
%--------------------------------------------------------------------------
b{5} = b{4};

% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{3},I{4},I{5}});
B{1} = spm_permute_kron(b,dim([1,2,3,4,5]),[1,2,3,4,5]);
clear b

% probabilistic transitions: B-cells
%==========================================================================

% marginal: B-cells {2} | T-cells {3}(1) naive
%--------------------------------------------------------------------------
%       naive        IgM plasma     inactive memory     active memory                      IgG plasma
%--------------------------------------------------------------------------
b{1} = [1            P.dpB          0                   P.dmB*(1-p{4}(2)*P.Bmp)            P.dpB;
        0            1 - P.dpB      0                   0                                  0;
        0            0              1-P.mba             0                                  0;
        0            0              P.mba               (1 - P.dmB)*(1 - p{4}(2)*P.Bmp)    0;
        0            0              0                   p{4}(2)*P.Bmp                      1 - P.dpB];

% marginal: B-cells {2} | T-cells {3}(2) CD4+
%--------------------------------------------------------------------------

%       naive        IgM plasma      inactive memory    active memory                      IgG plasma  
%--------------------------------------------------------------------------
b{2} = [1-P.BC          P.dpB        0                  P.dmB*(1-p{4}(2)*P.Bmp)            P.dpB; 
        P.BC*(1-P.BCm)  1 - P.dpB    0                  0                                  0;
        P.BC*P.BCm      0            1-P.mba            0                                  0;
        0               0            P.mba              (1 - P.dmB)*(1 - p{4}(2)*P.Bmp)    0;
        0               0            0                  p{4}(2)*P.Bmp                      1 - P.dpB];
    
% marginal: B-cells {2} | T-cells {3}(3) CD8+
%--------------------------------------------------------------------------
%       naive    plasma       memory  
%--------------------------------------------------------------------------
b{3} = b{1};
    
% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{1},I{4},I{5}});
B{2} = spm_permute_kron(b,dim([2,3,1,4,5]),[3,1,2,4,5]);
clear b

% probabilistic transitions: T-cells
%==========================================================================

% marginal: T-cells {3} | Pathogen {4}(1) absent
%--------------------------------------------------------------------------
%       naive        CD4+          CD8+  
%--------------------------------------------------------------------------
b{1} = [1            P.dT4         P.dT8; 
        0            1 - P.dT4     0;
        0            0             1 - P.dT8];

% marginal: T-cells {3} | Pathogen {4}(2) extracellular
%--------------------------------------------------------------------------
%       naive        CD4+          CD8+  
%--------------------------------------------------------------------------
b{2} = [1-P.CD4      P.dT4         P.dT8; 
        P.CD4        1 - P.dT4     0;
        0            0             1 - P.dT8];
    
% marginal: T-cells {3} | Pathogen {4}(3) intracellular
%--------------------------------------------------------------------------
%       naive        CD4+       CD8+  
%--------------------------------------------------------------------------
b{3} = [1-P.CD8      P.dT4         P.dT8; 
        0            1 - P.dT4     0;
        P.CD8        0             1 - P.dT8];
    
% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{1},I{2},I{5}});
B{3} = spm_permute_kron(b,dim([3,4,1,2,5]),[3,4,1,2,5]);
clear b

% probabilistic transitions: pathogen
%==========================================================================

% marginal: Pathogen {4} | IgM antibodies {1}(1) absent, T-cells {3}(1) naive 
%--------------------------------------------------------------------------
%       absent             extracellular                 intracellular  
%--------------------------------------------------------------------------
b{1,1} = [1-P.ext*p{4}(3)  P.dpe*(1-p{5}(3))+p{5}(3)             P.dpi; 
          P.ext*p{4}(3)    (1-P.dpe)*(1-P.int)*(1-p{5}(3))           0;
          0                (1-P.dpe)*P.int*(1-p{5}(3))        1-P.dpi];

% marginal: Pathogen {4} | IgM antibodies {1}(2) present, T-cells {3}(1) naive
%--------------------------------------------------------------------------
b{1,2} = b{1,1};
    
% marginal: Pathogen {4} | IgM antibodies {1}(3) neutralising, T-cells {3}(1) naive
%--------------------------------------------------------------------------
%       absent              extracellular    intracellular
%--------------------------------------------------------------------------
b{1,3} = [1-P.ext*p{4}(3)   1                        P.dpi; 
          P.ext*p{4}(3)     0                            0;
          0                 0                     1-P.dpi];

% marginal: Pathogen {4} | IgM antibodies {1}(1) absent, T-cells {3}(2) CD4+
%--------------------------------------------------------------------------
%         absent           extracellular                      intracellular
%--------------------------------------------------------------------------
b{2,1} = [1-P.ext*p{4}(3)  (1-p{5}(3))*(P.dpe + P.T4P)/2+p{5}(3)      P.dpi; 
          P.ext*p{4}(3)    (1-p{5}(3))*(2-P.dpe - P.T4P)*(1-P.int)/2      0;
          0                (1-p{5}(3))*(2-P.dpe - P.T4P)*P.int/2   1-P.dpi];
      
% marginal: Pathogen {4} | IgM antibodies {1}(1) present, T-cells {3}(2) CD4+
%--------------------------------------------------------------------------
b{2,2} = b{2,1};
% marginal: Pathogen {4} | IgM antibodies {1}(1) neutralising, T-cells {3}(2) CD4+
%--------------------------------------------------------------------------
%       absent              extracellular    intracellular
%--------------------------------------------------------------------------
b{2,3} = [1-P.ext*(p{4}(3)+1)/2  1                  P.dpi; 
          P.ext*(p{4}(3)+1)/2    0                      0;
          0                      0               1-P.dpi];
      
% marginal: Pathogen {4} | IgM antibodies {1}(1) absent, T-cells {3}(2) CD8+
%--------------------------------------------------------------------------
%         absent              extracellular                        intracellular
%--------------------------------------------------------------------------
b{3,1} = [1-P.ext*p{4}(3)         (1-p{5}(3))*P.dpe+p{5}(3)       (P.TCP + P.dpi)/2; 
          P.ext*p{4}(3)           (1-p{5}(3))*(1-P.dpe)*(1-P.int)                 0;
          0                       (1-p{5}(3))*(1-P.dpe)*P.int (2 - P.TCP -P.dpi)/2];
      
% marginal: Pathogen {4} | IgM antibodies {1}(1) present, T-cells {3}(2) CD8+
%--------------------------------------------------------------------------
b{3,2} = b{3,1};

% marginal: Pathogen {4} | IgM antibodies {1}(1) neutralising, T-cells {3}(2) CD8+
%--------------------------------------------------------------------------
%       absent               extracellular         intracellular
%--------------------------------------------------------------------------
b{3,3} = [1-P.ext*p{4}(3)         1                 (P.TCP+P.dpi)/2; 
          P.ext*p{4}(3)           0                               0;
          0                       0           (2 - P.TCP -P.dpi)/2];
      
% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b1    = spm_cat(spm_diag(b(1,:)));
b2    = spm_cat(spm_diag(b(2,:)));
b3    = spm_cat(spm_diag(b(3,:)));
b     = spm_cat(spm_diag({b1,b2,b3}));
b     = spm_kron({b,I{2},I{5}});
B{4}  = spm_permute_kron(b,dim([4,1,3,2,5]),[2,4,3,1,5]);
clear b


% probabilistic transitions: IgG antibodies
%==========================================================================

% marginal: IgG antibodies {1} | B-cells {2}(1) naive
%--------------------------------------------------------------------------
%       absent    present       neutralising  
%--------------------------------------------------------------------------
b{1} = [1         P.dAbG         P.dAbG; 
        0         1 - P.dAbG     0;
        0         0             1 - P.dAbG];

% marginal: IgG antibodies {1} | B-cells {2}(2) plasma (IgM)
%--------------------------------------------------------------------------
b{2} = b{1};
    
% marginal: IgG antibodies {1} | B-cells {2}(3) inactive memory
%--------------------------------------------------------------------------
%       absent    present       neutralising  
%--------------------------------------------------------------------------
b{3} = [1         P.dAbG         P.dAbG; 
        0         1 - P.dAbG     0;
        0         0             1 - P.dAbG];
    
% marginal: IgG antibodies {1} | B-cells {2}(4) active memory
%--------------------------------------------------------------------------
b{4} = b{3};

% marginal: IgG antibodies {1} | B-cells {2}(5) plasma (IgG)
%--------------------------------------------------------------------------
b{5} = [1-P.pAb          P.dAbG         P.dAbG; 
        P.pAb*(1-P.nAb)  1 - P.dAbG     0;
        P.pAb*P.nAb      0             1 - P.dAbG];

% kroneckor form (taking care to get the order of factors right)
%--------------------------------------------------------------------------
b    = spm_cat(spm_diag(b));
b    = spm_kron({b,I{3},I{4},I{1}});
B{5} = spm_permute_kron(b,dim([5,2,3,4,1]),[5,2,3,4,1]);
clear b

% probability transition matrix
%==========================================================================
T     = 1;
for i = 1:numel(B)
    T =  T*B{i};
end

