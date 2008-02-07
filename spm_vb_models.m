function [model] = spm_vb_models(SPM,factor)
% Specify models for ANOVAs
% FORMAT [model] = spm_vb_models(SPM,factor)
%
% SPM    - SPM structure
% factor - Structure specifying factors and levels
%          factor(i).name                  Name of ith factor
%          factor(i).levels                Number of levels
%          It is assumed that the levels of the first factor change
%          slowest with condition
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_vb_models.m 1143 2008-02-07 19:33:33Z spm $

%-Number of factors
%-----------------------------------------------------------------------
nf = length(factor);

%-Number of levels
%-----------------------------------------------------------------------
k1 = factor(1).levels;
if nf==2
    k2 = factor(2).levels;
else
    k2 = 1;
end

%-Number of conditions
%-----------------------------------------------------------------------
NC = length(SPM.Sess(1).U); 

if k1*k2 ~= NC
    error('Error in spm_vb_models: factors do not match conditions');
end

%-Get contingency table
% factor 1 levels on rows, factor 2 levels on columns
%-----------------------------------------------------------------------
ctable = reshape([1:1:NC],k2,k1);

dt = SPM.Sess(1).U(1).dt;
P  = SPM.Sess(1).U(1).P;

%-Create full model (same as original)
%-----------------------------------------------------------------------
model(6).name = 'Full';
model(6).U    = SPM.Sess(1).U;
model(6).X    = SPM.xX.X;

%-Create null model
%-----------------------------------------------------------------------
model(1).name = 'NULL';
null_cols     = [SPM.xX.iH,SPM.xX.iB,SPM.xX.iG];
model(1).X    = SPM.xX.X(:,null_cols);

%-Create average model
%-----------------------------------------------------------------------
model(2).name = 'Average';
model(2).U(1).name{1} = 'Average';
model(2).U(1).dt = dt;
ons = []; dur = [];
for i=1:NC
    %-Concatenate all conditions
    %-------------------------------------------------------------------
    ons = [ons; SPM.Sess(1).U(i).ons];
    dur = [dur; SPM.Sess(1).U(i).dur];
end
model(2).U(1).ons = ons;
model(2).U(1).dur = dur;
model(2).U(1).P   = P;

%-For two factor models, create model for factor A, B and both factors
%-----------------------------------------------------------------------
if nf == 2
    
    %-Create model for factor A
    %-------------------------------------------------------------------
    model(3).name = factor(1).name;
    nA = factor(1).levels;
    for i=1:nA
        %-Has as many inputs as levels of factor A
        %---------------------------------------------------------------
        ons = []; dur = [];
        nAi = factor(2).levels;
        for j=1:nAi
            %-Concatenate inputs from relevant conditions
            %-----------------------------------------------------------
            c   = ctable(i,j);
            ons = [ons; SPM.Sess(1).U(c).ons];
            dur = [dur; SPM.Sess(1).U(c).dur];
        end
        model(3).U(i).name{1} = ['Level ' int2str(i)];
        model(3).U(i).dt      = dt;
        model(3).U(i).ons     = ons;
        model(3).U(i).dur     = dur;
        model(3).U(i).P       = P;
    end
    
    %-Create model for factor B
    %-------------------------------------------------------------------
    model(4).name = factor(2).name;
    nB = factor(2).levels;
    for i=1:nB
        %-Has as many inputs as levels of factor B
        %---------------------------------------------------------------
        ons = []; dur = [];
        nBi = factor(1).levels;
        for j=1:nBi
            %-Concatenate inputs from relevant conditions
            %-----------------------------------------------------------
            c   = ctable(j,i);
            ons = [ons; SPM.Sess(1).U(c).ons];
            dur = [dur; SPM.Sess(1).U(c).dur];
        end
        model(4).U(i).name{1} = ['Level ',int2str(i)];
        model(4).U(i).dt      = dt;
        model(4).U(i).ons     = ons;
        model(4).U(i).dur     = dur;
        model(4).U(i).P       = P;
    end
    
    %-Create model for both factors
    %-------------------------------------------------------------------
    model(5).name = 'Both factors';
    nboth = nA + nB;
    i = 1;
    for k=1:nA
        model(5).U(i) = model(3).U(k);
        i = i + 1;
    end
    % Note: to avoid rank deficiency don't estimate 
    % average response to level nB of factor 2 as 
    % this is given by a linear combination of responses in
    % other levels eg. for 2by2: B2= A1+A2-B1
    for k=1:nB-1
        model(5).U(i) = model(4).U(k);
        i = i + 1;
    end
    
end
