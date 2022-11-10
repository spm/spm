function [P,C,str] = spm_CLIMATE_priors
% Prior expectation and covariance of parameters for a climate model
% FORMAT [P,C,str] = spm_CLIMATE_priors
%
% pE          - prior expectation (structure)
% pC          - prior covariances (structure)
% str.outcome - names
% str.states  - names
%
% This routine generates the prior density over model parameters in terms
% of a prior expectation and covariance structure. Crucially, there are
% three kinds of parameters. The first sets the initial values of the
% latent states. The second comprises the parameters of the equations of
% motion or flow of latent states correspond to the dynamic part of the
% model. The third kind of parameters map from the latent states to
% observable outcomes.
% 
% pE.x    - initial states 
% pE.P    - flow parameters
% pE.Y    - outcome parameters
%
% Because the flow parameters are (almost universally) rate or time
% constants, they are scale parameters. In other words, they are always
% greater than zero. This means that during estimation we will deal with
% log scale parameters that can take any value between plus and minus
% infinity. This allows one to place gaussian priors over nonnegative
% (scale) parameters. Practically, this means that this routine returns the
% logarithm of the flow parameters used to generate dynamics.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% names of outcome and latent states saved in a structure for potting
%==========================================================================
str.outcome = {...
    'Average temperature';
    'Extreme drought';
    'Exceptional drought';
    'Rainfall';
    'Crop production';
    'Irrigation';
    'Fertiliser use';
    'Milk production';
    'Food prices';
    'Income in exposed sector';
    'Childhood malnutrition';
    'Crop yield'};
   
str.states = {...
    'Meteorological (1)'                         ;
    'Meteorological (2)'                         ;
    'Meteorological (3)'                         ;
    'Anthropological activity'                   ;
    'Primary sector activity'                    ;
    'Yeild'                                      ;
    'Crop production'                            ;
    'Irrigation'                                 ;
    'Fertilisation'                              ;
    'Food production'                            ;
    'Food price'                                 ;
    'Malnutrition'                               };


% initial conditions (i.e., latent states at the start date): P.x
%==========================================================================
P.x(1)  = 0;               % 'Meteorological (fast)'                      ;
P.x(2)  = 1;               % 'Meteorological (first)'                     ;
P.x(3)  = 0;               % 'Meteorological (slow)'                      ;

P.x(4)  = 1;               % 'Anthropological activity'                   ;
P.x(5)  = 1;               % 'Primary sector activity'                    ;
P.x(6)  = 6;               % 'Yield'                                      ;
P.x(7)  = 1;               % 'Crop production'                            ;
P.x(8)  = 0;               % 'Irrigation'                                 ;
P.x(9)  = 1;               % 'Fertilisation'                              ;
P.x(10) = 0;               % 'Food production'                            ;
P.x(11) = 0;               % 'Food price'                                 ;
P.x(12) = 0;               % 'Malnutrition'                               ;

% flow parameters for updating latent states: P.P
%==========================================================================

% Meteorological
%--------------------------------------------------------------------------
P.P.L(1)  = 1/1700;          % (1) sensitivity to anthropomorphic activity
P.P.L(2)  = 1/6;             % (2) equilibrium point of slow state
P.P.L(3)  = 1/4098;          % (3) decay of slow state

% Anthropological activity
%--------------------------------------------------------------------------
P.P.R(1)  = 1/(12*72);       % (4)  P(dying)
P.P.R(2)  = (1.4)/(12*72);   % (5)  P(reproduction)
P.P.R(3)  = (1/8)/(12*72);   % (6)  sensitivity to food prices
P.P.R(4)  = (1/32)/(12*72);  % (7)  sensitivity to malnutrition

% Precipitation (yield)
%--------------------------------------------------------------------------
P.P.P(1)  = 8;               % (8)  bias
P.P.P(2)  = 8;               % (9)  sensitivity
P.P.P(3)  = 1/8;             % (10) decay

% Irrigation
%--------------------------------------------------------------------------
P.P.I(1)  = 4;               % (11) bias
P.P.I(2)  = 2;               % (12) sensitivity
P.P.I(3)  = 1/4;             % (13) decay

% Food production
%--------------------------------------------------------------------------
P.P.F(1)  = 8;               % (14) sensitivity (demand)
P.P.F(2)  = 1/2;             % (15) decay

% Crop production
%--------------------------------------------------------------------------
P.P.C(1)  = 8;               % (16) threshold
P.P.C(2)  = 1;               % (17) precision 
P.P.C(3)  = 8;               % (18) sensitivity to precipitation
P.P.C(4)  = 1/128;           % (19) sensitivity to irrigation
P.P.C(5)  = 1/2;             % (20) decay

% Resources (fertilsation)
%--------------------------------------------------------------------------
P.P.E(1)  = 1/4;             % (21) sensitivity to demand
P.P.E(2)  = 1/64;            % (22) decay

% Food price
%--------------------------------------------------------------------------
P.P.Z(1)  = 1/8;             % (23) sensitivity (fertilisation)
P.P.Z(2)  = 1/32;            % (24) sensitivity (food production)
P.P.Z(3)  = 1/8;             % (25) decay

% Malnutrition
%--------------------------------------------------------------------------
P.P.M(1)  = 1/32;            % (26) sensitivity to food price
P.P.M(2)  = 1;               % (27) sensitivity to food production
P.P.M(3)  = 1/8;             % (28) decay


% outcome parameters: prior expectations: P.Y
%==========================================================================

% temperature (U = 1)
%--------------------------------------------------------------------------
P.Y.t(1)  = 0;               % constant
P.Y.t(2)  = 1/128;           % coefficients
P.Y.t(3)  = 0;               % coefficients
P.Y.t(4)  = 0;               % coefficients

% Drought (U = 2,3)
%--------------------------------------------------------------------------
P.Y.d(1)  = 2;               % constant
P.Y.d(2)  = 1/2;             % coefficient (extreme)
P.Y.d(3)  = 3;               % constant
P.Y.d(4)  = 2/3;             % coefficient (exceptional)

% Rainfall (U = 4)
%--------------------------------------------------------------------------
P.Y.r(1)  = 0;               % constant
P.Y.r(2)  = 1;               % coefficient
P.Y.r(3)  = 0;               % coefficient
P.Y.r(4)  = 0;               % coefficient

% Crop production (U = 5)
%------------------------------------------------------------------
P.Y.c(1)  = 5;               % constant
P.Y.c(2)  = 4;               % coefficient

% Irrigation (U = 6)
%------------------------------------------------------------------
P.Y.i(1)  = 7;               % constant
P.Y.i(2)  = 1/2;             % coefficient

% Fertiliser use (U = 7)
%------------------------------------------------------------------
P.Y.u(1)  = 5;               % constant
P.Y.u(2)  = 1/4;             % coefficient

% Milk production (U = 8)
%------------------------------------------------------------------
P.Y.m(1)  = 6;               % constant
P.Y.m(2)  = 1;               % coefficient (crop production)
P.Y.m(3)  = 4;               % coefficient (food production)

% Food prices (U = 9)
%------------------------------------------------------------------
P.Y.f(1)  = -1/2;            % constant
P.Y.f(2)  = 1/4;             % coefficient

% Income in exposed sector (U = 10)
%------------------------------------------------------------------
P.Y.e(1)  = -2;              % constant
P.Y.e(2)  = 2;               % coefficient (crop production)
P.Y.e(3)  = 12;              % coefficient (food production)

% Childhood malnutrition (U = 11)
%------------------------------------------------------------------
P.Y.p(1)  = 12;              % constant
P.Y.p(2)  = 4;               % coefficient

% yield (U = 12)
%------------------------------------------------------------------
P.Y.y(1)  = -3;              % constant
P.Y.y(2)  = 1/2;             % coefficient


% Prior variances; specified as a structure 
%==========================================================================
C.x   = spm_vecfun(P.x,@(x)spm_zeros(x) + 0);     % very informative priors
C.P   = spm_vecfun(P.P,@(x)spm_zeros(x) + 1/64);  % mild shrinkage priors
C.Y   = spm_vecfun(P.Y,@(x)spm_zeros(x) + 1/64);  % uninformative priors

% fixed parameters
%--------------------------------------------------------------------------
% log transform flow parameters
%==========================================================================
P.P   = spm_vecfun(P.P,@log);

% The shrinkage priors on the flow parameters have a particular
% quantitative meaning because we are dealing with log parameters. This
% means that the prior variance controls the range of relative values that
% are plausibly expected. For example, if we think parameters can vary by
% less than an order of magnitude, then we would set their prior variance
% to about 1/16.


return