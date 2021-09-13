function [P,C,str] = spm_CLIMATE_priors
% Generate prior expectation and covariance log parameters
% FORMAT [P,C,str] = spm_CLIMATE_priors
% 
% pE          - prior expectation (structure)
% pC          - prior covariances (structure)
% str.outcome - names
% str.states  - names
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_CLIMATE_priors.m 8151 2021-09-13 09:12:37Z karl $

% names
%==========================================================================

str.outcome = {...
    'Average_temperature'                        ;
    'normalisedDrought_area_underSPEI_1_6'       ;
    'normalisedDrought_area_underSPEI_2'         ;
    'rainfall_mm'                                ;
    'yield_z_scores'                             ;
    'crop_production_ktons'                      ;
    'Irrigated_total_area'                       ;
    'Fertiliser_Consumption_tons'                ;
    'Milk__000T_'                                ;
    'Egg_millions_'                              ;
    'Meat__000T_'                                ;
    'price_z_scores'                             ;
    'PRIMARYSECTORCONSTANTPRICES_MillionsInRs_'  ;
    'TOTALCONSTANTPRICES_MillionsInRs_'          ;
    'PRIMARYSECTORCONSTANTSHARES_Percent_'       ;
    'STATEMALEAVERAGEFIELDLABOUR_RsPerDay_'      ;
    'growth_failure_5_14yo_dalys_per_100k'       ;
    'child_materNaNl_malnutrition_dalys_per_100k';
    'food_subsidy_million_R'                     };

str.states = {...
    'Meteorological (1)'                         ;
    'Meteorological (2)'                         ;
    'Meteorological (3)'                         ;
    'Anthropological activity'                   ;
    'Primary sector activity'                    ;
    'Precipitation'                              ;
    'crop production'                            ;
    'Irrigation'                                 ;
    'Fertilisation'                              ;
    'Food production'                            ;
    'Food price'                                 ;
    'Malnutrition'                               };

str.parameters = {...
    'Meteorological (1)'                         ;
    'Meteorological (2)'                         ;
    'Meteorological (3)'                         };



% initial conditions
%==========================================================================
P.x(1)  = 8;               % 'Meteorological (1)'                         ;
P.x(2)  = 8;               % 'Meteorological (2)'                         ;
P.x(3)  = 24;              % 'Meteorological (3)'                         ;
P.x(4)  = 1;               % 'Anthropological activity'                   ;
P.x(5)  = 1;               % 'Primary sector activity'                    ;
P.x(6)  = 0;               % 'Precipitation'                              ;
P.x(7)  = 0;               % 'Crop production'                            ;
P.x(8)  = 0;               % 'Irrigation'                                 ;
P.x(9)  = 0;               % 'Fertilisation'                              ;
P.x(10) = 0;               % 'Food production'                            ;
P.x(11) = 0;               % 'Food price'                                 ;
P.x(12) = 0;               % 'Malnutrition'                               ;

% Latent states
%==========================================================================

% Meteorological
%--------------------------------------------------------------------------
P.L(1)  = 10;              % Lorenz parameters [sig, beta, rho]
P.L(2)  = 8/3;             % Lorenz parameters [sig, beta, rho]
P.L(3)  = 24;              % Lorenz parameters [sig, beta, rho]
P.L(4)  = 17;              % Lorenz parameters (time)

% Anthropological activity
%--------------------------------------------------------------------------
P.R(1)  = 1/(12*72);       % P(dying)
P.R(2)  = (1.5)/(12*72);   % P(reproducing)
P.R(3)  = (1/32)/(12*72);   % Malnutrition
P.R(4)  = (1/4)/(12*72);     % Malnutrition (Primary)

P.LR    = 8;               % coupling

% Precipitation
%--------------------------------------------------------------------------
P.P(1)  = 32;              % bias
P.P(2)  = 1/128;           % sensitivity
P.P(3)  = 1/8;             % decay

% Irrigation
%--------------------------------------------------------------------------
P.I(1)  = 1/2;             % bias
P.I(2)  = 1/32;            % sensitivity
P.I(3)  = 1/8;             % decay

% Food production
%--------------------------------------------------------------------------
P.F(1)  = 1/2;             % sensitivity (anthropomorphic)
P.F(2)  = 1/16;            % decay
P.F(3)  = 1/512;           % sensitivity (climate)

P.FL(1) = 32;              % bias
P.FL(2) = 1/128;           % sensitivity

% Crop production
%--------------------------------------------------------------------------
P.C(1)  = 1/2;             % Precipitation
P.C(2)  = 1/2;             % Irrigation
P.C(3)  = 1;               % Fertilisation
P.C(4)  = 1/8;             % decay

% Fertilisation
%--------------------------------------------------------------------------
P.E(1)  = 1;               % demand
P.E(2)  = 1/16;            % decay

% Food price
%--------------------------------------------------------------------------
P.Z(1)  = 1/16;            % Crop production
P.Z(2)  = 1/16;            % Fertilisation
P.Z(3)  = 1/16;            % Food production
P.Z(4)  = 1/4;             % decay

% Malnutrition
%--------------------------------------------------------------------------
P.M(1)  = 1/32;            % sensitivityFood price
P.M(4)  = 1/16;            % decay

% outputs
%==========================================================================
P.y(1,1)  = 298;
P.y(2,1)  = 0;
P.y(3,1)  = 0;
P.y(4,1)  = 1;

P.y(5,1)  = -1/2;
P.y(6,1)  = 1;
P.y(7,1)  = 4;
P.y(8,1)  = -2/3;





% Variances (mildly informative priors)
%==========================================================================
C   = spm_vecfun(P,@(x)spm_zeros(x) + 1/16);


% fixed parameters
%--------------------------------------------------------------------------
C.L(1:2) = 0;                 % parameters [sig, beta]

% log transform
%==========================================================================
Q   = P;
Q   = spm_vecfun(P,@log);
Q.x = P.x;
Q.y = P.y;
P   = Q;


return