function [P,C] = spm_immune_priors
% Default priors for immune model
% FORMAT [P,C] = spm_immune_priors
% P - Prior expectations
% C - Prior covariances
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging
 
% Thomas Parr
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Log scale parameters for initial states
%--------------------------------------------------------------------------

P.n   = -2;            % Scaling for pre-existing (cross-reactive) immunity
P.m   =  0;            % Initial viral load

% Parameters of transition probablities (log probabilities)
%--------------------------------------------------------------------------
P.dAb = -6;              % Decay of antibodies over time (IgM)
P.dAbG= -8;              % Decay of antibodies over time (IgG)
P.pAb = -0;              % Production of antibodies by plasma cells
P.nAb = -1;              % Proportion of antibodies that are neutralising
P.dpB = -2;              % Decay of plasma cells over time
P.dmB = -32;             % Decay of memory cells over time
P.mba = -4;              % Activation of memory cells
P.Bmp = -1/2;            % Proportion of memory cells converting to plasma cells in presence of infection
P.BC  = -3;              % Specialisation of B-cells in presence of CD4+ T-cells
P.BCm = -1;              % Proportion of specialised B-cells becoming memory cells
P.dT4 = -4;              % Decay of CD4+ T-cells over time
P.dT8 = -5;              % Decay of CD8+ T-cells over time
P.CD4 = -3;              % Production of CD4+ T-cells with extracellular pathogen
P.CD8 = -3;              % Production of CD8+ T-cells with intracellular pathogen
P.TCP = -2;              % Neutralisation of intracellular pathogens by T-Cell mediated apoptosis
P.T4P = -8;              % Neutralisation of extracellular pathogen by antibody-independent CD4+ response
P.int = -2;              % Viral entry into cells
P.ext = -1;              % Viral shedding into extracellular space (absorbing replication rates)
P.dpe = -2;              % Decay of extracellular pathogen
P.dpi = -2;              % Decay of intracellular pathogen

% Measurement variables (log scale)
%--------------------------------------------------------------------------

P.Abm = 0;             % Scale for antibody titre
P.VLm = 0;             % Scale for viral load
P.INF = 0;             % Scaling for interferon measures

% Assume informative priors to preclude probabilities exceeding one
%--------------------------------------------------------------------------
C      = spm_unvec(ones(spm_length(P),1)/64,P);

% Relatively uninformative priors over initial states and measurement variables
%--------------------------------------------------------------------------
C.m    = 1/16;
C.Abm  = 1/16;
C.VLm  = 1/16;
C.INF  = 1/16;

% Very uninformative priors over variables of interest
%--------------------------------------------------------------------------
% The values here are chosen so that the alternative models tested in the
% DEM_Immune routine each involve a move of a single parameter
% approximately one standard deviation from the priors set out here. This
% is to ensure the accuracy is the key determinant of the model fit.

C.int  = 1;
C.n    = 16;
C.TCP  = 7/2;
