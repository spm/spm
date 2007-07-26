% Dynamic expectation maximisation (DEM) is a variational treatment of
% hierarchical, nonlinear dynamic or static models.  It uses a fixed-form
% Laplace assumption to approximate the conditional, variational or
% ensemble density of unknown states and parameters.  This is an
% approximation to the density that would obtain from Variational Filtering
% (VF) in generalized coordinates of motion. We start with a demonstration
% for VF using a simple convolution model and compare the results with a
% DEM.  We then demonstrate the inversion of increasingly complicated models;
% ranging from a simple General Linear Model to a Lorenz attractor.  It is
% anticipated that the reader will examine the routines called to fully
% understand the nature of the scheme.
%__________________________________________________________________________
 
clear str
 
% splash
%--------------------------------------------------------------------------
h      = help('DEM_demos');
str{1} = h(h ~= sprintf('\n'));
str{2} = ' ';
str{3} = '(GNU) Copyright (c) 2005 The Wellcome Trust Centre for Neuroimaging.';
if ~strcmp(questdlg(str),'Yes')
    return
end
 
file = {
    'DEM_demo_DFP',...
    'DEM_demo_GLM',...
    'DEM_demo_PEB',...
    'DEM_demo_factor_analysis',...
    'DEM_demo_convolution',...
    'DEM_demo_EM',...
    'DEM_demo_DEM',...
    'DEM_demo_Kalman_filtering',...
    'DEM_demo_particle_filtering',...
    'DEM_demo_Lorenz'};
 
for i = 1:length(file)
    h      = help(file{i});
    str{3} = h(h ~= sprintf('\n'));
    str{2} = ...
        '__________________________________________________________________________';
    str{1} = ['proceed with: ' file{i}];
    if ~strcmp(questdlg(str),'Yes')
        return
    else
        clear DEM M
    end
    eval(file{i})
end
