% runs DEMOs to check DEM scheme
%__________________________________________________________________________

clear functions
clear
clf
set(gcf,'name','Dynamic Expectation Maximisation')

file = {
    'DEM_demo_GLM',...
    'DEM_demo_PEB',...
    'DEM_demo_factor_analysis',...
    'DEM_demo_convolution',...
    'DEM_demo_EM',...
    'DEM_demo_Kalman_filtering',...
    'DEM_demo_particle_filtering',...
    'DEM_demo_Lorenz'};

for i = 1:length(file)
    h = help(file{i});
    str = {'proceed with' file{i} '_______________________________' ' ' h};
    if ~strcmp(questdlg(str),'Yes'),
        return
    else
        clear DEM M
    end
    eval(file{i})
end