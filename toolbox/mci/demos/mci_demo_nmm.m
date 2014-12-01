
clear all
close all

% Number of parameters to estimate
Np=2;
%Np=6;
%Np=21;

disp('Synthetic data from two region Neural Mass Model');
disp(sprintf('Model has %d parameters',Np));

% Backward connection ?
back=1;

% Observation noise SD
sd=0.01;

[M,U] = mci_nmm_struct (back,sd,Np);

% Parameters
switch Np
    case 2
        disp('Estimating extrinsic connections');
        P=[1,1]';
    case 6
        disp('Estimating intrinsic and extrinsic connections');
        P=[1,1,0,0,1,0]';
    case 21,
        disp('Estimating all connections');
        % All params set to prior mean except f/b
        P=M.pE;
        P.A{1}(2,1)=1; % Forward connection
        P.A{2}(1,2)=1; % Backward connection
end

Y = mci_nmm_gen (M,U,P);

mci_plot_outputs(M,Y);

%inference='mh';
%inference='vl';
inference='langevin';

verbose=1;
post = spm_mci_post (inference,M,U,Y,P,verbose);