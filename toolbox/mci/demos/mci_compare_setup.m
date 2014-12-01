function [P,M,U,Y,ind] = mci_compare_setup (model)
% Set up data structures for fwd/sens/grad comparisons
% FORMAT [P,M,U,Y,ind] = mci_compare_setup (model)
%
% model     'phase','fitz','nmm-r2p2','nmm-r2p21','nmm1','rlds'
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_compare_setup.m 6275 2014-12-01 08:41:18Z will $

switch model,
    case 'phase',
        d=7;
        disp(sprintf('Weakly coupled oscillator network : d=%d regions',d));
        [P,M,U,Y] = mci_phase_init(d);
        ind=[1:M.Np-M.n-1];
        
    case 'fitz',
        disp('Fitzhugh-Nagumo model');
        [P,M,U,Y] = mci_fitz_init();
        ind=[1:M.Np];
    
    case 'nmm-r2p2',
        disp('Two-region, two-parameter neural mass model');
        back=1;
        sd=0.01;
        Np=2;
        
        [M,U] = mci_nmm_struct(back,sd,Np);
        P=[1,1]';
        Y = mci_nmm_gen(M,U,P);
        ind=[1:M.Np];
        
    case 'nmm-r2p21',
        disp('Two-region, 21-parameter neural mass model');
        back=1;
        sd=0.01;
        Np=21;
        
        [M,U] = mci_nmm_struct(back,sd,Np);
        P=M.pE;
        P.A{1}(2,1)=1; % Forward connection
        P.A{2}(1,2)=1; % Backward connection
        Y = mci_nmm_gen(M,U,P);
        
        ind=[1:Np];
    
    case 'nmm1',
        disp('Single-region, ten-parameter neural mass model');
        sd=0.01;
        
        [P,M,U,Y] = nmm1_init(sd);
        ind=[1:M.Np];
        
    case 'rlds',
        disp('Linear Dynamical System with real modes');
        d=10;
        sd=0.01;
        R0=linspace(0.5,5,d)';
        P=linspace(-0.5,-3,d)';
        [M,U,Y] = irlds_init (d,sd,R0,P);
        %[P,M,U,Y] = rlds_init (d,sd);
        
        ind=[1:M.Np];
        
    otherwise
        disp('Unknown model');
        return
end
