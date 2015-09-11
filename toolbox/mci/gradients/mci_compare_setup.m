function [P,M,U,Y,ind] = mci_compare_setup (model)
% Set up data structures for fwd/sens/grad comparisons
% FORMAT [P,M,U,Y,ind] = mci_compare_setup (model)
%
% model     'phase', 'nmm-r2p2'
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_compare_setup.m 6548 2015-09-11 12:39:47Z will $

switch model,
    
    case 'phase',
        d=7;
        disp(sprintf('Weakly coupled oscillator network : d=%d regions',d));
        [P,M,U,Y] = mci_phase_init(d);
        ind=[1:M.Np-M.n-1];
    
    case 'nmm-r2p2',
        disp('Two-region, two-parameter neural mass model');
        back=1;
        sd=0.01;
        Np=2;
        
        [M,U] = mci_nmm_struct(back,sd,Np);
        P=[1,1]';
        Y = mci_nmm_gen(M,U,P);
        ind=[1:M.Np];
        
    otherwise
        disp('Unknown model');
        return
end
