function [P,M,U,Y,ind] = mci_compare_setup (model)
% Set up data structures for fwd/sens/grad comparisons
% FORMAT [P,M,U,Y,ind] = mci_compare_setup (model)
%
% model     'phase', 'nmm-r2p2'
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

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
