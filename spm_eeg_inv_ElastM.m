function ts = spm_eeg_inv_ElastM(ts);

%==========================================================================
% FORMAT ts = spm_eeg_inv_ElastM(ts);
%
% Modify the mesh in order to reduce overlong edges.
% The procedure uses an elastic model :
% At each vertex, the neighbouring triangles and vertices connected 
% directly are considered.
% Each edge is considered elastic and can be lengthened or shortened,
% depending on their length.
% Refs: G.Taubin, A signal processing approach to fair surface design, 1995
% This is a non-shrinking smoothing algorithm.
%
% Input : 
% ts         - tesselated surface
%
% Output :
% ts         - tesselated surface with corrected mesh
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips & Jeremie Mattout
% $Id: spm_eeg_inv_ElastM.m 719 2007-01-18 11:06:07Z christophe $

% Connection vertex-to-vertex
%--------------------------------------------------------------------------
M_con = sparse([ts.tri(1,:)';ts.tri(1,:)';ts.tri(2,:)';ts.tri(3,:)';ts.tri(2,:)';ts.tri(3,:)'], ...
               [ts.tri(2,:)';ts.tri(3,:)';ts.tri(1,:)';ts.tri(1,:)';ts.tri(3,:)';ts.tri(2,:)'], ...
               ones(ts.nr(2)*6,1),ts.nr(1),ts.nr(1));
                      
kpb   = .1;                       % Cutt-off frequency
lam   = .5; mu = lam/(lam*kpb-1); % Parameters for elasticity.
N     = 25;                       % Number of smoothing steps, the larger, the smoother
XYZmm = ts.XYZmm;

% smoothing iterations
%--------------------------------------------------------------------------
for j=1:N

    XYZmm_o = zeros(3,ts.nr(1)) ;
    XYZmm_o2 = zeros(3,ts.nr(1)) ;

    for i=1:ts.nr(1)
        ln = find(M_con(:,i));
        d_i = sqrt(sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).^2));
        w_i = d_i/sum(d_i);
        XYZmm_o(:,i) = XYZmm(:,i) + ...
            lam * sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).*(ones(3,1)*w_i),2);
	end

    for i=1:ts.nr(1)
        ln = find(M_con(:,i));
        d_i = sqrt(sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).^2));
        w_i = d_i/sum(d_i);
        XYZmm_o2(:,i) = XYZmm_o(:,i) + ...
            mu * sum((XYZmm_o(:,ln)-XYZmm_o(:,i)*ones(1,length(ln))).*(ones(3,1)*w_i),2);
    end

    XYZmm = XYZmm_o2;
    
end

ts.XYZmm = XYZmm;
