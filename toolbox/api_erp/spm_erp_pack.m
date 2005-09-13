function [E] = spm_erp_pack(P,m,n)
% [un]packs a parameter vector into a structure for erp models
% FORMAT [P] = spm_erp_pack(P,m,n)
%
%  P - free parameter vector
%  m - number of inputs
%  n - number of sources
%___________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% %W% Karl Friston %E%

global M

o = M.l; % number of outputs


% number of outputs
%---------------------------------------------------------------------------
u            = 8;                                        % noise components
s            = 0;

% get intrinic parameters
%---------------------------------------------------------------------------
E.T          = sparse(n,1);   l = length(E.T(:));
E.T(:)       = P([1:l] + s);  s = s + l;
E.H          = sparse(n,1);   l = length(E.H(:));
E.H(:)       = P([1:l] + s);  s = s + l;

% get observer parameters
%---------------------------------------------------------------------------
E.K          = sparse(n,1);   l = length(E.K(:));
E.K(:)       = P([1:l] + s);  s = s + l;

if isfield(M, 'dipfit')
    % parameterised leadfield
    E.Lpos          = sparse(3,n);   l = length(E.Lpos(:));
    E.Lpos(:)       = P([1:l] + s);  s = s + l;
    E.Lmom          = sparse(3,n);   l = length(E.Lmom(:));
    E.Lmom(:)       = P([1:l] + s);  s = s + l;
else
    % static leadfield
    E.L          = sparse(o,n);   l = length(E.L(:));
    E.L(:)       = P([1:l] + s);  s = s + l;
end

% get extrinsic connectivity
%---------------------------------------------------------------------------
for i = 1:3
    E.A{i}    = sparse(n,n);   l = length(E.A{i}(:));
    E.A{i}(:) = P([1:l] + s);  s = s + l;
end
for i = 1:m
    E.B{i}    = sparse(n,n);   l = length(E.B{i}(:));
    E.B{i}(:) = P([1:l] + s);  s = s + l;
end
E.C           = sparse(n,m);   l = length(E.C(:));
E.C(:)        = P([1:l] + s);  s = s + l;

% get conduction delays
%---------------------------------------------------------------------------
E.D           = sparse(n,n);   l = length(E.D(:));
E.D(:)        = P([1:l] + s);  s = s + l;

% get stimulus parameters
%---------------------------------------------------------------------------
E.R           = sparse(1,2);   l = length(E.R(:));
E.R(:)        = P([1:l] + s);  s = s + l;
E.N           = sparse(u,1);   l = length(E.N(:));
E.N(:)        = P([1:l] + s);  s = s + l;
E.U           = P(1 + s);      s = s + 1;
