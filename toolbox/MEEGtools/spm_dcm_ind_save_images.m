function spm_dcm_ind_save_images(S)
% Save DCM-IR original and reconstructed induced responses, coupling
% matrices (A) and effects (B) as images for performing statistical
% analysis
%
% FORMAT spm_dcm_ind_save_images(S)
%
% S         - struct (optional)
%(optional) fields of S:
% S.DCM - DCM struct or name of a DCM mat file
% S.predictedTF - 1 - export predicted time-frequency matrices
% S.observedTF  - 1 - export observed time-frequency matrices
% S.couplingA   - 1 - export baseline coupling (A) matrices
% S.couplingB   - 1 - export changes in coupling (B) matrices
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL

% Vladimir Litvak, based on spm_dcm_ind_results
% $Id: spm_dcm_ind_save_images.m 2617 2009-01-19 18:49:16Z vladimir $


if nargin == 0
    S = [];
end

if ~isfield(S, 'DCM')
    S.DCM =  spm_select(1, '\.mat$', 'Select DCM-IR mat file');
end

if ~isfield(S, 'predictedTF'),  S.predictedTF = 1;  end
if ~isfield(S, 'observedTF'),   S.observedTF = 1;   end
if ~isfield(S, 'couplingA'),    S.couplingA = 1;    end
if ~isfield(S, 'couplingB'),    S.couplingB = 1;    end
if ~isfield(S, 'outdir'),       S.outdir = pwd;     end


if isstruct(S.DCM)
    DCM = S.DCM;
else
    DCM = load(S.DCM);
    DCM = DCM.DCM;
end


xY     = DCM.xY;
nt     = length(xY.y);           % Nr of trial types
nr     = length(DCM.C);          % Nr of sources
nu     = length(DCM.B);          % Nr of experimental effects
nf     = size(xY.U,2);           % Nr of frequency modes
ns     = size(xY.y{1},1);        % Nr of time bins
pst    = xY.pst;                 % peri-stmulus time
Hz     = xY.Hz;                  % frequencies
name   = spm_str_manip(DCM.name, 't');


V=[];
V.dt=[spm_type('float64') 0];
V.mat = eye(4);
V.pinfo = [1 0 0]';

if S.predictedTF || S.observedTF
    V.dim = [length(Hz) length(pst)  1 ];

    % reconstitute time-frequency and get principle model over channels
    %----------------------------------------------------------------------
    nk    = length(Hz);
    for i = 1:nt
        for j = 1:nr
            TF{i,j} = sparse(ns,nk);
            RF{i,j} = sparse(ns,nk);
        end
    end
    for i = 1:nt
        for j = 1:nr
            for k = 1:nf
                TF{i,j} = TF{i,j} + DCM.H{i,k}(:,j)*xY.U(:,k)';
                RF{i,j} = RF{i,j} + DCM.R{i,k}(:,j)*xY.U(:,k)';
            end
        end
    end


    % loop over trials, sources (predicted and observed)
    %----------------------------------------------------------------------
    for i = 1:nt
        for j = 1:nr
            if S.observedTF

                V.fname = fullfile(S.outdir,  ['observedTF_' DCM.Sname{j} '_' xY.code{i} '_' name '.img']);
                spm_write_vol(V, TF{i,j}' + RF{i,j}');

            end

            if S.predictedTF

                V.fname =  fullfile(S.outdir, ['predictedTF_' DCM.Sname{j} '_' xY.code{i} '_' name '.img']);
                spm_write_vol(V, TF{i,j}');

            end
        end
    end
end

if S.couplingA
    V.dim = [length(Hz) length(Hz)  1 ];
    % reconstitute time-frequency coupling
    %----------------------------------------------------------------------
    for i = 1:nr
        for j = 1:nr
            ii = [1:nf]*nr - nr + i;
            jj = [1:nf]*nr - nr + j;
            A  = xY.U*DCM.Ep.A(ii,jj)*xY.U';
            V.fname =  fullfile(S.outdir, ['couplingA_from_' DCM.Sname{i} '_to_' DCM.Sname{j} '_' name '.img']);
            spm_write_vol(V, A);
        end
    end
end

if S.couplingB
    V.dim = [length(Hz) length(Hz)  1 ];

    for k = 1:nu
        % reconstitute time-frequency  coupling
        %----------------------------------------------------------------------
        for i = 1:nr
            for j = 1:nr

                ii = [1:nf]*nr - nr + i;
                jj = [1:nf]*nr - nr + j;
                B  = xY.U*DCM.Ep.B{k}(ii,jj)*xY.U';

                V.fname =  fullfile(S.outdir, ['couplingB_' DCM.xU.name{k} '_from_' DCM.Sname{i} '_to_' DCM.Sname{j} '_' name '.img']);
                spm_write_vol(V, B);
            end
        end
    end
end
