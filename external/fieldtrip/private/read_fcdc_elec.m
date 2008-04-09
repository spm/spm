function [sens] = read_fcdc_elec(varargin)

% READ_FCDC_ELEC read sensor positions from various manufacturers
% file formats. Currently supported are ASA, BESA, Polhemus and Matlab
% for EEG electrodes and CTF and Neuromag for MEG gradiometers.
%
% Use as
%   grad = read_fcdc_elec(filename)  % for gradiometers
%   elec = read_fcdc_elec(filename)  % for electrodes
%
% An electrode definition contain the following fields
%   elec.pnt     Nx3 matrix with carthesian (x,y,z) coordinates of each electrodes
%   elec.label   cell-array of length N with the label of each electrode
%
% A gradiometer definition generally consists of multiple coils per
% channel, e.g.two coils for a 1st order gradiometer in which the
% orientation of the coils is opposite. Each coil is described
% separately and a large "tra" matrix (can be sparse) has to be
% given that defines how the forward computed field is combined over
% the coils to generate the output of each channel. The gradiometer
% definition constsis of the following fields
%   grad.pnt     Mx3 matrix with the position of each coil
%   grad.ori     Mx3 matrix with the orientation of each coil
%   grad.tra     NxM matrix with the weight of each coil into each channel
%   grad.label   cell-array of length N with the label of each of the channels
%
% See also READ_FCDC_HEADER, READ_FCDC_DATA, READ_FCDC_EVENT

% Copyright (C) 2005-2008, Robert Oostenveld
%
% $Log: read_fcdc_elec.m,v $
% Revision 1.8  2008/03/05 10:46:36  roboos
% moved electrode reading functionality from read_fcdc_elec to read_sens, switched to the use of the new function
%
% Revision 1.7  2007/03/14 08:45:36  roboos
% also recognize ctf_ds and deal with it just as res4

% use the low-level reading function
[sens] = read_sens(varargin{:});