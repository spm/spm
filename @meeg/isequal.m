function res = isequal(this, that)
% Method to check if 2 MEEG objects are the same
% FORMAT res = isequal(this, that)
%__________________________________________________________________________

% Christophe Phillips
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = isequal(struct(this), struct(that));
