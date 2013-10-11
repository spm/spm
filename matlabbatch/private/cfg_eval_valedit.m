function [val, sts] = cfg_eval_valedit(varargin)
% Helper function to evaluate GUI inputs in MATLAB workspace
% FORMAT [val, sts] = cfg_eval_valedit(str)
% Evaluates GUI inputs in MATLAB 'base' workspace. Results are returned
% in val. Expressions in str can be either a pure rhs argument, or a set
% of commands that assign to a workspace variable named 'val'. If
% unsuccessful, sts is false and a message window is displayed. 
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_eval_valedit.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; 

val = [];
sts = false;
try
    % 1st, try to convert into numeric matrix without evaluation
    % This converts expressions like '1 -1' into [1 -1] instead of
    % evaluating them
    [val, sts] = str2num(varargin{1}); %#ok<ST2NM>
    if ~sts
        % try to evaluate str as rvalue
        val = evalin('base', varargin{1});
        sts = true;
    end
catch
    everr = lasterror;
    if strcmp(everr.identifier, 'MATLAB:m_invalid_lhs_of_assignment')
        try
            evalin('base', varargin{1});
            val = evalin('base','val');
            % test if val variable exists
            if ~exist('val','var')
                cfg_message('cfg_ui:local_eval_valedit:noval','No variable ''val'' assigned.');
            end;
            sts = true;
        catch
            sts = false;
            val = [];
            everr = lasterror;
            msgbox(everr.message,'Evaluation error','modal');
        end;
    end;
end;
