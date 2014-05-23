%        DYNAMIC CAUSAL MODELLING OF SMOOTH PURSUIT EYE MOVEMENTS
% 
% Overview:
% 
% This toolbox is a self-contained collection of routines and data that
% demonstrates the dynamic causal modelling of smooth pursuit eye movements
% (SPEM) - as measured using eye tracking. An illustration of model
% inversion is provided in the script spm_SEM_demo.m and the data (in
% DATA.mat) are described below.
%__________________________________________________________________________
%
% A full description of this approach and exemplar analysis can be found
% in:
% 
% Active inference and slow pursuit: the dynamic causal modelling of eye
% movements Rick A Adams and Karl J Friston
% 
% Keywords: oculomotor control, smooth pursuit, visual occlusion, active
% inference, schizophrenia, dynamic causal modelling, perception,
% precision
% 
% Abstract: This paper introduces a new paradigm that allows one to
% quantify the Bayesian beliefs evidenced by subjects during oculomotor
% pursuit. This paradigm uses non-invasive eye tracking responses to visual
% occlusion and rests on two innovations. The first is to treat eye
% tracking data in the same way that electrophysiological responses are
% averaged to form event related potentials. These response averages are
% then analysed using dynamic causal modelling (DCM). In DCM, observed
% responses are modelled using biologically plausible generative or forward
% models - usually biophysical models of neuronal activity. The second
% innovation is to use a generative model based on normative -
% Bayes-optimal - active inference. This allows us to model smooth pursuit
% eye movements in terms of a subject's beliefs about how visual targets
% move and how their oculomotor system should respond. Our aim here is to
% establish the face validity of the approach, using experimental
% manipulations of the content and precision of sensory information - and
% examining the ensuing changes in posterior beliefs. This combination of
% normative behavioural models and dynamic causal modelling features all
% the usual advantages of functionally grounded model comparison and
% quantitative parameter inference. In this application, the model
% parameters have an explicit interpretation in relation to beliefs about
% sensory exchanges with the world - and the confidence or expected
% precision associated with those beliefs. We hope to apply this paradigm
% to subjects with disorders like schizophrenia, to see if their responses
% to changes in the precision of sensory information differ from normal
% subjects.
%__________________________________________________________________________
%
% Data:
% 
% these data are the grand average over subjects of smooth pursuit eye
% movements - over one cycle of occluded pursuit. The averaged traces are
% in a structure called allsubj.
% 
% >> load DATA.mat
% 
% The structure's fields include:
% FN = Fast Noisy (i.e., 22 deg/sec)
% FS = Fast Smooth
% SN = Slow Noisy (i.e., 18 deg/sec)
% SS = Slow Smooth
% 
% These fields pertain to 4 conditions, under which data were acquired:
% these conditions conformed to a factorial design in which the target
% moved with a fast or slow speed, and moved in a noisy or smooth fashion.
% 
% Within these fields: '.one' contains the averaged eye position, measured
% in pixels (from 600 to -600) over several thousand milliseconds.
% 
% The target trajectories are in the 'target' structures target18.one and
% target22.one
% 
% The occluder was present between target18.plotocc.x(1,1) and (1,2), and
% between (2,1) and (2,2), in milliseconds.
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston and Rick Adams
% $Id: Contents.m 6014 2014-05-23 15:00:35Z guillaume $
