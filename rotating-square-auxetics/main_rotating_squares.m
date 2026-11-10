%% MAIN_ROTATING_SQUARES  Rotating-squares auxetic: analysis driver
% Run this file (F5), or run one section at a time (Ctrl+Enter).
% All inputs come from params_default.m - edit numbers there, not here.
%
% Stages:
%   1  Ideal kinematics (rigid squares, pin hinges)        -> stage1_run
%   2  Compliant hinges, effective stiffness               -> stage2_run
%   3  Own beam-element FE check                           -> stage3_run
%   4  Failure: hinge stress, yield, locking, buckling     -> stage4_run
%   5  Animation, slide figures, slider app                -> stage5_run,
%                                                             rs_explorer_app
%   6  CSV exports for SolidWorks / Ansys, print check     -> stage6_run,
%                                                             compare_with_ansys

%% 0. Setup
clear; close all; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(genpath(fullfile(here, 'src')));

P = params_default();
make_output_folders(P);

%% 1. Stage 1 - ideal rigid-square kinematics
S1 = stage1_run(P);

%% 2. Stage 2 - compliant hinges (energy method, sweeps)
S2 = stage2_run(P);

%% 3. Stage 3 - own frame FE solver (independent check)
S3 = stage3_run(P, S2);

%% 4. Stage 4 - failure: yield, locking, buckling, safe range
S4 = stage4_run(P, S2);

%% 5. Stage 5 - animation, slide stills
S5 = stage5_run(P, S2, S4);
% Interactive slider app (needs a display, so run it from the desktop):
%   rs_explorer_app(P, S4.edge_factor)

%% 6. Stage 6 - CAD/print exports, results CSV, Ansys template
S6 = stage6_run(P, S2, S3, S4);
% After you fill results/data/ansys_results.csv, run:
%   compare_with_ansys

%% Summary
fprintf('\n=== Summary: all stages passed = %d ===\n', ...
    all([S1.pass S2.pass S3.pass S4.pass S5.pass S6.pass]));
