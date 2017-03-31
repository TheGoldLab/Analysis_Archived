%% script for generating figures for LC Pupil
%
% created 10/10/14 by jig

%% Directory for dumping .eps figures
figdir = '../figs/eps';

%% FIG 1: MRI localization

%% FIG 2: Pupil stats
% Show basic properties of pupil diameter
figLCP_pupilStats(2)
print('-depsc',fullfile(figdir, 'fig02.eps'))

%% FIG 3: Tonic (across-trial) correlations
figLCP_tonic(3)
print('-depsc',fullfile(figdir, 'fig03.eps'))

%% FIG 4: Spike-triggered PD slope
figLCP_STA(4)
print('-depsc',fullfile(figdir, 'fig04.eps'))

%% FIG 5: PD event-locked PETH
figLCP_PETH(5)
print('-depsc',fullfile(figdir, 'fig05.eps'))

%% FIG 6: PD event-locked LFP

%% FIG 7: PD/spike responses to beeps
figLCP_pdBeep(7)
print('-depsc',fullfile(figdir, 'fig07.eps'))

%% FIG 8: PD/spike responses to ustim
figLCP_uStim(8)
print('-depsc',fullfile(figdir, 'fig08.eps'))



%%
%%%% OBSOLETE %%%%
% FIG 4: TIME/PHASE
%figLCP_phase(4)
%print('-depsc',fullfile(figdir, 'fig04.eps'))

