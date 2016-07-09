% This function which initialises a struct containing the preferences and
% details for a particular project; in this case the demo for the toolbox.
% It should be edited to taste before a project is run, and a new one
% created for each substantially different project (though the options
% struct will be saved each time the project is run under a new name, so
% all will not be lost if you don't do this).
%
% For a guide to how to fill out the fields in this file, consult the
% documentation pages.
%
% Cai Wingfield 11-2009, 6-2010, 7-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

function userOptions = projectOptions_DEMO1()

%% Project details
% This name identifies a collection of files which all belong to the same run of a project.
userOptions.analysisName = 'DEMO1';

% This is the root directory of the project.
userOptions.rootPath = [pwd,filesep,'DEMO1'];

%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL SETUP %%
%%%%%%%%%%%%%%%%%%%%%%%%

% The list of subjects to be included in the study.
% userOptions.subjectNames = { ...
% 	'BE', ...
%     'KO', ...
%     'SN', ...
%     'TI' ...
% 	};

% The default colour label for RDMs corresponding to RoI masks (as opposed to models).
userOptions.RoIColor = [0 0 1];
userOptions.ModelColor = [0 1 0];

% Should information about the experimental design be automatically acquired from SPM metadata?
% If this option is set to true, the entries in userOptions.conditionLabels MUST correspond to the names of the conditions as specified in SPM.
userOptions.getSPMData = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS PREFERENCES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First-order analysis

% Text lables which may be attached to the conditions for MDS plots.
[userOptions.conditionLabels{1:92}] = deal(' ');
userOptions.conditionColours = [repmat([1 0 0], 48,1); repmat([0 0 1], 44,1)];
% userOptions.alternativeConditionLabels = { ...
% 	'1', ...
% 	'2', ...
% 	'3', ...
% 	'4', ...
% 	'5' ...
% 	};
% userOptions.useAlternativeConditionLabels = false;

% Which distance measure to use when calculating first-order RDMs.
userOptions.distance = 'Correlation';
userOptions.dotSize = 8;

%% Second-order analysis

% Which similarity-measure is used for the second-order comparison.
userOptions.distanceMeasure = 'Spearman';

% How many permutations should be used to test the significance of the fits?  (10,000 highly recommended.)
userOptions.significanceTestPermutations = 10000;

% Should RDMs' entries be rank transformed into [0,1] before they're displayed?
userOptions.rankTransform = true;

% Should rubber bands be shown on the MDS plot?
userOptions.rubberbands = true;

% What criterion shoud be minimised in MDS display?
userOptions.criterion = 'metricstress';

userOptions.resampleSubjects = false;
userOptions.resampleConditions = true;

% How should figures be outputted?
userOptions.displayFigures = true;
userOptions.saveFiguresPDF = true;
userOptions.saveFiguresFig = false;
userOptions.saveFiguresPS = false;

userOptions.forcePromptReply = 'r';

end%function
