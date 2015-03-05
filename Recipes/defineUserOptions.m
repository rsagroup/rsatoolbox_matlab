function userOptions = defineUserOptions()
%
%  defineUserOptions is a nullary function which initialises a struct
%  containing the preferences and details for a particular project.
%  It should be edited to taste before a project is run, and a new
%  one created for each substantially different project (though the
%  options struct will be saved each time the project is run under
%  a new name, so all will not be lost if you don't do this).
%
%  For a guide to how to fill out the fields in this file, consult
%  the documentation folder (particularly the userOptions_guide.m)
%
%  Cai Wingfield 11-2009
%  Li Su updated 2-2012
%  Fawad updated 12-2013, 10-2014
%  Jana updated 10-2014
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

%% Project details

% This name identifies a collection of files which all belong to the same run of a project.
userOptions.analysisName = 'yourProjectName';

% This is the root directory of the project.
userOptions.rootPath = 'pathToRootDirectoryOfProject';

% The path leading to where the scans are stored (not including subject-specific identifiers).
% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
% "[[betaIdentifier]]" should be used as a placeholder to denote an output of betaCorrespondence.m if SPM is not being used; or an arbitrary filename if SPM is being used.
userOptions.betaPath = 'pathToYourSingleConditionResponses';% e.g. /imaging/mb01/lexpro/multivariate/ffx_simple/[[subjectName]]/[[betaIdentifier]]

% if set true, intermediate files will be saved
userOptions.debug = false;

% Regularization based on paper by Diedrichson et al. 2011.
userOptions.regularized = false;

%%%%%%%%%%%%%%%%%%%
%% Email Options %%
%%%%%%%%%%%%%%%%%%%
% Set to true to be informed when a script finishes.
userOptions.recieveEmail = true;
% Put your address to make 
userOptions.mailto = 'cw417@cam.ac.uk';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parallel Computing toolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run parallel locally set run_in_parallel *and*
% run_in_parallel_in_cluster to false. 
% To run on CBU adaptive queue set run_in_parallel_in_cluster to true. 
% Do NOT set this true for fixed effect analysis (searchlight and
% permutation.
userOptions.run_in_parallel = true;
userOptions.run_in_parallel_in_cluster= true;
% Sometimes the performance will drop if there are large number of tiny
% jobs due to the communication and setting-up overhead.
userOptions.jobSize = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adaptive Computing Cluster Queueing OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set true to delete jobs from the queue otherwise set to false.
% If you do not want to delete all jobs but only sepecific one do not set
% this variable to true.
userOptions.flush_Queue = false; 
% Used only when using CBU cluster.
% i.e. when run_in_parallel_in_cluster = true;
userOptions.wallTime = '12:00:00';
% Cluster machines requested.
userOptions.nodesReq = 8;
% Processors requested per processor machine.
userOptions.proPNode = 1;
% The product of nodesReq and proPNode should be greater or equal to the
% number of workers requested.
userOptions.nWorkers = 8;
% In gigabytes, to be distributed amongst all nodes.
userOptions.memReq = 400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEG SENSOR-LEVEL ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set this true for both searchlight and fixed time window RoI analyses
userOptions.sensorLevelAnalysis = true;

% ===== use searchlight settings for doing sensor level searchlight ==== %

% == use following settings for sensor level RoI analysis only========== %

% specify a name for your mask, used in saving your results
userOptions.maskSpec.maskName = {'allSensors'};

% remove/add sensors in your mask by setting them true/false
userOptions.maskSpec.MEGSensors.Gradiometers = true;
userOptions.maskSpec.MEGSensors.Magnetometers = true;
userOptions.maskSpec.MEGSensorSites = (1:102); % for each mag and grad: create above mentioned 
            % mask by putting sensor numbers in arrays put any values as [90 98 100 1 4];

userOptions.maskSpec.MEGSensors.EEG = false;
userOptions.maskSpec.EEGSensorSites = (1:70); % eeg

% time window for RoI analysis
userOptions.maskSpec.timeWindow = [-200 100];

% pattern of analysis 'spatial', 'temporal' or 'spatiotemporal'
userOptions.maskSpec.patternType =  'Spatial';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEATUERS OF INTEREST SELECTION OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
userOptions.maskingFlag = true;
userOptions.slidingTimeWindow = true; % will use temporalSearchlightWidth & resolution

	%% %% %% %% %%
    %% fMRI/MEG %% Use these next three options if you're working in fMRI native space or MEG source space:
	%% %% %% %% %%
	
	% The path to a stereotypical mask data file is stored (not including subject-specific identifiers).
	% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
	% "[[maskName]]" should be used as a placeholder to denote an entry in userOptions.maskNames
	userOptions.maskPath = 'pathToWhereYourMasksAreStored';%'/imaging/mb01/lexpro/multivariate/ffx_simple/[[subjectName]]/[[maskName]].img';
		
		% The list of mask filenames (minus .hdr extension) to be used.
        % For MEG, names should be in pairs, such as maskName-lh, maskName-rh
		userOptions.maskNames = { ...
			'mask',...
			};
			
		% time windows are specified for each regions. For ROI analysis, different 
		% time windows can be set for different ROIs, but it can also be used for 
		% searchlight analysis to reduce the search space, in searchlight analysis,
		% time window should be same for all ROIs.
		userOptions.maskTimeWindows = {
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500], ...
		    [0 500],[0 500]
		    };

%%%%%%%%%%%%%%%%%%%%%%%%%
%% SEARCHLIGHT OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% set this true for sensor and source level searchlight
userOptions.searchlight = false; 
% == Source Level Note: if mask set to true, will use 
% temporalsearchlightLimits rather than maskTimeWindows for searchlight== %
% == No masking available for sensor level searchlight. Works only on GRADIOMETERS == %

	%% %% %% %% %%
	%% fMRI  %% Use these next three options if you're working in fMRI native space:
	%% %% %% %% %%

		% What is the path to the anatomical (structural) fMRI scans for each subject?
		% "[[subjectName]]" should be used to denote an entry in userOptions.subjectNames
		userOptions.structuralsPath = 'paathToWhereYourSubject''s structuralImagesAreStored ';% e.g. /imaging/mb01/lexpro/[[subjectName]]/structurals/
	
		% What are the dimensions (in mm) of the voxels in the scans?
		userOptions.voxelSize = [3 3 3.75];
	
		% What radius of searchlight should be used (mm)?
		userOptions.searchlightRadius = 15;
		
		% Correlate over space ('spatial') or regularized ('regularized')
		userOptions.searchlightPatterns = 'spatial';

%% %% %% %% %%
%%  MEG  %% Use these next four options if you're working in MEG:
%% %% %% %% %%

% ================ common options for sensor and source ==================%
% all time values in ms

% The width of the sliding window (ms)
userOptions.temporalSearchlightWidth = 20; %20;

% Correlate over space ('spatial'),  time ('temporal') or
% spatiotemporal ('spatiotemporal') (or regularized ('regularized') for
% source level)
userOptions.searchlightPatterns = 'spatial';

% The timestep for sliding window (ms)
userOptions.temporalSearchlightResolution = 10; %10; % (data point equivalent = total_in_ms/total_dataPoints)

% The overall window of interest (ms)
userOptions.temporalSearchlightLimits = [-200 800];
    
    % ===========================  (sensor space) ========================%

% The radius of the sensor searchlight (in adjacent sensors, excluding the centre: 0 => 1 sensor, 1 => ~9)
userOptions.sensorSearchlightRadius = 1;

    % ========================= OR (source space) ========================%

% The radius of the source-space searchlight (in mm)
userOptions.sourceSearchlightRadius = 20;

% The average surface files
userOptions.averageSurfaceFile = '/imaging/cw03/decom2/subjects/average/surf/lh.inflated';


%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL SETUP %%
%%%%%%%%%%%%%%%%%%%%%%%%

% The list of subjects to be included in the study.
userOptions.subjectNames = { ...
	'subject1','subject2',...
	};% eg CBUXXXXX

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

% temporal downsampling
userOptions.temporalDownsampleRate = 1;
userOptions.temporalSearchlightWidth = userOptions.temporalSearchlightWidth / userOptions.temporalDownsampleRate;
userOptions.temporalSearchlightResolution = userOptions.temporalSearchlightResolution / userOptions.temporalDownsampleRate;
userOptions.temporalSearchlightLimits = userOptions.temporalSearchlightLimits / userOptions.temporalDownsampleRate;

% Text lables which may be attached to the conditions for MDS plots.
[userOptions.conditionLabels{1:92}] = deal(' ');
% userOptions.alternativeConditionLabels = { ...
% 	' ', ...
% 	' ', ...
% 	' ', ...
% 	' ', ...
% 	' ' ...
% 	};
% userOptions.useAlternativeConditionLabels = false;

% What colours should be given to the conditions?
userOptions.conditionColours = [repmat([1 0 0], 48,1); repmat([0 0 1], 44,1)];

% Which distance measure to use when calculating first-order RDMs.
userOptions.distance = 'Correlation';

%% Second-order analysis

% % Which model RDM to test? (This number corresponds to the order in variable Models, which is specified in ModelRDMs.m)
userOptions.modelNumber = 1; % this is default value overwritten by which_model variable when Recipe is used
userOptions.partial_correlation = false;
userOptions.partial_modelNumber = {5,7}; % all models listed here will be partialed out from the original model

% Which similarity-measure is used for the second-order comparison.
userOptions.distanceMeasure = 'Spearman';

% How many permutations should be used to test the significance of the fits?  (10,000 highly recommended.)
userOptions.significanceTestPermutations = 10000;

% Bootstrap options
userOptions.nResamplings = 1000;
userOptions.resampleSubjects = true;
userOptions.resampleConditions = false;

userOptions.fisher = true;

% Group statistics options: random effect ('RFX') or fixed effect ('FFX')
userOptions.groupStats = 'RFX';
% if set true, group stats will be computed using t maps, if false then r
% maps will be used.
userOptions.tmap = true;

% Clustering analysis: primary cluster-forming threshold in terms of the
% top x% for all vertexes across space and time.
% It means that primary threshold is set to p<x ,one tailed for positive
% clusters. This doesnot need to be adjusted according to the degree of
% freedom of data. If the requirement is to select top ONLY x% of all data, set
% this value to nan.
% If not using tmaps for RFX, this value only will define the primary threshold.
% DoF does not help compute threshold for that case.
userOptions.primaryThreshold = 0.05;

% spatial smoothing and upsampling parameters
userOptions.targetResolution = 10242;
userOptions.minDist = 5; % 5mm is the smallest distance between two adjacent vertex in 10242 resolution.
% 10mm is the smallest distance between two adjacent vertex in 2562 resolution.
userOptions.smoothingWidth = 10; %mm

% Should RDMs' entries be rank transformed into [0,1] before they're displayed?
userOptions.rankTransform = true;

% Should rubber bands be shown on the MDS plot?
userOptions.rubberbands = true;

% What criterion shoud be minimised in MDS display?
userOptions.criterion = 'metricstress';

% What is the colourscheme for the RDMs?
userOptions.colourScheme = jet(128); 

% Set any of the following to true to delete files saved as work-in-progress
% throughout the analysis. This will save space.
userOptions.deleteTMaps_Dir= false;
userOptions.deleteImageData_Dir= false;
userOptions.deletePerm= true;

% How should figures be outputted?
userOptions.displayFigures = true;
userOptions.saveFiguresPDF = false;
userOptions.saveFiguresFig = false;
userOptions.saveFiguresPS = false;
% Which dots per inch resolution do we output?
userOptions.dpi = 300;
% Remove whitespace from PDF/PS files?
% Bad if you just want to print the figures since they'll
% no longer be A4 size, good if you want to put the figure
% in a manuscript or presentation.
userOptions.tightInset = false;

userOptions.forcePromptReply = 'r';

% Present user with graphical feedback?
userOptions.dialogueBox = false;


end%function
