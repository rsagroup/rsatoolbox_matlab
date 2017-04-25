function userOptions = defineUserOptions()
%
%  projectOptions is a nullary function which initialises a struct
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
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

%% Project details

% This name identifies a collection of files which all belong to the same run of a project.
userOptions.projectName  = 'iRSA_fMRI';

% This name identifies a collection of files which all belong to the same analysis within a project.
userOptions.analysisName = 'vanilla';

% This is the root directory of the project.
userOptions.rootPath = '/home/adf/charesti/Documents/iRSA_fMRI/mridata2015';

% The path leading to where the scans are stored (not including subject-specific identifiers).
% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
% "[[betaIdentifier]]" should be used as a placeholder to denote an output of betaCorrespondence.m if SPM is not being used; or an arbitrary filename if SPM is being used.
userOptions.betaPath = fullfile(userOptions.rootPath,'[[subjectName]]','spmmodel','[[betaIdentifier]]');% e.g. /imaging/mb01/lexpro/multivariate/ffx_simple/[[subjectName]]/[[betaIdentifier]]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEATUERS OF INTEREST SELECTION OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%% %% %% %% %%
	%% fMRI  %% Use these next three options if you're working in fMRI native space:
	%% %% %% %% %%
	
	% The path to a stereotypical mask data file is stored (not including subject-specific identifiers).
	% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
	% "[[maskName]]" should be used as a placeholder to denote an entry in userOptions.maskNames
	userOptions.maskPath = fullfile(userOptions.rootPath,'[[subjectName]]','masks','[[maskName]].nii');%'/imaging/mb01/lexpro/multivariate/ffx_simple/[[subjectName]]/[[maskName]].img';
		
		% The list of mask filenames (minus .hdr extension) to be used.
		userOptions.maskNames = { ...
			'hit',...
            'evc',...
			};

%%%%%%%%%%%%%%%%%%%%%%%%%
%% SEARCHLIGHT OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%

	%% %% %% %% %%
	%% fMRI  %% Use these next three options if you're working in fMRI native space:
	%% %% %% %% %%

		% What is the path to the anatomical (structural) fMRI scans for each subject?
		% "[[subjectName]]" should be used to denote an entry in userOptions.subjectNames
		userOptions.structuralsPath = fullfile(userOptions.rootPath,'[[subjectName]]','structurals');% e.g. /imaging/mb01/lexpro/[[subjectName]]/structurals/
	
		% What are the dimensions (in mm) of the voxels in the scans?
		userOptions.voxelSize = [3 3 3.75];
	
		% What radius of searchlight should be used (mm)?
		userOptions.searchlightRadius = 15;
	
%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL SETUP %%
%%%%%%%%%%%%%%%%%%%%%%%%

% The list of subjects to be included in the study.
userOptions.subjectNames = { ...
	'CBU101295',...
	};% eg CBUXXXXX

% The default colour label for RDMs corresponding to RoI masks (as opposed to models).
userOptions.RoIColor = [0 0 1];
userOptions.ModelColor = [0 1 0];

% Should information about the experimental design be automatically acquired from SPM metadata?
% If this option is set to true, the entries in userOptions.conditionLabels MUST correspond to the names of the conditions as specified in SPM.
userOptions.getSPMData = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS PREFERENCES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First-order analysis

userOptions.nConditions = 72;
% Text lables which may be attached to the conditions for MDS plots.
[userOptions.conditionLabels{1:userOptions.nConditions}] = deal(' ');
% userOptions.alternativeConditionLabels = { ...
% 	' ', ...
% 	' ', ...
% 	' ', ...
% 	' ', ...
% 	' ' ...
% 	};
% userOptions.useAlternativeConditionLabels = false;

% What colours should be given to the conditions?
OwnBodyParts = 1:3; OwnFaces = 4:8; OwnPet = 9; OwnPlaces = 10:15; OwnObjects=16:18;
OtherBodyParts = 19:21; OtherFaces = 22:26; OtherPet = 27; OtherPlaces = 28:33; OtherObjects=34:36;
GeneralBodyParts = 37:44; GeneralFaces = 45:52; GeneralPets = [53 54]; GeneralPlaces=55:66; GeneralObjects=67:72;

animates = [OwnBodyParts OwnFaces OwnPet OtherBodyParts OtherFaces OtherPet GeneralBodyParts GeneralFaces GeneralPets];
inanimates = [OwnPlaces OwnObjects OtherPlaces OtherObjects GeneralPlaces GeneralObjects];%setdiff(1:nconditions,animates)

userOptions.conditionColours = zeros(userOptions.nConditions,3);
userOptions.conditionColours(animates,:) = repmat([1 0 0],length(animates),1);
userOptions.conditionColours(inanimates,:) = repmat([0 0 1],length(inanimates),1);

% Which distance measure to use when calculating first-order RDMs.
userOptions.distance = 'Correlation';

%% Second-order analysis

% Which similarity-measure is used for the second-order comparison.
userOptions.distanceMeasure = 'Spearman';

% How many permutations should be used to test the significance of the fits?  (10,000 highly recommended.)
userOptions.significanceTestPermutations = 10000;

% Bootstrap options
userOptions.nResamplings = 1000;
userOptions.resampleSubjects = true;
userOptions.resampleConditions = false;

% Should RDMs' entries be rank transformed into [0,1] before they're displayed?
userOptions.rankTransform = true;

% Should rubber bands be shown on the MDS plot?
userOptions.rubberbands = true;

% What criterion shoud be minimised in MDS display?
userOptions.criterion = 'metricstress';

% What is the colourscheme for the RDMs?
userOptions.colourScheme = bone(128);

% How should figures be outputted?
userOptions.displayFigures = true;
userOptions.saveFiguresPDF = true;
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

end%function
